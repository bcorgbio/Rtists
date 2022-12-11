library(rgbif)
library(tidyverse)
library(MuMIn)
library(rnoaa)
library(data.table)
library(ggmap)
library(usmap)
library(magick)#for examples
library(cowplot)#for examples
library(lme4) #for linear mixed effect models
library(car) #for LME anova testing

raven <- occ_data(scientificName = "Corvus corax", stateProvince="Maine", limit=200,year=2018)
ME<- map_data('state', 'maine')
raven.p <- ggplot(ME, aes(long,lat,group=subregion) )+
  geom_polygon(colour = "gray",fill="gray90")+geom_point(data=raven[[2]],aes(x=decimalLongitude,y=decimalLatitude,size=individualCount),alpha=0.3,inherit.aes = F)+ coord_quickmap()+theme_void()

#and add image to the map with cowplot
raven.p2 <- ggdraw() +
  draw_image("https://www.allaboutbirds.org/guide/assets/photo/63739491-480px.jpg",scale = 0.3,halign=0,valign=1) +
  draw_plot(raven.p)
print(raven.p2)

names(raven)
lapply(raven,head)

species <- c("Megaceryle alcyon","Sayornis phoebe","Vireo philadelphicus","Catharus ustulatus","Setophaga caerulescens","Setophaga dominica")

y <- paste0("1990",",","2019")
m <- paste0("3",",","6")

?occ_data


dat.l <- list()


for (i in species){
  n.obs <- occ_data(scientificName = i,year = y,month = m,limit=0,country = "US",basisOfRecord = "HUMAN_OBSERVATION", stateProvince = "Massachusetts")$meta$count
  print(n.obs)
  dat.l[[paste0(i)]] <- occ_data(scientificName = i, year=y, month = m, limit = n.obs,country = "US",basisOfRecord = "HUMAN OBSERVATION",stateProvince = "Massachusetts")[[2]]
}

data = rbindlist(dat.l,fill = T)
head(data)

saveRDS(data,"massbird.data.RDS")
dat<- readRDS("massbird.data.RDS")


dat%>%
  group_by(year,species)%>%
  summarise(count=sum(individualCount,na.rm = T))%>%
  ggplot(aes(x=year,y=count,col=species))+geom_point()


options(noaakey="OOZIuxBSmRgSCBNQPuYyWVUZGaWoDIBo")

sts <- c(
  "GHCND:USW00013894", #Mobible, AL 2k away about 10 days away @200 km/day
  "GHCND:USW00013881", #Charlotte, NC 1000 km away about 6 days away @200 km/day
  "GHCND:USW00014739" #Boston
)

bos <- ncdc_stations(stationid = "GHCND:USW00014739")

print(bos)

sta.d <- bind_rows( #bind the rows
  lapply(sts,function(x) ncdc_stations(stationid = x)$data ) #use lapply to run through stations
)%>%
  mutate(usmap_transform(.,input_names = c("longitude","latitude"),output_names = c("longitude.1", "latitude.1")))%>% #join transformation of lat/long for projection with usmap
  mutate(name=str_sub(name, -5,-4))%>%#simplify the name column, grab just the state
  mutate(migr.day=c(10,5,0))%>% #so we can look at wind speed 0, 5 or 10 days before arrive in boston
  separate(id,into = c("station.type","id"))%>%#need to cut station type out from station id number
  print()

plot_usmap(include = c(.northeast_region,.south_region,.east_north_central))+geom_point(data=sta.d,aes(x=longitude.1,y=latitude.1,col=name),size=5)+geom_label(data=sta.d,aes(x=longitude.1,y=latitude.1,col=name,label=name),size=5,nudge_x = 1e6*0.25)+theme(legend.position = "none")

weather.d <- meteo_pull_monitors(sta.d$id,date_min = "2000-01-01")

head(weather.d)

#all data 
all.dat<- dat%>%
  group_by(year)%>%
  mutate(date=as.Date(paste0(year,"-",month,"-",day)),
         j.day=julian(date,origin=as.Date(paste0(unique(year),"-01-01")))
  )%>%
  group_by(species,year,j.day,date)%>%
  summarise(day.tot=sum(individualCount,na.rm=T))%>%
  group_by(species,year)%>%
  mutate(prop=cumsum(day.tot/sum(day.tot,na.rm = T)))%>%
  filter(year>1999)
all.dat%>%
  ggplot(aes(j.day,prop,color=species))+geom_point()+facet_wrap(year~.)

#all data log prediction 
all.dat.pred <- all.dat%>%
  group_by(year)%>%
  summarize(
    pred=predict(nls(prop~SSlogis(j.day,Asym, xmid, scal)),newdata=data.frame(j.day=min(j.day):max(j.day))),#predict the logistic curve for each species
    j.day=min(j.day):max(j.day),
  )%>%
  left_join(all.dat%>%dplyr::select(j.day,date)) ## add date back to tibble
all.dat%>%
  ggplot(aes(j.day,prop,color = species))+geom_point(aes=0.3)+geom_line(data=all.dat.pred,aes(x=j.day,y=pred),col="black",size=2)+facet_wrap(year~.)


#all data arrival date 
alldat.arrive.date <-all.dat.pred%>%
  group_by(year)%>%
  filter(j.day==j.day[which.min(abs(pred-0.25))])

alldat.arrive.date%>%
  ggplot(aes(year,j.day,color=species))+geom_point()

#all data with unaveraged weather 
alldat.arr.weath <- alldat.arrive.date%>%
  left_join(weather.d)%>%
  left_join(all.dat%>%dplyr::select(year,date,j.day))

#all data with 2 week weather average 
alldat.arr.weath2 <- alldat.arrive.date%>%
  left_join(weather.wk)


#all data linear models 
#linear model 1 -- 
alldat.lmer <- lmer(j.day~tmin*tmax*wvec+(1|name)+(1|species),alldat.arr.weath,na.action = "na.fail")
Anova(alldat.lmer)

#0Mean two week weather preceding arrival
alldat.lmer2 <- lmer(j.day~wk.tmin*wk.tmax*wk.wvec+(1|name)+(1|species),alldat.arr.weath2,na.action = "na.fail")
Anova(alldat.lmer2)

#secondary model 
alldat.arr.aic <- dredge(alldat.lmer2,fixed = c("wk.tmin","wk.tmax","wk.wvec"),)
pv.kb <- kable(pv.arr.aic[1:4,],caption = "Fit values for nested models of the most complicated lme model")
kable_styling(pv.kb)

#all data best model 
best.alldat.lmer <-  lmer(j.day~wk.tmin+wk.tmax+wk.wvec+(1|name)+(1|species),alldat.arr.weath2,na.action = "na.fail")
Anova(best.alldat.lmer)

#Megaceryle alcyon - belted kingfisher
mea<- dat%>%
  filter(species=="Megaceryle alcyon")%>%
  group_by(year)%>%
  mutate(date=as.Date(paste0(year,"-",month,"-",day)),
         j.day=julian(date,origin=as.Date(paste0(unique(year),"-01-01")))
  )%>%
  group_by(species,year,j.day,date)%>%
  summarise(day.tot=sum(individualCount,na.rm=T))%>%
  group_by(species,year)%>%
  mutate(prop=cumsum(day.tot/sum(day.tot,na.rm = T)))%>%
  filter(year>1999)

mea%>%
  ggplot(aes(j.day,prop))+geom_point()+facet_wrap(year~.)

#Sayornis phoebe - Eastern Phoebe 
ep<- dat%>%
  filter(species=="Sayornis phoebe")%>%
  group_by(year)%>%
  mutate(date=as.Date(paste0(year,"-",month,"-",day)),
         j.day=julian(date,origin=as.Date(paste0(unique(year),"-01-01")))
  )%>%
  group_by(species,year,j.day,date)%>%
  summarise(day.tot=sum(individualCount,na.rm=T))%>%
  group_by(species,year)%>%
  mutate(prop=cumsum(day.tot/sum(day.tot,na.rm = T)))%>%
  filter(year>1999)


ep%>%
  ggplot(aes(j.day,prop))+geom_point()+facet_wrap(year~.)

#Vireo philadelphicus - Philadelphia Virelo 
pv<- dat%>%
  filter(species=="Vireo philadelphicus")%>%
  group_by(year)%>%
  mutate(date=as.Date(paste0(year,"-",month,"-",day)),
         j.day=julian(date,origin=as.Date(paste0(unique(year),"-01-01")))
  )%>%
  group_by(species,year,j.day,date)%>%
  summarise(day.tot=sum(individualCount,na.rm=T))%>%
  group_by(species,year)%>%
  mutate(prop=cumsum(day.tot/sum(day.tot,na.rm = T)))%>%
  filter(year>1999)


pv%>%
  ggplot(aes(j.day,prop))+geom_point()+facet_wrap(year~.)

#Catharus ustulatus - Swainson's Thrush 
st <-  dat%>%
  filter(species=="Catharus ustulatus")%>%
  group_by(year)%>%
  mutate(date=as.Date(paste0(year,"-",month,"-",day)),
         j.day=julian(date,origin=as.Date(paste0(unique(year),"-01-01")))
  )%>%
  group_by(species,year,j.day,date)%>%
  summarise(day.tot=sum(individualCount,na.rm=T))%>%
  group_by(species,year)%>%
  mutate(prop=cumsum(day.tot/sum(day.tot,na.rm = T)))%>%
  filter(year>1999)


st%>%
  ggplot(aes(j.day,prop))+geom_point()+facet_wrap(year~.)


#Setophaga caerulescens - Black-throated Blue Warbler 
btw <-  dat%>%
  filter(species=="Setophaga caerulescens")%>%
  group_by(year)%>%
  mutate(date=as.Date(paste0(year,"-",month,"-",day)),
         j.day=julian(date,origin=as.Date(paste0(unique(year),"-01-01")))
  )%>%
  group_by(species,year,j.day,date)%>%
  summarise(day.tot=sum(individualCount,na.rm=T))%>%
  group_by(species,year)%>%
  mutate(prop=cumsum(day.tot/sum(day.tot,na.rm = T)))%>%
  filter(year>1999)


btw%>%
  ggplot(aes(j.day,prop))+geom_point()+facet_wrap(year~.)

#Setophaga dominica - Yellow-throated  Warbler 
yw <-  dat%>%
  filter(species=="Setophaga dominica")%>%
  group_by(year)%>%
  mutate(date=as.Date(paste0(year,"-",month,"-",day)),
         j.day=julian(date,origin=as.Date(paste0(unique(year),"-01-01")))
  )%>%
  group_by(species,year,j.day,date)%>%
  summarise(day.tot=sum(individualCount,na.rm=T))%>%
  group_by(species,year)%>%
  mutate(prop=cumsum(day.tot/sum(day.tot,na.rm = T)))%>%
  filter(year>1999)


yw%>%
  ggplot(aes(j.day,prop))+geom_point()+facet_wrap(year~.)

#logistic predictions for each species 

#Megaceryle alcyon - belted kingfisher
mea.pred <- mea%>%
  group_by(year)%>%
  summarize(
    pred=predict(nls(prop~SSlogis(j.day,Asym, xmid, scal)),newdata=data.frame(j.day=min(j.day):max(j.day))),#predict the logistic curve for each species
    j.day=min(j.day):max(j.day),
  )%>%
  left_join(mea%>%dplyr::select(j.day,date)) ## add date back to tibble

mea%>%
  ggplot(aes(j.day,prop))+geom_point(aes=0.3)+geom_line(data=mea.pred,aes(x=j.day,y=pred),col="steelblue",size=2)+facet_wrap(year~.)


#Sayornis phoebe - Eastern Phoebe 
ep.pred <- ep%>%
  group_by(year)%>%
  summarize(
    pred=predict(nls(prop~SSlogis(j.day,Asym, xmid, scal)),newdata=data.frame(j.day=min(j.day):max(j.day))),#predict the logistic curve for each species
    j.day=min(j.day):max(j.day),
  )%>%
  left_join(ep%>%dplyr::select(j.day,date)) ## add date back to tibble

ep%>%
  ggplot(aes(j.day,prop))+geom_point(aes=0.3)+geom_line(data=ep.pred,aes(x=j.day,y=pred),col="navajowhite4",size=2)+facet_wrap(year~.)

#Vireo philadelphicus - Philadelphia Virelo 
pv.forlog <- pv %>%  filter(year==2001 | year == 2002 | year == 2005|year ==2007|year==2011|year==2017|year==2018|year==2019)
pv.pred <- pv.forlog%>%
  group_by(year)%>%
  summarize(
    pred=predict(nls(prop~SSlogis(j.day,Asym, xmid, scal)),newdata=data.frame(j.day=min(j.day):max(j.day))),#predict the logistic curve for each species
    j.day=min(j.day):max(j.day),
  )%>%
  left_join(pv.forlog%>%dplyr::select(j.day,date)) ## add date back to tibble


pv.forlog%>%
  ggplot(aes(j.day,prop))+geom_point(aes=0.3)+geom_line(data=pv.pred,aes(x=j.day,y=pred),col="yellow4",size=2)+facet_wrap(year~.)

#Catharus ustulatus - Swainson's Thrush 
st.pred <- st%>%
  group_by(year)%>%
  summarize(
    pred=predict(nls(prop~SSlogis(j.day,Asym, xmid, scal)),newdata=data.frame(j.day=min(j.day):max(j.day))),#predict the logistic curve for each species
    j.day=min(j.day):max(j.day),
  )%>%
  left_join(st%>%dplyr::select(j.day,date)) ## add date back to tibble


st%>%
  ggplot(aes(j.day,prop))+geom_point(aes=0.3)+geom_line(data=st.pred,aes(x=j.day,y=pred),col="khaki4",size=2)+facet_wrap(year~.)

#Setophaga caerulescens - Black-throated Blue Warbler 
btw.pred <- btw%>%
  group_by(year)%>%
  summarize(
    pred=predict(nls(prop~SSlogis(j.day,Asym, xmid, scal)),newdata=data.frame(j.day=min(j.day):max(j.day))),#predict the logistic curve for each species
    j.day=min(j.day):max(j.day),
  )%>%
  left_join(btw%>%dplyr::select(j.day,date)) ## add date back to tibble


btw%>%
  ggplot(aes(j.day,prop))+geom_point(aes=0.3)+geom_line(data=btw.pred,aes(x=j.day,y=pred),col="royalblue3",size=2)+facet_wrap(year~.)

#Setophaga dominica - Yellow-throated  Warbler
yw.forlog <- yw%>%filter(year==2004 |year==2009 |year== 2011|year==2013| year==2014| year==2016| year==2017| year==2019)
yw.pred <- yw.forlog%>%
  group_by(year)%>%
  summarize(
    pred=predict(nls(prop~SSlogis(j.day,Asym, xmid, scal)),newdata=data.frame(j.day=min(j.day):max(j.day))),#predict the logistic curve for each species
    j.day=min(j.day):max(j.day),
  )%>%
  left_join(yw.forlog%>%dplyr::select(j.day,date)) ## add date back to tibble


yw.forlog%>%
  ggplot(aes(j.day,prop))+geom_point(aes=0.3)+geom_line(data=yw.pred,aes(x=j.day,y=pred),col="gold1",size=2)+facet_wrap(year~.)

#julian day closest to .25 population arriving 

#Megaceryle alcyon - belted kingfisher
mea.arrive.date <-mea.pred%>%
  group_by(year)%>%
  filter(j.day==j.day[which.min(abs(pred-0.25))])

mea.arrive.date%>%
  ggplot(aes(year,j.day))+geom_point()

#Sayornis phoebe - Eastern Phoebe
ep.arrive.date <-ep.pred%>%
  group_by(year)%>%
  filter(j.day==j.day[which.min(abs(pred-0.25))])

ep.arrive.date%>%
  ggplot(aes(year,j.day))+geom_point()

#Vireo philadelphicus - Philadelphia Virelo 
pv.arrive.date <-pv.pred%>%
  group_by(year)%>%
  filter(j.day==j.day[which.min(abs(pred-0.25))])

pv.arrive.date%>%
  ggplot(aes(year,j.day))+geom_point()

#Catharus ustulatus - Swainson's Thrush 
st.arrive.date <-st.pred%>%
  group_by(year)%>%
  filter(j.day==j.day[which.min(abs(pred-0.25))])

st.arrive.date%>%
  ggplot(aes(year,j.day))+geom_point()

#Setophaga caerulescens - Black-throated Blue Warbler 
btw.arrive.date <-btw.pred%>%
  group_by(year)%>%
  filter(j.day==j.day[which.min(abs(pred-0.25))])

btw.arrive.date%>%
  ggplot(aes(year,j.day))+geom_point()

#Setophaga dominica - Yellow-throated  Warbler
yw.arrive.date <-yw.pred%>%
  group_by(year)%>%
  filter(j.day==j.day[which.min(abs(pred-0.25))])

yw.arrive.date%>%
  ggplot(aes(year,j.day))+geom_point()

weather.d <- weather.d%>%
  mutate(year=as.integer(str_sub(date,1,4)), #add year
         date=as.Date(date))%>%
  group_by(year)%>% #group by year so we can compute julian day
  mutate(j.day=julian(date,origin=as.Date(paste0(unique(year),"-01-01"))), #add julian day
         date2=date,
         wdir.rad=(180-abs(wdf2-180))*pi/180, #radians so we can use a trig function to compute wind vector, scale degrees first to 180 scale to 2x pi and subtract from 180 (wind comes out of a direction)
         wvec=cos(wdir.rad)*-1*awnd # we want a negative value for positive value for 2x pi
  )%>% #store day in new column
  dplyr::select(id,year,date2,j.day,tmin,tmax,wvec)%>% #select the rows we need
  left_join(sta.d%>%select(id,name,migr.day))%>% #add the station id info (ie. name)
  mutate(j.day=j.day+migr.day)#make j.day ahead of BOS according to the migration days away so we can join weather along path

#first go at combining weather data + arrival date data 
#Megaceryle alcyon - belted kingfisher
mea.arr.weath <- mea.arrive.date%>%
  left_join(weather.d)%>%
  left_join(mea%>%dplyr::select(year,date,j.day))

#eastern phoebe 
ep.arr.weath <- ep.arrive.date%>%
  left_join(weather.d)%>%
  left_join(ep%>%dplyr::select(year,date,j.day))

#Vireo philadelphicus - Philadelphia Virelo 
pv.arr.weath <- pv.arrive.date%>%
  left_join(weather.d)%>%
  left_join(pv%>%dplyr::select(year,date,j.day))

#Catharus ustulatus - Swainson's Thrush 
st.arr.weath <- st.arrive.date%>%
  left_join(weather.d)%>%
  left_join(st%>%dplyr::select(year,date,j.day))

#Setophaga caerulescens - Black-throated Blue Warbler 
btw.arr.weath <- btw.arrive.date%>%
  left_join(weather.d)%>%
  left_join(btw%>%dplyr::select(year,date,j.day))

#Setophaga dominica - Yellow-throated  Warbler
yw.arr.weath <- yw.arrive.date%>%
  left_join(weather.d)%>%
  left_join(yw%>%dplyr::select(year,date,j.day))

weather.wk <-weather.d %>% 
  group_by(year,name) %>% 
  mutate(wk.tmin = frollmean(tmin, n=14,align="right"),
         wk.tmax = frollmean(tmax, n=14,align="right"),
         wk.wvec = frollmean(wvec, n=14,align="right")
  )%>%
  dplyr::select(j.day,date2,name,wk.tmin,wk.tmax,wk.wvec)

#following weather table modification (computing mean weather at each location) added to bird data

#Megaceryle alcyon - belted kingfisher
mea.arr.weath2 <- mea.arrive.date%>%
  left_join(weather.wk)

#linear model 1 -- 
mea.lmer <- lmer(j.day~tmin*tmax*wvec+(1|name),mea.arr.weath,na.action = "na.fail")
Anova(mea.lmer)

#0Mean two week weather preceding arrival
mea.lmer2 <- lmer(j.day~wk.tmin*wk.tmax*wk.wvec+(1|name),mea.arr.weath2,na.action = "na.fail")
Anova(mea.lmer2)

#secondary model 
mea.arr.aic <- dredge(mea.lmer2,fixed = c("wk.tmin","wk.tmax","wk.wvec"),)
mea.kb <- kable(mea.arr.aic[1:4,],caption = "Fit values for nested models of the most complicated lme model")
kable_styling(mea.kb)


#eastern phoebe 
ep.arr.weath2 <- ep.arrive.date%>%
  left_join(weather.wk)

#linear model 1 -- 
ep.lmer <- lmer(j.day~tmin*tmax*wvec+(1|name),ep.arr.weath,na.action = "na.fail")
Anova(ep.lmer)

#0Mean two week weather preceding arrival
ep.lmer2 <- lmer(j.day~wk.tmin*wk.tmax*wk.wvec+(1|name),ep.arr.weath2,na.action = "na.fail")
Anova(ep.lmer2)

#secondary model 
ep.arr.aic <- dredge(ep.lmer2,fixed = c("wk.tmin","wk.tmax","wk.wvec"),)
ep.kb <- kable(ep.arr.aic[1:4,],caption = "Fit values for nested models of the most complicated lme model")
kable_styling(ep.kb)

#Vireo philadelphicus - Philadelphia Virelo 
pv.arr.weath2 <- pv.arrive.date%>%
  left_join(weather.wk)

#linear model 1 -- 
pv.lmer <- lmer(j.day~tmin*tmax*wvec+(1|name),pv.arr.weath,na.action = "na.fail")
Anova(pv.lmer)

#0Mean two week weather preceding arrival
pv.lmer2 <- lmer(j.day~wk.tmin*wk.tmax*wk.wvec+(1|name),pv.arr.weath2,na.action = "na.fail")
Anova(pv.lmer2)

#secondary model 
pv.arr.aic <- dredge(pv.lmer2,fixed = c("wk.tmin","wk.tmax","wk.wvec"),)
pv.kb <- kable(pv.arr.aic[1:4,],caption = "Fit values for nested models of the most complicated lme model")
kable_styling(pv.kb)

#Catharus ustulatus - Swainson's Thrush 
st.arr.weath2 <- st.arrive.date%>%
  left_join(weather.wk)

#linear model 1 -- 
st.lmer <- lmer(j.day~tmin*tmax*wvec+(1|name),st.arr.weath,na.action = "na.fail")
Anova(st.lmer)

#0Mean two week weather preceding arrival
st.lmer2 <- lmer(j.day~wk.tmin*wk.tmax*wk.wvec+(1|name),st.arr.weath2,na.action = "na.fail")
Anova(st.lmer2)

#secondary model 
st.arr.aic <- dredge(st.lmer2,fixed = c("wk.tmin","wk.tmax","wk.wvec"),)
st.kb <- kable(st.arr.aic[1:4,],caption = "Fit values for nested models of the most complicated lme model")
kable_styling(st.kb)

#Setophaga caerulescens - Black-throated Blue Warbler 
btw.arr.weath2 <- btw.arrive.date%>%
  left_join(weather.wk)

#linear model 1 -- 
btw.lmer <- lmer(j.day~tmin*tmax*wvec+(1|name),btw.arr.weath,na.action = "na.fail")
Anova(btw.lmer)

#0Mean two week weather preceding arrival
btw.lmer2 <- lmer(j.day~wk.tmin*wk.tmax*wk.wvec+(1|name),btw.arr.weath2,na.action = "na.fail")
Anova(btw.lmer2)

#secondary model 
btw.arr.aic <- dredge(btw.lmer2,fixed = c("wk.tmin","wk.tmax","wk.wvec"),)
btw.kb <- kable(btw.arr.aic[1:4,],caption = "Fit values for nested models of the most complicated lme model")
kable_styling(btw.kb)

#Setophaga dominica - Yellow-throated  Warbler
yw.arr.weath2 <- yw.arrive.date%>%
  left_join(weather.wk)

#linear model 1 -- 
yw.lmer <- lmer(j.day~tmin*tmax*wvec+(1|name),yw.arr.weath,na.action = "na.fail")
Anova(yw.lmer)

#0Mean two week weather preceding arrival
yw.lmer2 <- lmer(j.day~wk.tmin*wk.tmax*wk.wvec+(1|name),yw.arr.weath2,na.action = "na.fail")
Anova(yw.lmer2)

#secondary model 
yw.arr.aic <- dredge(yw.lmer2,fixed = c("wk.tmin","wk.tmax","wk.wvec"),)
yw.kb <- kable(yw.arr.aic[1:4,],caption = "Fit values for nested models of the most complicated lme model")
kable_styling(yw.kb)


#best linear models (single day) 
best.mea.lmer1 <-  lmer(j.day~tmin+tmax+wvec+(1|name),mea.arr.weath,na.action = "na.fail")
bestmea1 <- Anova(best.mea.lmer1)
print(bestmea1)
best.ep.lmer1 <-  lmer(j.day~tmin+tmax+wvec+(1|name),ep.arr.weath,na.action = "na.fail")
bestep1 <- Anova(best.ep.lmer1)
print(bestep1)
best.pv.lmer1 <-  lmer(j.day~tmin+tmax+wvec+(1|name),pv.arr.weath,na.action = "na.fail")
bestpv1 <- Anova(best.pv.lmer1)
print(bestpv1)
best.st.lmer1 <-  lmer(j.day~tmin+tmax+wvec+(1|name),st.arr.weath,na.action = "na.fail")
bestst1 <- Anova(best.st.lmer1)
print(bestst1)
best.btw.lmer1 <-  lmer(j.day~tmin+tmax+wvec+(1|name),btw.arr.weath,na.action = "na.fail")
bestbtw1 <- Anova(best.btw.lmer1)
print(bestbtw1)
best.yw.lmer1 <-  lmer(j.day~tmin+tmax+wvec+(1|name),yw.arr.weath,na.action = "na.fail")
bestbtw1 <- Anova(best.yw.lmer1)
print(bestyw1)


#best linear models
best.mea.lmer <-  lmer(j.day~wk.tmin+wk.tmax+wk.wvec+(1|name),mea.arr.weath2,na.action = "na.fail")
bestmea <- Anova(best.mea.lmer)
print(bestmea)
best.ep.lmer <-  lmer(j.day~wk.tmin+wk.tmax+wk.wvec+(1|name),ep.arr.weath2,na.action = "na.fail")
bestep <- Anova(best.ep.lmer)
print(bestep)
best.pv.lmer <-  lmer(j.day~wk.tmin+wk.tmax+wk.wvec+(1|name),pv.arr.weath2,na.action = "na.fail")
bestpv <- Anova(best.pv.lmer)
print(bestpv)
best.st.lmer <-  lmer(j.day~wk.tmin+wk.tmax+wk.wvec+(1|name),st.arr.weath2,na.action = "na.fail")
bestst <- Anova(best.st.lmer)
print(bestst)
best.btw.lmer <-  lmer(j.day~wk.tmin+wk.tmax+wk.wvec+(1|name),btw.arr.weath2,na.action = "na.fail")
bestbtw <- Anova(best.btw.lmer)
print(bestbtw)
best.yw.lmer <-  lmer(j.day~wk.tmin+wk.tmax+wk.wvec+(1|name),yw.arr.weath2,na.action = "na.fail")
bestyw <- Anova(best.yw.lmer)
print(bestyw)
