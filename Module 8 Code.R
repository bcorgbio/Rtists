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

data%>%
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
  lapply(sts,function(x) ncdc_stations(stationid = x)$data )
  )%>%
    left_join(usmap_transform(.[,c("longitude","latitude")]))%>%
    mutate(name=str_sub(name, -5,-4))%>%#simplify the name column, grab just the state
    mutate(migr.day=c(10,5,0))%>% #so we can look at wind speed 0, 5 or 10 days before arrive in boston
    separate(id,into = c("station.type","id"))%>%#need to cut station type out from station id number
    print()
