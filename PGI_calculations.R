library(ncdf4)
library(RNetCDF)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(sp)
library(raster)
library(akima)
library(ggExtra)
library(spam)
library(lubridate)


################################
#  Data
################################

# Date to calculate
dts<-"2025-01-01"
(doy<-yday(as.Date(dts)))

R_clear<-readRDS("R_clear.rds")
mrsos_p5<-readRDS("mrsos_p5.rds")
mrsos_p95<-readRDS("mrsos_p95.rds")


################################
#  Functions
################################

type=c("C3", "C4", "SH", "TR", "CR")
Tlow<-data.frame(type=type,val=c(5, 10, 5, 5, 5))
Topt<-data.frame(type=type,val=c(20, 30, 25, 23, 20))
Tupp<-data.frame(type=type,val=c(32, 42, 40, 35, 35))



KtoC<-function(K){
  K - 273.15
}

BriereF<- function(a, clss){
  ((pmax(a, Tlow$val[Tlow$type==clss]) - Tlow$val[Tlow$type==clss])*
     sqrt(Tupp$val[Tupp$type==clss] -
            pmin(a, Tupp$val[Tupp$type==clss]))) / 
    ((Topt$val[Topt$type==clss] - Tlow$val[Tlow$type==clss]) *
       sqrt(Tupp$val[Tupp$type==clss] - Topt$val[Topt$type==clss]))
}

TT<-seq(-10, 50, 0.1)
BF<-rbind(data.frame(Type="C3", Temp=TT, Value=BriereF(TT, "C3")),
          data.frame(Type="C4", Temp=TT, Value=BriereF(TT, "C4")),
          data.frame(Type="SH", Temp=TT, Value=BriereF(TT, "SH")),
          data.frame(Type="TR", Temp=TT, Value=BriereF(TT, "TR")),
          data.frame(Type="CR", Temp=TT, Value=BriereF(TT, "CR")))

ggplot(BF)+
  geom_path(aes(x=Temp, y=Value, col=Type))+
  theme_bw()

################################
#  Map etc
################################

world <- ne_countries(scale = "medium", returnclass = "sf")
coord.lim<-c(113, 154, -44, -10.1)
fig<-ggplot(data = world) +
  geom_sf(fill="white") +
  coord_sf(xlim = c(coord.lim[1], coord.lim[2]), 
           ylim = c(coord.lim[3], coord.lim[4]), expand = FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

fig

################################################################
#  Analysis
################################################################

flnm<-"PFT_EVI_C3_C4_2020_ESA.rds"
grd<-readRDS(flnm)

grd$C3n<-0
grd$C3n<-(grd$C3/(grd$C3+grd$C4))* pmin(grd$C3+grd$C4, grd$GRASSNAT/100)
grd$C4n<-0
grd$C4n<-(grd$C4/(grd$C3+grd$C4))* pmin(grd$C3+grd$C4, grd$GRASSNAT/100)
grd$C3n[which(is.na(grd$C3n))]<-0
grd$C4n[which(is.na(grd$C4n))]<-0

 
fig+
  geom_tile(data=grd[grd$C3n>0,], aes(x=lon, y=lat, fill=C3n))+
  scale_fill_continuous(name="C3%", type="viridis", direction=-1, limits=c(0,1))+
  theme_bw()+xlab("")+ylab("")

fig+
  geom_tile(data=grd[grd$C4n>0,], aes(x=lon, y=lat, fill=C4n))+
  scale_fill_continuous(name="%C4", type="viridis", direction=-1, limits=c(0,1))+
  theme_bw()+xlab("")+ylab("")

fig+
  geom_tile(data=grd[grd$GRASSMAN>0,], aes(x=lon, y=lat, fill=GRASSMAN))+
  scale_fill_continuous(name="%Crop", type="viridis", direction=-1, limits=c(0,100))+
  theme_bw()+xlab("")+ylab("")

grd$SHRB<-(grd$SHRUBSBD+grd$SHRUBSBE+grd$SHRUBSNE)/100
fig+
  geom_tile(data=grd[grd$SHRB>0,], aes(x=lon, y=lat, fill=SHRB*100))+
  scale_fill_continuous(name="%", type="viridis", direction=-1, limits=c(0,100))+
  theme_bw()+xlab("")+ylab("")

grd$TREES<-(grd$TREESBD+ grd$TREESBE+ grd$TREESND+ grd$TREESNE)/100
fig+
  geom_tile(data=grd[grd$TREES>0,], aes(x=lon, y=lat, fill=TREES*100))+
  scale_fill_continuous(name="%", type="viridis", direction=-1, limits=c(0,100))+
  theme_bw()+xlab("")+ylab("")

flnm<-paste0("rsds/rsds_daily_", dts, ".rds")
xx<-readRDS(flnm)
grd$swsf<-xx$swsf
grd$LI<-pmin(pmax(grd$swsf/R_clear[, doy], 0),1)

fig+
  geom_tile(data=grd, aes(x=lon, y=lat, fill=LI))+
  scale_fill_continuous(name="", type="viridis", direction=-1, limits=c(0,1))+
  theme_bw()+xlab("")+ylab("")+ theme(legend.position ="bottom")

fig+
  geom_tile(data=grd[grd$LI>0,], aes(x=lon, y=lat, fill=LI))+
  scale_fill_continuous(name="%", type="viridis", direction=-1, limits=c(0,1))+
  theme_bw()+xlab("")+ylab("")

fl<-paste0("mrsos/mrsos_daily_", dts, ".rds")
xx<-readRDS(fl)
grd$MI<-pmax(pmin((xx$mrsos-mrsos_p5[, doy])/(mrsos_p95[, doy]-mrsos_p5[, doy]),1),0)


fig+
  geom_tile(data=grd[grd$MI>0,], aes(x=lon, y=lat, fill=MI))+
  scale_fill_continuous(name="%", type="viridis", direction=-1, limits=c(0,1))+
  theme_bw()+xlab("")+ylab("")

fl<-paste0("tas/tas_daily_", dts, ".rds")
xx<-readRDS(fl)
xx$cels<-KtoC(xx$tas)

aa<-as.numeric(xx$cels > Tlow$val[Tlow$type=="C3"]) * 
  as.numeric(xx$cels < Tupp$val[Tupp$type=="C3"]) * 
  (grd$C3n)*BriereF(xx$cels, "C3") +
  as.numeric(xx$cels > Tlow$val[Tlow$type=="C4"]) * 
  as.numeric(xx$cels < Tupp$val[Tupp$type=="C4"]) * 
  (grd$C4n)*BriereF(xx$cels, "C4") +
  as.numeric(xx$cels > Tlow$val[Tlow$type=="SH"]) * 
  as.numeric(xx$cels < Tupp$val[Tupp$type=="SH"]) * 
  (grd$SHRB)*BriereF(xx$cels, "SH") +
  as.numeric(xx$cels > Tlow$val[Tlow$type=="TR"]) * 
  as.numeric(xx$cels < Tupp$val[Tupp$type=="TR"]) * 
  (grd$TREES)*BriereF(xx$cels, "TR")+
  as.numeric(xx$cels > Tlow$val[Tlow$type=="CR"]) * 
  as.numeric(xx$cels < Tupp$val[Tupp$type=="CR"]) * 
  (grd$GRASSMAN/100)*BriereF(xx$cels, "CR")

grd$CR<-  as.numeric(xx$cels > Tlow$val[Tlow$type=="CR"]) * 
  as.numeric(xx$cels < Tupp$val[Tupp$type=="CR"]) * 
  (grd$GRASSMAN/100)*BriereF(xx$cels, "CR")

fig+
  geom_tile(data=grd[grd$CR>0,], aes(x=lon, y=lat, fill=CR))+
  scale_fill_continuous(name="%", type="viridis", direction=-1)+
  theme_bw()+xlab("")+ylab("")


grd$TI<-pmin(pmax(aa,0),1)
fig+
  geom_tile(data=grd[grd$TI>0,], aes(x=lon, y=lat, fill=TI))+
  scale_fill_continuous(name="%", type="viridis", direction=-1, limits=c(0,1))+
  theme_bw()+xlab("")+ylab("")

grd$PGI<-grd$LI*grd$MI*grd$TI

fig+
  geom_tile(data=grd[grd$PGI>0,], aes(x=lon, y=lat, fill=PGI))+
  scale_fill_continuous(name="%", type="viridis", direction=-1, limits=c(0,1))+
  theme_bw()+xlab("")+ylab("")

