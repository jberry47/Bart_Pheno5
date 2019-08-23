library(plyr)
library(ggplot2)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(car)
library(reshape2)
library(grid)
library(stringr)
library(RColorBrewer)
library(scales)
library(multcomp)
library(lme4)
library(lsmeans)
library(gridExtra)
library(ggdendro)
library(dendextend)
library(dendextendRcpp)
library(survminer)
library(ggridges)
library(gstat)
library(sp)

setwd("/media/jberry/Extra Drive 1/Danforth/Sorghum/Pheno5")

img_to_barcode <- read.csv("SnapshotInfo.csv",header = T,stringsAsFactors = F)
img_to_barcode <- img_to_barcode[img_to_barcode$tiles != "",]
colnames(img_to_barcode)[3] <- "Barcodes"
img_to_barcode <- img_to_barcode[,c("id","Barcodes","timestamp")]
sv_shapes <- read.table("pheno5_shapes.txt",header = F,stringsAsFactors = F,sep = " ")
colnames(sv_shapes) <- c("meta","area","hull_area","solidity","perimeter","width","height","cmx","cmy","hull_verticies","ex","ey","emajor","eminor","angle","eccen","circ","round","ar","fd","oof", "det")
sv_shapes$id <- substring(as.character(sapply(sv_shapes$meta,function(i) strsplit(i,"/")[[1]][2])),9)
sv_shapes$imgname <- as.character(sapply(sv_shapes$meta,function(i) strsplit(strsplit(i,"/")[[1]][3],"[.]")[[1]][1]))
sv_shapes <- join(sv_shapes,img_to_barcode[,c("id","Barcodes","timestamp")],by="id")
assoc <- read.csv("pheno5_design.csv",header=T,stringsAsFactors = F)
sv_shapes <- join(sv_shapes,assoc,by="Barcodes")
sv_shapes <- sv_shapes[!is.na(sv_shapes$Drought),]
sv_shapes$timestamp <- as.POSIXct(strptime(sv_shapes$timestamp,format = "%Y-%m-%d %H:%M:%S"))
beg <- as.POSIXct(strptime(paste(
  as.character(min(strptime(unlist(lapply(strsplit(as.character(sv_shapes$timestamp)," "),function(i)i[1])),format = "%Y-%m-%d"))),
  strsplit(as.character(min(strptime(unlist(lapply(strsplit(as.character(sv_shapes$timestamp)," "),function(i)i[2])),format = "%H:%M:%S")))," ")[[1]][2]
),format = "%Y-%m-%d %H:%M:%S"))
sv_shapes$DAP <- floor(as.numeric(difftime(sv_shapes$timestamp,beg,units = "days")))+6
sv_shapes$Drought <- ordered(sv_shapes$Drought, levels=c("WS_25","WW_80"))
sv_shapes$Genotype <- "BTx623"
sv_shapes$hour <- lubridate::hour(sv_shapes$timestamp)
sv_shapes <- sv_shapes[sv_shapes$Microbes != "Blank",]
sv_shapes <- sv_shapes[(read.csv("pheno5_outliers_bool.csv",header = T,stringsAsFactors = F)$x),]
area_convert <- 13.2*3.7/46856
tail(sv_shapes)

aggregate(data=sv_shapes,Barcodes~DAP,FUN=function(i)length(unique(i)))
#*************************************************************************************************
# Outlier Detection - Done above, don't rerun
#*************************************************************************************************
library(gputools)
chooseGpu(1)
cooksd <- cooks.distance(gpuGlm(data=sv_shapes,area~Microbes:Drought:as.factor(DAP)))
write.csv(cooksd < 3*mean(cooksd),"pheno5_outliers_bool.csv",row.names = F,quote = F)


#*************************************************************************************************
# Positional effects
#*************************************************************************************************
statefile <- NULL
for(my_dir in list.dirs("Pheno5_statefile")[-1]){
  for(i in list.files(my_dir)){
    temp <- read.table(paste0(substr(my_dir,0,24),"/",i),header = F,stringsAsFactors = F,sep = ";",skip = 1)
    num <- as.numeric(read.table(paste0(substr(my_dir,0,204),"/",i),header = F,stringsAsFactors = F,sep = ";",nrows = 1))
    meta <- substr(strsplit(i,"[.]")[[1]][3],5,7)
    statefile <- rbind(statefile,data.frame("Date"=substr(my_dir,18,24),"GH"=strsplit(meta,"_")[[1]][1],"Lane"=as.numeric(strsplit(meta,"_")[[1]][2]),"Position"=1:num,"Barcodes"=temp$V3,stringsAsFactors = F))
  }
}
statefile$GH <- ordered(statefile$GH,levels=5:1)
statefile$timestamp <- as.POSIXlt(strptime(statefile$Date,format = "%m%d%y"))
statefile$Date <- paste0(lubridate::month(statefile$timestamp),lubridate::day(statefile$timestamp),lubridate::year(statefile$timestamp))
tail(statefile)

dat <- aggregate(data=sv_shapes,area~Barcodes+timestamp+Drought+Microbes+DAP,"mean")
dat$Date <- paste0(lubridate::month(dat$timestamp),lubridate::day(dat$timestamp),lubridate::year(dat$timestamp))
dat <- join(statefile,dat,by=c("Barcodes","Date"))
dat <- na.omit(dat[,-6])
tail(dat)

trait <- "area"
all_spat <- data.frame("x"=rep(1:38,each=30),"y"=rep(1:30,38))
for(day in 9:25){
  sub <- dat[dat$DAP == day,]
  design_effects <- aggregate(data=sub,as.formula(paste0(trait,"~Drought+Microbes")),function(i) mean(as.numeric(unlist(i))))
  design_sd <- aggregate(data=sub,as.formula(paste0(trait,"~Drought+Microbes")),function(i) sd(as.numeric(unlist(i))))
  sub$trait_n <- apply(sub[,c("Drought","Microbes",trait)],1,function(i){(as.numeric(i[3])-design_effects[design_effects$Drought == i[1] & design_effects$Microbes ==i[2],trait])/design_sd[design_sd$Drought == i[1] & design_sd$Microbes ==i[2],trait]})
  sub$Lane1 <- sub$Lane+(as.numeric(as.character(sub$GH))-1)*8
  coordinates(sub) <- ~Lane1+Position
  vgm1 <- variogram(trait_n~1,sub)
  fit <- fit.variogram(vgm1,model=vgm(0,"Sph"))
  sub.grid <- data.frame("x"=rep(1:38,each=30),"y"=rep(1:30,38))
  coordinates(sub.grid) <- ~x+y
  kriged <- krige(trait_n~1,sub,sub.grid,fit,maxdist=3)
  out <- setNames(data.frame(kriged)[,1:3],c("x","y",paste0("pred_",day)))
  all_spat <- join(all_spat,out,by=c("x","y"))
}
all_spat$avg <- rowMeans(all_spat[,3:(ncol(all_spat))],na.rm = T)
all_spat$var <- apply(all_spat[,3:(ncol(all_spat)-1)],1,function(i) var(i,na.rm = T))
all_spat$cv <- all_spat$var/all_spat$avg
all_spat$gh <- sapply(all_spat$x,function(i) if(i %in% 1:8){"GH1"}else if(i %in% 9:16){"GH2"}else if(i %in% 17:24){"GH3"}else{"GH4"})
all_spat$gh <- ordered(all_spat$gh,levels=c("GH4","GH3","GH2","GH1"))
head(all_spat)

p <- ggplot(all_spat,aes(x,y))+
  facet_grid(~gh,space = "free_x",scales = "free_x")+
  geom_tile(aes(fill=avg))+
  scale_x_reverse(expand = c(0, 0),breaks=c(seq(2,38,2)))+
  scale_y_continuous(expand = c(0, 0),breaks=c(seq(2,30,2)))+
  scale_fill_gradient2(mid = "gray95",high="darkgreen",limits=c(-2.5,2.5))+
  theme_light()+
  theme(panel.spacing = unit(0.03, "lines"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=12,color="white"),
    strip.text.y=element_text(size=12,color="white"))
p  
ggsave("pheno5_bothww-ws_allDAP_areaNorm_spatial.png",width=6.9,height=4.9,plot = p, dpi = 300)

my_list <- do.call(rbind,lapply(unique(dat$DAP),function(i){
  sub <- dat[dat$DAP == i & dat$Drought == "WS_25",]
  sub$Lane1 <- sub$Lane+(as.numeric(as.character(sub$GH))-1)*8
  colnames(sub)[4] <- "y"
  colnames(sub)[11] <- "x"
  df <- join(sub,all_spat[,c("x","y","avg")],by=c("x","y"))
  df$area_c <- df$area-df$area*df$avg*0.5
  df
}))

ggplot(my_list,aes(area*area_convert,area_c*area_convert))+
  geom_point(aes(color=avg),size=3)+
  scale_color_gradient2(mid = "gray95",high="darkgreen")+
  geom_abline(slope=1, intercept=0,color="gray20")+
  theme_light()

ggplot(my_list[my_list$Microbes %in% c("Control","SynCom A","SynCom B","SynCom C"),],aes(DAP,area_c*area_convert))+
  geom_smooth(aes(color=Microbes),method = "loess")

sub <- dat[dat$DAP == 20 & dat$Drought == "WS_25",]
sub$Lane1 <- sub$Lane+(as.numeric(as.character(sub$GH))-1)*8
colnames(sub)[4] <- "y"
colnames(sub)[11] <- "x"

df <- join(sub,all_spat[,c("x","y","avg")],by=c("x","y"))
df$status <- "Uncalibrated"
df1 <- df
df1$area <- df1$area-df1$area*df1$avg*0.5
df1$status <- "Calibrated"
df <- rbind(df,df1)
df$status <- ordered(df$status, levels=c("Uncalibrated","Calibrated"))

ggplot(df[df$Microbes %in% c("Control","SynCom A","SynCom B","SynCom C"),],aes(Microbes,area*area_convert))+
  facet_grid(~status)+
  geom_boxplot()

#*************************************************************************************************
# Quick plots
#*************************************************************************************************
sub <- sv_shapes
p <- ggplot(sub[sub$Microbes %in% c("Control","SynCom A","SynCom B","SynCom C"),],aes(DAP,area*area_convert))+
  facet_wrap(~Drought)+
  geom_smooth(aes(color=Microbes),method = "loess")+
  ylab(~~Area~(cm^2))+
  scale_color_manual(values=c("gray20",muted("red",60,100),muted("green",60,100),muted("cyan",60,100)))+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p
ggsave("pheno5_area_syncom_trends.png",width=8.03,height=3.9,plot = p, dpi = 300)


#*************************************************************************************************
# Comparing to unloading images (make sure outliers ARE NOT removed)
#*************************************************************************************************
ul_shapes <- read.table("unloading_shapes.txt",header = T,stringsAsFactors = F,sep = " ")
colnames(ul_shapes) <- paste0("ul_",colnames(ul_shapes))
ul_shapes$Barcodes <- as.character(sapply(ul_shapes$ul_meta,function(i) str_sub(strsplit(i,"/")[[1]][2],1,13)))
ul_shapes <- join(aggregate(.~Barcodes,ul_shapes[,c(2:22)],mean),assoc,by="Barcodes")

endpoint <- sv_shapes[sv_shapes$Barcodes %in% ul_shapes$Barcodes,]
endpoint <- do.call(rbind,lapply(split(endpoint,endpoint$Barcodes),function(i) i[i$DAP == max(i$DAP),]))
endpoint <- aggregate(.~Barcodes,endpoint[,c(2:22,25)],mean)

merged <- join(endpoint,ul_shapes,by=c("Barcodes"))

shapes <- c("area","hull_area","solidity","perimeter","width","height","cmx","cmy","hull_verticies","ex","ey","emajor","eminor","angle","eccen","circ","round","ar","fd")
df <- data.frame("shape"=shapes,"r2"=sapply(shapes,function(i) cor(merged[,i],merged[,paste0("ul_",i)])^2),row.names = NULL)
df$shape <- ordered(df$shape,df$shape[order(df$r2,decreasing = T)])
p <- ggplot(df,aes(shape,r2))+
  geom_bar(stat = "identity")+
  theme_bw()+
  scale_y_continuous(limits=c(0,1))+
  ylab(~~R^2)+
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))+
  theme(axis.text = element_text(size = 14),
    axis.title.y= element_text(size = 18),
    axis.title.x = element_blank())+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(panel.border = element_rect(colour = "gray60", fill=NA, size=1,linetype = 1))+
  theme(legend.position = "top")+
  guides(fill = guide_legend(title = ""))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave("pheno5_oof_agreement.png",width=6.03,height=3.9,plot = p, dpi = 300)


p <- ggplot(merged,aes(ul_area,area))+
  geom_point(aes(color=as.factor(oof)))+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p 
ggsave("pheno5_cmy_oof-agreement.png",width=6.35,height=4.69,plot = p, dpi = 300)


#*************************************************************************************************
# Looking at leaf data
#*************************************************************************************************
leaves <- setNames(read.table("pheno5_leaves.txt",stringsAsFactors = F,header=F,sep=" "),c("meta","leaf_num","area","eucD","path_length","tort"))
leaves$id <- substring(as.character(sapply(leaves$meta,function(i) strsplit(i,"/")[[1]][2])),9)
leaves$imgname <- as.character(sapply(leaves$meta,function(i) strsplit(strsplit(i,"/")[[1]][3],"[.]")[[1]][1]))
leaves <- join(leaves,img_to_barcode[,c("id","Barcodes","timestamp")],by="id")
leaves <- join(leaves,assoc,by="Barcodes")
leaves$timestamp <- as.POSIXct(strptime(leaves$timestamp,format = "%Y-%m-%d %H:%M:%S"))
beg <- as.POSIXct(strptime(paste(
  as.character(min(strptime(unlist(lapply(strsplit(as.character(leaves$timestamp)," "),function(i)i[1])),format = "%Y-%m-%d"))),
  strsplit(as.character(min(strptime(unlist(lapply(strsplit(as.character(leaves$timestamp)," "),function(i)i[2])),format = "%H:%M:%S")))," ")[[1]][2]
),format = "%Y-%m-%d %H:%M:%S"))
leaves$DAP <- floor(as.numeric(difftime(leaves$timestamp,beg,units = "days")))+6
leaves <- leaves[!is.na(leaves$Drought),]
#leaves <- leaves[leaves$tort < 10,]
leaves <- leaves[leaves$imgname %in% sv_shapes$imgname,]

back <- data.frame(with(leaves[leaves$DAP==17,],table(leaf_num,Microbes,Drought))[,,1])
back_list <- split(back,back$Microbes)
df <- data.frame(do.call("rbind",lapply(back_list,function(i) {i$prob <- i$Freq/i$Freq[1];i})))

p <- ggplot(df[df$Microbes %in% c("Control","SynCom A","SynCom B","SynCom C"),],aes(leaf_num,prob))+
  geom_line(aes(color=Microbes,group=Microbes))+
  xlab("Number of leaves")+
  ylab("Proportion of Plants")+
  theme_light()+
  theme(axis.text = element_text(size = 14),
    axis.title= element_text(size = 18))+ 
  theme(strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p
ggsave("pheno5_WS_dap15_numLeavesProp_syncoms.png",width=6.01,height=4.82,plot = p, dpi = 300)



str(back)
library(gputools)
chooseGpu(1)
df <- aggregate(data=leaves,tort~Barcodes+Microbes+Drought+DAP,FUN="mean")
cooksd <- cooks.distance(gpuGlm(data=df,tort~Microbes:Drought:as.factor(DAP)))
df <- df[cooksd < 3*mean(cooksd),]
leaves_list <- split(leaves,leaves$imgname)

leaves_list[[10000]]


ggplot(df[df$DAP == 10,],aes(Microbes,tort))+
  facet_wrap(~Drought)+
  geom_boxplot(aes(fill=Microbes))

ggplot(df[df$Microbes %in% c("Control","SynCom A","SynCom B"),],aes(DAP,tort))+
  facet_wrap(~Drought)+
  geom_smooth(aes(color=Microbes))+
  theme_light()+
  theme(axis.text = element_text(size = 12),
    axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
    strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
  


#*************************************************************************************************
# Helper functions
#*************************************************************************************************
get_color <- function(file_name,start,stop){
  color_data <- read.table(file_name,header = F,stringsAsFactors = F,sep = " ")[,-257]
  color_data$id <- as.character(sapply(color_data$V1,function(i) strsplit(strsplit(i,"/")[[1]][2],"snapshot")[[1]][2]))
  color_data$imgname <- as.character(sapply(color_data$V1,function(i) strsplit(strsplit(i,"/")[[1]][3],"[.]")[[1]][1]))
  color_data <- join(color_data,img_to_barcode[,c("id","Barcodes","timestamp")],by="id")
  color_data <- join(color_data,assoc,by="Barcodes")
  color_data$timestamp <- strptime(color_data$timestamp,format = "%Y-%m-%d %H:%M:%S")
  color_data$DAP <- floor(as.numeric(difftime(color_data$timestamp,beg,units = "days")))+2
  color_data[,start:stop] <- t(apply(color_data[,start:stop],1,function(i){i/(sum(i,na.rm = T)+1)}))*100
  color_data$hr <- as.POSIXlt(color_data$timestamp)$hour
  #color_data$status <- sapply(color_data$hr,function(i) if(i %in% 7:21){"Daytime"}else{"Nighttime"})
  return(color_data)
}

hist_avg <- function(data,start,stop){
  sub <- data
  test <- data.frame(do.call("rbind",lapply(split(sub,sub$Drought),function(t){
    data.frame(do.call("rbind",lapply(split(t,t$Microbes),function(g){
      data.frame(do.call("rbind",lapply(split(g,g$DAP),function(m){
        colMeans(m[,start:stop],na.rm = T)
      }
      )))
    }
    )))
  })))
  return(test)
}

hist_sd <- function(data,start,stop){
  sub <- data
  test <- data.frame(do.call("rbind",lapply(split(sub,sub$Drought),function(t){
    data.frame(do.call("rbind",lapply(split(t,t$Microbes),function(g){
      data.frame(do.call("rbind",lapply(split(g,g$DAP),function(m){
        apply(m[,start:stop],2,function(i){sd(i,na.rm = T)})
      }
      )))
    }
    )))
  })))
  return(test)
}

#*************************************************************************************************
# NIR data
#*************************************************************************************************
nir <- get_color("pheno5_nir.txt",2,256)
nir$intensityAVG <- apply(nir[,3:255],1,function(i){sum((i/100)*(2:255),na.rm = T)})
nir <- nir[!is.na(nir$Microbes),]
wd <- read.csv("pheno5_dryweight.csv",stringsAsFactors = F)
wf <- rbind(read.csv("freshweight1.csv",stringsAsFactors = F),read.csv("freshweight2.csv",stringsAsFactors = F))

nir_sub <- setNames(data.frame(do.call(rbind,lapply(split(nir[nir$DAP == 21,],nir$Barcodes[nir$DAP == 21]),function(i)mean(i$intensityAVG)))),c("nir_avg"))
nir_sub$Barcodes <- rownames(nir_sub); rownames(nir_sub) <- NULL
test <- join(join(join(wd,wf),nir_sub),assoc)
test$weightDIFF <- test$freshweight-test$Dryweight
head(test)

p <- ggplot(test,aes(weightDIFF,nir_avg))+
  facet_wrap(~Microbes,nrow = 3)+
  geom_point(aes(color=Drought))+
  xlab("Water Content (mg)")+
  ylab("NIR Intensity AVG")+
  geom_smooth()+
  theme_light()+
  theme(axis.text = element_text(size = 12),
    axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
    strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p
ggsave("pheno5_nir-waterContent.png",plot=p,dpi = 300,width = 6.91,height=5.49)


p <- ggplot(test,aes(weightDIFF,nir_avg))+
  geom_point(aes(color=Microbes))+
  xlab("Water Content (mg)")+
  ylab("NIR Intensity AVG")+
  geom_smooth(aes(color=Microbes),se=F)+
  theme_light()+
  theme(axis.text = element_text(size = 12),
    axis.title= element_text(size = 18))+
  theme(plot.title = element_text(hjust = 0.5),
    strip.background=element_rect(fill="gray50"),
    strip.text.x=element_text(size=14,color="white"),
    strip.text.y=element_text(size=14,color="white"))
p


p <- ggplot(nir[nir$Microbes != "Blank",],aes(DAP,intensityAVG))+
  facet_grid(~Drought)+
  geom_smooth(aes(color=Microbes),method = "loess")
p
