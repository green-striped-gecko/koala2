library(plyr)
library(tidyverse)
library(lubridate)
library(RODBC)

params <- list()
params$database <- "database/GDA94z55Database2020June.accdb"
params$study_periods <-  "periods/periods_coastal.csv"
params$study_site <- "sites/sites_coastal_"


newdb <- params$database

# retrieve data from database  
channel_new <- odbcConnectAccess2007(newdb) #open channel
locations<- sqlFetch(channel_new,"Locations") #load location data
subs <- sqlFetch(channel_new, "SUB_AREAS")

trees <- sqlFetch(channel_new,"TREE_MASTER") 
treelist<-sqlFetch(channel_new,"Tree_Species")
sitevisits<-sqlFetch(channel_new,"SITE_VISITS")
odbcClose(channel_new)
#### core sites and periods need to be loaded!!!!

defperiods <- read.csv(params$study_periods, stringsAsFactors = TRUE)
nperiods <- nrow(defperiods)



locations<-locations[,c("LOCATION_ID","PLANNED_EASTING","PLANNED_NORTHING","LATITUDE","LONGITUDE","GRID")]

#round to no decimals!!
locations$PLANNED_EASTING <- round(locations$PLANNED_EASTING)
locations$PLANNED_NORTHING <- round(locations$PLANNED_NORTHING)


keep<-c("SITE_VISIT_ID","LOCATION_ID","S_EASTINGS","S_NORTHINGS","A_LEVEL","VISITED_SITE","DATE_SURVEYED","STUDY_AREA","GEOLOGY","SOIL_TYPE")
sitevisits<-sitevisits[,(names(sitevisits) %in% keep)]


# now clean up the tree data...
treedrop<-c("ID","GGFP_CONFIRMED","COMMENTS","SCAT AGE","BARK CHEWING")
trees<-trees[,!(names(trees) %in% treedrop)]

#clean up site visits (adding xys and remove missing data)
combloc<-plyr::join(sitevisits,locations, by="LOCATION_ID")
combloc<-combloc[!is.na(combloc$VISITED_SITE),]



# now combine everything...

loctrees<- plyr::join(combloc,trees, by="SITE_VISIT_ID", type="left")
#loctrees<-merge(combloc,trees,all.x=TRUE)
loctrees$counter<-1
loctrees$subarea<- NA

# Clean up the treelist data set
names(treelist)[grep("ID$",names(treelist))]<-"TREE_SPECIES_ID"
treelist<-treelist[,c("SPECIES","TREE_SPECIES_ID","SCIENTIFIC_NAME","COMMON_NAME")]
treelist<-treelist[treelist$SPECIES!="Blank", ] #remove species ID blank


#add subarea to the data set
#load each subarea file xys
#find all sites within 500m of for each subareas
#add subarea number to loctrees


#unique sites visited in the database
unlocsites <- data.frame(gridcoords=unique(paste(loctrees$PLANNED_EASTING, loctrees$PLANNED_NORTHING, sep="-")))
unlocsites$x <- as.numeric(lapply(strsplit(unlocsites$gridcoords,"-"), function(x) x[1]))
unlocsites$y <- as.numeric(lapply(strsplit(unlocsites$gridcoords,"-"), function(x) x[2]))


loctrees$gridcoords <- paste(loctrees$PLANNED_EASTING, loctrees$PLANNED_NORTHING, sep="-")
loctrees$subarea <- NA


for (i in 1:5)  #5 subareas
{
coresites <- read.csv(paste0(params$study_site,i,".csv"), stringsAsFactors = TRUE)

coresites$gridcoords <- paste(as.character(coresites$x), as.character(coresites$y), sep="-")

#find for subareas the sites within 500 meters
allsites <- rbind(coresites[,c("x","y")], unlocsites[,c("x","y")])

distall <- as.matrix(dist(allsites))

#check only distances between coresites and visited sites
distsites <- distall[-(1:nrow(coresites)),]
distsites <- distsites[,1:nrow(coresites)]

index <- which((apply(distsites , 1, min))<500 )


plot(coresites$x, coresites$y)
points(unlocsites$x, unlocsites$y, cex=0.8, pch=16, col="blue")
points(unlocsites$x[index], unlocsites$y[index], cex=0.5, pch=17, col="red")

#now select all the sites 
within <- unlocsites[index,]


  index2<- which(loctrees$gridcoords %in% within$gridcoords)

loctrees$subarea[index2] <- i

}


loctrees_sub <- loctrees[!is.na(loctrees$subarea),]

#define sampling period as well

loctrees_sub$year <- year(loctrees_sub$DATE_SURVEYED)
loctrees_sub$period <- cut(loctrees_sub$year,breaks = c(2000,2009,2014,2019), labels = 1:3)


#join tree species
loctrees_sub_trees <- join(loctrees_sub, treelist, by="TREE_SPECIES_ID", type="left")

#remove multiple visits in a site at one period.....
#just rename tomake live easier
koalas <- loctrees_sub_trees

k1 <- filter(koalas, !is.na(KFP &!is.na(period)) & year>=2007 & year<=2019)

#now get rid of sites within a period that were visited more than once

#filter out sites that have been visited more than once in a period (and subarea)
kmt1<-  ddply( k1,.(LOCATION_ID, period), summarise,  len=length(KFP))

#multiple visits consistent?
mv_locid <- kmt1$LOCATION_ID[kmt1$len>30]
tt <- list()
#select only the first visit and remove the rest

takeoutvid <- NULL

for (i in 1:length(mv_locid))
{
  sel <- mv_locid[i]
  for (ii in 1:3) {
    take <- k1$LOCATION_ID==sel & k1$period==ii
    if (sum(take)>30) #now take the first ones
    {
      if (length(takeoutvid)==0) takeoutvid <- unique(k1$SITE_VISIT_ID[take])[-1] else
        takeoutvid<- c(takeoutvid,unique(k1$SITE_VISIT_ID[take])[-1])
    }
  }
  
}

k2 <- k1[!k1$SITE_VISIT_ID %in% takeoutvid,]


#write.csv(k2, "loctrees_sub_trees_08072021.csv", row.names = FALSE)

#k2 <- read.csv("loctrees_sub_trees_29062021.csv")

#write.csv(k2, "loctrees_sub_trees_29062021b.csv", row.names = FALSE)


#write.csv(treelist, "tree_species.csv",row.names=FALSE)

par(mfrow=c(3,2))

for (i in 1:5)  #5 subareas
{
  coresites <- read.csv(paste0(params$study_site,i,".csv"), stringsAsFactors = TRUE)
  

plot(coresites$x, coresites$y)
points(loctrees_sub_trees$PLANNED_EASTING[loctrees_sub_trees$subarea==i],loctrees_sub_trees$PLANNED_NORTHING[loctrees_sub_trees$subarea==i], cex=0.5, pch=17, col="red")

}




#summarise occurrence accross 
kocc<-  plyr::ddply( k2,.(SITE_VISIT_ID,LOCATION_ID, PLANNED_EASTING, PLANNED_NORTHING, subarea, period), summarise, occ=as.numeric(sum(KFP)>0), act=sum(KFP), len=length(KFP))

#number of visits per subarea and perdod
table( kocc$subarea,kocc$period)
table( kocc$subarea,kocc$period, kocc$occ)



library(tidyverse)
library(plyr)
fn <- paste0("d:/bernd/r/ikoala/sites/sites_coastal_",1:5,".csv")

ff <- list()
for (i in 1:5)
{
  ff[[i]]<- read.csv(fn[i], stringsAsFactors = TRUE)
  ff[[i]]$subarea <- i
}


sites <- do.call(rbind,ff)





#plot sites
ggplot()+geom_point(data = sites, aes(x=x, y=y))+geom_point(data=kocc, aes(x=PLANNED_EASTING, y=PLANNED_NORTHING, col=factor(occ)))+facet_wrap(~subarea, scales = "free")


ggplot()+geom_point(data = sites, aes(x=x, y=y))+geom_point(data=kocc, aes(x=PLANNED_EASTING, y=PLANNED_NORTHING, col=factor(occ)))+facet_wrap(~period)+scale_color_discrete(name="Occupancy")


ggplot()+geom_point(data = sites, aes(x=x, y=y))+geom_point(data=kocc, aes(x=PLANNED_EASTING, y=PLANNED_NORTHING, col=factor(occ)))+ coord_fixed()+scale_color_discrete(name="Occupancy")







