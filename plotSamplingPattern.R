library(raster) # for georeferencing grid data
library(tidyverse) # for ggplot and dplyr
library(viridis) # for colour scale
library(ggthemes) # for theme_map()

##### Specifying directories -----
# location of the directory containing runDir and saDir folders in relation to current working directory
# line 9 needs to be altered before running the script
baseDir <- "" 

runDir  <- "epidemicRuns" # name of the folder containing the epidemic simulations
saDir   <- "samplingPattern" # name of the folder containing the optimisation results
baseMap <- "propFull.txt" # name of the basemap to be plotted at the end (include extension)

crs <- "+proj=utm +zone=17 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" # specify the CRS used for spatial data (UTM zone 17N in the case of Florida)

##### Extracting host and spatial information -----
setwd(baseDir)

nRuns   <- length(grep("endReason",dir(runDir))) # identifies number of simulation runs

hostInfoFile <- paste(runDir,"\\activeLandscape.txt",sep="") # find activeLandscape file
hostData <- read.csv(hostInfoFile,sep=" ",header=F) # read in activeLandscape file
names(hostData) <- c("hostCol","hostRow","propFull","hostID") # name the variables

basemapDat <- read.table(baseMap, header=F, nrows=6) # Pulls out the spatial parameters
basemapDatNames <- basemapDat$V1 # takes the variable names
basemapDat <- data.frame(t(sapply(basemapDat$V2, as.numeric))) # takes transpose
colnames(basemapDat) <- basemapDatNames # reattaches names

basemapCells <- read.table(baseMap, header=F, skip=6) # extracts grid of basemap
basemapCells[basemapCells==basemapDat$NODATA_value]<-NA # converts NODATA values to NA

#rm(list=c("basemapDat","basemapDatNames"))

##### Creating a georeferenced dataframe -----

# raster map of basemap grid
basemapRaster <- 
  raster(as.matrix(basemapCells), 
         xmn=basemapDat$xllcorner,
         xmx=(basemapDat$xllcorner + basemapDat$ncols*basemapDat$cellsize),
         ymn=basemapDat$yllcorner,
         ymx=(basemapDat$yllcorner + basemapDat$nrows*basemapDat$cellsize),
         crs=crs) # convert the grid into a georeferenced raster map

# raster map of hostID
hostIDRaster <- 
  raster(nrow=basemapDat$nrows,
         ncol=basemapDat$ncols, 
         xmn=basemapDat$xllcorner,
         xmx=(basemapDat$xllcorner + basemapDat$ncols*basemapDat$cellsize),
         ymn=basemapDat$yllcorner,
         ymx=(basemapDat$yllcorner + basemapDat$nrows*basemapDat$cellsize),
         crs=crs) # create an empty georeferenced raster map
cells <- cellFromRowCol(hostIDRaster, row=(hostData$hostRow+1), col=(hostData$hostCol+1)) # extracts hostID numbers (note that C counts from 0 whilst R counts from 1)
hostIDRaster[cells] <- hostData$hostID # attaches hostIDs

rm(list=c("basemapCells","cells"))

allHostLocsDF <- as.data.frame(rasterToPoints(hostIDRaster)) # converting to a dataframe with coordinates
names(allHostLocsDF) <- c("easting","northing","hostID") # attach variable names
allHostLocsDF <- merge(allHostLocsDF,hostData,by="hostID",all.x=TRUE) # attaches column, row, and citrus density data

##### including simulation model summary data -----

modelSimulations = list() # to store outputs
for (runID in 0:(nRuns-1)){ # for all simulation runs
  sFile <- paste(runDir,"\\",runDir,"_",runID,".txt",sep="") # identifying epi_Final file
  runData <- read.csv(sFile,sep=" ",header=F)  # load in file
  if (nrow(runData)>0){ # only for realisations where infection occurred:
    names(runData)<-c("hostCol","hostRow","tInf","infType","sourceCol","sourceRow","propFull","relInf","relSus","relPri","j","numCells","hostID","totalFull","finalProp") # attach names
    runData$simNum <- (runID+1) # attach the number of the simulation
  } 
  modelSimulations[[(runID+1)]] <- runData # stores in list
} # end of loop through all epi_Final files

epiFinalRuns <- do.call(rbind.data.frame, modelSimulations) # create a single dataframe with all information
rm(modelSimulations)

# extracting site-specific information:
allLocsRiskSum <- epiFinalRuns %>%
  group_by(hostID,hostCol,hostRow) %>% 
  summarise(infProb=n()/nRuns,
            meanRisk=sum(finalProp/nRuns)) %>%
  ungroup()
allLocsRiskSum <- data.frame(allLocsRiskSum)

allHostLocsDF <- merge(allHostLocsDF,allLocsRiskSum, all.x=TRUE) # merge with master dataframe

rm(allLocsRiskSum)

##### information on final sampling pattern -----

patternData <- read.csv(paste(saDir,"\\objectiveFunction.txt",sep=""), sep=" ",header=F)
selectedSites <- data.frame(hostID=as.numeric(patternData[nrow(patternData),-(1:2)]),
                       selected=1)

allHostLocsDF <- merge(allHostLocsDF,selectedSites, all.x=TRUE)

allHostLocsDF <- allHostLocsDF %>%
  mutate(infProb=replace_na(infProb,0),
         meanRisk=replace_na(meanRisk,0),
         selected=replace_na(selected,0))

rm(list=c("patternData","selectedSites"))
   
##### Saving the output -----

dir.create("output") # create an output folder

ggplot(data=allHostLocsDF, aes(x=easting, y=northing, fill=propFull)) +
  geom_tile() +
  scale_fill_viridis(name="Citrus density", direction=1, discrete=FALSE, option="E", begin=0.5, end=1.0) +
  geom_point(data=filter(allHostLocsDF,selected==1)) + # comment out if optimised locations not wanted
  coord_equal() +
  theme_map() +
  theme(legend.position = c(0.2, 0.2),
        legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(size = 15)) 
ggsave("output/propFull.tiff") # plot of citrus density

ggplot(data=allHostLocsDF, aes(x=easting, y=northing, fill=infProb)) +
  geom_tile() +
  scale_fill_viridis(name="Probability of infection", direction=1, discrete=FALSE, option="E", begin=0.5, end=1.0) +
  geom_point(data=filter(allHostLocsDF,selected==1)) + # comment out if optimised locations not wanted
  coord_equal() +
  theme_map() +
  theme(legend.position = c(0.2, 0.2),
        legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(size = 15)) 
ggsave("output/infProb.tiff") # plot of infection probability

ggplot(data=allHostLocsDF, aes(x=easting, y=northing, fill=meanRisk)) +
  geom_tile() +
  scale_fill_viridis(name="Mean end prevalence", direction=1, discrete=FALSE, option="E", begin=0.5, end=1.0) +
  geom_point(data=filter(allHostLocsDF,selected==1)) + # comment out if optimised locations not wanted
  coord_equal() +
  theme_map() +
  theme(legend.position = c(0.2, 0.2),
        legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(size = 15)) 
ggsave("output/meanRisk.tiff") # plot of mean end prevalence

write_csv(allHostLocsDF,path="output/allHostLocsDF.csv") # save the data frame of sites