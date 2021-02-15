library(tidyverse) #For generic tools 
library(sdm) #Species distribution modeling
library(sp)#Point data
library(sf) #raster compatibility
library(raster) #Raster data
library(rgbif) #Calling GBIF API
library(usdm) #variable selection


#Data gathering and preparation-------------------------

#Fetch data from GBIF
gbif_query <- occ_search(scientificName = "Plebejus argus", country = "GB", limit = 30000, 
                                                        fields = c("scientificName", "key", "decimalLatitude", "decimalLongitude", 
                                                                   "issues", "basisOfRecord", "coordinateUncertaintyInMeters"), 
                                                        hasCoordinate = TRUE)

#Query returns a list, subsetting for dataframe
occ_data <- as.data.frame(gbif_query[3])


#Change colnames
colnames(occ_data) <- c("key", "basis", "species", "long", "lat", "variation_m", "issues")


#Cleaning
clean_species <- occ_data %>% 
  as_tibble() %>% 
  filter(species == "Plebejus argus (Linnaeus, 1758)",
         variation_m < 1500 | is.na(variation_m)) %>% 
  mutate(key = as.numeric(key),
         species = str_extract(species, "Plebejus argus"),
         species = recode(species, "Plebejus argus" = "Plebejus_argus"))


#Create df to include just coordinates
point_data <- dplyr::select(clean_species, long, lat) %>% distinct()

#Create Spatial Points Dataframe
#Add 1 to represent species
coordinates(point_data) <- c("long", "lat")
point_data$species <- 1


#Get Climate Data 
bioclim <- raster::getData("worldclim", var = 'bio', res = 10)

#CC Data
future_climate <- raster::getData('CMIP5', var = 'bio', res = 10, year = 70, model = 'CN', rcp = 85)

#Standardize names
names(future_climate) <- names(bioclim)

#Read UK shapefile 
#Data from https://geoportal.statistics.gov.uk/datasets/48b6b85bb7ea43699ee85f4ecd12fd36_0/data?orderBy=st_area(shape)&orderByAsc=false

UK <- st_read("UK_shapefile")
UK <- UK %>% st_union()
UK <- as_Spatial(UK)
UK <- spTransform(UK, crs(bioclim))

#Create extent of the UK shapefile, then to sp object for use with raster 
extent <- UK %>% st_bbox() %>% st_as_sfc() %>% as_Spatial()

#Cropping World Raster Stack to UK BBox
UK_Climate_Data <- crop(bioclim, extent)
future_climate <- crop(future_climate, extent)

#Variable selection------------------------------

#Checking for multicollinearity 
#Strong collinerarity between environmental variables can increase model uncertainty and reduce statistical power 

#Create matrix for vifcor to evaluate problematic variables
vifcor_matrix <- raster::extract(UK_Climate_Data, point_data)

#6 variables with problematic collinerarity 
#1:4 relate to precipitation & dry/wetness, 5:6 relate to min temp in coldest month and annual temp range. No unexpected issues. 
#Default PCC threshold = 0.9, max kept is 0.87.

vifcor_results <- usdm::vifcor(vifcor_matrix)

#Exclude from bioclim raster to proceed
bioclim <- exclude(bioclim, vifcor_results)

#Do same to Future_climate 
future_climate <- exclude(future_climate, vifcor_results)

#Train SDM-------------------------------

#Create a 'recipie' for the SDM. 
#bg denotes parameters to determine the 'absence' to increase model accuracy

#Parameters:
#SVM, Random Forest and GLMnet models produced
#Bootsrap resampling
#75/25 train/test split
#3 subsampling repeats 

recipe <- sdmData(species~., point_data, predictors = bioclim, bg = list(method = "gRandom"), n = 1250)



#Train and save model
model <- sdm(species~., recipe, methods = c("glm", "rf", "brt", "fda"), replication = "boot", test.p = 25, n = 3, 
             parallelSettings = list(ncore=4, method = "parallel"))


#Predictions and Ensembles----------------------
#Use trained model and climatic data to predict probability of occurrence 
#Compare our different models, each has been produced 3 times so need to subset ID's
prediction <- predict(model, UK_Climate_Data)
plot(prediction[[c(1,4,8,12)]])

#Ensemble our different models together 
ensemble <- sdm::ensemble(model, prediction, setting = list(method = 'weighted', stat='tss', opt = 2))
plot(ensemble)


#Repeat for future data 
prediction2 <- predict(model, future_climate)
ensemble2 <- sdm::ensemble(model, prediction2, setting = list(method = 'weighted', stat='tss', opt = 2))
plot(ensemble2)

#Analysis--------------------------

#Compare the difference 
#Uses raster math, indicates increased or decreased probability of occurrence levels 
diff <- ensemble2 - ensemble
plot(diff)


#Convert rasters to dataframes, returns XY and probability of occurence values 
ens1df <-  as.data.frame(ensemble, xy = TRUE, na.rm = T) %>% as_tibble() %>% mutate(model = "Current")
ens2df <-  as.data.frame(ensemble2, xy = TRUE, na.rm = T) %>% as_tibble() %>% mutate(model="CC")

#Join for comparison 
full <- full_join(ens1df, ens2df)

#Compare distribution in current and future climate scenarios 
full %>% 
  ggplot(aes(layer, fill = model))+
  geom_histogram(bins = 20, alpha = 0.6)



#Could even set threshold of 70% for example 
#https://www.r-graph-gallery.com/histogram_several_group.html


#Number of cells suitable increases, but not all currently suitable cells are occupied 
#Current distribution probabbility declines but still remains quite high 
#https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.00825

#Likely to have less affect/be less accurate on small scale species affected by microclimates 


#to github: https://gist.github.com/JoshuaTPierce/b919168421b40e06481080eb53c3fb2f
#video: https://www.youtube.com/watch?v=83dMS3bcjJM&ab_channel=Biogeoinformatics
#Ideally should use the name number of pesudoabsences as presence
#https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2011.00172.x
#https://kevintshoemaker.github.io/NRES-746/SDM_v6.html
#extent object: https://www.rdocumentation.org/packages/raster/versions/3.4-5/topics/extent
#Manually defining extent: https://datacarpentry.org/r-raster-vector-geospatial/11-vector-raster-integration/index.html
#https://www.youtube.com/watch?v=83dMS3bcjJM&ab_channel=Biogeoinformatics
#https://gis.stackexchange.com/questions/285548/selecting-a-subset-of-a-raster-based-on-aspect-and-slope-as-input-to-calculate