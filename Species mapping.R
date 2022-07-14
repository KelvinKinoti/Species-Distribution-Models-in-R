install.packages("dismo")
install.packages("maptools")
install.packages("rgdal")
install.packages("raster")
install.packages("sp")

library("sp")
library("raster")
library("maptools")
library("rgdal")
library("dismo")

dir.create(path = "data")
dir.create(path = "output")

bioclim.data<-raster::getData(name = "worldclim",
                        var = "bio",
                        res = 10,
                        path = "data/")
head(bioclim)

obs.data <- read.csv(file = "data/Athene_cunicularia_gbif.csv")
summary(obs.data)# No NA's so proceed
summary(bioclim.data)
bioclim.data[is.na(bioclim.data)] <- 0
# Determine the limits of the data
max.lat <- ceiling(max(obs.data$Latitude))
min.lat <- floor(min(obs.data$Latitude))
max.lon <- ceiling(max(obs.data$Longitude))
min.lon <- floor(min(obs.data$Longitude))
geographic.extent <- extent(x = c(min.lon, max.lon, min.lat, max.lat))

# Load the data to use for our base map
data(wrld_simpl)

plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")

points(x = obs.data$Longitude, 
       y = obs.data$Latitude, 
       col = "olivedrab", 
       pch = 20, 
       cex = 0.75)
box()

##MAP 2
#Building and visualizing a model
# Crop bioclim data 
bioclim.data <- crop(x = bioclim.data, y = geographic.extent)
#Replae NA's with zero
bioclim.data[is.na(bioclim.data)] <- 0
is.na(bioclim.data)

# Drop unused column and interchange the columns
obs.data <- obs.data[,2:3]
obs.data <- obs.data[, c("Longitude", "Latitude")]
head(obs.data)
# Build species distribution model
bc.model <- bioclim(x = bioclim.data, p = obs.data)
#Prediction from the model"s probability using the dismo package
predict.presence <- dismo::predict(object = bc.model, 
                                   x = bioclim.data, 
                                   ext = geographic.extent,na.rm=TRUE)
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")
# Add model probabilities
plot(predict.presence, add = TRUE)
# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
# Add original observations
points(obs.data$Longitude, obs.data$Latitude, col = "olivedrab", pch = 20, cex = 0.75)
box()

##MAP 3
bil.files <- list.files(path = "data/wc2-5", 
                        pattern = "*.bil$", 
                        full.names = TRUE)
# We only need one file, so use the first one in the list of .bil files
mask <- raster(bil.files[1])
# Set the seed for the random-number generator to ensure results are similar
set.seed(20210707)
# Randomly sample points (same number as our observed points)
background <- randomPoints(mask = mask, # Provides resolution of sampling points
                           n = nrow(obs.data), # Number of random points
                           ext = geographic.extent, # Spatially restricts sampling
                           extf = 1.25) # Expands sampling a little bit##Step 2
# Arbitrarily assign group 1 as the testing data group
testing.group <- 1
# Create vector of group memberships
group.presence <- kfold(x = obs.data, k = 5) # kfold is in dismo package
# Should see even representation in each group
table(group.presence)

# Separate observations into training and testing groups
presence.train <- obs.data[group.presence != testing.group, ]
presence.test <- obs.data[group.presence == testing.group, ]
# Repeat the process for pseudo-absence points
group.background <- kfold(x = background, k = 5)
background.train <- background[group.background != testing.group, ]
background.test <- background[group.background == testing.group, ]

#Train and test model
#Build model
bc.model <- bioclim(x = bioclim.data, p = presence.train)
# Predict presence from model (same as previously, but with the update model)
predict.presence <- dismo::predict(object = bc.model, 
                                   x = bioclim.data, 
                                   ext = geographic.extent)

# Use testing data for model evaluation
bc.eval <- evaluate(p = presence.test, # The presence testing data
                    a = background.test, # The absence testing data
                    model = bc.model, # The model we are evaluating
                    x = bioclim.data) # Climatic variables for use by model
# Determine minimum threshold for "presence"
bc.threshold <- threshold(x = bc.eval, stat = "spec_sens")

# Plot base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")
# Only plot areas where probability of occurrence is greater than the threshold
plot(predict.presence > bc.threshold, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "olivedrab"))
# And add those observations
points(x = obs.data$Longitude, 
       y = obs.data$Latitude, 
       col = "black",
       pch = "+", 
       cex = 0.75)
# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()

## Map 4
forecast.data<-getData(name = "CMIP5", 
                         var = "bio",
                         res = 2.5, 
                         path = "data/", 
                         model = "GD", 
                         rcp = "45", 
                         year = 70) 

names(forecast.data) <- names(bioclim.data)
# Predict presence from model with forecast data
forecast.presence <- dismo::predict(object = bc.model, 
                                    x = forecast.data, 
                                    ext = geographic.extent)

# Plot base map
plot(wrld_simpl, 
     xlim = c(min.lon, max.lon),
     ylim = c(min.lat, max.lat),
     axes = TRUE, 
     col = "grey95")
# Only plot areas where probability of occurrence is greater than the threshold
plot(forecast.presence > bc.threshold, 
     add = TRUE, 
     legend = FALSE, 
     col = c(NA, "olivedrab"))
# And add those observations
points(x = obs.data$longitude, 
       y = obs.data$latitude, 
       col = "black",
       pch = "+", 
       cex = 0.6)
# Redraw those country borders
plot(wrld_simpl, add = TRUE, border = "grey5")
box()