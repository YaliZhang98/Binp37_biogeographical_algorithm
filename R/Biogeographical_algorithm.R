# This script builds a algorithm than can localize origin of ancient individuals according to thier allele frequencies patterns.
# It also contain visualization of results and accuracy detection part.

library(nominatim)
library(mapdataNE)
library(plyr)
library(dplyr)
library(sp)
library(rworldmap) 
library(caret)  
library(rpart)
library(maps)
library(MASS)
library(geosphere)
library(doParallel)
library(forcats)
library(DMwR)
library(png)

# Dataset preparation ----------------------------------------------------------

Dataset <- read.csv(file="34_gene_pool.csv",header=TRUE)

dim(Dataset)
sapply(Dataset,class)

Dataset$Country <- make.names(Dataset$Country)
Dataset$Continent_detail <- make.names(Dataset$Continent_detail)
Dataset$Continent <- make.names(Dataset$Continent)

Dataset$Continent_detail <- factor(Dataset$Continent_detail)
Dataset$Continent <- factor(Dataset$Continent)
Dataset$Country <- factor(Dataset$Country)

# longitude adjustment : adjust longitude near the boundary especially Oceania and east Asia
for (i in (1:length(Dataset$Long.))){
  if (Dataset[i,"Continent"] == "Asia" & Dataset[i,"Long."] < 0 ){
    Dataset[i,"adj_Long."] <- 360 + Dataset[i,"Long."]
   }else{
     Dataset[i,"adj_Long."] <- Dataset[i,"Long."]
   }
  if (Dataset[i,"Continent"] == "Oceania" & Dataset[i,"Long."] < 0){
      Dataset[i,"adj_Long."] <- 360 + Dataset[i,"Long."]
  }
}


# Extract headers of components
optimumVars <- names(Dataset)[7:40] # should contain all the components columns


# model construction -----------------------------------------------------------
ML_Location <-  function(training,testing,classTarget,variables){
  
  set.seed(1234) # create a same random numeric list
  
  # create 5 folds for 5 folds cross-validation
  folds <- createFolds(training[,classTarget], k = 5, returnTrain = T) # used in index item in parameters below (trainControl)
  
  # used in train model
  trControlClass <-  trainControl(
    method = "cv",
    number = 5,  
    verboseIter = FALSE,
    returnData = FALSE,
    search = "grid",
    savePredictions = "final",
    classProbs = T, # different with trControl
    allowParallel = T,
    index = folds )
  
  # used in train model
  trControl <-  trainControl(
    method = "cv",
    number = 5,  
    verboseIter = FALSE,
    returnData = FALSE,
    search = "grid",
    savePredictions = "final",
    allowParallel = T,
    index = folds)
  
  # used in train model
  tune_grid <- expand.grid(
    nrounds = c(400,600),
    eta = c( 0.05, 0.1),
    max_depth = c(3,6,9),
    gamma = 0,
    colsample_bytree = c(0.6,0.8),
    min_child_weight = c(1),
    subsample = (0.7)
  )
  
  
  ##### model training part ____________________
  
  # prediction model for continent
  Xgb_region <- train(x = training[,variables],y = training[,"Continent"], 
                      method = "xgbTree", # train model selection
                      trControl = trControlClass,
                      tuneGrid = tune_grid,
                      nthread = 1)
 
  print('continent training model has been done')
  
  l1_train <- data.frame(training[,c(variables)],Xgb_region[["pred"]][order(Xgb_region$pred$rowIndex),levels(training[,"Continent"]) ])
  
  # prediction model for continent subregions
  Xgb_class <- train(x = l1_train,y = training[,classTarget],
                     method = "xgbTree",
                     trControl = trControlClass,
                     tuneGrid = tune_grid,
                     nthread = 1)
  
  print('continent detail training model has been done')
  
  l2_train <- data.frame(l1_train,Xgb_class[["pred"]][order(Xgb_class$pred$rowIndex),levels(training[,classTarget]) ])

  ##### country prediction: country label can be added into chain prediction
  # Xgb_country <- train(x = l2_train,y = training[,'Country'],
  #                      method = "xgbTree",
  #                      trControl = trControlClass,
  #                      tuneGrid = tune_grid,
  #                      nthread = 1)
  # 
  # print('country training model has been done')
  # 
  # l22_train <- data.frame(l2_train,Xgb_country[["pred"]][order(Xgb_country$pred$rowIndex),levels(training[,'Country']) ])
  # 
  #
  ##### date prediction: age label (date mean in BP) can be added into chain prediction
  # Xgb_date <- train(x = l2_train,y = training[,"Date.mean.in.BP"],
  #                   method = "xgbTree",
  #                   trControl = trControl, #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #                   tuneGrid = tune_grid,
  #                   nthread = 1)
  # 
  # print('date training model has been done')
  # 
  # l_date_train <- data.frame(l2_train, "datePred" = Xgb_date[["pred"]][order(Xgb_date$pred$rowIndex),"pred" ])
 
  # prediction model for latitude
  Xgb_latitude <- train(x = l2_train ,y = training[,"Lat."], 
                        method = "xgbTree",
                        trControl = trControl,
                        tuneGrid = tune_grid,
                        nthread = 1)
  
  print('latitude training model has been done')
  
  l3_train <- data.frame(l2_train, "latPred" = Xgb_latitude[["pred"]][order(Xgb_latitude$pred$rowIndex),"pred" ])
  
  # prediction model for longitude
  Xgb_longitude <- train(x = l3_train ,y = training[,"adj_Long."],
                         method = "xgbTree",
                         trControl = trControl,
                         tuneGrid = tune_grid,
                         nthread = 1)
  
  print('longtitude training model has been done')
  
  
  ##### prediction part ___________________
  
  # predict continent
  regProbs <- predict(Xgb_region, newdata = testing[,variables],type ="prob")
  print('continent prediction has been done')
  
  l1_test <- data.frame(testing[,variables], regProbs)
  
  # predict continent subregion
  classPred <- predict(Xgb_class, newdata = l1_test)
  classProbs <- predict(Xgb_class, newdata = l1_test,type ="prob")
  print('continent detail prediction has been done')
  
  l2_test <-  data.frame(l1_test, classProbs) 

  ##### predict country
  # countryPred <- predict(Xgb_country, newdata = l2_test)
  # countryProbs <- predict(Xgb_country, newdata = l2_test,type ="prob")
  # print('country prediction has been done')
  # 
  # l22_test <-  data.frame(l2_test, countryProbs)
  # 
  ##### predict age (date mean in BP)
  # datePred <- predict(Xgb_date, newdata = l2_test)
  # print('date prediction has been done')
  # date_test <- data.frame(l2_test, datePred)
  
  # predict latitude
  latPred <- predict(Xgb_latitude, newdata = l2_test)
  print('latitude prediction has been done')
  
  l3_test <- data.frame(date_test, latPred)
  
  # predict longitude
  longPred <- predict(Xgb_longitude, newdata = l3_test)
  print('longtitude prediction has been done')
  
  return(list(classPred, latPred, longPred))
  
}


# Prediction of original coordinates for ancient individuals -------------------
set.seed(18)
trainFolds <-  createFolds(Dataset$Continent_detail, k = 5, returnTrain = T)
GeoPreds <- list()
registerDoParallel(7) 

for (i in 1:5){ 
  
  print(i)
  print('calculating ..............................................')
  
  train <- Dataset[trainFolds[[i]],] # select 4/5 individuals to model training part
  test <- Dataset[-trainFolds[[i]],] # exclude 1/5 individuals form model training part and predict thier origin
  
  # SMOTE: increase the number of individuals in region with a small sample size
  for (j in (1:18)){
    
    name <- as.data.frame(table(train$Continent_detail))
    
    if (min(name$Freq) < 350){
      min_continent <- name[name$Freq == min(name$Freq),]$Var1
      N <- sum(train$Continent_detail == min_continent)
      
      if (min(name$Freq) < 100){
        a <- 600 %/% N * 100
      }
      if (min(name$Freq) > 100 & min(name$Freq) < 300){
        a <- 500 %/% N * 100
      }
      if (min(name$Freq) > 300 & min(name$Freq) < 380){
        a <- 600 %/% N * 100
      }
      
      b <- 10000 * (length(train$Continent_detail) - N) / a / N 
      train2 <- SMOTE(Continent_detail~.,train,perc.over=a, perc.under=b)
      new_ind <- train2[train2$Continent_detail == min_continent, ]
      train_other <- train[train$Continent_detail != min_continent, ]
      train <- rbind(train_other,new_ind)
    }
  }
  print('SMOTE have been done !')
  
  # predict origin of individuals through constructed model 
  testPreds <- ML_Location(training = train, testing = test, classTarget = "Continent_detail",variables = optimumVars)
  GeoPreds[[i]] <- testPreds
  
  print(i) 
  print("finished !!!!!! -----------------------------------------")
}

# merge the prediction results of all individuals
add_preds <- list()
for (i in 1:5){
  add_preds[[i]] <- cbind(Dataset[-trainFolds[[i]],] , 
                          "regionPred"= GeoPreds[[i]][[1]],
                          # "datePred" = GeoPreds[[i]][[2]],
                          "latPred" = GeoPreds[[i]][[2]], 
                          "longPred" = GeoPreds[[i]][[3]] )
}

DataPreds <- rbind.fill(add_preds)

# adjust longitude to right expression
for (i in (1:length(DataPreds$longPred))){
  if (DataPreds[i,"longPred"] > 180){
    DataPreds[i,"adj_longPred"] <- DataPreds[i,"longPred"] - 360
  }else{
    DataPreds[i,"adj_longPred"] <- DataPreds[i,"longPred"]
  }
}

write.csv(DataPreds,"36_DataPreds.csv",row.names = T)

# coast adjustment
# get world coastlines
coastlines <- cbind("x"  = SpatialLines2map(coastsCoarse)$x ,"y" =SpatialLines2map(coastsCoarse)$y)
coastlines <- coastlines[complete.cases(coastlines),]
coastlines <- coastlines[coastlines[,1] < 180 ,]

# Function that can find the nearest land point for individuals whose predicted location is in sea
find_coast <- function(long,lat){
  distances_from_coastline <-  spDistsN1(coastlines , c(long,lat), longlat = TRUE)
  closest_point <-  which.min(distances_from_coastline)
  new_coords <- coastlines[closest_point,]
  return(new_coords)
}

# locate predicted coordinates of individuals into country
locate_country <- map.where(database = "world", DataPreds$adj_longPred, DataPreds$latPred)
adjust_land <- DataPreds[which(is.na(locate_country)),]

# adjust location in sea to nearest land
adjust_coordinate <- mapply(find_coast, long = adjust_land$adj_longPred, lat = adjust_land$latPred )
DataPreds[which(is.na(locate_country)), "latPred"] <- adjust_coordinate[2,]
DataPreds[which(is.na(locate_country)), "adj_longPred"] <- adjust_coordinate[1,]

write.csv(DataPreds,"36_DataPreds_adjustcoast.csv",row.names = T)


# result visualization ---------------------------------------------------------

# plot excavation location of individuals on world map

palette <-c( "brown","red3","maroon1","slateblue1","olivedrab1","darkorchid4","orange","cyan4",
             "cyan2","mediumorchid1","seashell2","dodgerblue3","violetred3","palevioletred1",
             "rosybrown1","mediumspringgreen","gold2","orangered")

map <- getMap(resolution = "coarse")

size1 <- c()
for (i in 1: length(levels(DataPreds$Continent_detail))){
  
  this_continent <- levels(DataPreds$Continent_detail)[i]
  size1[i] <- length(which(DataPreds$Continent_detail == this_continent))
}

png("origin_global_map.png", width = 13,height = 8, units = 'in', res = 600)
plot(map,xlim = c(-160,160), col = "grey", border = "darkgrey", bg = "lightskyblue1", xlab = "", ylab = "")
title(ylab="Latitude", mgp=c(2,1,0),cex.lab=1.2)

for (i in 1:length(levels(DataPreds$Continent_detail))){
  this_continent <- levels(DataPreds$Continent_detail)[i]
  find_lats <- DataPreds[DataPreds[,"Continent_detail"] == this_continent,]$Lat.
  find_longs <- DataPreds[DataPreds[,"Continent_detail"] == this_continent,]$Long.
  points(find_longs, find_lats, col = palette[i], pch = "+", cex = 1.1)
}

legend(-165,20,legend=c(legend=c(paste0(levels(DataPreds$Continent_detail),"  (",size1,")"))),
       col=palette,pch = "+",cex=1, bg = "lightskyblue1")  

map.axes()
dev.off()

# plot prediction location of all individuals on world map

png("prediction_global_map.png", width = 13,height = 8, units = 'in', res = 600)
plot(map,xlim = c(-160,160), col = "grey", border = "darkgrey", bg = "lightskyblue1", xlab = "", ylab = "")
title(ylab="Latitude",xlab = "Longitude", mgp=c(2,1,0),cex.lab=1.2)

for (i in 1:length(levels(DataPreds$Continent_detail))){
  this_continent <- levels(DataPreds$Continent_detail)[i]
  find_lats <- DataPreds[DataPreds[,"Continent_detail"] == this_continent,]$latPred
  find_longs <- DataPreds[DataPreds[,"Continent_detail"] == this_continent,]$adj_longPred
  points(find_longs, find_lats, col = palette[i], pch = "+", cex = 1.1)
}

legend(-165,20,legend=c(paste0(levels(DataPreds$Continent_detail),"  (",size1,")")),
       col=palette,pch = "+",cex=1, bg = "lightskyblue1")  

map.axes()
dev.off()


# distance between excavation coordinates and prediction coordinates

# continent subregions as unit

for (i in 1:length((DataPreds$Continent_detail))){ 
  region_lats <- DataPreds[i,]$Lat.
  region_longs <- DataPreds[i,]$Long.
  pred_lat <- DataPreds[i,]$latPred
  pred_long  <- DataPreds[i,]$adj_longPred
  distance <- c()
  for (n in 1:length(region_lats)){ 
    distance[n] <- distHaversine(c(pred_long ,pred_lat ),c(region_longs[n],region_lats[n]))/1000
  }
  DataPreds[i,"distance_from_continent"] <- min(distance, na.rm = TRUE)
}

bar_df <- data.frame(row.names = c( "Overall",levels(DataPreds$Continent_detail)))

for (i in 1: length(levels(DataPreds$Continent_detail))){
  overall_prop <- mean(DataPreds[,"distance_from_continent"] < 200)
  bar_df[1,"0 - 200km"] <- overall_prop
  
  this_continent <- levels(DataPreds$Continent_detail)[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] < 200)
  bar_df[i+1,"0 - 200km"] <- prop
}

for (i in 1: length(levels(DataPreds$Continent_detail))){
  this_continent <- levels(DataPreds$Continent_detail)[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] > 200 & DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] < 500)
  bar_df[i+1,"200 - 500km"] <- prop
  
  overall_prop <- mean(DataPreds[,"distance_from_continent"] > 200 & DataPreds[,"distance_from_continent"] < 500)
  bar_df[ 1,"200 - 500km"] <- overall_prop
}

for (i in 1: length(levels(DataPreds$Continent_detail))){
  this_continent <- levels(DataPreds$Continent_detail)[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] > 500 & DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] < 1000)
  bar_df[i+1,"500 - 1000km"] <- prop
  
  overall_prop <- mean(DataPreds[,"distance_from_continent"] > 500 & DataPreds[,"distance_from_continent"] < 1000)
  bar_df[ 1,"500 - 1000km"] <- overall_prop
}

for (i in 1: length(levels(DataPreds$Continent_detail))){
  this_continent<- levels(DataPreds$Continent_detail)[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] > 1000 & DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] < 2000)
  bar_df[i+1,"1000 - 2000km"] <- prop
  
  overall_prop <- mean(DataPreds[,"distance_from_continent"] > 1000 & DataPreds[,"distance_from_continent"] < 2000)
  bar_df[ 1,"1000 - 2000km"] <- overall_prop
}

for (i in 1: length(levels(DataPreds$Continent_detail))){
  this_continent <- levels(DataPreds$Continent_detail)[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] > 2000 & DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] < 3000)
  bar_df[i+1,"2000 - 3000km"] <- prop
  
  overall_prop <- mean(DataPreds[,"distance_from_continent"] > 2000 & DataPreds[,"distance_from_continent"] < 3000)
  bar_df[ 1,"2000 - 3000km"] <- overall_prop
}

for (i in 1: length(levels(DataPreds$Continent_detail))){
  this_continent <- levels(DataPreds$Continent_detail)[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"distance_from_continent"] > 3000 )
  bar_df[i+1,"> 3000km"] <- prop
  
  overall_prop <- mean(DataPreds[,"distance_from_continent"] > 3000)
  bar_df[ 1,"> 3000km"] <- overall_prop
}

  ## change the order of columns

bar_df2 <- bar_df[c("Overall","Central.Africa","North.Africa","East.Africa","South.Africa",
                     "Middle.America","North.America","South.America",
                     "Central.Asia","North.Asia","South.Asia","East.Asia","Southeast.Asia","West.Asia",
                     "North.Europe","South.Europe","East.Europe","West.Europe","Oceania"),]

order <- c("Overall","Central.Africa","North.Africa","East.Africa","South.Africa",
           "Middle.America","North.America","South.America",
           "Central.Asia","North.Asia","South.Asia","East.Asia","Southeast.Asia","West.Asia",
           "North.Europe","South.Europe","East.Europe","West.Europe","Oceania")

size1 <- c()
size1[1] <- length(DataPreds$Continent)
for (i in 2: length(order)){
  this_continent <- order[i]
  size1[i] <- length(which(DataPreds$Continent_detail == this_continent))
}

par(xpd = T, mar = par()$mar + c(1,0,0,7))
bp <- barplot(t(bar_df2*100), col=c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), 
              names.arg=c(paste0(order,"  (",size1,")")) ,args.legend = list(x = "topright", inset=c(-0.5,0)), las =2, cex.names=.7,ylab = "Proportion of sample predictions %")
legend("topright",inset = c(-0.2,0.4), rev(c(colnames(bar_df2))), fill = rev(c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) , bty = 1, cex = 0.8)
par(mar=c(5, 4, 4, 2) + 0.1)

png("bar_distance_region.png", width = 13,height = 8, units = 'in', res = 600)

par(xpd = T, mar = par()$mar + c(1,0,0,7), mgp = c(0,0.7,0), las=2)
bp <- barplot(t(bar_df2*100), col=c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"),
              names.arg=c(paste0(order,"  (",size1,")")) ,
              args.legend = list(x = "topright", inset=c(-0.5,0)), las =2,
              cex.names=.6,ylab = "", axisnames = F,axes = F, space =0)
axis(side =2, pos = 0)
mtext(text = c(paste0(order,"  (",size1,")")) , side = 1, at = bp, line = 0, padj = 1, cex = 0.7)
title(ylab="Proportion of sample predictions %", mgp=c(0,0,0),cex.lab=1.2)
legend("topright",inset = c(-0.12,0.4), rev(c(colnames(bar_df2))), fill = rev(c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) , bty = 1, cex = 1)

par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()


# country as unit

for (i in 1:length((DataPreds$Country))){

  country_lats <- DataPreds[i,]$Lat.
  country_longs <- DataPreds[i,]$Long.
  pred_lat <- DataPreds[i,]$latPred
  pred_long  <- DataPreds[i,]$adj_longPred
  distance <- c()
  for (n in 1:length(country_lats)){
    distance[n] <- distHaversine(c(pred_long ,pred_lat ),c(country_longs[n],country_lats[n]))/1000
  }
  DataPreds[i,"distance_from_country"] <- min(distance, na.rm = TRUE)
}

bar_df <- data.frame(row.names = c( "Overall_country",levels(DataPreds$Country)))

for (i in 1: length(levels(DataPreds$Country))){
  overall_prop <- mean(DataPreds[,"distance_from_country"] < 200)
  bar_df[1,"0 - 200km"] <- overall_prop

  this_country <- levels(DataPreds$Country)[i]
  prop <- mean(DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] < 200)
  bar_df[i+1,"0 - 200km"] <- prop
}

for (i in 1: length(levels(DataPreds$Country))){
  this_country <- levels(DataPreds$Country)[i]
  prop <- mean(DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] > 200 & DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] < 500)
  bar_df[i+1,"200 - 500km"] <- prop

  overall_prop <- mean(DataPreds[,"distance_from_country"] > 200 & DataPreds[,"distance_from_country"] < 500)
  bar_df[ 1,"200 - 500km"] <- overall_prop
}


for (i in 1: length(levels(DataPreds$Country))){
  this_country <- levels(DataPreds$Country)[i]
  prop <- mean(DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] > 500 & DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] < 1000)
  bar_df[i+1,"500 - 1000km"] <- prop

  overall_prop <- mean(DataPreds[,"distance_from_country"] > 500 & DataPreds[,"distance_from_country"] < 1000)
  bar_df[ 1,"500 - 1000km"] <- overall_prop
}

for (i in 1: length(levels(DataPreds$Country))){
  this_country<- levels(DataPreds$Country)[i]
  prop <- mean(DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] > 1000 & DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] < 2000)
  bar_df[i+1,"1000 - 2000km"] <- prop

  overall_prop <- mean(DataPreds[,"distance_from_country"] > 1000 & DataPreds[,"distance_from_country"] < 2000)
  bar_df[ 1,"1000 - 2000km"] <- overall_prop
}

for (i in 1: length(levels(DataPreds$Country))){
  this_country <- levels(DataPreds$Country)[i]
  prop <- mean(DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] > 2000 & DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] < 3000)
  bar_df[i+1,"2000 - 3000km"] <- prop

  overall_prop <- mean(DataPreds[,"distance_from_country"] > 2000 & DataPreds[,"distance_from_country"] < 3000)
  bar_df[ 1,"2000 - 3000km"] <- overall_prop
}

for (i in 1: length(levels(DataPreds$Country))){
  this_country <- levels(DataPreds$Country)[i]
  prop <- mean(DataPreds[DataPreds$Country == this_country,][,"distance_from_country"] > 3000 )
  bar_df[i+1,"> 3000km"] <- prop

  overall_prop <- mean(DataPreds[,"distance_from_country"] > 3000)
  bar_df[ 1,"> 3000km"] <- overall_prop
}

size1 <- c()
for (i in 1: length(levels(DataPreds$Country))){

  this_country <- levels(DataPreds$Country)[i]
  size1[i] <- length(which(DataPreds$Country == this_country))
}

par(xpd = T, mar = par()$mar + c(1,0,0,7))
bp <- barplot(t(bar_df*100), col=c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), names.arg=c("Overall",paste0(levels(DataPreds$Country),"  (",size1,")")) ,args.legend = list(x = "topright", inset=c(-0.5,0)), las =2, cex.names=.6,ylab = "Proportion of sample predictions %")
legend("topright",inset = c(-0.15,0.4), rev(c(colnames(bar_df))), fill = rev(c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) , bty = 1, cex = 0.8)

par(mar=c(5, 4, 4, 2) + 0.1)

png("36_date_SMOTE_bar_distance_country.png", width = 13,height = 8, units = 'in', res = 600)

par(xpd = T, mar = par()$mar + c(1,0,0,7), mgp = c(0,0.7,0), las=2)
bp <- barplot(t(bar_df*100), col=c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"),
              names.arg=c("Overall",paste0(levels(DataPreds$Country),"  (",size1,")")) ,
              args.legend = list(x = "topright", inset=c(-0.5,0)), las =2,
              cex.names=.6,ylab = "", axisnames = F,axes = F, space =0)

axis(side =2, pos = 0)
mtext(text = c("Overall",paste0(levels(DataPreds$Country),"  (",size1,")")) , side = 1, at = bp, line = 0, padj = 1, cex = 0.35)
title(ylab="Proportion of sample predictions %", mgp=c(0,0,0),cex.lab=1.2)
legend("topright",inset = c(-0.12,0.4), rev(c(colnames(bar_df))), fill = rev(c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) , bty = 1, cex = 1)

par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()


# If age (date mean in BP) label has been added in prediction model

# region age distance

for (i in 1:length(DataPreds$Date.mean.in.BP)){
  real_date <- DataPreds[i,]$Date.mean.in.BP
  pred_date <- DataPreds[i,]$datePred
  DataPreds[i,"date_distance"] <- abs(pred_date - real_date)
}

bar_df <- data.frame(row.names = c(order))

for (i in 2: length(order)){
  overall_prop <- mean(DataPreds[,"date_distance"] < 200)
  bar_df[1,"0 - 200 BP"] <- overall_prop

  this_continent <- order[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] < 200)
  bar_df[i,"0 - 200 BP"] <- prop
}

for (i in 2: length(order)){
  this_continent <- order[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] > 200 & DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] < 500)
  bar_df[i,"200 - 500 BP"] <- prop

  overall_prop <- mean(DataPreds[,"date_distance"] > 200 & DataPreds[,"date_distance"] < 500)
  bar_df[ 1,"200 - 500 BP"] <- overall_prop
}

for (i in 2: length(order)){
  this_continent <- order[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] > 500 & DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] < 1000)
  bar_df[i,"500 - 1000 BP"] <- prop

  overall_prop <- mean(DataPreds[,"date_distance"] > 500 & DataPreds[,"date_distance"] < 1000)
  bar_df[ 1,"500 - 1000 BP"] <- overall_prop
}

for (i in 2: length(order)){
  this_continent<- order[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] > 1000 & DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] < 2000)
  bar_df[i,"1000 - 2000 BP"] <- prop

  overall_prop <- mean(DataPreds[,"date_distance"] > 1000 & DataPreds[,"date_distance"] < 2000)
  bar_df[ 1,"1000 - 2000 BP"] <- overall_prop
}

for (i in 2: length(order)){
  this_continent <- order[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] > 2000 & DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] < 3000)
  bar_df[i,"2000 - 3000 BP"] <- prop

  overall_prop <- mean(DataPreds[,"date_distance"] > 2000 & DataPreds[,"date_distance"] < 3000)
  bar_df[ 1,"2000 - 3000 BP"] <- overall_prop
}

for (i in 2: length(order)){
  this_continent <- order[i]
  prop <- mean(DataPreds[DataPreds$Continent_detail == this_continent,][,"date_distance"] > 3000 )
  bar_df[i,"> 3000 BP"] <- prop

  overall_prop <- mean(DataPreds[,"date_distance"] > 3000)
  bar_df[ 1,"> 3000 BP"] <- overall_prop
}

bar_df2 <- bar_df[c("Overall","Central.Africa","North.Africa","East.Africa","South.Africa",
                    "Middle.America","North.America","South.America",
                    "Central.Asia","North.Asia","South.Asia","East.Asia","Southeast.Asia","West.Asia",
                    "North.Europe","South.Europe","East.Europe","West.Europe","Oceania"),]

size1 <- c()
size1[1] <- length(DataPreds$Continent)
for (i in 2: length(order)){
  this_continent <- order[i]
  size1[i] <- length(which(DataPreds$Continent_detail == this_continent))
}

par(xpd = T, mar = par()$mar + c(1,0,0,7))

bp <- barplot(t(bar_df2*100), col=c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"), names.arg=c(paste0(order,"  (",size1,")")) ,args.legend = list(x = "topright", inset=c(-0.5,0)), las =2, cex.names=.7,ylab = "Proportion of sample predictions %")
legend("topright",inset = c(-0.2,0.4), rev(c(colnames(bar_df2))), fill = rev(c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) , bty = 1, cex = 0.8)

par(mar=c(5, 4, 4, 2) + 0.1)

png("age_bar_region.png", width = 13,height = 8, units = 'in', res = 600)

par(xpd = T, mar = par()$mar + c(1,0,0,7), mgp = c(0,0.7,0), las=2)
bp <- barplot(t(bar_df2*100), col=c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue"),
              names.arg=c(paste0(order,"  (",size1,")")) ,
              args.legend = list(x = "topright", inset=c(-0.5,0)), las =2,
              cex.names=.6,ylab = "", axisnames = F,axes = F, space =0)

axis(side =2, pos = 0)
mtext(text = c(paste0(order,"  (",size1,")")) , side = 1, at = bp, line = 0, padj = 1, cex = 0.7)
title(ylab="Proportion of sample predictions %", mgp=c(0,0,0),cex.lab=1.2)
legend("topright",inset = c(-0.15,0.4), rev(c(colnames(bar_df2))), fill = rev(c("lightcyan","slategray1","lightblue", "skyblue", "royalblue3", "darkblue")) , bty = 1, cex = 1)

par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()