library(raster)
library(dplyr)
library(sf)
library(readr)
library(ggplot2)

setwd("D:/TESI_GIADA/Giada/RIPROIETTATO/RHONE")

#Statistical correlation between Kriging and cost map
##importing raster grids
Ccost_obs <- raster("Cost_map/cost_map_river_obstacle/cost_cumu.tif")
Ccost_aid <- raster("Cost_map/river/river_cumulative.tif")
CCost_aid_ani <- raster("Cost_map/river_anisotropic/Rhone_basin_cumulative.tif")
no_riv <- raster("Cost_map/no_river/no_riv_cumulative_cost_map.tif")
Krig <- raster("Kriging/Rhone_UniKriging.tif")
Krig_res <- resample(Krig, Ccost_aid, method = "ngb")
CCost_no_riv <- resample(no_riv, Ccost_aid, method = "ngb")

##calculating correlation between Co-Kriging and Cumulative cost map with river as obstacle
obs_cor <- cor(values(Ccost_obs),
    values(Krig_res),
    use = "na.or.complete")
print(obs_cor)
##calculating correlation between Co-Kriging and Cumulative cost map with river as aid
aid_cor <- cor(values(Ccost_aid),
               values(Krig_res),
               use = "na.or.complete")
print(aid_cor)
##calculating correlation between Co-Kriging and Cumulative cost map with river as aid (anisotropic)
aid_ani_cor <- cor(values(CCost_aid_ani),
                   values(Krig_res),
                   use = "na.or.complete")
print(aid_ani_cor)
##calculating correlation between Co-Kriging and Cumulative cost map with no river
no_riv_cor <- cor(values(CCost_no_riv),
                  values(Krig_res),
                  use = "na.or.complete")
print(no_riv_cor)

#Linear regression
## Krig ~ Cost obstacle
data.r1 <- data.frame(x = values(Krig_res), y = values(Ccost_obs))
lin_reg1 <- lm(x ~ y, data.r1)
summary(lin_reg1)

## Krig ~ Cost aid
data.r2 <- data.frame(x = values(Krig_res), y = values(Ccost_aid))
lin_reg2 <- lm(x ~ y, data.r2)
summary(lin_reg2)

## Krig ~ Cost aid anisotropic
data.r3 <- data.frame(x = values(Krig_res), y = values(CCost_aid_ani))
lin_reg3 <- lm(x ~ y, data.r3)
summary(lin_reg3)

## Krig ~ Cost no river
data.r4 <- data.frame(x = values(Krig_res), y = values(CCost_no_riv))
lin_reg4 <- lm(x ~ y, data.r4)
summary(lin_reg4)

 
 #BIC ANALYSIS
 library(MASS)
 ## Bayesian Information Criterion (BIC): stepwise variable selection
 ### Model 0
 n1 <- nrow(data.r1)
 Regr_model_0BIC<-stepAIC(lin_reg1, k=log(length(n1)))
 summary(Regr_model_0BIC)
 
 ### model 1
 n2 <- nrow(data.r2)
 Regr_model_1BIC<-stepAIC(lin_reg2, k=log(length(n2)))
 summary(Regr_model_1BIC)
 
 ### model 2
 n3 <- nrow(data.r3)
 Regr_model_2BIC<-stepAIC(lin_reg3, k=log(length(n3)))
 summary(Regr_model_2BIC)
 
 ### model 3
 n4 <- nrow(data.r4)
 Regr_model_3BIC<-stepAIC(lin_reg4, k=log(length(n4)))
 summary(Regr_model_3BIC)
 
 ## Create model lists with homogeneous point process model5=Regr_model_5BIC,
 Regr_modelstot<-list(model0=Regr_model_0BIC, model1=Regr_model_1BIC, model2=Regr_model_2BIC, model3=Regr_model_3BIC)
 
 ##Export model summary with homogeneous point process
 write.table(capture.output(print(Regr_modelstot)), file = file.path("Cost_map/statistic_correlation_R", "Regr_modelstot_def1.txt"))
 
 #Compare BIC scores for the different models homogeneous point process
 Regr_modelstot_BICtot<-lapply(Regr_modelstot,BIC)
 write.table(Regr_modelstot_BICtot,file = file.path("Cost_map/statistic_correlation_R", "Regr_modelstot_BIC_scorestot_def1.txt"))
 
 AIC.BIC.weight<-function(x){
   for(i in 1:length(x)){
     x.vect<-as.numeric(c(x[1:i]))}
   delta<-x.vect-min(x.vect)
   L<-exp(-0.5*delta)
   result<-round(L/sum(L),digits=7)
   return(result)
 }
 
 write.table(AIC.BIC.weight(Regr_modelstot_BICtot[1:4]),file = file.path("Cost_map/statistic_correlation_R", "Regr_modelstot_BICtot_weights_def1.txt"))
 
 ### Chi- square Model
 null_data <- data.frame(Krig_res = values(Krig_res))
 null_model <- lm(Krig_res ~ 1, data = null_data)
 null_dev <- deviance(null_model)
 dev <- deviance(lin_reg1)
 
 with(lin_reg1, 1 - pchisq(null_dev - dev, null_model$df.residual - lin_reg1$df.residual)) 
 with(lin_reg1, 1 - pchisq(dev, lin_reg1$df.residual))
 
 ## stepAIC between kriging and covariates
 DEM <- raster("DEM_100x100.tif", proj4string = crs.reg)
 cos_aspect <- raster("cos_ASPECT.tif", proj4string = crs.reg)
 sin_aspect <- raster("Sin_ASPECT.tif", proj4string = crs.reg)
 Slope <- raster("SLOPE_100x100.tif", proj4string = crs.reg)
 Stream_dist4 <- raster("distance_matrix_Rhone_Strahler_4.tif", proj4string = crs.reg)
 Stream_dist6 <- raster("distance_matrix_Rhone_Strahler_6.tif", proj4string = crs.reg)
 Flint_dist <- raster("distance_matrix_flint.tif", proj4string = crs.reg)
 Coast_dist <- raster("distance_matrix_coast_line.tif", proj4string = crs.reg)
 Inland_dist <- raster("distance_matrix_inland.tif", proj4string = crs.reg)
 
 data.r5 <- data.frame(x = values(Krig_res), y = values(DEM), y1 = values(cos_aspect), y2 = values(sin_aspect), y3 = values(Slope), y4 = values(Stream_dist4), y5 = values(Stream_dist6), y6 = values(Coast_dist), y7 = values(Flint_dist))
 data.r5.na <- na.omit(data.r5)
 lin_reg5 <- lm(x ~ y + y1 + y2 + y3 + y4 + y5 + y6 + y7, data.r5.na)
 stepAIC(lin_reg5, k = 2)
 
 #focal correlation
 ##river obstacle
 Krig_Cost_obs <- stack("Cost_map/statistic_correlation_R/Krig_Cost_obs_clipped.tif") #import 2-band raster with Kriging and cost map obs
 
 KC_pos<- raster(Krig_Cost_obs,1) #grid defining positions
 values(KC_pos) <- 1:ncell(Krig_Cost_obs)
 
 matrix <- values(Krig_Cost_obs) # stack as raster
 focal_cor_obs <- focal(
   x = KC_pos,
   w = matrix(1, 11, 11), #defining a window
   fun = function(x, y = matrix){ 
     cor(y[x, 1], y[x, 2], 
         use = "na.or.complete")
   },
   filename = file.path("Cost_map/statistic_correlation_R", "focal_cor_obs.tif"),
   overwrite = TRUE
 )
 plot(focal_cor_obs, main = "FOCAL CORRELATION OBS")
 
 
 #SPEED RATES 
 speed_rate_obs <- terrain(Ccost_obs, opt = 'slope', unit = 'degrees', neighbors = 8,
                           filename = file.path("Cost_map/statistic_correlation_R", "speed_rate_obs.tif")) 
 plot(speed_rate_obs)
 
 
 speed_rate_aid <- terrain(Ccost_aid, opt = 'slope', unit = 'degrees', neighbors = 8,
                           filename = file.path("Cost_map/statistic_correlation_R", "speed_rate_aid.tif")) 
 plot(speed_rate_aid)
 
 
 speed_rate_krig <- terrain(Krig_res, opt = 'slope', unit = 'degrees', neighbors = 8,
                            filename = file.path("Cost_map/statistic_correlation_R", "speed_rate_krig.tif"))
 plot(speed_rate_krig)
 
 sl_KRsubOB <- speed_rate_krig - speed_rate_obs
 plot(sl_KRsubOB)
 
 sl_KRsubAI <- speed_rate_krig - speed_rate_aid
 plot(sl_KRsubAI)
 