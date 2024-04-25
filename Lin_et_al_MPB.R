#' ---
#' title: "Environmental influence on coastal benthic communities in a complex setting"
#' author: "Yuting Vicky Lin"
#' date: "2024-04-23"
#' output: 
#'   html_document:
#'     toc_float: true
#'     toc: true
#'     toc_depth: 3
#'     number_sections: true
#' ---
#' 
#' Data and script to replicate analyses in Lin et al. (XXX). If used in full or in part, please cite the original publication: 
#' 
#' Lin YV, Château P-A, Nozawa Y, Wei C-L, Wunderlich RF, Denis V (Under Revision in Marine Pollution Bulletin) Environmental influence on coastal benthic communities in a complex setting. DOI: XXX
#' 
#' Raw data also available on [Dryad](XXX)
#' 
#' A study from [FRElab](https://www.dipintothereef.com/)
#' 
#' 
#' # **Generate R Script**
#' 

#' 
#' # **Packages**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
library(vegan)
library(NbClust)
library(bootcluster)
library(indicspecies)
library(ape)
library(ggplot2)
library(partykit)
library(party)
library(rgdal)
library(ggrepel)
library(scatterpie)

#' 
#' # **Dataset**
#' 
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
rm(list=ls())
# setwd('') # setup the working directory
benthic.fct <- read.csv('Data/Linetal_dataset_Benthic.csv',header = T) # import the benthic matrix
rownames(benthic.fct) <- benthic.fct[,1] # assign transects' information to the row names of the benthic matrix
benthic.fct <- benthic.fct[,-1] # delete the first column, which contains transects' information
benthic.fct.hell <- as.data.frame(decostand(benthic.fct, "hell")) # hellinger transformation on the benthic matrix
env.info <- read.csv('Data/Linetal_dataset_Env.csv',header = T )# import the environmental driver matrix
coord.df <- read.csv('Data/Linetal_dataset_Cor.csv' , head = T , sep = ",")

#' 
#' 
#' # **Benthic communities**
#' ## **Partition of benthic communities**
#' 
#' The outcomes from the k-means clustering, stability of partitions, conditional random forest models, and the variable importance from models are slightly different every time when running, so we pre-save the outcomes of the current study and used it here to make sure the results showing here are consistent with the study.
#' 
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
#k.means <- NbClust(benthic.fct.hell, diss = NULL, method = "kmeans" ,min.nc=2, max.nc=15, index = "all") # k-means - Appendix C

k.means <- readRDS("Data/k.means.RData") # read the outcome of the current study  

#' 
#' 
#' ## **Stability of community delineation**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
#set.seed(1)
#stability1 <- stability(benthic.fct.hell, 5, B = 20, r = 5, scheme_2 = T) # estimate k-means partition stability with bootstrapping 

stability1 <- readRDS('Data/stability1.RData') # read the outcome of the current study  

cluster.stability1 <- data.frame(rownames(benthic.fct.hell), stability1$obs_wise)  # save the stability information
colnames(cluster.stability1) <- c('transect','stability') # rename the column name of the stability information matrix
cluster.stability1[,2] <- round(cluster.stability1[,2],2) # show the value of stability until two decimal places
stability.clu <- data.frame(as.factor(as.vector(stability1$membership)), as.vector(stability1$obs_wise)) # add the partition information
colnames(stability.clu) <- c('BC','stab') # add column names
rownames(stability.clu) <- rownames(benthic.fct.hell) # add row names
head(stability.clu) # Appendix D

#' 
#' ## **Significance of k-means**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
mrpp(benthic.fct.hell, stability.clu$BC, permutations = 9999)# permutation test to check the significance of k-means

#' 
#' ## **Moran I's test**
#'  
#' Moran's I test on community occurrence to see if the spatial autocorrelation exist for each community
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
dist <- as.matrix(dist(cbind(env.info$Longitude,env.info$Latitude)))
dist.inv <- 1/dist
diag(dist.inv) <- 0
dist.inv[is.infinite(dist.inv)] <- 20000

BC1 <- ifelse(stability.clu$BC == 1,stability.clu$stab*1,stability.clu$stab*0)
BC2 <- ifelse(stability.clu$BC == 2,stability.clu$stab*1,stability.clu$stab*0)
BC3 <- ifelse(stability.clu$BC == 3,stability.clu$stab*1,stability.clu$stab*0)
BC4 <- ifelse(stability.clu$BC == 4,stability.clu$stab*1,stability.clu$stab*0)
BC5 <- ifelse(stability.clu$BC == 5,stability.clu$stab*1,stability.clu$stab*0)

moran.BC1 <- Moran.I(BC1,dist.inv,na.rm=T)
moran.BC2 <- Moran.I(BC2,dist.inv,na.rm=T)
moran.BC3 <- Moran.I(BC3,dist.inv,na.rm=T)
moran.BC4 <- Moran.I(BC4,dist.inv,na.rm=T)
moran.BC5 <- Moran.I(BC5,dist.inv,na.rm=T)

moran.BC1 # Appendix E
moran.BC2 # Appendix E
moran.BC3 # Appendix E
moran.BC4 # Appendix E
moran.BC5 # Appendix E

#' 
#' ## **Indicator species analysis**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
ind <-  multipatt(benthic.fct.hell, stability.clu$BC, duleg = TRUE, control = how(nperm=999)) # find out indicator groups of each community
summary(ind) # summarize the result
ind.df <- round(ind$sign,3) # round to three decimal places
head(ind.df) # Appendix F

#' 
#' ## **Ordination - nMDS**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
set.seed(123)
nmds.fct <- metaMDS(benthic.fct.hell, distance = "euclidean", try = 10, #nMDS
                    maxit = 2000, sratmax = 1.2, sfgrmin = 1e-19,  k= 2)  
nmds.fct$stress # check the stress value of nMDS
stressplot(nmds.fct) # Appendix G

#' 
#' extract the scores of nMDS
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
site.scores.fct <- as.data.frame(scores(nmds.fct)$sites) # extract site scores
site.scores.fct$Region <- env.info$Region # add region information to site score matrix
site.scores.fct$Site <- env.info$Site # add site information to site score matrix
site.scores.fct$Depth <- env.info$Depth # add depth information to site score matrix
site.scores.fct$BC <- as.factor(stability.clu$BC) # add partition of communities information to site score matrix

species.scores.fct <- as.data.frame(scores(nmds.fct, "species")) # extract species scores
species.scores.fct$species <- rownames(species.scores.fct) # give the row names to species score matrix

# calculate the centroid of each identified BC
cent <- aggregate(cbind(NMDS1, NMDS2) ~ BC, data = site.scores.fct, FUN = mean)
segs <- merge(site.scores.fct, setNames(cent, c('BC','oNMDS1','oNMDS2')),
              by = 'BC', sort = FALSE)

#' 
#' 
#' # **Environmental drivers**
#' 
#' ## **Visually check the distribution of continuous environmental factors**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
hist(env.info$Mean_SST)
hist(env.info$SD_SST)
hist(env.info$Light_intensity)
hist(env.info$Wave_exposure) # skewed distribution
hist(env.info$Wave_height)
hist(env.info$Mean_chl_a)
hist(env.info$SD_chl_a)
hist(env.info$Nitrate)
hist(env.info$Nitrite)
hist(env.info$Phosphate)
hist(env.info$DHW)
hist(env.info$DHW_recovery)
hist(env.info$Typhoon_disturbance)
hist(env.info$Typhoon_recovery)
hist(env.info$Typhoon_frequency)
hist(env.info$Anthropogenic_land_use) # skewed distribution
hist(env.info$Forest_land_use) # skewed distribution
hist(env.info$Population_density) # skewed distribution
hist(env.info$Tourist_visitors) # skewed distribution
hist(env.info$Unstable_substrate_cover) # skewed distribution

#' 
#' ## **Log transformation for the factors with the skewed distribution**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
env.info$Anthropogenic_land_use_log <- log(env.info$Anthropogenic_land_use+1)
env.info$Forest_land_use_log <- log(env.info$Forest_land_use+1)
env.info$Tourist_visitors_log <- log(env.info$Tourist_visitors+1)
env.info$Population_density_log <- log(env.info$Population_density+1)
env.info$Wave_exposure_log <- log(env.info$Wave_exposure+1)
env.info$Unstable_substrate_cover_log <- log(env.info$Unstable_substrate_cover+1)

env.info$Anthropogenic_land_use <- NULL
env.info$Forest_land_use <- NULL
env.info$Tourist_visitors <- NULL
env.info$Population_density <- NULL
env.info$Wave_exposure <- NULL
env.info$Unstable_substrate_cover <- NULL

Management_status <- as.factor(env.info$Management_status) # set Management_status as factor
env.info$Management_status <- Management_status

#' 
#' 
#' 
#' # **Conditional random forest models**
#' ## **Best hyper-parameters to build models**
#' **Number of splits**
#' The ideal number of splits (mtry) is suggested as the square root of the number of explained factors. Accordingly, the ideal number of splits for this study is the square root of 21, which is between 4 to 5. Therefore, we decide mtry = 5
#' 
#' **Number of trees**
#' Functions to calculate the models with the smallest MSE by adjusting the number of trees (mtree)
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
cforest.best.mtree <- function(data,mtry, mtree_min, mtree_max){
  
  mse <- matrix(NA, ncol=1, nrow=mtree_max)
  
  for (i in mtree_min : mtree_max){
    c.rf <- cforest(data[,ncol(data)] ~ . , data = data, controls = cforest_unbiased(mtry = mtry, ntree = i)) # build the random forest model
    pred.c.rf <- predict(c.rf, OOB = TRUE) #get predicted (fitted) values
    c.rf.resid <- pred.c.rf - data[,ncol(data)] # the residuals = #cross-tabulate the predicted and observed values
    mse[i,] <- mean(c.rf.resid^2) # mean squared error
  }
  mtree.plot <- plot(mse, type = "b", col = "blue")
  mtree.min <- c(which.min(mse),mse[which.min(mse)])
  c.mtree <- list(mtree.plot,mtree.min)
  return(c.mtree)
  
}

#' 
#' ## **General conditional random forest (CRF)**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
env.tran.all <- env.info[,c(6:26)] # a matrix for only environmental factors
BC <- as.factor(stability.clu$BC) # set community type as a f actor
c.forest.benthic.env <- data.frame(env.tran.all, BC) # combine the matrix of environmental factors with the matrix of benthic type

# c.rf.mintree <- cforest.best.mtree(c.forest.benthic.env,5,50,3000) # MSEs become stable when mtree > 1000, so setup mtree = 2000
# don't run since it takes for a whilw

#c.rf.near.tran <- cforest(BC ~ Mean_SST + SD_SST + 
#                            Light_intensity + Wave_exposure_log + Wave_height + 
#                            Mean_chl_a + SD_chl_a + Unstable_substrate_cover_log + 
#                            Nitrate + Nitrite + Phosphate +
#                            DHW + DHW_recovery + 
#                            Typhoon_disturbance + Typhoon_recovery + Typhoon_frequency +
#                            Anthropogenic_land_use_log + Forest_land_use_log +
#                            Population_density_log + Tourist_visitors_log + Management_status , 
#                          data = c.forest.benthic.env, 
#                          controls = cforest_unbiased(mtry = 5, ntree = 2000)) # add management later on
c.rf.near.tran <- readRDS('Data/c.rf.near.tran.BC.RData') # read the outcome of the current study  

#' 
#' 
#' **Variable importance**
#' 
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
#c.importance <- varimp(c.rf.near.tran, conditional = T) # variable importance
c.importance <- readRDS('Data/c.importance.test.RData') # read the outcome of the current study  
imp.mean <-mean(c.importance)
par(mar = c(5,18.5, 4.1, 2.1))
barplot(main = 'General CRF', sort(c.importance)[10:21],horiz = T,las=1, xlim = c(0,0.02),cex.main=1, cex.lab=1.5,
        col='white', alpha = 0.5, xlab="Mean decrease in accuracy")
abline(v = imp.mean ,col="red4",lty=2) # Appendix I

#' 
#' **Model accuracy - Appendix H**
#' 
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
pred.c.rf.tran.oob <- predict(c.rf.near.tran, OOB = TRUE) #  the predict values of general CRF 
pred.sum <- sum(table(pred.c.rf.tran.oob, c.forest.benthic.env$BC)) # check if the prediction fits the observation
pred.err.OOB <- (8+7+3+13+22+5+15+22+1+4+1+8)/ pred.sum #  out-of-bag error of general CRF
pred.err.OOB

#' 
#' ## **Individual CRFs - Community 1** 
#' 
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
c.df.near.tran.BC1 <- data.frame(env.tran.all,BC1)

# don't run since it takes for a while
#c.rf.BC1.mintree <- cforest.best.mtree(c.df.near.tran.BC1,5,50,3000) # MSEs become stable when mtree > 1000, so setup mtree = 2000

#c.rf.near.tran.BC1 <- cforest(BC1 ~ Mean_SST + SD_SST + 
#                                Light_intensity + Wave_exposure_log + Wave_height + 
#                                Mean_chl_a + SD_chl_a + Unstable_substrate_cover_log + 
#                                Nitrate + Nitrite + Phosphate +
#                                DHW + DHW_recovery + 
#                                Typhoon_disturbance + Typhoon_recovery + Typhoon_frequency +
#                                Anthropogenic_land_use_log + Forest_land_use_log +
#                                Population_density_log + Tourist_visitors_log + Management_status ,
#                               data = c.df.near.tran.BC1,
#                               controls = cforest_unbiased(mtry = 5, ntree = 2000)) # build the random forest model
c.rf.near.tran.BC1 <- readRDS('Data/c.rf.near.tran.BC1.RData') # read the outcome of the current study  

#' 
#' **Variable importance**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
#c.importance.BC1 <- varimp(c.rf.near.tran.BC1, conditional = T)
c.importance.BC1 <- readRDS('Data/c.importance.BC1.RData') # read the outcome of the current study  
imp.mean.BC1 <- mean(c.importance.BC1)
par(mar = c(5.1,18.5, 4.1, 2.1))
barplot(main = 'Community 1', sort(c.importance.BC1)[10:21],horiz = T,las=1, xlim = c(0,0.025),
        cex.main=1,cex.lab=1.5,col='#FF66FF', alpha = 0.5, xlab="Mean decrease in accuracy")
abline(v = imp.mean.BC1 ,col="red4",lty=2)

#' 
#' **Model accuracy - Appendix H**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
pred.c.rf.tran.oob.BC1 <- predict(c.rf.near.tran.BC1, OOB = TRUE) #get predicted (fitted) values
resid.BC1 <- c.df.near.tran.BC1$BC1 -pred.c.rf.tran.oob.BC1
mse.BC1 <- mean(resid.BC1^2)
mse.BC1 # mean square error of the model
var.exp.BC1 <- 1 - sum((resid.BC1)^2)/sum((c.df.near.tran.BC1$BC1-mean(c.df.near.tran.BC1$BC1))^2)
var.exp.BC1 # variance explained of the model

#' 
#' **Partial dependence plots for each environmental factor - Appendix J**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
load("Data/c.rf.BC1.partial.RData") # read the outcome of the current study  
par(mfrow=c(6, 4), mar = c(4, 2, 1, 2))
plot(c.rf.BC1.light, type = 'l', xlab = 'Light intensity')
plot(c.rf.BC1.Typhoon_frequency, type = 'l', xlab = 'Typhoon frequency')
plot(c.rf.BC1.Wave_exposure, type = 'l', xlab = 'Wave expsure')
plot(c.rf.BC1.DHW_recovery, type = 'l', xlab = 'DHW recovery')
plot(c.rf.BC1.Management_status, type = 'l', xlab = 'Management status')
plot(c.rf.BC1.Forest_land_use_log, type = 'l', xlab = 'Forest land use')
plot(c.rf.BC1.Wave_height, type = 'l', xlab = 'Wave height')
plot(c.rf.BC1.Tourist_visitors_log, type = 'l', xlab = 'Tourist visitors')
plot(c.rf.BC1.Phosphate, type = 'l', xlab = 'Phosphate')
plot(c.rf.BC1.SD_SST, type = 'l', xlab = 'SD SST')
plot(c.rf.BC1.Nitrite, type = 'l', xlab = 'Nitrite')
plot(c.rf.BC1.SD_chl_a, type = 'l', xlab = 'SD Chl a')
plot(c.rf.BC1.Typhoon_recovery, type = 'l', xlab = 'Typhoon recovery')
plot(c.rf.BC1.Nitrate, type = 'l', xlab = 'Nitrate')
plot(c.rf.BC1.Mean_chl_a, type = 'l', xlab = 'Mean Chl a')
plot(c.rf.BC1.DHW, type = 'l', xlab = 'DHW')
plot(c.rf.BC1.Mean_SST, type = 'l', xlab = 'Mean SST')
plot(c.rf.BC1.Population_density_log, type = 'l', xlab = 'Population density')
plot(c.rf.BC1.Typhoon_disturbance, type = 'l', xlab = 'Typhoon disturbance')
plot(c.rf.BC1.Anthropogenic_land_use_log, type = 'l', xlab = 'Anthropogenic land use')
plot(c.rf.BC1.Unstable_substrate_cover, type = 'l', xlab = 'Unstable substrate cover')

#' 
#' 
#' ## **Individual CRFs - Community 2** 
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
c.df.near.tran.BC2 <- data.frame(env.tran.all,BC2)

# don't run since it takes for a whilw
#c.rf.BC2.mintree <- cforest.best.mtree(c.df.near.tran.BC2,5,20,3500) # MSEs become stable when mtree > 1000, so setup mtree = 2000

#c.rf.near.tran.BC2 <- cforest(BC2 ~ Mean_SST + SD_SST + 
#                                Light_intensity + Wave_exposure_log + Wave_height + 
#                                Mean_chl_a + SD_chl_a + Unstable_substrate_cover_log + 
#                                Nitrate + Nitrite + Phosphate +
#                                DHW + DHW_recovery + 
#                                Typhoon_disturbance + Typhoon_recovery + Typhoon_frequency +
#                                Anthropogenic_land_use_log + Forest_land_use_log +
#                                Population_density_log + Tourist_visitors_log + Management_status ,
#                               data = c.df.near.tran.BC2,
#                               controls = cforest_unbiased(mtry = 5, ntree = 2000)) # build the random forest model
c.rf.near.tran.BC2 <- readRDS('Data/c.rf.near.tran.BC2.RData') # read the outcome of the current study  

#' 
#' 
#' **Variable importance**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
#c.importance.BC2 <- varimp(c.rf.near.tran.BC2, conditional = T)
c.importance.BC2 <- readRDS('Data/c.importance.BC2.RData') # read the outcome of the current study  
imp.mean.BC2 <-mean(c.importance.BC2)
par(mar = c(5.1,18.5, 4.1, 2.1))
barplot(main = 'Community 2', sort(c.importance.BC2)[10:21],horiz = T,las=1, xlim = c(0,0.005),cex.main=1,
        col='#99CCFF', alpha = 0.5, xlab="Mean decrease in accuracy",cex.lab=1.5)
abline(v = imp.mean.BC2,col="red4",lty=2)

#' 
#' **Model accuracy - Appendix H**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
pred.c.rf.tran.oob.BC2 <- predict(c.rf.near.tran.BC2, OOB = TRUE) #get predicted (fitted) values
resid.BC2 <- pred.c.rf.tran.oob.BC2- c.df.near.tran.BC2$BC2
mse.BC2 <- mean(resid.BC2^2)
mse.BC2 # mean square error of the model
var.exp.BC2 <- 1 - sum((resid.BC2)^2)/sum((c.df.near.tran.BC2$BC2-mean(c.df.near.tran.BC2$BC2))^2)
var.exp.BC2 # variance explained of the model

#' 
#' **Partial dependence plots for each environmental factor - Appendix J**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
load("Data/c.rf.BC2.partial.RData") # read the outcome of the current study  
par(mfrow=c(6,4), mar = c(4,2, 1, 2))
plot(c.rf.BC2.Nitrite, type = 'l', xlab = 'Nitrite')
plot(c.rf.BC2.light, type = 'l', xlab = 'Light intensity')
plot(c.rf.BC2.SD_chl_a, type = 'l', xlab = 'SD chl a')
plot(c.rf.BC2.Phosphate, type = 'l', xlab = 'Phosphate')
plot(c.rf.BC2.Anthropogenic_land_use_log, type = 'l', xlab = 'Anthropogenic land use')
plot(c.rf.BC2.Tourist_visitors_log, type = 'l', xlab = 'Tourist visitors')
plot(c.rf.BC2.Mean_chl_a, type = 'l', xlab = 'Mean chl a')
plot(c.rf.BC2.Forest_land_use_log, type = 'l', xlab = 'Forest land use')
plot(c.rf.BC2.Nitrate, type = 'l', xlab = 'Nitrate')
plot(c.rf.BC2.Wave_exposure, type = 'l', xlab = 'Wave exposure')
plot(c.rf.BC2.Unstable_substrate_cover, type = 'l', xlab = 'Unstable substrate cover')
plot(c.rf.BC2.Typhoon_disturbance, type = 'l', xlab = 'Typhoon disturbance')
plot(c.rf.BC2.Population_density_log, type = 'l', xlab = 'Population density')
plot(c.rf.BC2.Mean_SST, type = 'l', xlab = 'Mean SST')
plot(c.rf.BC2.Wave_height, type = 'l', xlab = 'Wave height')
plot(c.rf.BC2.SD_SST, type = 'l', xlab = 'SD SST')
plot(c.rf.BC2.DHW, type = 'l', xlab = 'DHW')
plot(c.rf.BC2.Typhoon_frequency, type = 'l', xlab = 'Typhoon frequency')
plot(c.rf.BC2.DHW_recovery, type = 'l', xlab = 'DHW recovery')
plot(c.rf.BC2.Typhoon_recovery, type = 'l', xlab = 'Typhoon recovery')
plot(c.rf.BC2.Management_status, type = 'l', xlab = 'Management status')

#' 
#' ## **Individual CRFs - Community 3** 
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
c.df.near.tran.BC3 <- data.frame(env.tran.all,BC3)

# don't run since it takes for a while
#c.rf.BC3.mintree <- cforest.best.mtree(c.df.near.tran.BC3,5,20,3500) # MSEs become stable when mtree > 1000, so setup mtree = 2000

#c.rf.near.tran.BC3 <- cforest(BC3 ~ Mean_SST + SD_SST + 
#                                Light_intensity + Wave_exposure_log + Wave_height + 
##                                Mean_chl_a + SD_chl_a + Unstable_substrate_cover_log + 
##                                Nitrate + Nitrite + Phosphate +
#                                DHW + DHW_recovery + 
#                                Typhoon_disturbance + Typhoon_recovery + Typhoon_frequency +
#                                Anthropogenic_land_use_log + Forest_land_use_log +
#                                Population_density_log + Tourist_visitors_log + Management_status ,
#                              data = c.df.near.tran.BC3,
#                              controls = cforest_unbiased(mtry = 5, ntree = 2000)) # build the random forest model
c.rf.near.tran.BC3 <- readRDS('Data/c.rf.near.tran.BC3.RData') # read the outcome of the current study  

#' 
#' **Variable importance**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
#c.importance.BC3 <- varimp(c.rf.near.tran.BC3, conditional = T)
c.importance.BC3 <- readRDS('Data/c.importance.BC3.RData') # read the outcome of the current study  
imp.mean.BC3 <-mean(c.importance.BC3)
par(mar = c(5.1,18.5, 4.1, 2.1))
barplot(main = 'Community 3', sort(c.importance.BC3)[10:21],horiz = T,las=1, xlim = c(0,0.006),cex.main=1,
        col='grey', alpha = 0.5, xlab="Mean decrease in accuracy",cex.lab=1.5)
abline(v = imp.mean.BC3,col="red4",lty=2)

#' 
#' **Model accuracy - Appendix H**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
pred.c.rf.tran.oob.BC3 <- predict(c.rf.near.tran.BC3, OOB = TRUE) #get predicted (fitted) values
resid.BC3 <- pred.c.rf.tran.oob.BC3- c.df.near.tran.BC3$BC3
mse.BC3 <-mean(resid.BC3^2)
mse.BC3 # mean square error of the model
var.exp.BC3 <- 1 - sum((resid.BC3)^2)/sum((c.df.near.tran.BC3$BC3-mean(c.df.near.tran.BC3$BC3))^2)
var.exp.BC3 # variance explained of the model

#' 
#' **Partial dependence plots for each environmental factor - Appendix J**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
load("Data/c.rf.BC3.partial.RData") # read the outcome of the current study  
par(mfrow=c(6,4), mar = c(4,2, 1, 2))
plot(c.rf.BC3.Unstable_substrate_cover, type = 'l', xlab = 'Unstable substrate cover')
plot(c.rf.BC3.light, type = 'l', xlab = 'Light intensity')
plot(c.rf.BC3.Management_status, type = 'l', xlab = 'Management status')
plot(c.rf.BC3.Anthropogenic_land_use_log, type = 'l', xlab = 'Anthropogenic land use')
plot(c.rf.BC3.Tourist_visitors_log, type = 'l', xlab = 'Tourist visitors')
plot(c.rf.BC3.Wave_exposure, type = 'l', xlab = 'Wave exposure')
plot(c.rf.BC3.Nitrite, type = 'l', xlab = 'Nitrite')
plot(c.rf.BC3.Phosphate, type = 'l', xlab = 'Phosphate')
plot(c.rf.BC3.Population_density_log, type = 'l', xlab = 'Population density')
plot(c.rf.BC3.Nitrate, type = 'l', xlab = 'Nitrate')
plot(c.rf.BC3.Forest_land_use_log, type = 'l', xlab = 'Forest land use')
plot(c.rf.BC3.Wave_height, type = 'l', xlab = 'Wave height')
plot(c.rf.BC3.Typhoon_disturbance, type = 'l', xlab = 'Typhoon disturbance')
plot(c.rf.BC3.Mean_chl_a, type = 'l', xlab = 'Mean chl a')
plot(c.rf.BC3.Typhoon_frequency, type = 'l', xlab = 'Typhoon frequency')
plot(c.rf.BC3.SD_chl_a, type = 'l', xlab = 'SD chl a')
plot(c.rf.BC3.DHW, type = 'l', xlab = 'DHW')
plot(c.rf.BC3.Typhoon_recovery, type = 'l', xlab = 'Typhoon recovery')
plot(c.rf.BC3.Mean_SST, type = 'l', xlab = 'Mean SST')
plot(c.rf.BC3.SD_SST, type = 'l', xlab = 'SD SST')
plot(c.rf.BC3.DHW_recovery, type = 'l', xlab = 'DHW recovery')

#' 
#' ## **Individual CRFs - Community 4** 
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
c.df.near.tran.BC4 <- data.frame(env.tran.all,BC4)

# don't run since it takes for a while
#c.rf.BC4.mintree <- cforest.best.mtree(c.df.near.tran.BC2,5,20,3500) # MSEs become stable when mtree > 1000, so setup mtree = 2500

#c.rf.near.tran.BC4 <- cforest(BC4 ~ Mean_SST + SD_SST + 
#                                Light_intensity + Wave_exposure_log + Wave_height + 
#                                Mean_chl_a + SD_chl_a + Unstable_substrate_cover_log + 
#                                Nitrate + Nitrite + Phosphate +
#                                DHW + DHW_recovery + 
#                                Typhoon_disturbance + Typhoon_recovery + Typhoon_frequency +
#                                Anthropogenic_land_use_log + Forest_land_use_log +
#                                Population_density_log + Tourist_visitors_log + Management_status ,
#                              data = c.df.near.tran.BC4,
#                              controls = cforest_unbiased(mtry = 5, ntree = 2500)) # build the random forest model
c.rf.near.tran.BC4 <- readRDS('Data/c.rf.near.tran.BC4.RData') # read the outcome of the current study  

#' 
#' **Variable importance**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
#c.importance.BC4 <- varimp(c.rf.near.tran.BC4, conditional = T)
c.importance.BC4 <- readRDS('Data/c.importance.BC4.RData') # read the outcome of the current study  
imp.mean.BC4 <- mean(c.importance.BC4)
par(mar = c(5.1,18.5, 4.1, 2.1))
barplot(main = 'Community 4', sort(c.importance.BC4)[10:21],horiz = T,las=1, xlim = c(0,0.005),cex.main=1,
        col='#B8E186', alpha = 0.5, xlab="Mean decrease in accuracy",cex.lab=1.5)
abline(v = imp.mean.BC4 ,col="red4",lty=2)

#' 
#' **Model accuracy - Appendix H**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
pred.c.rf.tran.oob.BC4 <- predict(c.rf.near.tran.BC4, OOB = TRUE) #get predicted (fitted) values
resid.BC4 <- pred.c.rf.tran.oob.BC4- c.df.near.tran.BC4$BC4
mse.BC4 <- mean(resid.BC4^2)
mse.BC4 # mean square error of the model
var.exp.BC4 <- 1 - sum((resid.BC4)^2)/sum((c.df.near.tran.BC4$BC4-mean(c.df.near.tran.BC4$BC4))^2)
var.exp.BC4 # variance explained of the model

#' 
#' **Partial dependence plots for each environmental factor - Appendix J**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
load("Data/c.rf.BC4.partial.RData") # read the outcome of the current study  
par(mfrow=c(6,4), mar = c(4,2, 1, 2))
plot(c.rf.BC4.SD_chl_a, type = 'l', xlab = 'SD chl a')
plot(c.rf.BC4.Population_density_log, type = 'l', xlab = 'Population density')
plot(c.rf.BC4.DHW, type = 'l', xlab = 'DHW')
plot(c.rf.BC4.Nitrate, type = 'l', xlab = 'Nitrate')
plot(c.rf.BC4.Nitrite, type = 'l', xlab = 'Nitrite')
plot(c.rf.BC4.Wave_exposure, type = 'l', xlab = 'Wave exposure')
plot(c.rf.BC4.SD_SST, type = 'l', xlab = 'SD SST')
plot(c.rf.BC4.light, type = 'l', xlab = 'Light intensity')
plot(c.rf.BC4.Wave_height, type = 'l', xlab = 'Wave height')
plot(c.rf.BC4.Tourist_visitors_log, type = 'l', xlab = 'Tourist visitors')
plot(c.rf.BC4.Mean_chl_a, type = 'l', xlab = 'Mean Chl a')
plot(c.rf.BC4.Mean_SST, type = 'l', xlab = 'Mean SST')
plot(c.rf.BC4.Typhoon_disturbance, type = 'l', xlab = 'Typhoon disturbance')
plot(c.rf.BC4.Phosphate, type = 'l', xlab = 'Phosphate')
plot(c.rf.BC4.Forest_land_use_log, type = 'l', xlab = 'Forest land use')
plot(c.rf.BC4.Unstable_substrate_cover, type = 'l', xlab = 'Unstable substrate cover')
plot(c.rf.BC4.Management_status, type = 'l', xlab = 'Management status')
plot(c.rf.BC4.Anthropogenic_land_use_log, type = 'l', xlab = 'Anthropogenic land use')
plot(c.rf.BC4.DHW_recovery, type = 'l', xlab = 'DHW recovery')
plot(c.rf.BC4.Typhoon_recovery, type = 'l', xlab = 'Typhoon recovery')
plot(c.rf.BC4.Typhoon_frequency, type = 'l', xlab = 'Typhoon frequency')

#' 
#' ## **Individual CRFs - Community 5** 
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
c.df.near.tran.BC5 <- data.frame(env.tran.all,BC5)

# don't run since it takes for a while
#c.rf.BC5.mintree <- cforest.best.mtree(c.df.near.tran.BC5,5,20,3500)  # MSEs become stable when mtree > 1000, so setup mtree = 2500

#c.rf.near.tran.BC5 <- cforest(BC5 ~ Mean_SST + SD_SST + 
#                                Light_intensity + Wave_exposure_log + Wave_height + 
#                                Mean_chl_a + SD_chl_a + Unstable_substrate_cover_log + 
#                                Nitrate + Nitrite + Phosphate +
#                                DHW + DHW_recovery + 
#                                Typhoon_disturbance + Typhoon_recovery + Typhoon_frequency +
#                                Anthropogenic_land_use_log + Forest_land_use_log +
#                                Population_density_log + Tourist_visitors_log + Management_status ,
#                              data = c.df.near.tran.BC5,
#                              controls = cforest_unbiased(mtry = 5, ntree = 2500)) # build the random forest model
c.rf.near.tran.BC5 <- readRDS('Data/c.rf.near.tran.BC5.RData') # read the outcome of the current study  

#' 
#' **Variable importance**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
#c.importance.BC5 <- varimp(c.rf.near.tran.BC5, conditional = T)
c.importance.BC5 <- readRDS('Data/c.importance.BC5.RData') # read the outcome of the current study  
imp.mean.BC5 <- mean(c.importance.BC5)
par(mar = c(5.1,18.5, 4.1, 2.1))
barplot(main = 'Community 5', sort(c.importance.BC5)[10:21],horiz = T,las=1, xlim = c(0,0.02),cex.main=1,
        col='#FFED6F', alpha = 0.5, xlab="Mean decrease in accuracy",cex.lab=1.5)
abline(v = imp.mean.BC5 ,col="red4",lty=2)

#' 
#' **Model accuracy - Appendix H**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
pred.c.rf.tran.oob.BC5 <- predict(c.rf.near.tran.BC5, OOB = TRUE) #get predicted (fitted) values
resid.BC5 <- pred.c.rf.tran.oob.BC5- c.df.near.tran.BC5$BC5
mse.BC5 <- mean(resid.BC5^2)
mse.BC5 # mean square error of the model
var.exp.BC5 <- 1 - sum((resid.BC5)^2)/sum((c.df.near.tran.BC5$BC5-mean(c.df.near.tran.BC5$BC5))^2)
var.exp.BC5 # variance explained of the model

#' 
#' **Partial dependence plots for each environmental factor - Appendix J**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
load("Data/c.rf.BC5.partial.RData") # read the outcome of the current study  
par(mfrow=c(6,4), mar = c(4,2, 1, 2))
plot(c.rf.BC5.Typhoon_recovery, type = 'l', xlab = 'Typhoon recovery')
plot(c.rf.BC5.Typhoon_frequency, type = 'l', xlab = 'Typhoon frequency')
plot(c.rf.BC5.SD_chl_a, type = 'l', xlab = 'SD Chl a')
plot(c.rf.BC5.SD_SST, type = 'l', xlab = 'SD SST')
plot(c.rf.BC5.Mean_SST, type = 'l', xlab = 'Mean SST')
plot(c.rf.BC5.Mean_chl_a, type = 'l', xlab = 'Mean Chl a')
plot(c.rf.BC5.Forest_land_use_log, type = 'l', xlab = 'Forest land use')
plot(c.rf.BC5.Phosphate, type = 'l', xlab = 'Phosphate')
plot(c.rf.BC5.Wave_exposure, type = 'l', xlab = 'Wave exposure')
plot(c.rf.BC5.Anthropogenic_land_use_log, type = 'l', xlab = 'Anthropogenic land use')
plot(c.rf.BC5.Wave_height, type = 'l', xlab = 'Wave height')
plot(c.rf.BC5.Tourist_visitors_log, type = 'l', xlab = 'Tourist visitors')
plot(c.rf.BC5.Nitrite, type = 'l', xlab = 'Nitrite')
plot(c.rf.BC5.Management_status, type = 'l', xlab = 'Management status')
plot(c.rf.BC5.Population_density_log, type = 'l', xlab = 'Population density')
plot(c.rf.BC5.DHW, type = 'l', xlab = 'DHW')
plot(c.rf.BC5.Nitrate, type = 'l', xlab = 'Nitrate')
plot(c.rf.BC5.DHW_recovery, type = 'l', xlab = 'DHW recovery')
plot(c.rf.BC5.light, type = 'l', xlab = 'Light intensity')
plot(c.rf.BC5.Typhoon_disturbance, type = 'l', xlab = 'Typhoon disturbance')
plot(c.rf.BC5.Unstable_substrate_cover, type = 'l', xlab = 'Unstable substrate cover')

#' 
#' # **Figures** 
#' ## **Figure 1a**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
taiwan <- readOGR('Data/Linetal_dataset_TW.shp') # import map data
taiwan_map <- fortify(taiwan)  # convert .shp file into dataframe

ggplot() +
  geom_hline(yintercept = 23.5,linetype=5,col='black',alpha=.3)+
  geom_polygon(data= taiwan_map, aes(x= long, y= lat, group= group), color = "black", fill = 'antiquewhite') +
  geom_point(data= coord.df, aes(x = Longitude, y = Latitude), size = 1.5,col='red')+
  geom_text_repel(data= coord.df, aes(x = Longitude, y = Latitude,label=Label1),max.overlaps=Inf,
                  min.segment.length = 0)+
  scale_x_continuous(limits = c(116.2, 123)) +                                                          
  scale_y_continuous(limits = c(20.5, 25.5)) +
  annotate("text", x = 120.8, y = 23.6, label = "Taiwan", size=8) +
  theme_minimal() +                                             
  theme(legend.position="none", panel.background = element_rect(fill = 'aliceblue'))+
  coord_map(projection ="mollweide")

#' 
#' ## **Figure 1b**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
col <- c('#FF66FF','#99CCFF','grey','#B8E186','#FFED6F')
ggplot(coord.df, aes(x = 1, y= site)) +  
  geom_scatterpie(aes(x=x, y=y, r=.5), data=coord.df,cols=c('X1','X2','X3','X4','X5')) +
  geom_hline(yintercept = 40.5,linetype=5,col='black',alpha=.8)+
  scale_fill_manual(values = col)+
  scale_x_continuous("",expand = c(0,0),breaks = 1:5,labels = c("5","10","15","20","40"))+
  scale_y_continuous("",expand = c(0,0),breaks = 1:65,labels = unique(coord.df$Label2)[-2])+
  coord_fixed(ratio = 0.8)+ 
  labs(x = 'Depth (m)')

#' 
#' ## **Figure 2**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
ggplot() +
  geom_point(data = site.scores.fct, aes(x = NMDS1, y = NMDS2, color = BC),size=2) +
  annotate(size=5, geom = "label", x = -0.45, y = 0.5, 
           label = paste("Stress: ", round(nmds.fct$stress, digits = 2))) +
  geom_segment(data = segs,  mapping = aes(x = NMDS1, y = NMDS2,xend = oNMDS1, yend = oNMDS2, color = BC, alpha = 0.4)) + # spiders
  geom_point(data = cent, aes(x = NMDS1, y = NMDS2, color = BC), size = 5) + # centroid
  scale_color_manual(values = col, name = "BC")+
  geom_text(data = species.scores.fct, mapping = aes(x = NMDS1, y = NMDS2, label = species), alpha = 0.5, size = 5) +
  theme_minimal() +
  xlim (-0.5,0.4) +
  ylim (-0.5,0.5) +
  theme(legend.position = "right",text = element_text(size = 24))

#'  
#' ## **Figure 3**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
par(mfrow=c(3,2), mar = c(3,7, 3, 1))
barplot(main = 'General CRF', sort(c.importance)[11:21],horiz = T,las=1, xlim = c(0,0.02),cex.names=0.6,cex.main=1,col='white', alpha = 0.5)
abline(v = imp.mean ,col="red4",lty=2, xlab="Mean decrease in accuracy")
barplot(main = 'Community 1', sort(c.importance.BC1)[12:21],horiz = T,las=1, xlim = c(0,0.02),cex.names=0.7,cex.main=1,col='#FF66FF', alpha = 0.5)
abline(v = imp.mean.BC1 ,col="red4",lty=2, xlab="Mean decrease in accuracy")
barplot(main = 'Community 2', sort(c.importance.BC2)[12:21],horiz = T,las=1, xlim = c(0,0.005),cex.names=0.7,cex.main=1,col='#99CCFF', alpha = 0.5)
abline(v = imp.mean.BC2 ,col="red4",lty=2, xlab="Mean decrease in accuracy")
barplot(main = 'Community 3', sort(c.importance.BC3)[12:21],horiz = T,las=1, xlim = c(0,0.005),cex.names=0.7,cex.main=1,col='grey', alpha = 0.5)
abline(v = imp.mean.BC3 ,col="red4",lty=2, xlab="Mean decrease in accuracy")
barplot(main = 'Community 4', sort(c.importance.BC4)[12:21],horiz = T,las=1, xlim = c(0,0.005),cex.names=0.7,cex.main=1,col='#B8E186', alpha = 0.5)
abline(v = imp.mean.BC4 ,col="red4",lty=2, xlab="Mean decrease in accuracy")
barplot(main = 'Community 5', sort(c.importance.BC5)[12:21],horiz = T,las=1, xlim = c(0,0.02),cex.names=0.7,cex.main=1,col='#FFED6F', alpha = 0.5)
abline(v = imp.mean.BC5 ,col="red4",lty=2, xlab="Mean decrease in accuracy")

#' 
#' ## **Figure 4 - Partial dependence plots**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
par(mfrow=c(5,5), mar = c(4,2, 1, 2))
plot(c.rf.BC1.light, type = 'l')
plot(c.rf.BC1.Typhoon_frequency, type = 'l')
plot(c.rf.BC1.Wave_exposure, type = 'l')
plot(c.rf.BC1.DHW_recovery, type = 'l')
plot(c.rf.BC1.Management_status, type = 'l')

plot(c.rf.BC2.Nitrite, type = 'l')
plot(c.rf.BC2.light, type = 'l')
plot(c.rf.BC2.SD_chl_a, type = 'l')
plot(c.rf.BC2.Phosphate, type = 'l')
plot(c.rf.BC2.Anthropogenic_land_use_log, type = 'l')

plot(c.rf.BC3.Unstable_substrate_cover, type = 'l')
plot(c.rf.BC3.light, type = 'l')
plot(c.rf.BC3.Management_status, type = 'l')
plot(c.rf.BC3.Anthropogenic_land_use_log, type = 'l')
plot(c.rf.BC3.Tourist_visitors_log, type = 'l')

plot(c.rf.BC4.SD_chl_a, type = 'l')
plot(c.rf.BC4.Population_density_log, type = 'l')
plot(c.rf.BC4.DHW, type = 'l')
plot(c.rf.BC4.Nitrate, type = 'l')
plot(c.rf.BC4.Nitrite, type = 'l')

plot(c.rf.BC5.Typhoon_recovery, type = 'l')
plot(c.rf.BC5.Typhoon_frequency, type = 'l')
plot(c.rf.BC5.SD_chl_a, type = 'l')
plot(c.rf.BC5.SD_SST, type = 'l')
plot(c.rf.BC5.Mean_SST, type = 'l')

#' 
#' # **Table**
#' ## **Table 1**
## ----eval = TRUE, echo=T, message=F, warning=F--------------------------------
# major benthic composition of each identified benthic community
major <- c('ascidian','algae','algae','algae','cca','cyanobacterians',
           'hard_corals','hard_corals','hard_corals','hard_corals',
           'hard_corals','hard_corals','hard_corals','octocorals',
           'octocorals','octocorals','octocorals','octocorals',
           'octocorals','octocorals','seagrass','sponges','sponges',
           'sponges','turf','zoanthids','zoanthids') # create a vector with major benthic categories corresponding to each morpho-functional group
benthic.major <- data.frame(major,t(benthic.fct)) # combine the vector with benthic matrix
benthic.major <- aggregate(.~major,benthic.major,sum) # aggregate the covers of morpho-functional groups by major benthic categories
row.names(benthic.major) <- benthic.major[,1] # give row names to the major benthic composition matrix
benthic.major[,1] <- NULL # delete the row containing the major benthic composition information
benthic.major <- t(benthic.major) # transport the matrix 

major.mean <- list() # create an empty list for mean covers of each major benthic category
major.sd <- list() # create an empty list for covers of standard deviation of each major benthic category
major.per <- list() # create an empty list for relative covers of each major benthic category
BC.numeric <- as.numeric(stability.clu$BC)

for (i in 1:max(BC.numeric)) {
  major.per[[i]] <- benthic.major[BC.numeric==i,]/rowSums(benthic.fct[BC.numeric==i,]) # calculating the relative covers of each major benthic category
  major.mean[[i]] <- 100*round(apply(major.per[[i]],2, mean) ,3) # calculating the mean covers of each major benthic category
  major.sd[[i]] <- 100*round(apply(major.per[[i]],2, sd) ,3) # calculating covers of standard deviation of each major benthic category
}
major.mean.df <- matrix(unlist(major.mean),byrow = T,nrow=5) # mean covers of each major benthic category
colnames(major.mean.df) <-  names(major.mean[[1]]) # give the column name to the mean cover matrix
major.sd.df <- matrix(unlist(major.sd),byrow = T,nrow=5) # covers of standard deviation of each major benthic category
colnames(major.sd.df) <-  names(major.mean[[1]]) # give the column name to the SD cover matrix 
major.mean.df # check the mean cover matrix 
major.sd.df # check the SD cover matrix 

# morpho-functional composition of each identified benthic community
par.mean <- list() # create an empty list for mean covers of each morpho-funtional group
par.sd <- list() # create an empty list for covers of standard deviation of each morpho-funtional group
par.per <- list() # create an empty list for relative covers of each morpho-funtional group 

for (i in 1:max(BC.numeric)) {
  par.per[[i]] <- benthic.fct[BC.numeric==i,]/rowSums(benthic.fct[BC.numeric==i,])  # calculating the relative covers of each morpho-funtional group 
  par.mean[[i]] <- 100*round(apply(par.per[[i]],2, mean) ,3) # calculating the mean covers of each morpho-funtional group
  par.sd[[i]] <- 100*round(apply(par.per[[i]],2, sd) ,3) # calculating covers of standard deviation of each morpho-funtional group
}
par.mean.df <- matrix(unlist(par.mean),byrow = T,nrow=5) # mean covers of each morpho-funtional group
colnames(par.mean.df) <-  names(par.mean[[1]]) # give the column name to the mean cover matrix
par.sd.df <- matrix(unlist(par.sd),byrow = T,nrow=5) # covers of standard deviation of each morpho-funtional group
colnames(par.sd.df) <-  names(par.mean[[1]]) # give the column name to the SD cover matrix 
par.mean.df # check the mean cover matrix 
par.sd.df # check the SD cover matrix 

