rm(list = ls())

list.of.packages <- c("ggplot2", "dplyr","tidyr","stringr","parallel","doParallel","data.table","here")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

require(ggplot2)
library(dplyr)
library(tidyr)
require(stringr)
require(parallel)
require(doParallel)
require(data.table)
require(here)

# set your working directory to the folder 'Output'
setwd("Output")

parallel::detectCores()

# define number of course you want to use, make sure to leave at least two cores for background processes
numCores <- parallel::detectCores() - 2

mycluster <- parallel::makeCluster(
  numCores, 
  type = "PSOCK",
  outfile=""
)

doParallel::registerDoParallel(cl = mycluster)

# list all the output files from the Movement Tool
# double check if the list contains the files you want to analyse
fls <- list.files(getwd(),"csv")

# choose group you want analyse (VHF- very high freqency, HF - high frequency, LF - low frequency, PCW - phocids in water)
groups <- c("LF", "VHF", "HF", "PCW")

# define thresholds for transition (as a step function) from impulsive to non-impulsive [m]
# you have to choose all 5 thresholds. If you want to test only one
# you can just define the same threshold for the 5 values

tr1 <- 1*1000
tr2 <- 2*1000
tr3 <- 3*1000
tr4 <- 7*1000
tr5 <- 10*1000 

# Cum_SEL growth rates depending on impulsiveness [dB], same for all hearing groups
tsI <- 2.3 #for impulsive
tsNI <- 1.6 # for non-impulsive

#SEL corrections [dB] for all hearing groups based on weighing by Southall
SELweighVHF <- 39.9 
SELweighHF <- 35.8 
SELweighLF <- 5.6 
SELweighPCW <- 18.8 

# TTS thresholds for each group
ttsVHF <- 140 
ttsHF <- 170 
ttsLF <- 168 
ttsPCW <- 170 

# PTS thresholds is the same for all groups
pts <- 40 # see description way below

##################################################################
# this is the main loop:
# don't change any code below this line 


st <- Sys.time()
for (k in 1:length(fls)){ 
ff <- fls[k]
sm1 <- fread(ff) # fread is much faster than read.csv

cn <- c("time", #seconds from the start
        "anim_id",
        "dist", #distance [m] from the source
        "exper_SL", # sound level experienced by an animal in a given time step
        "blast") # is there a blast in a given time step

colnames(sm1) <- cn

# extracting settings of your simulations from the file name

ee <- str_split_1(ff,"_|.c")
sp <- as.numeric(ee[6])
br <- as.numeric(str_split_1(ee[4],"seq")[2])
sl <- as.numeric(ee[10])
tl <- as.numeric(ee[12])

# setting exper_SL to zero when there is no blow
sm1$exper_SL[sm1$blast==0] <- 0

# because we assume that in between blows animals's accumulation of sound = 0
# we can delete all zeros
sm1 <- sm1[sm1$blast==1,]

# 1. splitting file for each animal

sm1PerAnim <- split(sm1, sm1$anim_id)
dists <- list()

for (gr in 1:length(groups)){
group <- groups[gr]

# groups defines Southall weighting and TTS threshold
# SEL correction [dB] for groups based on weighing by Southall
if (group == "VHF") SELweigh <- SELweighVHF
if (group == "HF") SELweigh <- SELweighHF
if (group == "LF") SELweigh <- SELweighLF 
if (group == "PCW") SELweigh <- SELweighPCW 

# TTS thresholds for each group
if (group == "VHF") tts <- ttsVHF 
if (group == "HF") tts <- ttsHF 
if (group == "LF") tts <- ttsLF 
if (group == "PCW") tts <- ttsPCW 

# calculate cum SEL for each animal 

# 2. calculating cumulative SEL (weighted) per animal per time step
#    calculating distance of each animal at the beginning of simulation

Cum_SELs_weighted <- foreach (j = 1:length(sm1PerAnim)) %dopar% { #length(sm1PerAnim)
  temp <- sm1PerAnim[[j]]$exper_SL - SELweigh 
  temp[temp<=0] <- 0.00001 #it cannot be zero
  Cum_SEL_weighted <- c(rep(NA,length(temp)))
  for (i in 1:length(temp)){
    Cum_SEL_weighted[i] <- 10*log10(sum(10^(temp[1:i]/10)))
  }
  
  Cum_SEL_weighted
}

for (j in 1:length(sm1PerAnim)){
sm1PerAnim[[j]]$Cum_SEL_weighted <- Cum_SELs_weighted[[j]]
sm1PerAnim[[j]]$start_dist <- sm1PerAnim[[j]]$dist[1]
}

sm <- do.call("rbind",sm1PerAnim)

# recalculating category of impulsiveness at each step based on threshold
sm$imp_tr1 <- "impulsive"
sm$imp_tr2 <- "impulsive"
sm$imp_tr3 <- "impulsive"
sm$imp_tr4 <- "impulsive"
sm$imp_tr5 <- "impulsive"

sm$imp_tr1[sm$dist>=tr1] <- "non-impulsive"
sm$imp_tr2[sm$dist>=tr2] <- "non-impulsive"
sm$imp_tr3[sm$dist>=tr3] <- "non-impulsive"
sm$imp_tr4[sm$dist>=tr4] <- "non-impulsive"
sm$imp_tr5[sm$dist>=tr5] <- "non-impulsive"


sm$imp_tr1 <- factor(sm$imp_tr1, levels=c("impulsive","non-impulsive"))
sm$imp_tr2 <- factor(sm$imp_tr2, levels=c("impulsive","non-impulsive"))
sm$imp_tr3 <- factor(sm$imp_tr3, levels=c("impulsive","non-impulsive"))
sm$imp_tr4 <- factor(sm$imp_tr4, levels=c("impulsive","non-impulsive"))
sm$imp_tr5 <- factor(sm$imp_tr5, levels=c("impulsive","non-impulsive"))


# calculating hearing thresholds
# we first calculate when animals (at which distance) experience TTS (using tts) based on Cum_sel_weighted
# and this is irrespective whether animals are in impulsive or non-impulsive zone
# once they reach TTS, how calculation from getting TTS to PTS differs depending whether animals are in the impulsive or non-impulsive zone
# and how it differs depends on this two growth rates (tsI, tsNI)

# step 1 - calculating distance when animals experienced TTS

sm$TTS <- "N"
sm$TTS[sm$Cum_SEL_weighted>=tts] <- "Y"

perId <- split(sm,sm$anim_id)
tts_perID <- list()
for (i in 1:length(perId)){
  temp <- perId[[i]][perId[[i]]$TTS=="Y",]
  tts_perID[[i]] <- temp[1,] # we only want the distance when TTS happened.
}
tts_perID <- do.call("rbind",tts_perID)
# the NAs are these animals which never experienced TTS

tts_perID <- tts_perID[!is.na(tts_perID$anim_id),]

#save(tts_perID, file=paste("tts_perID_",group,"_BlowRate_",ee[7],"_",'_Speed_',ee[5],'_sl_',ee[9],'_TL_',ee[11],".Rdata",sep=""))

# for this given speed and given blow rate
# all animals which started at:
# max(tts_perID$start_dist) experienced TTS within the model duration


# step 2: calculating PTS
# we first calculate difference between Cum_sel between each steps for which tts where reached
# and multiply by the two growth rates (tsI, tsNI) depending if this is impulsive or not
# this step is only calculated for the animals which experiences TTS for a given simulation
# so we first filter for this animals which did experience tts
anim_tts <- tts_perID$anim_id
perId <- do.call("rbind",perId)
pts_perID <- perId[perId$anim_id %in% c(anim_tts),]
pts_perID <- split(pts_perID,pts_perID$anim_id)

# there may be simulations when none of the simulated animat received tts, in such a case below loop
# wont apply
if (length(pts_perID) != 0)
{
  for (i in 1:length(pts_perID)){
    pts_perID[[i]] <- pts_perID[[i]][pts_perID[[i]]$TTS=="Y",]
    
    # in cases where animals received TTS at the very last step of the model, you cant really calculate delta
    # in the previous step
    if (nrow(pts_perID[[i]])==1){
      pts_perID[[i]]$deltaCumSEL <- pts_perID[[i]]$Cum_SEL_weighted - tts
    }
    else {
      pts_perID[[i]]$deltaCumSEL <- c(pts_perID[[i]]$Cum_SEL_weighted[1]-tts,pts_perID[[i]]$Cum_SEL_weighted[2:nrow(pts_perID[[i]])]-pts_perID[[i]]$Cum_SEL_weighted[1:(nrow(pts_perID[[i]])-1)])
      
    }
    pts_perID[[i]]$SEL_pts_tr1 <- NA
    pts_perID[[i]]$SEL_pts_tr1[pts_perID[[i]]$imp_tr1=="impulsive"] <- pts_perID[[i]]$deltaCumSEL[pts_perID[[i]]$imp_tr1=="impulsive"] * tsI
    pts_perID[[i]]$SEL_pts_tr1[pts_perID[[i]]$imp_tr1=="non-impulsive"] <- pts_perID[[i]]$deltaCumSEL[pts_perID[[i]]$imp_tr1=="non-impulsive"] * tsNI
    
    pts_perID[[i]]$SEL_pts_tr2 <- NA
    pts_perID[[i]]$SEL_pts_tr2[pts_perID[[i]]$imp_tr2=="impulsive"] <- pts_perID[[i]]$deltaCumSEL[pts_perID[[i]]$imp_tr2=="impulsive"] * tsI
    pts_perID[[i]]$SEL_pts_tr2[pts_perID[[i]]$imp_tr2=="non-impulsive"] <- pts_perID[[i]]$deltaCumSEL[pts_perID[[i]]$imp_tr2=="non-impulsive"] * tsNI
    
    pts_perID[[i]]$SEL_pts_tr3 <- NA
    pts_perID[[i]]$SEL_pts_tr3[pts_perID[[i]]$imp_tr3=="impulsive"] <- pts_perID[[i]]$deltaCumSEL[pts_perID[[i]]$imp_tr3=="impulsive"] * tsI
    pts_perID[[i]]$SEL_pts_tr3[pts_perID[[i]]$imp_tr3=="non-impulsive"] <- pts_perID[[i]]$deltaCumSEL[pts_perID[[i]]$imp_tr3=="non-impulsive"] * tsNI
    
  }
  
  
  # now we calculate cumulative SEL_pts_trx. It is linear so we just sum them up
  
  Cum_SEL_pts_tr1_res <- foreach (j = 1:length(pts_perID)) %dopar% { 
    Cum_SEL_pts_tr1 <- c(rep(NA,nrow(pts_perID[[j]])))
    for (i in 1:nrow(pts_perID[[j]])){
      Cum_SEL_pts_tr1[i] <- sum(pts_perID[[j]]$SEL_pts_tr1[1:i], na.rm=T)
    }
    Cum_SEL_pts_tr1
  }
  
  
  Cum_SEL_pts_tr2_res <- foreach (j = 1:length(pts_perID)) %dopar% { 
    Cum_SEL_pts_tr2 <- c(rep(NA,nrow(pts_perID[[j]])))
    for (i in 1:nrow(pts_perID[[j]])){
      Cum_SEL_pts_tr2[i] <- sum(pts_perID[[j]]$SEL_pts_tr2[1:i])
    }
    Cum_SEL_pts_tr2
  }
  
  
  Cum_SEL_pts_tr3_res <- foreach (j = 1:length(pts_perID)) %dopar% { 
    Cum_SEL_pts_tr3 <- c(rep(NA,nrow(pts_perID[[j]])))
    for (i in 1:nrow(pts_perID[[j]])){
      Cum_SEL_pts_tr3[i] <- sum(pts_perID[[j]]$SEL_pts_tr3[1:i])
    }
    Cum_SEL_pts_tr3
  }
  
  
  # merging all results for pts and finding distance at which animals experiences pts
  
  for (i in 1:length(pts_perID)){
    pts_perID[[i]]$Cum_SEL_pts_tr1 <- Cum_SEL_pts_tr1_res[[i]]
    pts_perID[[i]]$Cum_SEL_pts_tr2 <- Cum_SEL_pts_tr2_res[[i]]
    pts_perID[[i]]$Cum_SEL_pts_tr3 <- Cum_SEL_pts_tr3_res[[i]]
    
    pts_perID[[i]]$pts_tr1 <- "N"
    pts_perID[[i]]$pts_tr1[pts_perID[[i]]$Cum_SEL_pts_tr1>=pts] <- "Y"
    
    pts_perID[[i]]$pts_tr2 <- "N"
    pts_perID[[i]]$pts_tr2[pts_perID[[i]]$Cum_SEL_pts_tr2>=pts] <- "Y"
    
    pts_perID[[i]]$pts_tr3 <- "N"
    pts_perID[[i]]$pts_tr3[pts_perID[[i]]$Cum_SEL_pts_tr3>=pts] <- "Y"
    
    
  }
  
  
  #save(pts_perID, file=paste("pts_perID_",group,"_BlowRate_",ee[7],"_",'_Speed_',ee[5],'_sl_',ee[9],'_TL_',ee[11],".Rdata",sep=""))
  
  
  # finally finding the max start distance at which animals experiences pts
  # note that for very slow animals the threshold of impulsiveness may not matter as they never reach the threshold distance
  
  
  pts_max_dist <- list()
  for (i in 1:length(pts_perID)){
    temp <- pts_perID[[i]]$start_dist[pts_perID[[i]]$pts_tr1=="Y"]
    ifelse(length(temp)==0,dist_tr1 <- NA, dist_tr1<- max(pts_perID[[i]]$start_dist[pts_perID[[i]]$pts_tr1=="Y"]))
    
    temp <- pts_perID[[i]]$start_dist[pts_perID[[i]]$pts_tr2=="Y"]
    ifelse(length(temp)==0,dist_tr2 <- NA, dist_tr2<- max(pts_perID[[i]]$start_dist[pts_perID[[i]]$pts_tr2=="Y"]))
    
    temp <- pts_perID[[i]]$start_dist[pts_perID[[i]]$pts_tr3=="Y"]
    ifelse(length(temp)==0,dist_tr3 <- NA, dist_tr3<- max(pts_perID[[i]]$start_dist[pts_perID[[i]]$pts_tr3=="Y"]))
    
    pts_max_dist[[i]] <- c(dist_tr1,dist_tr2,dist_tr3)
    
  }
  pts_max_dist <- do.call("rbind",pts_max_dist)
  pts_max_dist <- as.data.frame(pts_max_dist)
  colnames(pts_max_dist) <- c("pts_tr1","pts_tr2","pts_tr3")
  
  max_tr1 <- max(pts_max_dist$pts_tr1, na.rm=T)
  max_tr2 <- max(pts_max_dist$pts_tr2, na.rm=T)
  max_tr3 <- max(pts_max_dist$pts_tr3, na.rm=T)
  
  h_group <- groups[gr]
  dists[[gr]] <- as.data.frame(cbind(h_group,max_tr1,max_tr2,max_tr3,br,sp, sl, tl))
  
}

else
{
  h_group <- groups[gr]
  dists[[gr]] <- as.data.frame(cbind(h_group,NA,NA,NA,br,sp, sl, tl))
}

#save(pts_max_dist, file=paste("pts_max_dist_",group,"_BlowRate_",ee[7],"_",'_Speed_',ee[5],'_sl_',ee[9],'_TL_',ee[11],".Rdata",sep=""))
print(paste("Group:", groups[gr], "of file:", k))

} #  this is the end of gr loop
distss <- do.call("rbind",dists)
write.csv(distss,file=paste("pts_max_dist_BlowRate_",ee[7],"_",'_Speed_',ee[5],'_sl_',ee[9],'_TL_',ee[11],".csv",sep=""), row.names=FALSE)

} # this is the end of k loop
parallel::stopCluster(cl = mycluster)

en <- Sys.time()
en-st


