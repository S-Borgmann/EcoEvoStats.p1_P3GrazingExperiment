#############################Grazing rate simulation#############################
library(dplyr)

## Set the random seed ##
set.seed(43)

## Parameter
# -> Could also include a linear regression model for the mesocosm and only use parameter for normal distirbution of species, would then leed to 2 normal dis. for species and 
# 4 uniform dis. for beta parameter of regression for intra and inter Meso per species
u1= 1.4 # mean of Normal Distribution of Dm, intraspecific Meso
u2= 1.4 # mean of Normal Distribution of Dm, interspecific Meso
u3= 1.2 # mean of Normal Distribution of Dp, intraspecific Meso
u4= 1.2 # mean of Normal Distribution of Dp, interspecific Meso

o1= 0.15 # variance of Normal Distribution of Dm, intraspecific Meso
o2= 0.15 # variance of Normal Distribution of Dm, interspecific Meso
o3= 0.15 # variance of Normal Distribution of Dp, intraspecific Meso
o4= 0.15 # variance of Normal Distribution of Dp, interspecific Meso

o5= 0.1 # interclonal variance

o6= 0.05 # intraclonal variance

b2=runif(1, 0, 0.1)
b3=runif(1, -0.1, 0.001)
b4=runif(1, -0.1, 0.001)
b5=runif(1, -0.1, 0.001)
b6=runif(1, -0.1, 0.001)
b7=runif(1, -0.1, 0.001)
b8=runif(1, -0.1, 0.001)

##Data
#fulldata= read.csv("/Users/sebastianborgmann/Nextcloud/Stuff/Projects/EcoEvoStats.p1/Project3_Grazing_experiment/grazing_analysis_example/HierarchyDataTest.csv")
fulldata= read.csv("C:/Users/Test/Nextcloud/Projects/EcoEvoStats.p1/Project3_Grazing_experiment/grazing_analysis_example/HierarchyDataTest.csv")
clonedata= distinct(fulldata, clone, .keep_all= TRUE)

##Hyperparameter for now

a1= rnorm(1,u1,o1) # Normal Distribution of Dm, intraspecific Meso
a2= rnorm(1,u2,o2) # Normal Distribution of Dm, interspecific Meso
a3= rnorm(1,u3,o3) # Normal Distribution of Dp, intraspecific Meso
a4= rnorm(1,u4,o4) # Normal Distribution of Dp, interspecific Meso

#interclonal variation in mesocosms
for (i in 1:nrow(clonedata)){
  if (clonedata[i,2] == 0 & clonedata[i,3] == 0 ){
    a5= rnorm(1,a1,o5)
  }
  else if (clonedata[i,2] == 0 & clonedata[i,3] == 1 ){
    a5= rnorm(1,a2,o5) 
  }
  else if (clonedata[i,2] == 1 & clonedata[i,3] == 0 ){ 
    a5= rnorm(1,a3,o5)
  }
  else if (clonedata[i,2] == 1 & clonedata[i,3] == 1 ){
    a5= rnorm(1,a4,o5)
  }
  clonedata[i,14] = a5
}

##Treatment model
fulldata <- merge(fulldata,clonedata,by = 'clone',all.x=T)
fulldata <- fulldata[ -c(14:25) ]

for (i in 1:nrow(fulldata)){
  fulldata[i,13]= fulldata[i,14] + b2*fulldata[i,5] + b3*fulldata[i,6] + b4*fulldata[i,7] + b5*fulldata[i,8] + b6*fulldata[i,9] + b7*fulldata[i,10] + b8*fulldata[i,11]
}

## intraclonal variation in treatment repetitions
for (i in 1:nrow(fulldata)){
  fulldata[i,15:17]= rnorm(3,fulldata[i,13],o6)
}

gr= fulldata[-c (4:14)]

##plotting
par(mfrow=c(2,2))
# Plot 1: dm intra
boxplot(fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_a.x=="1"],fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_b.x=="1"],fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_c.x=="1"],fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_d.x=="1"],fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_e.x=="1"],fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_f.x=="1"],fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_g.x=="1"],fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_h.x=="1"],ylab="grazing rate (micrograms C)",names=c("a","b","c","d","e","f","g","h"), main="D. magna, intraspecific",ylim = c(0.5, 2))
# Plot 2: dm inter
boxplot(fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_a.x=="1"],fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_b.x=="1"],fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_c.x=="1"],fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_d.x=="1"],fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_e.x=="1"],fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_f.x=="1"],fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_g.x=="1"],fulldata$V15[fulldata$species.x=="0" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_h.x=="1"],ylab="grazing rate (micrograms C)",names=c("a","b","c","d","e","f","g","h"), main="D. magna, interspecific",ylim = c(0.5, 2))
# Plot 3: dp intra
boxplot(fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_a.x=="1"],fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_b.x=="1"],fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_c.x=="1"],fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_d.x=="1"],fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_e.x=="1"],fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_f.x=="1"],fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_g.x=="1"],fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="0" & fulldata$graze_trt_h.x=="1"],ylab="grazing rate (micrograms C)",names=c("a","b","c","d","e","f","g","h"), main="D. pulex, intraspecific",ylim = c(0.5, 2))
# Plot 4: dp inter
boxplot(fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_a.x=="1"],fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_b.x=="1"],fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_c.x=="1"],fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_d.x=="1"],fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_e.x=="1"],fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_f.x=="1"],fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_g.x=="1"],fulldata$V15[fulldata$species.x=="1" & fulldata$meso_spec.x=="1" & fulldata$graze_trt_h.x=="1"],ylab="grazing rate (micrograms C)",names=c("a","b","c","d","e","f","g","h"), main="D. pulex, interspecific",ylim = c(0.5, 2))

