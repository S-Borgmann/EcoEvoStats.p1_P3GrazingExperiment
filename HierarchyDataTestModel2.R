#############################Grazing rate simulation#############################
library(dplyr)

## Set the random seed ##
set.seed(43)

#Create Dataframe: 35 clones = 8 Dm (4 intra + 4 inter) + 27 Dp (13 intra + 14 inter) and 8 reps

fulldata=data.frame(clone = c(rep(1:35, each=8)), Species = c(rep(0,280)), MesoDm= c(rep(0,280)), MesoDp= c(rep(0,280))
                , Trta=rep(c(1,0,0,0,0,0,0,0),times=35), Trtb=rep(c(0,1,0,0,0,0,0,0),times=35)
                , Trtc=rep(c(0,0,1,0,0,0,0,0),times=35), Trtd=rep(c(0,0,0,1,0,0,0,0),times=35)
                , Trte=rep(c(0,0,0,0,1,0,0,0),times=35), Trtf=rep(c(0,0,0,0,0,1,0,0),times=35)
                , Trtg=rep(c(0,0,0,0,0,0,1,0),times=35), Trth=rep(c(0,0,0,0,0,0,0,1),times=35))

fulldata[65:280,2]=1
fulldata[33:64,3]=1
fulldata[169:280,4]=1
clonedata= distinct(fulldata, clone, .keep_all= TRUE)

##Parameter

#Species
a1= 1.4 #Mean Dm
a2= 1.4 #Mean Dp

o1= 0.2 #var Dm (? Add an exponetial prior for this and other variance?)
o2= 0.2 #var Dp

#Mesocosm
b1= runif(1, -0.1, 0.05) # influence of Mesocoms on Dm intraspecific
b2= runif(1, -0.1, 0.05) # influence of Mesocoms on Dm interspecific
b3= runif(1, -0.1, 0.05) # influence of Mesocoms on Dp intraspecific
b4= runif(1, -0.1, 0.05) # influence of Mesocoms on Dp interspecific

#interclonal variance
o3= 0.1

#Treatment
b5= runif(1, -0.001, 0.001) # influence of intraC treatment
b6= runif(1, -0.001, 0.001) # influence of intraS treatment
b7= runif(1, -0.001, 0.001) # influence of intraS treatment
b8= runif(1, -0.001, 0.001) # influence of intraS treatment
b9= runif(1, -0.001, 0.001) # influence of interS treatment
b10=runif(1, -0.001, 0.001) # influence of interS treatment 
b11=runif(1, -0.001, 0.001) # influence of interS treatment

#intraclonal variance
o4= 0.01


##Model for Meso-clones grazing rate
for (i in 1:nrow(clonedata)){
  if (clonedata[i,2] == 0){
    ax= rnorm(1,a1,o1) # Normal Distribution of Dm
    y= ax + b1 * (1 - clonedata[i,3])+ b2 * clonedata[i,3] #Evolution in mesocosms
    clonedata[i,13]= rnorm(1, y, o3) # Inter-Clone variance, is it needed?
  } else {
    ax= rnorm(1,a2,o2) # Normal Distribution of Dp
    y= ax + b3 * (1 - clonedata[i,4])+ b4 * clonedata[i,4] #Effect of mesocosm (evolution)
    clonedata[i,13]=rnorm(1, y, o3) # Inter-Clone variance, is it needed?
  }
}

##Model for treatment grazing rate (? Problem is, beta parameter similar for each mesocosm treatment? Example:Maybe inter-meso-clones are better adapted to inter-treatment?)
v = merge(clonedata[c("clone","V13")], fulldata, by="clone")

for (j in 1:nrow(fulldata)){
  z= v[j,2] + b5*fulldata[j,6] + b6*fulldata[j,7] + b7*fulldata[j,8] + b8*fulldata[j,9] + b9*fulldata[j,10] + b10*fulldata[j,11] + b11*fulldata[j,12] # Effect of treatments
  fulldata[j,13:15]= rnorm(3,z,o4)
}


##plotting
par(mfrow=c(2,2))
# Plot 1: dm intra
boxplot(fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="0" & fulldata$Trta=="1"],fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="0" & fulldata$Trtb=="1"],fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="0" & fulldata$Trtc=="1"],fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="0" & fulldata$Trtd=="1"],fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="0" & fulldata$Trte=="1"],fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="0" & fulldata$Trtf=="1"],fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="0" & fulldata$Trtg=="1"],fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="0" & fulldata$Trth=="1"],ylab="grazing rate (micrograms C)",names=c("a","b","c","d","e","f","g","h"), main="D. magna, intraspecific, (n=4)",ylim = c(0.5, 2))
# Plot 2: dm inter
boxplot(fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="1" & fulldata$Trta=="1"],fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="1" & fulldata$Trtb=="1"],fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="1" & fulldata$Trtc=="1"],fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="1" & fulldata$Trtd=="1"],fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="1" & fulldata$Trte=="1"],fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="1" & fulldata$Trtf=="1"],fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="1" & fulldata$Trtg=="1"],fulldata$V15[fulldata$Species=="0" & fulldata$MesoDm=="1" & fulldata$Trth=="1"],ylab="grazing rate (micrograms C)",names=c("a","b","c","d","e","f","g","h"), main="D. magna, interspecific, (n=4)",ylim = c(0.5, 2))
# Plot 3: dp intra
boxplot(fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="0" & fulldata$Trta=="1"],fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="0" & fulldata$Trtb=="1"],fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="0" & fulldata$Trtc=="1"],fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="0" & fulldata$Trtd=="1"],fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="0" & fulldata$Trte=="1"],fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="0" & fulldata$Trtf=="1"],fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="0" & fulldata$Trtg=="1"],fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="0" & fulldata$Trth=="1"],ylab="grazing rate (micrograms C)",names=c("a","b","c","d","e","f","g","h"), main="D. pulex, intraspecific, (n=13)",ylim = c(0.5, 2))
# Plot 4: dp inter
boxplot(fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="1" & fulldata$Trta=="1"],fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="1" & fulldata$Trtb=="1"],fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="1" & fulldata$Trtc=="1"],fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="1" & fulldata$Trtd=="1"],fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="1" & fulldata$Trte=="1"],fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="1" & fulldata$Trtf=="1"],fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="1" & fulldata$Trtg=="1"],fulldata$V15[fulldata$Species=="1" & fulldata$MesoDp=="1" & fulldata$Trth=="1"],ylab="grazing rate (micrograms C)",names=c("a","b","c","d","e","f","g","h"), main="D. pulex, interspecific, (n=14)",ylim = c(0.5, 2))

