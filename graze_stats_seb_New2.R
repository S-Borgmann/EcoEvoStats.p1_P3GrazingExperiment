#############################Grazing rate simulation#############################

## Set the random seed ##
#set.seed(43)

## constants ##
rep <- 3

## Parameter

#Species / Ancestral clones
g1 <- 1.5 # Fixed grazing rate, D. pulex
g2 <- 1.7 # Fixed grazing rate, D. magna 
vg1 <- 0.02 # Fixed variance in grazing rate, D. pulex 
vg2 <- 0.02 # Fixed variance in grazing rate, D. magna 

#Mesocosm
a1 <- 1.01 # change in mean grazing rate %, intraspecific mesocosm, Dp
a2 <- 1.05 # change in mean grazing rate %, intraspecific mesocosm, Dm
va1 <- 0.01 # grazing variance, intraspecific mesocosm, Dp
va2 <- 0.01 # grazing variance, intraspecific mesocosm, Dm

b1 <- 0.9 # change in mean grazing rate %, interspecific mesocosm, Dp
b2 <- 0.98 # change in mean grazing rate %, interspecific mesocosm, Dm
vb1 <- 0.01 # grazing variance, interspecific mesocosm, Dp
vb2 <- 0.01 # grazing variance, interspecific mesocosm, Dm

# a1 <- 0.15 # change in mean grazing rate, intraspecific mesocosm, Dp
# a2 <- 0.2 # change in mean grazing rate, intraspecific mesocosm, Dm
# va1 <- 0.1 # grazing variance, intraspecific mesocosm, Dp
# va2 <- 0.1 # grazing variance, intraspecific mesocosm, Dm
# 
# b1 <- -0.15 # change in mean grazing rate, interspecific mesocosm, Dp
# b2 <- -0.2 # change in mean grazing rate, interspecific mesocosm, Dm
# vb1 <- 0.1 # grazing variance, interspecific mesocosm, Dp
# vb2 <- 0.1 # grazing variance, interspecific mesocosm, Dm

#Setup
c1 <- 0.9 #change in mean grazing rate %, intraclonal setup, Dp 
c2 <- 0.9 #change in mean grazing rate %, intraclonal setup, Dm

d1 <- 0.9 #change in mean grazing rate %, intraspec setup, Dp
d2 <- 0.9 #change in mean grazing rate %, intraspec setup, Dm

e1 <- 0.88 #change in mean grazing rate %, interspec setup, Dp
e2 <- 0.98 #change in mean grazing rate %, interspec setup, Dm

# c1 <- -0.2 #change in mean grazing rate, intraclonal setup, Dp 
# c2 <- -0.2 #change in mean grazing rate, intraclonal setup, Dm
# 
# d1 <- -0.1 #change in mean grazing rate, intraspec setup, Dp
# d2 <- -0.1 #change in mean grazing rate, intraspec setup, Dm
# 
# e1 <- -0.1 #change in mean grazing rate, interspec setup, Dp
# e2 <- -0.02 #change in mean grazing rate, interspec setup, Dm

#clonal variation
vf1= 0.02 #clonal variance Dp
vf2= 0.02 #clonal variance Dm

#binaries
# Spec= 0 (Dp)/ 1 (Dm)
# Meso = 0 (intra), 1 (inter)
# Setup = 0 (alone), 1 (intraC), 2 (intraS), 3 (interS)


## Model function
g_gen <- function(rep,g1,g2,vg1,vg2,a1,a2,va1,va2,b1,b2,vb1,vb2,c1,c2,d1,d2,e1,e2,vf1,vf2,Spec,Meso,Setup){

  if (Setup == 0){
    S= rnorm(3, M, (vf1*(1-Spec)+vf2*(Spec)))
  
  } else if (Setup == 1) {
    S = rnorm(3, M * (c1*(1-Spec)+c2*(Spec)), (vf1*(1-Spec)+vf2*(Spec)))
  
  } else if (Setup == 2) {
    S = rnorm(3, M * (d1*(1-Spec)+d2*(Spec)), (vf1*(1-Spec)+vf2*(Spec)))
    
  } else if (Setup == 3) {
    S = rnorm(3, M * (e1*(1-Spec)+e2*(Spec)), (vf1*(1-Spec)+vf2*(Spec)))  }
}   

#Test=g_gen(rep,g1,g2,vg1,vg2,a1,a2,va1,va2,b1,b2,vb1,vb2,c1,c2,d1,d2,e1,e2,vf1,vf2,Spec,Meso,Setup)


##grazing rate sim function##
g_sim <- function(focal_species,Meso,treatment){
  # Generates a grazing rate value given the treatment of the grazing rate experiment and the species of the focal clone
  if(focal_species=="Dp"){
    Spec <- 0
    if(treatment=="A"){
      Setup <- 0
    } else if(treatment=="B"){
      Setup <- 1
    } else if(treatment == "C" | treatment == "D" | treatment == "E"){
      Setup <- 2
    } else if(treatment == "F" | treatment == "G" | treatment == "H"){
      Setup <- 3
    }
  }
  else if (focal_species=="Dm"){
    Spec <- 1
    if(treatment=="A"){
      Setup <- 0
    } else if(treatment=="B"){
      Setup <- 1
    } else if(treatment == "C" | treatment == "D" | treatment == "E"){
      Setup <- 3
    } else if(treatment == "F" | treatment == "G" | treatment == "H"){
      Setup <- 2
    }
  }
  g_gen(rep,g1,g2,vg1,vg2,a1,a2,va1,va2,b1,b2,vb1,vb2,c1,c2,d1,d2,e1,e2,vf1,vf2,Spec=Spec,Meso=Meso,Setup=Setup)
}


#############################testing#############################
treatments=c(LETTERS[1:8])
infos <- read.delim("Nextcloud/Stuff/Projects/EcoEvoStats.p1/Project3_Grazing_experiment/grazing_analysis_example/GrazingExperimentInfos.txt")
infos["species"][infos["species"] == "m"] <- 1
infos["species"][infos["species"] == "p"] <- 0
infos["Coexisting.before."][infos["Coexisting.before."] == "n"] <- 0
infos["Coexisting.before."][infos["Coexisting.before."] == "y"] <- 1
all_gr= data.frame()

for (n in 1:nrow(infos)) {
  
  if (infos[n,3] == 0 & infos[n,4] == 0 ){
    ##Dp - Intra##
    all_gr_intra_Dp=data.frame()
    Spec=0
    Meso=0
    
    Sp = rnorm(1,(g1*(1-Spec) + g2*Spec),(vg1*(1-Spec) + vg2*Spec))
    
    M = rnorm(1, (Sp* ((a1*(1-Spec)+a2*(Spec))*(1-Meso) + (b1*(1-Spec)+b2*(Spec))*(Meso)))
              ,((va1*(1-Spec)+va2*(Spec))*(1-Meso) + (vb1*(1-Spec)+vb2*(Spec))*(Meso)))
    
    for(i in 1:8){
      all_gr_intra_Dp[i,1:3]=g_sim(focal_species="Dp",Meso=0,treatment=treatments[i])
    }
    colnames(all_gr_intra_Dp) <- c('rep1','rep2','rep3') 
    all_gr <- rbind(all_gr,all_gr_intra_Dp) }
  else if (infos[n,3] == 0 & infos[n,4] == 1 ){
    ##Dp - Inter##
    all_gr_inter_Dp=data.frame()
    Spec=0
    Meso=1
    
    Sp = rnorm(1,(g1*(1-Spec) + g2*Spec),(vg1*(1-Spec) + vg2*Spec))
    
    M = rnorm(1, (Sp* ((a1*(1-Spec)+a2*(Spec))*(1-Meso) + (b1*(1-Spec)+b2*(Spec))*(Meso)))
              ,((va1*(1-Spec)+va2*(Spec))*(1-Meso) + (vb1*(1-Spec)+vb2*(Spec))*(Meso)))
    
    for(i in 1:8){
      all_gr_inter_Dp[i,1:3]=g_sim(focal_species="Dp",Meso=1,treatment=treatments[i])
    }
    colnames(all_gr_inter_Dp) <- c('rep1','rep2','rep3') 
    all_gr <- rbind(all_gr,all_gr_inter_Dp)}
    else if (infos[n,3] == 1 & infos[n,4] == 0 ){
      ##Dm - Intra##
      all_gr_intra_Dm=data.frame()
      Spec=1
      Meso=0
      
      Sp = rnorm(1,(g1*(1-Spec) + g2*Spec),(vg1*(1-Spec) + vg2*Spec))
      
      M = rnorm(1, (Sp* ((a1*(1-Spec)+a2*(Spec))*(1-Meso) + (b1*(1-Spec)+b2*(Spec))*(Meso)))
                ,((va1*(1-Spec)+va2*(Spec))*(1-Meso) + (vb1*(1-Spec)+vb2*(Spec))*(Meso)))
      
      for(i in 1:8){
        all_gr_intra_Dm[i,1:3]=g_sim(focal_species="Dm",Meso=0,treatment=treatments[i])
      }
      colnames(all_gr_intra_Dm) <- c('rep1','rep2','rep3') 
      all_gr <- rbind(all_gr,all_gr_intra_Dm)}
    else if (infos[n,3] == 1 & infos[n,4] == 1 ){
      ##Dm - Inter''
      all_gr_inter_Dm=data.frame()
      Spec=1
      Meso=1
      
      Sp = rnorm(1,(g1*(1-Spec) + g2*Spec),(vg1*(1-Spec) + vg2*Spec))
      
      M = rnorm(1, (Sp* ((a1*(1-Spec)+a2*(Spec))*(1-Meso) + (b1*(1-Spec)+b2*(Spec))*(Meso)))
                ,((va1*(1-Spec)+va2*(Spec))*(1-Meso) + (vb1*(1-Spec)+vb2*(Spec))*(Meso)))
      
      for(i in 1:8){
        all_gr_inter_Dm[i,1:3]=g_sim(focal_species="Dm",Meso=1,treatment=treatments[i])
      }
      colnames(all_gr_inter_Dm) <- c('rep1','rep2','rep3') 
      all_gr <- rbind(all_gr,all_gr_inter_Dm)}
}

# ##Dp - Intra##
# all_gr_intra_Dp=data.frame()
# Spec=0
# Meso=0
# 
# Sp = rnorm(1,(g1*(1-Spec) + g2*Spec),(vg1*(1-Spec) + vg2*Spec))
# 
# M = rnorm(1, (Sp* ((a1*(1-Spec)+a2*(Spec))*(1-Meso) + (b1*(1-Spec)+b2*(Spec))*(Meso)))
#           ,((va1*(1-Spec)+va2*(Spec))*(1-Meso) + (vb1*(1-Spec)+vb2*(Spec))*(Meso)))
# 
# for(i in 1:8){
#   all_gr_intra_Dp[i,1:3]=g_sim(focal_species="Dp",Meso=0,treatment=treatments[i])
# }
# colnames(all_gr_intra_Dp) <- c('rep1','rep2','rep3') 
# 
# ##Dp - Inter##
# all_gr_inter_Dp=data.frame()
# Spec=0
# Meso=1
# 
# Sp = rnorm(1,(g1*(1-Spec) + g2*Spec),(vg1*(1-Spec) + vg2*Spec))
# 
# M = rnorm(1, (Sp* ((a1*(1-Spec)+a2*(Spec))*(1-Meso) + (b1*(1-Spec)+b2*(Spec))*(Meso)))
#           ,((va1*(1-Spec)+va2*(Spec))*(1-Meso) + (vb1*(1-Spec)+vb2*(Spec))*(Meso)))
# 
# for(i in 1:8){
#   all_gr_inter_Dp[i,1:3]=g_sim(focal_species="Dp",Meso=1,treatment=treatments[i])
# }
# colnames(all_gr_inter_Dp) <- c('rep1','rep2','rep3') 
# 
# ##Dm - Intra##
# all_gr_intra_Dm=data.frame()
# Spec=1
# Meso=0
# 
# Sp = rnorm(1,(g1*(1-Spec) + g2*Spec),(vg1*(1-Spec) + vg2*Spec))
# 
# M = rnorm(1, (Sp* ((a1*(1-Spec)+a2*(Spec))*(1-Meso) + (b1*(1-Spec)+b2*(Spec))*(Meso)))
#           ,((va1*(1-Spec)+va2*(Spec))*(1-Meso) + (vb1*(1-Spec)+vb2*(Spec))*(Meso)))
# 
# for(i in 1:8){
#   all_gr_intra_Dm[i,1:3]=g_sim(focal_species="Dm",Meso=0,treatment=treatments[i])
# }
# colnames(all_gr_intra_Dm) <- c('rep1','rep2','rep3') 
# 
# ##Dm - Inter''
# all_gr_inter_Dm=data.frame()
# Spec=1
# Meso=1
# 
# Sp = rnorm(1,(g1*(1-Spec) + g2*Spec),(vg1*(1-Spec) + vg2*Spec))
# 
# M = rnorm(1, (Sp* ((a1*(1-Spec)+a2*(Spec))*(1-Meso) + (b1*(1-Spec)+b2*(Spec))*(Meso)))
#           ,((va1*(1-Spec)+va2*(Spec))*(1-Meso) + (vb1*(1-Spec)+vb2*(Spec))*(Meso)))
# 
# for(i in 1:8){
#   all_gr_inter_Dm[i,1:3]=g_sim(focal_species="Dm",Meso=1,treatment=treatments[i])
# }
# colnames(all_gr_inter_Dm) <- c('rep1','rep2','rep3') 
# 




#for(i in 1:8){
 # all_gr_inter[i,1:3]=g_sim(focal_species="Dm",Meso=1,treatment=treatments[i])
  #all_gr_inter[i,4:6]=g_sim(focal_species="Dp",Meso=1,treatment=treatments[i])
#}


#colnames(all_gr_inter) <- c('Dm_rep1','Dm_rep2','Dm_rep3','Dp_rep1','Dp_rep2','Dp_rep3') 

# attach(mtcars)
# par(mfrow=c(2,2))
# 
# plot(c(1:8), all_gr_intra_Dp$rep1, type = "o", col = 1, ylim = c(0, 3))
# lines(c(1:8), all_gr_intra_Dp$rep2, type = "o", col = 2, ylim = c(0, 3)) 
# lines(c(1:8), all_gr_intra_Dp$rep3, type = "o", col = 3, ylim = c(0, 3)) 
# 
# plot(c(1:8), all_gr_inter_Dp$rep1, type = "o", col = 1, ylim = c(0, 3))
# lines(c(1:8), all_gr_inter_Dp$rep2, type = "o", col = 2, ylim = c(0, 3)) 
# lines(c(1:8), all_gr_inter_Dp$rep3, type = "o", col = 3, ylim = c(0, 3)) 
# 
# plot(c(1:8), all_gr_intra_Dm$rep1, type = "o", col = 1, ylim = c(0, 3))
# lines(c(1:8), all_gr_intra_Dm$rep2, type = "o", col = 2, ylim = c(0, 3)) 
# lines(c(1:8), all_gr_intra_Dm$rep3, type = "o", col = 3, ylim = c(0, 3)) 
# 
# plot(c(1:8), all_gr_inter_Dm$rep1, type = "o", col = 1, ylim = c(0, 3))
# lines(c(1:8), all_gr_inter_Dm$rep2, type = "o", col = 2, ylim = c(0, 3)) 
# lines(c(1:8), all_gr_inter_Dm$rep3, type = "o", col = 3, ylim = c(0, 3)) 


