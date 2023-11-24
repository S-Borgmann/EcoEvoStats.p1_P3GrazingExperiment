#############################Grazing rate simulation#############################

## Set the random seed ##
set.seed(43)

## constants ##
rep <- 3

## Parameter

#Species / Ancestral clones
g1 <- 1.5 # Fixed grazing rate, D. pulex
g2 <- 1.7 # Fixed grazing rate, D. magna 
vg1 <- 0.2 # Fixed variance in grazing rate, D. pulex 
vg2 <- 0.2 # Fixed variance in grazing rate, D. magna 

#Mesocosm
a1 <- 0.15 # change in mean grazing rate, intraspecific mesocosm, Dp 
a2 <- 0.2 # change in mean grazing rate, intraspecific mesocosm, Dm
va1 <- 0.1 # intraclonal grazing variance, intraspecific mesocosm, Dp
va2 <- 0.1 # intraclonal grazing variance, intraspecific mesocosm, Dm

b1 <- -0.15 # change in mean grazing rate, interspecific mesocosm, Dp
b2 <- -0.2 # change in mean grazing rate, interspecific mesocosm, Dm
vb1 <- 0.1 # intraclonal grazing variance, interspecific mesocosm, Dp
vb2 <- 0.1 # intraclonal grazing variance, interspecific mesocosm, Dm

#Setup
c1 <- -0.2 #change in mean grazing rate, intraclonal setup, Dp 
c2 <- -0.2 #change in mean grazing rate, intraclonal setup, Dm

d1 <- -0.1 #change in mean grazing rate, intraspec setup, Dp
d2 <- -0.1 #change in mean grazing rate, intraspec setup, Dm

e1 <- -0.5 #change in mean grazing rate, interspec setup, Dp
e2 <- -0.02 #change in mean grazing rate, interspec setup, Dm

#clonal variation
vf1= 0.2 #clonal variance Dp
vf2= 0.2 #clonal variance Dm

#binaries
# Spec= 0 (Dp)/ 1 (Dm)
# Meso = 0 (intra), 1 (inter)
# Setup = 0 (alone), 1 (intraC), 2 (intraS), 3 (interS)


## Model function
g_gen <- function(rep,g1,g2,vg1,vg2,a1,a2,va1,va2,b1,b2,vb1,vb2,c1,c2,d1,d2,e1,e2,vf1,vf2,Spec,Meso,Setup){
  
  Sp = rnorm(1,(g1*(1-Spec) + g2*Spec),(vg1*(1-Spec) + vg2*Spec))
  
  M = rnorm(1, (Sp+ ((a1*(1-Spec)+a2*(Spec))*(1-Meso) + (b1*(1-Spec)+b2*(Spec))*(Meso)))
            ,((va1*(1-Spec)+va2*(Spec))*(1-Meso) + (vb1*(1-Spec)+vb2*(Spec))*(Meso)))
  
  if (Setup == 0){
    S= rnorm(3, M, (vf1*(1-Spec)+vf2*(Spec)))
  
  } else if (Setup == 1) {
    S = rnorm(3, M + (c1*(1-Spec)+c2*(Spec)), (vf1*(1-Spec)+vf2*(Spec)))
  
  } else if (Setup == 2) {
    S = rnorm(3, M + (d1*(1-Spec)+d2*(Spec)), (vf1*(1-Spec)+vf2*(Spec)))
    
  } else if (Setup == 3) {
    S = rnorm(3, M + (e1*(1-Spec)+e2*(Spec)), (vf1*(1-Spec)+vf2*(Spec)))  }
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
all_gr_intra=data.frame()
all_gr_inter=data.frame()

for(i in 1:8){
  all_gr_intra[i,1:3]=g_sim(focal_species="Dm",Meso=0,treatment=treatments[i])
  all_gr_intra[i,4:6]=g_sim(focal_species="Dp",Meso=0,treatment=treatments[i])
}

for(i in 1:8){
  all_gr_inter[i,1:3]=g_sim(focal_species="Dm",Meso=1,treatment=treatments[i])
  all_gr_inter[i,4:6]=g_sim(focal_species="Dp",Meso=1,treatment=treatments[i])
}

colnames(all_gr_intra) <- c('Dm_rep1','Dm_rep2','Dm_rep3','Dp_rep1','Dp_rep2','Dp_rep3') 
colnames(all_gr_inter) <- c('Dm_rep1','Dm_rep2','Dm_rep3','Dp_rep1','Dp_rep2','Dp_rep3') 


plot(c(1:8), all_gr_intra$Dm_rep1, type = "o", col = 1, ylim = c(0, 3))
lines(c(1:8), all_gr_intra$Dm_rep2, type = "o", col = 1, ylim = c(0, 3)) 
lines(c(1:8), all_gr_intra$Dm_rep3, type = "o", col = 1, ylim = c(0, 3)) 
lines(c(1:8), all_gr_intra$Dp_rep1, type = "o", col = 2, ylim = c(0, 3)) 
lines(c(1:8), all_gr_intra$Dp_rep2, type = "o", col = 2, ylim = c(0, 3)) 
lines(c(1:8), all_gr_intra$Dp_rep3, type = "o", col = 2, ylim = c(0, 3)) 

plot(c(1:8), all_gr_inter$Dm_rep1, type = "o", col = 1, ylim = c(0, 3))
lines(c(1:8), all_gr_inter$Dm_rep2, type = "o", col = 1, ylim = c(0, 3)) 
lines(c(1:8), all_gr_inter$Dm_rep3, type = "o", col = 1, ylim = c(0, 3)) 
lines(c(1:8), all_gr_inter$Dp_rep1, type = "o", col = 2, ylim = c(0, 3)) 
lines(c(1:8), all_gr_inter$Dp_rep2, type = "o", col = 2, ylim = c(0, 3)) 
lines(c(1:8), all_gr_inter$Dp_rep3, type = "o", col = 2, ylim = c(0, 3)) 
     