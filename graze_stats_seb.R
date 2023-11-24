#############################Constants#############################

## Set the random seed ##
set.seed(43)

## constants ##
rep <- 3

## parameters ##

# Level 3 - Species #
B_dp <- 1.3 # Fixed D. pulex grazing rate
B_dm <- 1.7 # Fixed D. magna grazing rate

# Level 4 # Treatment
B_dp_intraC <- 1.1 # Fixed D. pulex intra-clonal grazing rate - I should probably modify that this depends on clone.
B_dm_intraC <- 1.5 # Fixed D. magna intra-clonal grazing rate - I should probably modify that this depends on clone.

d1 <- 0.2 # Fixed D. pulex among-clonal variance in grazing rate
d2 <- 0.2 # Fixed D. magna among-clonal variance in grazing rate

dp_intraS <- 1.2 # Center of D. pulex intraspecifc grazing rate
d3 <- 0.1 # Fixed D. pulex variance in intraspecific grazing

dp_interS <- 0.8 # Center of D. pulex interspecific grazing rate with D. magna
d4 <- 0.1 # Fixed D. pulex variance in interspecific grazing with D. magna

dm_intraS <- 1.6 # Center of D. magna intraspecifc grazing rate
d5 <- 0.1 # Fixed D. magna variance in intraspecific grazing

dm_interS <- 1.68 # Center of D. magna interspecific grazing rate with D. pulex
d6 <- 0.1 # Fixed D. magna variance in interspecific grazing with D. pulex


#############################functions#############################
##grazing rate generating function ##
g_gen <- function(rep=3,B_dp=1.3,B_dp_intraC=1.1,B_dm=1.7,B_dm_intraC=1.5,d1=0.2,d2=0.2,dp_intraS=1.2,d3=0.1,dp_interS=0.8,d4=0.1,dm_intraS=1.6,d5=0.1,dm_interS=1.68,d6=0.1,Spec,Alone,IntraC,Comp){
  
  ##Geting 3 values out of normal distribution of 'daphnia species' 'intra/interspecific' grazing rates
  B_dp_intraS <- rnorm(3,dp_intraS,d3)
  B_dp_interS <- rnorm(3,dp_interS,d4)
  B_dm_intraS <- rnorm(3,dm_intraS,d5)
  B_dm_interS <- rnorm(3,dm_interS,d6)
  
  ## Setting the grazing rates to if there is Comp 0(same species setup) or 1(different species setup) = If they would come out of coexistence or no coexistence mesocosms 
  B3 <- B_dp_intraS*(1-Comp) + B_dp_interS*Comp # B_dp_intraS is the center of the intraspecific grazing rate for D. pulex, and Comp takes a value of 0 for Treatments with the same species. B_dp_interS is the center of the interspecific grazing rate for D. pulex against D. magna, and Comp takes a value of 1 for treatments with a seperate species.
  B4 <- B_dm_intraS*(1-Comp) + B_dm_interS*Comp # B_dm_intraS is the center of the intraspecific grazing rate for D. magna, and Comp takes a value of 0 for Treatments with the same species. B_dm_interS is the center of the interspecific grazing rate for D. magna against D. pulex, and Comp takes a value of 1 for treatments with a seperate species.
  
  B_dp_comp <- B_dp_intraC*IntraC + B3*(1-IntraC) # B_dp_intraC is the intra-clonal grazing rate (competition should be most intense when another individual of the same clone is present) and IntraC takes a value of 1 for Treatment B (when the Daphnia is in competition with an individual of the same clone) and a value of 0 for the other treatments.
  B_dm_comp <- B_dm_intraC*IntraC + B4*(1-IntraC) # B_dm_intraC is the intra-clonal grazing rate (competition should be most intense when another individual of the same clone is present) and IntraC takes a value of 1 for Treatment B (when the Daphnia is in competition with an individual of the same clone) and a value of 0 for the other treatments.
  
  B1 <- B_dp*Alone + B_dp_comp*(1-Alone) # B_dp is the fixed D. pulex grazing rate center, Alone takes a value of 1 for Treatment A (when the Daphnia is in isolation) and 0 for any other treatment (when the Daphnia is in competition).
  B2 <- B_dm*Alone + B_dm_comp*(1-Alone) # B_dm is the fixed D. magna grazing rate center, Alone takes a value of 1 for Treatment A (when the Daphnia is in isolation) and 0 for any other treatment (when the Daphnia is in competition).
  
  g_s <- B1*Spec + B2*(1-Spec) # B1 is the D. pulex grazing rate center, B2 is the D. magna grazing rate center. Spec takes a value of 1 for Dp and 0 for Dm.
  delta_c <- d1*Spec + d2*(1-Spec) # d1 is the Dp among-clone variance, d2 is the Dm among-clone variance.
  
  g <- rnorm(rep,g_s,delta_c) # rep is the number of replicates per treatment
  return(g)
}

##grazing rate sim function##

g_sim <- function(focal_species,treatment){
  # Generates a grazing rate value given the treatment of the grazing rate experiment and the species of the focal clone
  if(focal_species=="Dm"){
    Spec <- 0
    if(treatment=="A"){
      Alone <- 1
      IntraC <- 0
      Comp <- 0
    } else if(treatment=="B"){
      Alone <- 0
      IntraC <- 1
      Comp <- 0
    } else if(treatment == "C" | treatment == "D" | treatment == "E"){
      Alone <- 0
      IntraC <- 0
      Comp <- 1
    } else if(treatment == "F" | treatment == "G" | treatment == "H"){
      Alone <- 0
      IntraC <- 0
      Comp <- 0
      }
    }
  else if (focal_species=="Dp"){
    Spec <- 1
    if(treatment=="A"){
      Alone <- 1
      IntraC <- 0
      Comp <- 0
    } else if(treatment=="B"){
      Alone <- 0
      IntraC <- 1
      Comp <- 0
    } else if(treatment == "C" | treatment == "D" | treatment == "E"){
      Alone <- 0
      IntraC <- 0
      Comp <- 0
    } else if(treatment == "F" | treatment == "G" | treatment == "H"){
      Alone <- 0
      IntraC <- 0
      Comp <- 1
      }
    }
  g_gen(Spec=Spec,Alone=Alone,IntraC=IntraC,Comp=Comp)
}


#############################testing#############################
treatments=c(LETTERS[1:8])
all_gr=data.frame()

for(i in 1:8){
  all_gr[i,1:3]=g_sim(focal_species="Dm",treatment=treatments[i])
  all_gr[i,4:6]=g_sim(focal_species="Dp",treatment=treatments[i])
}

colnames(all_gr) <- c('Dm_rep1','Dm_rep2','Dm_rep3','Dp_rep1','Dp_rep2','Dp_rep3') 


hist(all_gr)


