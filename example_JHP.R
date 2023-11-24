## example analysis for Daphnia grazing rate experiment
## author: Lutz Becks
## This script takes features extracted from images of the grazing rate experiment samples (seen in line #20 below) and performs k-means clustering on them. In the datasheet, each row is a unique cell - you can see that there are multiple rows per well (from the 96-well plate), and so this is information on features extracted from each cell that was identified in the image. The script performs k-means clustering to classify each row (cell) (per well) into one of 5 categories. The number of cells in each cluster is then calculated. This is all done in the loop in lines 17-25. You can see that line 24 gives the cell counts per "species" (per k-means cluster). We put 5 algal species in each experimental bottle (so that is the reason for 5 clusters).

setwd("/Users/jhpantel/Nextcloud/AG_Pantel/phd_Borgmann_EcoEvoStats/code/grazing_analysis_example")
#setwd("/Users/jhpantel/Library/CloudStorage/OneDrive-Personal/Documents/AUP/Research/Konstanz/Compex/code/phytoplankton")



test=read.delim("ID_5910.txt",dec='.',header = TRUE)
str(test)

treat=c('A','B','C','D','E','F','G','H')
number=c('one','two','two','two','two','two','two','two')
comb=c('none', 'focal','anc1_DP','anc2_DP','anc3_DP','anc1_DM','anc2_DM','anc3_DM')
#run loop for with subset for each well and get frequency table 
freqs_5910=data.frame()
wells=unique(test$Well.Name)
for(i in 1:length(wells)){
  re=wells[i]
  data=subset(test,Well.Name == re)
  kl=kmeans(data[,10:21],5)
  data$cluster=kl$cluster
  tb=table(data$cluster)
  temp=data.frame(well.name=re,cells=sum(tb[1],tb[2],tb[3],tb[4],tb[5]),G1=tb[1],G3=tb[3],G2=tb[2],G4=tb[4],G5=tb[5])
  freqs_5910=rbind(freqs_5910,temp)
}
# Here we organize the counts according to the treatments (1-8), the focal clone identity (Dp_25.2, Dp_18.8, Dp_26.8. Dp_27.9), and the mesocosm treatment that the Daphnia came from.
freqs_5910$treat=rep(c(treat,rev(treat)),c((dim(freqs_5910)[1])/16))
freqs_5910$number=rep(c(number,rev(number)),c((dim(freqs_5910)[1])/16))
freqs_5910$comb=rep(c(comb,rev(comb)),c((dim(freqs_5910)[1])/16))
freqs_5910$rep=c(rep("1",8),rep("2",8),rep("3",8))
freqs_5910$spec=c(rep('Dp',c((dim(freqs_5910)[1]))))
freqs_5910$clone=c(rep('Dp_25.2',24),rep('Dp_18.8',24),rep('Dp_26.8',24),rep('Dp_27.9',24))
freqs_5910$nut=c(rep('low',24),rep('low',24),rep('low',24),rep('low',24))
freqs_5910$Daph=c(rep('both',24),rep('single',24),rep('single',24),rep('single',24))

# calculate ingestion rates in ug C using the start values 
# Lutz must have counts of cells in each cluster from separate images - that is why numbers are entered directly in the formula (1610,1498,1845,2225,1633). These numbers give the cell count at the beginning of the experiment (hour 0). So we subtract to get the change in cell number, and divide that by 24 to get the per-hour rate. This cell count value is multiplied by the total mass of Carbon (in micrograms) associated with each species (again, these values are taken from another source: 210,122,31,22,23). The formula calculated is thus:
# [(cell_count_h0 - cell_count_h24) / 24] * cell_mass_ug_C
# the values thus give the grazing rate of Daphnia in treatment 1-8 per hour, for a given species, in units of micrograms Carbon / hour.
freqs_5910$IngG1=210*((1610-freqs_5910$G1)/24)
freqs_5910$IngG2=122*((1498-freqs_5910$G2)/24)
freqs_5910$IngG3=31*((1845-freqs_5910$G3)/24)
freqs_5910$IngG4=22*((2225-freqs_5910$G4)/24)
freqs_5910$IngG5=23*((1633-freqs_5910$G5)/24)

### Here I make a plot with grazing rate all the treatments, pooling the pairings with ancestral Dp and ancestral Dm.
### JHP 29-06-2022
boxplot(freqs_5910$IngG1[freqs_5910$treat=="A"],freqs_5910$IngG1[freqs_5910$treat=="B"],freqs_5910$IngG1[freqs_5910$treat=="C" | freqs_5910$treat=="D" | freqs_5910$treat=="E" ],freqs_5910$IngG1[freqs_5910$treat=="F" | freqs_5910$treat=="G" | freqs_5910$treat=="H" ],ylab="grazing rate (micrograms C)",names=c("alone","Dp+Dp","Dp+ref_Dp","Dp+ref_Dm"))

boxplot(freqs_5910$IngG2[freqs_5910$treat=="A"],freqs_5910$IngG2[freqs_5910$treat=="B"],freqs_5910$IngG2[freqs_5910$treat=="C" | freqs_5910$treat=="D" | freqs_5910$treat=="E" ],freqs_5910$IngG2[freqs_5910$treat=="F" | freqs_5910$treat=="G" | freqs_5910$treat=="H" ],ylab="grazing rate (micrograms C)",names=c("alone","Dp+Dp","Dp+ref_Dp","Dp+ref_Dm"))

boxplot(freqs_5910$IngG3[freqs_5910$treat=="A"],freqs_5910$IngG3[freqs_5910$treat=="B"],freqs_5910$IngG3[freqs_5910$treat=="C" | freqs_5910$treat=="D" | freqs_5910$treat=="E" ],freqs_5910$IngG3[freqs_5910$treat=="F" | freqs_5910$treat=="G" | freqs_5910$treat=="H" ],ylab="grazing rate (micrograms C)",names=c("alone","Dp+Dp","Dp+ref_Dp","Dp+ref_Dm"))

boxplot(freqs_5910$IngG4[freqs_5910$treat=="A"],freqs_5910$IngG4[freqs_5910$treat=="B"],freqs_5910$IngG4[freqs_5910$treat=="C" | freqs_5910$treat=="D" | freqs_5910$treat=="E" ],freqs_5910$IngG4[freqs_5910$treat=="F" | freqs_5910$treat=="G" | freqs_5910$treat=="H" ],ylab="grazing rate (micrograms C)",names=c("alone","Dp+Dp","Dp+ref_Dp","Dp+ref_Dm"))

boxplot(freqs_5910$IngG5[freqs_5910$treat=="A"],freqs_5910$IngG5[freqs_5910$treat=="B"],freqs_5910$IngG5[freqs_5910$treat=="C" | freqs_5910$treat=="D" | freqs_5910$treat=="E" ],freqs_5910$IngG5[freqs_5910$treat=="F" | freqs_5910$treat=="G" | freqs_5910$treat=="H" ],ylab="grazing rate (micrograms C)",names=c("alone","Dp+Dp","Dp+ref_Dp","Dp+ref_Dm"))

### Here I make a plot with grazing rate all of the 8 treatments, per algal species.
boxplot(IngG1 ~ comb,data=freqs_5910,ylab="grazing rate (micrograms C)")
mod1 <- lm(IngG1 ~ comb,data=freqs_5910)

boxplot(IngG2 ~ comb,data=freqs_5910,ylab="grazing rate (micrograms C)")
mod2 <- lm(IngG2 ~ comb,data=freqs_5910)


