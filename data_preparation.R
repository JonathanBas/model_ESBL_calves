rm(list=ls(all=TRUE))

library(readxl)
library(rjags)
library(reshape2)

# Data preparation:

# Importation of carriage data with 11/12 samplings:
portage_old = as.data.frame(read.csv("~/carriage_calves.csv"))

# Table with every day:
portage_real = as.data.frame(matrix(NA, nrow=45, ncol=164))
portage_real[,(1:3)] = portage_old[,(1:3)]
colnames(portage_real) = c(colnames(portage_old)[1:3], as.character(1:161))

# NA between sampling dates:
for (i in 1:45){
  insert_na = cbind(t(portage_old[i,-(1:3)]), matrix(as.factor(NA),length(portage_old[i,-(1:3)]),13))
  insert_na = as.vector(t(insert_na))[-(156:168)]
  portage_real[i,-(1:3)] = c(rep(as.factor(NA),6), insert_na)
  if(portage_real$ID[i] != portage_old$ID[i]){print("Error")}
}

# Because one of the sampling date is shifted by one day:
replace = portage_real$"105"
portage_real$"105" = portage_real$"106"
portage_real$"106" = replace

sampling_days = which(!is.na(portage_real[1,-(1:3)]))
N_calves = nrow(portage_real)
t_end = ncol(portage_real[,-(1:3)])
from_56 = c((1:31),(33:N_calves))

# Antibiotics exposure:
expo = array(data=0, dim=c(N_calves,161,5))
dimnames(expo)[[1]] = portage_real$ID
dimnames(expo)[[2]] = seq(1,161)
dimnames(expo)[[3]] = c("Pe","Po","T","M","ST")

expo[which(portage_real$Farm == "A"), (1:10), c("Po","ST")] = 1
expo[which(portage_real$Farm == "A"), (11:20), "T"] = 1
expo[which(portage_real$Farm == "A"), 53, "T"] = 1
expo[which(portage_real$Farm == "A"), 136, "Pe"] = 1
expo[which(portage_real$Farm == "B"), 26, c("T","M")] = 1
expo[which(portage_real$Farm == "B"), 90, "T"] = 1
expo[which(portage_real$Farm == "B"), 101, "Pe"] = 1
expo[which(portage_real$Farm == "C"), (3:8), "ST"] = 1
expo[which(portage_real$Farm == "C"), (10:16), "T"] = 1
expo[which(portage_real$Farm == "C"), (20:24), "T"] = 1
expo[which(portage_real$Farm == "C"), (25:27), "M"] = 1
expo[which(portage_real$Farm == "C"), (80:82), "T"] = 1

expo_tot = rowSums(expo, dims = 2)
expo_tot = (expo_tot > 0)
expo_tot = 1*expo_tot

# Matrices with farm and pen of each calf:

farm_bel = as.data.frame(matrix(0, N_calves, 4))
colnames(farm_bel) = c("A","B","C","num_farm")

box_bel = as.data.frame(matrix(0, N_calves, 10))
colnames(box_bel) = c("A1","A2","A3","B1","B2","B3","C1","C2","C3","num_box")

rownames(farm_bel) = rownames(box_bel) = portage_real$ID

for(i in 1:N_calves){
  which_farm = portage_real$Farm[i]
  which_box = portage_real$Box[i]
  
  farm_bel$num_farm[i] = which(colnames(farm_bel) == which_farm)
  box_bel$num_box[i] = which(colnames(box_bel) == which_box)
  
  farm_bel[i, which_farm] = 1
  box_bel[i, which_box] = 1
}
box_bel[which(is.na(box_bel), arr.ind=T)] = 0

indiv_farm_1 = 1:15
indiv_farm_2 = 16:30
indiv_farm_3 = 31:45

indiv_farm_1_from_56 = indiv_farm_2_from_56 = 1:15
indiv_farm_3_from_56 = c(1, 3:15)

rm(replace, insert_na, i, j, which_box, which_farm)

# Transformation of carriage matrix according to modelling hypothesis:
portage_real = 1*(portage_real[,-c(1,2,3)] !=0)
portage_supp = as.data.frame(portage_real)
for (i in c(1:N_calves)){
  for (j in 1:(length(sampling_days)-1)){
    
    samp_day = as.numeric(sampling_days[j])
    next_samp_day = as.numeric(sampling_days[j+1])
    
    if (!(is.na(portage_real[i, samp_day]) | is.na(portage_real[i, next_samp_day]))){
      if (portage_real[i, samp_day] == portage_real[i, next_samp_day]){
        portage_supp[i, (samp_day:next_samp_day)] = portage_real[i, samp_day]
      }
    }
  }
}

# Positive and negative calves on first day:
init_pos = c(1:15, 21:30, 35:45)
init_neg = c(16:20, 31:34)

# 3D matrix with, for each calf and each day, the delay since last exposure to each AB class:
# "last_expo_class"
# - Calf 1 to 45
# - Week 1 to t_end
# - AB class 1 to 5
last_expo_class = expo

for(c in 1:5){
  for(i in 1:N_calves){
    last_date = NA
    for(j in 1:t_end){
      if(expo[i, j, c] == 1){
        last_date = j
      }
      last_expo_class[i, j, c] = j-last_date
    }
  }
}
last_expo_class[which(is.na(last_expo_class))] = 10000

# Function for conversion from format mcarray (output from jags.samples) to mcmc.list (output from coda.samples):
to_mcmc.list = function(mcarray_obj, to_burn, which_par = "all"){
  nchains = length(as.mcmc.list(mcarray_obj[["deviance"]]))
  
  if(which_par == "all"){
    var_mod = names(mcarray_obj)[which(!(names(mcarray_obj) %in% c("pD")))]
  }else{
    var_mod = names(mcarray_obj)[which(!(names(mcarray_obj) %in% c("pD", "N_intro_plus1_f2", "N_intro_plus1_f3", "p1_f3[3]", "p2_f3[3]")))]
  }
  
  Mch = as.list(rep(NA, nchains))
  for (chain_i in 1:nchains){
    for (na in var_mod){
      post_var_chain_i = as.mcmc.list(mcarray_obj[[na]])[[chain_i]]
      Mch[[chain_i]] = cbind(Mch[[chain_i]], post_var_chain_i)
    }
    included_in_chain = (1:nrow(Mch[[chain_i]]))[! (1:nrow(Mch[[chain_i]])) %in% to_burn]
    Mch[[chain_i]] = mcmc(Mch[[chain_i]][included_in_chain,-1])
  }
  output_Mch = as.mcmc.list(Mch)
  output_combined_Mch = do.call("rbind", Mch)
  
  list(output_var_mod = var_mod, output_Mch=output_Mch, output_combined_Mch=output_combined_Mch)
}

portage_obs = portage_supp
