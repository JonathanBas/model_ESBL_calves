# Estimation, chains diagnostic and DIC calculation

# Needs data and functions from data_preparation.R

library(rjags)
library(pracma)
library(parallel)
library(ggplot2)
library(reshape2)
library(bayesplot)

# 1) Estimation:

model_allFarms = function(numb_mod, saving_path = "~/"){
  
  load.module("dic")
  
  n_chains_all = 3
  n_adapt_all = 2000
  length_chain_all = 8000
  thin_all = 40
  
  if(numb_mod == 0){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'from_56' = from_56),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('beta_0_f1','beta_0_f2','beta_0_f3','nu_0','mu','N_intro_plus1_f1','N_intro_plus1_f2','N_intro_plus1_f3','p1_f1[2]','p1_f1[3]','p2_f1[3]','p1_f2[2]','p1_f2[3]','p2_f2[3]','p1_f3[2]','p1_f3[3]','p2_f3[3]','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 1){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'from_56' = from_56),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('beta_f_f1','beta_f_f2','beta_f_f3','nu_0','mu','N_intro_plus1_f1','N_intro_plus1_f2','N_intro_plus1_f3','p1_f1[2]','p1_f1[3]','p2_f1[3]','p1_f2[2]','p1_f2[3]','p2_f2[3]','p1_f3[2]','p1_f3[3]','p2_f3[3]','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 2){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'box_bel' = box_bel,
                                   'from_56' = from_56),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('beta_w_f1','beta_w_f2','beta_w_f3','beta_b_f1','beta_b_f2','beta_b_f3','nu_0','mu','N_intro_plus1_f1','N_intro_plus1_f2','N_intro_plus1_f3','p1_f1[2]','p1_f1[3]','p2_f1[3]','p1_f2[2]','p1_f2[3]','p2_f2[3]','p1_f3[2]','p1_f3[3]','p2_f3[3]','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 3){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'from_56' = from_56,
                                   'last_expo_class' = last_expo_class),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('tau','alpha_d_Pe','alpha_d_Po','alpha_d_T','alpha_d_M','alpha_d_ST','beta_f_f1','beta_f_f2','beta_f_f3','nu_0','mu','N_intro_plus1_f1','N_intro_plus1_f2','N_intro_plus1_f3','p1_f1[2]','p1_f1[3]','p2_f1[3]','p1_f2[2]','p1_f2[3]','p2_f2[3]','p1_f3[2]','p1_f3[3]','p2_f3[3]','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 4){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'box_bel' = box_bel,
                                   'from_56' = from_56,
                                   'last_expo_class' = last_expo_class),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('tau','alpha_d_Pe','alpha_d_Po','alpha_d_T','alpha_d_M','alpha_d_ST','beta_w_f1','beta_w_f2','beta_w_f3','beta_b_f1','beta_b_f2','beta_b_f3','nu_0','mu','N_intro_plus1_f1','N_intro_plus1_f2','N_intro_plus1_f3','p1_f1[2]','p1_f1[3]','p2_f1[3]','p1_f2[2]','p1_f2[3]','p2_f2[3]','p1_f3[2]','p1_f3[3]','p2_f3[3]','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 5){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'from_56' = from_56,
                                   'last_expo_class' = last_expo_class),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('tau','alpha_c_Pe','alpha_c_Po','alpha_c_T','alpha_c_M','alpha_c_ST','beta_f_f1','beta_f_f2','beta_f_f3','nu_0','mu','N_intro_plus1_f1','N_intro_plus1_f2','N_intro_plus1_f3','p1_f1[2]','p1_f1[3]','p2_f1[3]','p1_f2[2]','p1_f2[3]','p2_f2[3]','p1_f3[2]','p1_f3[3]','p2_f3[3]','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 6){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'box_bel' = box_bel,
                                   'from_56' = from_56,
                                   'last_expo_class' = last_expo_class),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('tau','alpha_c_Pe','alpha_c_Po','alpha_c_T','alpha_c_M','alpha_c_ST','beta_w_f1','beta_w_f2','beta_w_f3','beta_b_f1','beta_b_f2','beta_b_f3','nu_0','mu','N_intro_plus1_f1','N_intro_plus1_f2','N_intro_plus1_f3','p1_f1[2]','p1_f1[3]','p2_f1[3]','p1_f2[2]','p1_f2[3]','p2_f2[3]','p1_f3[2]','p1_f3[3]','p2_f3[3]','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 7){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'from_56' = from_56,
                                   'last_expo_class' = last_expo_class),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('tau','alpha_d_Pe','alpha_d_Po','alpha_d_T','alpha_d_M','alpha_d_ST','beta_0_f1','beta_0_f2','beta_0_f3','nu_0','mu','N_intro_plus1_f1','N_intro_plus1_f2','N_intro_plus1_f3','p1_f1[2]','p1_f1[3]','p2_f1[3]','p1_f2[2]','p1_f2[3]','p2_f2[3]','p1_f3[2]','p1_f3[3]','p2_f3[3]','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 8){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'from_56' = from_56,
                                   'last_expo_class' = last_expo_class),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('tau','alpha_c_Pe','alpha_c_Po','alpha_c_T','alpha_c_M','alpha_c_ST','beta_0_f1','beta_0_f2','beta_0_f3','nu_0','mu','N_intro_plus1_f1','N_intro_plus1_f2','N_intro_plus1_f3','p1_f1[2]','p1_f1[3]','p2_f1[3]','p1_f2[2]','p1_f2[3]','p2_f2[3]','p1_f3[2]','p1_f3[3]','p2_f3[3]','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 10){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'from_56' = from_56),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('beta_0_f1','beta_0_f2','beta_0_f3','nu_0','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 11){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'from_56' = from_56),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('beta_f_f1','beta_f_f2','beta_f_f3','nu_0','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 12){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'box_bel' = box_bel,
                                   'from_56' = from_56),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('beta_w_f1','beta_w_f2','beta_w_f3','beta_b_f1','beta_b_f2','beta_b_f3','nu_0','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 13){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'from_56' = from_56,
                                   'last_expo_class' = last_expo_class),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('tau','alpha_d_Pe','alpha_d_Po','alpha_d_T','alpha_d_M','alpha_d_ST','beta_f_f1','beta_f_f2','beta_f_f3','nu_0','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 14){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'box_bel' = box_bel,
                                   'from_56' = from_56,
                                   'last_expo_class' = last_expo_class),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('tau','alpha_d_Pe','alpha_d_Po','alpha_d_T','alpha_d_M','alpha_d_ST','beta_w_f1','beta_w_f2','beta_w_f3','beta_b_f1','beta_b_f2','beta_b_f3','nu_0','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 15){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'from_56' = from_56,
                                   'last_expo_class' = last_expo_class),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('tau','alpha_c_Pe','alpha_c_Po','alpha_c_T','alpha_c_M','alpha_c_ST','beta_f_f1','beta_f_f2','beta_f_f3','nu_0','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 16){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'box_bel' = box_bel,
                                   'from_56' = from_56,
                                   'last_expo_class' = last_expo_class),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('tau','alpha_c_Pe','alpha_c_Po','alpha_c_T','alpha_c_M','alpha_c_ST','beta_w_f1','beta_w_f2','beta_w_f3','beta_b_f1','beta_b_f2','beta_b_f3','nu_0','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 17){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'from_56' = from_56,
                                   'last_expo_class' = last_expo_class),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('tau','alpha_d_Pe','alpha_d_Po','alpha_d_T','alpha_d_M','alpha_d_ST','beta_0_f1','beta_0_f2','beta_0_f3','nu_0','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == 18){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'from_56' = from_56,
                                   'last_expo_class' = last_expo_class),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('tau','alpha_c_Pe','alpha_c_Po','alpha_c_T','alpha_c_M','alpha_c_ST','beta_0_f1','beta_0_f2','beta_0_f3','nu_0','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
  
  if(numb_mod == "5withbeta0"){
    jags <- jags.model(paste0('~/BUGS/Model_', numb_mod, '_AllFarms.bug'),
                       data = list('portage_obs' = portage_obs,
                                   'N_calves' = N_calves,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel,
                                   'from_56' = from_56,
                                   'last_expo_class' = last_expo_class),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('tau','alpha_c_Pe','alpha_c_Po','alpha_c_T','alpha_c_M','alpha_c_ST','beta_0_f1','beta_0_f2','beta_0_f3','beta_f_f1','beta_f_f2','beta_f_f3','nu_0','mu','N_intro_plus1_f1','N_intro_plus1_f2','N_intro_plus1_f3','p1_f1[2]','p1_f1[3]','p2_f1[3]','p1_f2[2]','p1_f2[3]','p2_f2[3]','p1_f3[2]','p1_f3[3]','p2_f3[3]','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_", numb_mod, "_AllFarms.rdata"))
  }
}

model_perFarm = function(numb_mod, numb_farm, saving_path = "~/"){
  
  load.module("dic")
  
  n_chains_all = 3
  n_adapt_all = 4000
  length_chain_all = 10000
  thin_all = 40
  
  N_calves_mod = 15
  if(numb_farm == 1){
    indiv_farm = indiv_farm_1
    indiv_farm_from_56 = indiv_farm_1_from_56
    init_pos = 2:15
    init_neg = 1
  }else if(numb_farm == 2){
    indiv_farm = indiv_farm_2
    indiv_farm_from_56 = indiv_farm_2_from_56
    init_pos = 6:15
    init_neg = 1:5
  }else if(numb_farm == 3){
    indiv_farm = indiv_farm_3
    indiv_farm_from_56 = indiv_farm_3_from_56
    init_pos = 5:15
    init_neg = 1:4
  }
  
  if(numb_mod == 5){
    jags <- jags.model('~/BUGS/Model_5_perFarm.bug',
                       data = list('portage_obs' = portage_obs[indiv_farm,],
                                   'N_calves' = N_calves_mod,
                                   't_end' = t_end,
                                   'numb_farm' = numb_farm,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel[indiv_farm,],
                                   'from_56' = indiv_farm_from_56,
                                   'last_expo_class' = last_expo_class),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('tau','alpha_c_Pe','alpha_c_Po','alpha_c_T','alpha_c_M','alpha_c_ST','beta_f','nu_0','mu','N_intro_plus1','p1[2]','p1[3]','p2[3]','deviance','pD'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_5_Farm", numb_farm, ".rdata"))
  }
  
  if(numb_mod == 11){
    jags <- jags.model('~/BUGS/Model_11_perFarm.bug',
                       data = list('portage_obs' = portage_obs[indiv_farm,],
                                   'N_calves' = N_calves_mod,
                                   't_end' = t_end,
                                   'init_pos' = init_pos,
                                   'init_neg' = init_neg,
                                   'farm_bel' = farm_bel[indiv_farm,],
                                   'from_56' = indiv_farm_from_56),
                       n.chains = n_chains_all,
                       n.adapt = n_adapt_all)
    samples <- jags.samples(jags,
                            c('beta_f','nu_0','deviance'),
                            length_chain_all,
                            thin = thin_all)
    save(samples, file = paste0(saving_path, "Model_11_Farm", numb_farm, ".rdata"))
  }
}

cl<-makeCluster(7)
clusterEvalQ(cl, library(rjags))
clusterExport(cl,c("portage_obs","N_calves","box_bel","farm_bel","last_expo_class","t_end","from_56","init_pos","init_neg"))
t = as.numeric(Sys.time())
parLapply(cl, X = c(0:8,10:18), fun = model_allFarms)

stopCluster(cl)
as.numeric(Sys.time()) - t
rm(portage_obs)


######################################################################

# 2) Diagnostic and posterior distributions:

load("~/Model_5_AllFarms.rdata")

transformed_chain = to_mcmc.list(mcarray_obj = samples, to_burn = 0)
var_mod = transformed_chain[["output_var_mod"]]
Mch = transformed_chain[["output_Mch"]]
combined_Mch = transformed_chain[["output_combined_Mch"]]

# Chains convergence:
plot(Mch)

# Autocorrelation (acf):
par(mfrow = c(3, 3)) ; for(i in 1:length(var_mod)){acf(combined_Mch[,i], lag.max = 30, ylim = c(-1,1), main = colnames(combined_Mch)[i])} ; par(mfrow = c(1, 1))

# Correlations:
cor(combined_Mch)
ggplot(data = melt(cor(combined_Mch)), aes(x = Var1, y = Var2, fill = value, col = (abs(value) >0.4))) + geom_tile(size = 1.5) + scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) + scale_colour_manual(values = c("TRUE" = "forestgreen", "FALSE" = "white"))

# Effective sample size:
effectiveSize(Mch)

# Posterior distributions:
par(mfrow = c(3, 3)) ; for(i in 1:length(var_mod)){hist(combined_Mch[,i], nclass=20, main = colnames(combined_Mch)[i], col = "skyblue")} ; par(mfrow = c(1, 1))
summary(Mch)
HPDinterval(as.mcmc(combined_Mch))
table(combined_Mch[,"N_intro_plus1_f1"])*100/length(combined_Mch[,"N_intro_plus1_f1"])
table(combined_Mch[,"N_intro_plus1_f2"])*100/length(combined_Mch[,"N_intro_plus1_f2"])
table(combined_Mch[,"N_intro_plus1_f3"])*100/length(combined_Mch[,"N_intro_plus1_f3"])

# Figure in Supplementary:
#####
distrib_in_sm = combined_Mch[,which(! colnames(combined_Mch) %in% c("p1_f2[2]", "p1_f2[3]", "p2_f2[3]", "p1_f3[2]"))]
colnames(distrib_in_sm)[colnames(distrib_in_sm) == "p1_f1[2]"] = "D_fA(N_fA = 1)"
colnames(distrib_in_sm)[colnames(distrib_in_sm) == "p1_f1[3]"] = "D1_fA(N_fA = 2)"
colnames(distrib_in_sm)[colnames(distrib_in_sm) == "p2_f1[3]"] = "D2_fA(N_fA = 2)"
colnames(distrib_in_sm)[colnames(distrib_in_sm) == "p1_f3[3]"] = "D1_fC(N_fC = 2)"
colnames(distrib_in_sm)[colnames(distrib_in_sm) == "p2_f3[3]"] = "D2_fC(N_fC = 2)"
colnames(distrib_in_sm)[colnames(distrib_in_sm) == "N_intro_plus1_f1"] = "N_fA"
distrib_in_sm[,"N_fA"] = distrib_in_sm[,"N_fA"] -1 # as we estimated the number of contaminations +1
colnames(distrib_in_sm)[colnames(distrib_in_sm) == "N_intro_plus1_f2"] = "N_fB"
distrib_in_sm[,"N_fB"] = distrib_in_sm[,"N_fB"] -1 # as we estimated the number of contaminations +1
colnames(distrib_in_sm)[colnames(distrib_in_sm) == "N_intro_plus1_f3"] = "N_fC"
distrib_in_sm[,"N_fC"] = distrib_in_sm[,"N_fC"] -1 # as we estimated the number of contaminations +1
colnames(distrib_in_sm) = gsub("f1", "fA", colnames(distrib_in_sm))
colnames(distrib_in_sm) = gsub("f2", "fB", colnames(distrib_in_sm))
colnames(distrib_in_sm) = gsub("f3", "fC", colnames(distrib_in_sm))
par(mfrow = c(4, 5), mai = c(0.3, 0.3, 0.3, 0.2)) ; for(i in 1:ncol(distrib_in_sm)){hist(distrib_in_sm[,i], nclass=20, main = colnames(distrib_in_sm)[i], col = "skyblue", xlab="", ylab="")} ; par(mfrow = c(1, 1))
#####

# Gelman-Rubin statistic:
transformed_chain = to_mcmc.list(mcarray_obj = samples, to_burn = 0, which_par = 1)
var_mod = transformed_chain[["output_var_mod"]]
Mch = transformed_chain[["output_Mch"]]
gelman.diag(Mch)

######################################################################

# 3) DIC:

# Function calculating DIC from jags.samples output (mcarray):
cal_dic = function(samples){
  return(summary(samples$deviance, mean)$stat + summary(samples$pD, mean)$stat)
}

DIC_all_farms_all_mod = function(which_mod = c(0:8, 10:18), path){
  DIC_mod = rep(NA, length(which_mod))
  names(DIC_mod) = as.character(which_mod)
  
  load(paste0(path, which_mod[1], "_AllFarms.rdata"))
  DIC_mod_0 = cal_dic(samples)
  
  for (n_mod in which_mod){
    load(paste0(path, n_mod, "_AllFarms.rdata"))
    DIC_mod[as.character(n_mod)] = cal_dic(samples)
  }
  
  DIC_mod
}
DIC_all_farms_all_mod(path = "~/Model_")
