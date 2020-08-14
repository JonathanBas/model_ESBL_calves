# Simulations of models

# Needs data and functions from data_preparation.R

library(ggplot2)

simul_mod_perFarm = function(mod, numb_farm, axis.title = T, nsimu=100, path="~/Model_", plot_all_simu=F, burnin = 0, pred_val = "Median", pred_int = 95, last_expo_class = last_expo_class, deter_param = c(), for_save = F){
  
  load(paste0(path, mod, "_Farm", numb_farm, ".rdata"))
  
  transformed_chain = to_mcmc.list(mcarray_obj = samples, to_burn = burnin)
  Mch = transformed_chain[["output_Mch"]]
  combined_Mch = transformed_chain[["output_combined_Mch"]]
  
  # Population for this farm:
  
  if(numb_farm == 1){
    indiv_farm = indiv_farm_1
    indiv_farm_from_56 = indiv_farm_1_from_56
    init_pos = 1:15
    init_neg = NULL
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
  N_calves = 15
  farm_bel = farm_bel[indiv_farm,]
  box_bel = box_bel[indiv_farm,]
  from_56_general = from_56
  from_56 = indiv_farm_from_56
  
  # Simulations:
  portage_all_simu = array(data = NA, dim = c(N_calves, t_end, nsimu))
  
  if(mod == 5){
    
    time_iteration_2to55 = function(indiv){
      beta_tot = 1-exp(-(p_intro + beta_f * sum(portage_simu[,t-1] * farm_bel[, farm_bel[indiv,4]]) / sum(farm_bel[, farm_bel[indiv,4]])))
      nu_tot = 1-exp(-(nu_0 * alpha_c_Pe^(max(0, 1 - last_expo_class[indiv, t, 1]*tau)) * alpha_c_Po^(max(0, 1 - last_expo_class[indiv, t, 2]*tau)) * alpha_c_T^(max(0, 1 - last_expo_class[indiv, t, 3]*tau)) * alpha_c_M^(max(0, 1 - last_expo_class[indiv, t, 4]*tau)) * alpha_c_ST^(max(0, 1 - last_expo_class[indiv, t, 5]*tau))))
      
      prob_equal_1 = ifelse(portage_simu[indiv,t-1] == 0, beta_tot, 1-nu_tot)
      if(is.na(prob_equal_1)){print(paste0("ERROR: prob_equal_1 is NA for t=", t, " and indiv=", indiv))}
      
      rbinom(n = 1, size = 1, prob = min(0.999, max(0.001, prob_equal_1)))
    }
    
    time_iteration_56toEnd = function(indiv){
      beta_tot = 1-exp(-(p_intro + beta_f * sum(portage_simu[from_56,t-1] * farm_bel[from_56, farm_bel[indiv,4]]) / sum(farm_bel[from_56, farm_bel[indiv,4]])))
      nu_tot = 1-exp(-(nu_0 * alpha_c_Pe^(max(0, 1 - last_expo_class[indiv, t, 1]*tau)) * alpha_c_Po^(max(0, 1 - last_expo_class[indiv, t, 2]*tau)) * alpha_c_T^(max(0, 1 - last_expo_class[indiv, t, 3]*tau)) * alpha_c_M^(max(0, 1 - last_expo_class[indiv, t, 4]*tau)) * alpha_c_ST^(max(0, 1 - last_expo_class[indiv, t, 5]*tau))))
      
      prob_equal_1 = ifelse(portage_simu[indiv,t-1] == 0, beta_tot, 1-nu_tot)
      if(is.na(prob_equal_1)){print(paste0("ERROR: prob_equal_1 is NA for t=", t, " and indiv=", indiv))}
      
      rbinom(n = 1, size = 1, prob = min(0.999, max(0.001, prob_equal_1)))
    }
    
    for (simu in 1:nsimu){
      cat("\r", paste0("Repetition ", simu, "/", nsimu))
      portage_simu = matrix(NA, N_calves, t_end)
      
      # Draw set of parameters:
      
      set_param = base::sample(x = 1:nrow(combined_Mch), size = 1)
      
      list_param = c("nu_0", "beta_f", "mu", "N_intro_plus1", "p1[2]", "p1[3]", "p2[3]", "tau", "alpha_c_Pe", "alpha_c_Po", "alpha_c_T", "alpha_c_M", "alpha_c_ST")
      val_param = rep(NA, length(list_param))
      names(val_param) = list_param
      
      for(param_i in list_param){
        if(param_i %in% names(deter_param)){
          val_param[param_i] = deter_param[param_i]
        }else{
          val_param[param_i] = combined_Mch[set_param, param_i]
        }
      }
      
      nu_0 = val_param["nu_0"]
      beta_f = val_param["beta_f"]
      
      mu = val_param["mu"]
      N_intro_plus1 = val_param["N_intro_plus1"]
      p1 = c(0, val_param["p1[2]"], val_param["p1[3]"])
      p2 = c(0, 0, val_param["p2[3]"])
      
      tau = val_param["tau"]
      alpha_c_Pe = val_param["alpha_c_Pe"]
      alpha_c_Po = val_param["alpha_c_Po"]
      alpha_c_T = val_param["alpha_c_T"]
      alpha_c_M = val_param["alpha_c_M"]
      alpha_c_ST = val_param["alpha_c_ST"]
      
      portage_simu[init_pos,1] = 1
      portage_simu[init_neg,1] = 0
      
      for (t in 2:55){
        p_intro = mu * (p1[N_intro_plus1] == t) + mu * (p2[N_intro_plus1] == t)
        portage_simu[1:N_calves, t] = sapply(X = 1:N_calves, FUN = time_iteration_2to55)
      }
      for (t in 56:t_end){
        p_intro = mu * (p1[N_intro_plus1] == t) + mu * (p2[N_intro_plus1] == t)
        portage_simu[from_56, t] = sapply(X = from_56, FUN = time_iteration_56toEnd)
      }
      
      portage_all_simu[,,simu] = portage_simu
    }
    cat("\n")
  }
  
  if(mod == 11){
    
    time_iteration_2to55 = function(indiv){
      beta_tot = 1-exp(-(beta_f * sum(portage_simu[,t-1] * farm_bel[, farm_bel[indiv,4]]) / sum(farm_bel[, farm_bel[indiv,4]])))
      nu_tot = 1-exp(-nu_0)
      
      prob_equal_1 = ifelse(portage_simu[indiv,t-1] == 0, beta_tot, 1-nu_tot)
      if(is.na(prob_equal_1)){print(paste0("ERROR: prob_equal_1 is NA for t=", t, " and indiv=", indiv))}
      
      rbinom(n = 1, size = 1, prob = min(0.999, max(0.001, prob_equal_1)))
    }
    
    time_iteration_56toEnd = function(indiv){
      beta_tot = 1-exp(-(beta_f * sum(portage_simu[from_56,t-1] * farm_bel[from_56, farm_bel[indiv,4]]) / sum(farm_bel[from_56, farm_bel[indiv,4]])))
      nu_tot = 1-exp(-nu_0)
      
      prob_equal_1 = ifelse(portage_simu[indiv,t-1] == 0, beta_tot, 1-nu_tot)
      if(is.na(prob_equal_1)){print(paste0("ERROR: prob_equal_1 is NA for t=", t, " and indiv=", indiv))}
      
      rbinom(n = 1, size = 1, prob = min(0.999, max(0.001, prob_equal_1)))
    }
    
    for (simu in 1:nsimu){
      cat("\r", paste0("Repetition ", simu, "/", nsimu))
      portage_simu = matrix(NA, N_calves, t_end)
      
      # Draw set of parameters:
      
      set_param = base::sample(x = 1:nrow(combined_Mch), size = 1)
      
      list_param = c("nu_0", "beta_f")
      val_param = rep(NA, length(list_param))
      names(val_param) = list_param
      
      for(param_i in list_param){
        if(param_i %in% names(deter_param)){
          val_param[param_i] = deter_param[param_i]
        }else{
          val_param[param_i] = combined_Mch[set_param, param_i]
        }
      }
      
      nu_0 = val_param["nu_0"]
      beta_f = val_param["beta_f"]
      
      if("eta" %in% names(deter_param)){
        portage_simu[,1] = rbinom(n = N_calves, size = 1, prob = deter_param["eta"])
      }else{
        portage_simu[init_pos,1] = 1
        portage_simu[init_neg,1] = 0
      }
      
      for (t in 2:55){
        portage_simu[1:N_calves, t] = sapply(X = 1:N_calves, FUN = time_iteration_2to55)
      }
      for (t in 56:t_end){
        portage_simu[from_56, t] = sapply(X = from_56, FUN = time_iteration_56toEnd)
      }
      
      portage_all_simu[,,simu] = portage_simu
    }
    cat("\n")
  }
  
  # Plot:
  
  summa_prev = as.data.frame(matrix(NA, t_end, 10))
  colnames(summa_prev) = c("day", "val", "inf", "sup", "obs", "Pe", "Po", "T", "M", "ST")
  summa_prev$day = 1:t_end
  
  for(t in 1:55){
    iter = which(summa_prev$day == t)
    prev_t = colSums(portage_all_simu[,t,])/N_calves
    summa_prev[iter,"val"] = summary(100 * prev_t)[pred_val]
    summa_prev[iter,c("inf", "sup")] = 100 * quantile(prev_t, probs=c((100-pred_int)/200, 1-(100-pred_int)/200))
    if(t %in% sampling_days){
      summa_prev$obs[iter] = 100 * sum(portage_real[indiv_farm,t+3] != 0) /N_calves
    }
    
    # Pour afficher les jours d'exposition aux antibiotiques en haut du graphe:
    for (class_i in c("Pe", "Po", "T", "M", "ST")){
      if(expo[1+15*(numb_farm-1), t, class_i]){
        summa_prev[iter, class_i] = 1
      }
    }
  }
  for(t in 56:t_end){
    iter = which(summa_prev$day == t)
    prev_t = colSums(portage_all_simu[from_56,t,])/length(from_56)
    summa_prev[iter,"val"] = summary(100 * prev_t)[pred_val]
    summa_prev[iter,c("inf", "sup")] = 100 * quantile(prev_t, probs=c((100-pred_int)/200, 1-(100-pred_int)/200))
    if(t %in% sampling_days){
      summa_prev$obs[iter] = 100 * sum(portage_real[intersect(indiv_farm, from_56_general),t+3] != 0) /length(from_56)
    }
    
    # Pour afficher les jours d'exposition aux antibiotiques en haut du graphe:
    
    for (class_i in c("Pe", "Po", "T", "M", "ST")){
      if(expo[1+15*(numb_farm-1), t, class_i]){
        summa_prev[iter, class_i] = 1
      }
    }
  }
  
  library(ggplot2)
  
  p = ggplot()
  
  if(plot_all_simu == F){
    p = p + geom_ribbon(data = summa_prev, aes(x=day, ymin=inf, ymax=sup), alpha=0.5, fill="#E69F00")
  }else{
    prev_all_simu = melt(colMeans(portage_all_simu[,,], dims=1, na.rm=T)*100)
    colnames(prev_all_simu) = c("time", "simu", "val")
    
    p = p + geom_line(data = prev_all_simu, aes(x = time, y = val, group = simu), col = "grey")
    p = p + geom_line(data = summa_prev, aes(y=inf, x=day), col="indianred3")
    p = p + geom_line(data = summa_prev, aes(y=sup, x=day), col="indianred3")
  }
  
  p = p + geom_line(data = summa_prev, aes(y=val, x=day), col="black", size=1)
  if(mod %in% c(2, 4, 6, 12, 14, 16)){
    p = p + geom_vline(xintercept = 49, linetype = "dotted")
  }
  p = p + geom_point(data = summa_prev, aes(y = obs, x=day), size = 3, shape = 23,  fill = "red", color = "black")
  p = p + geom_point(data = summa_prev[! is.na(summa_prev$obs),], aes(y = -10, x=day), size = 3, shape = 4, color = "black")
  p = p + geom_rect(data = summa_prev, aes(xmin=day*Pe-0.5, xmax=day*Pe+0.5, ymin=105, ymax=107, fill="Penicillin"))
  p = p + geom_rect(data = summa_prev, aes(xmin=day*Po-0.5, xmax=day*Po+0.5, ymin=107, ymax=109, fill="Colistin"))
  p = p + geom_rect(data = summa_prev, aes(xmin=day*T-0.5, xmax=day*T+0.5, ymin=109, ymax=111, fill="Tetracycline"))
  p = p + geom_rect(data = summa_prev, aes(xmin=day*M-0.5, xmax=day*M+0.5, ymin=111, ymax=113, fill="Macrolide"))
  p = p + geom_rect(data = summa_prev, aes(xmin=day*ST-0.5, xmax=day*ST+0.5, ymin=113, ymax=115, fill="Sulfonamide-Trimethoprim"))
  p = p + scale_y_continuous(breaks = seq(0,100,25), limits = c(-12,115))
  p = p + scale_fill_manual(name = "Antibiotic class used:",
                            values = c("Penicillin" = "darkviolet", "Colistin" = "darkgreen", "Tetracycline" = "deepskyblue3", "Macrolide" = "darkorange", "Sulfonamide-Trimethoprim" = "brown3"))
  
  if(axis.title){
    p = p + xlab("Time (days)") + ylab("Prevalence of ESBL-EC in calves (%)")
    p = p + theme_bw()
    if(! for_save){
      p = p + ggtitle(paste("Observed and predicted carriage in Farm", c("A","B","C")[numb_farm]))
      # p = p + theme(legend.position="bottom")
    }else{
      p = p + theme(legend.position="none")
    }
  }else{
    p = p + guides(fill=FALSE)
    p = p + theme_void()
  }
  
  p
  
  return(list(portage_all_simu, p, summa_prev))
}

simul_mod_allFarms = function(mod, axis.title = T, nsimu=100, path="~/Model_", plot_all_simu=F, burnin = 0, pred_val = "Median", pred_int = 95, last_expo_class = last_expo_class, deter_param = c(), save_fig = F){
  
  load(paste0(path, mod, "_AllFarms.rdata"))
  
  transformed_chain = to_mcmc.list(mcarray_obj = samples, to_burn = burnin)
  Mch = transformed_chain[["output_Mch"]]
  combined_Mch = transformed_chain[["output_combined_Mch"]]
  
  # Simulations:
  portage_all_simu = array(data = NA, dim = c(N_calves, t_end, nsimu))

  if(mod == 5){
    
    time_iteration_2to55 = function(indiv){
      beta_tot = 1-exp(-(p_intro[farm_bel[indiv,4]] + beta_f[farm_bel[indiv,4]] * sum(portage_simu[,t-1] * farm_bel[, farm_bel[indiv,4]]) / sum(farm_bel[, farm_bel[indiv,4]])))
      nu_tot = 1-exp(-(nu_0 * alpha_c_Pe^(max(0, 1 - last_expo_class[indiv, t, 1]*tau)) * alpha_c_Po^(max(0, 1 - last_expo_class[indiv, t, 2]*tau)) * alpha_c_T^(max(0, 1 - last_expo_class[indiv, t, 3]*tau)) * alpha_c_M^(max(0, 1 - last_expo_class[indiv, t, 4]*tau)) * alpha_c_ST^(max(0, 1 - last_expo_class[indiv, t, 5]*tau))))
      
      prob_equal_1 = ifelse(portage_simu[indiv,t-1] == 0, beta_tot, 1-nu_tot)
      if(is.na(prob_equal_1)){print(paste0("ERROR: prob_equal_1 is NA for t=", t, " and indiv=", indiv))}
      
      rbinom(n = 1, size = 1, prob = prob_equal_1)
    }
    
    time_iteration_56toEnd = function(indiv){
      beta_tot = 1-exp(-(p_intro[farm_bel[indiv,4]] + beta_f[farm_bel[indiv,4]] * sum(portage_simu[from_56,t-1] * farm_bel[from_56, farm_bel[indiv,4]]) / sum(farm_bel[from_56, farm_bel[indiv,4]])))
      nu_tot = 1-exp(-(nu_0 * alpha_c_Pe^(max(0, 1 - last_expo_class[indiv, t, 1]*tau)) * alpha_c_Po^(max(0, 1 - last_expo_class[indiv, t, 2]*tau)) * alpha_c_T^(max(0, 1 - last_expo_class[indiv, t, 3]*tau)) * alpha_c_M^(max(0, 1 - last_expo_class[indiv, t, 4]*tau)) * alpha_c_ST^(max(0, 1 - last_expo_class[indiv, t, 5]*tau))))
      
      prob_equal_1 = ifelse(portage_simu[indiv,t-1] == 0, beta_tot, 1-nu_tot)
      if(is.na(prob_equal_1)){print(paste0("ERROR: prob_equal_1 is NA for t=", t, " and indiv=", indiv))}
      
      rbinom(n = 1, size = 1, prob = prob_equal_1)
    }
    
    for (simu in 1:nsimu){
      cat("\r", paste0("Repetition ", simu, "/", nsimu))
      portage_simu = matrix(NA, N_calves, t_end)
      
      # Draw set of parameters:
      
      set_param = base::sample(x = 1:nrow(combined_Mch), size = 1)
      
      list_param = c('tau','alpha_c_Pe','alpha_c_Po','alpha_c_T','alpha_c_M','alpha_c_ST','beta_f_f1','beta_f_f2','beta_f_f3','nu_0','mu','N_intro_plus1_f1','N_intro_plus1_f2','N_intro_plus1_f3','p1_f1[2]','p1_f1[3]','p2_f1[3]','p1_f2[2]','p1_f2[3]','p2_f2[3]','p1_f3[2]','p1_f3[3]','p2_f3[3]')
      val_param = rep(NA, length(list_param))
      names(val_param) = list_param
      
      for(param_i in list_param){
        if(param_i %in% names(deter_param)){
          val_param[param_i] = deter_param[param_i]
        }else{
          val_param[param_i] = combined_Mch[set_param, param_i]
        }
      }
      
      tau = val_param["tau"]
      alpha_c_Pe = val_param["alpha_c_Pe"]
      alpha_c_Po = val_param["alpha_c_Po"]
      alpha_c_T = val_param["alpha_c_T"]
      alpha_c_M = val_param["alpha_c_M"]
      alpha_c_ST = val_param["alpha_c_ST"]
      
      beta_f = c(val_param["beta_f_f1"], val_param["beta_f_f2"], val_param["beta_f_f3"])
      
      nu_0 = val_param["nu_0"]
      mu = val_param["mu"]
      
      N_intro_plus1_f1 = val_param["N_intro_plus1_f1"]
      N_intro_plus1_f2 = val_param["N_intro_plus1_f2"]
      N_intro_plus1_f3 = val_param["N_intro_plus1_f3"]
      
      p1_f1 = c(0, val_param["p1_f1[2]"], val_param["p1_f1[3]"])
      p2_f1 = c(0, 0, val_param["p2_f1[3]"])
      p1_f2 = c(0, val_param["p1_f2[2]"], val_param["p1_f2[3]"])
      p2_f2 = c(0, 0, val_param["p2_f2[3]"])
      p1_f3 = c(0, val_param["p1_f3[2]"], val_param["p1_f3[3]"])
      p2_f3 = c(0, 0, val_param["p2_f3[3]"])
      
      if("eta" %in% names(deter_param)){
        portage_simu[,1] = rbinom(n = N_calves, size = 1, prob = deter_param["eta"])
      }else{
        portage_simu[init_pos,1] = 1
        portage_simu[init_neg,1] = 0
      }
      
      for (t in 2:55){
        p_intro = c(mu * (p1_f1[N_intro_plus1_f1] == t) + mu * (p2_f1[N_intro_plus1_f1] == t),
                    mu * (p1_f2[N_intro_plus1_f2] == t) + mu * (p2_f2[N_intro_plus1_f2] == t),
                    mu * (p1_f3[N_intro_plus1_f3] == t) + mu * (p2_f3[N_intro_plus1_f3] == t))
        portage_simu[1:N_calves, t] = sapply(X = 1:N_calves, FUN = time_iteration_2to55)
      }
      for (t in 56:t_end){
        p_intro = c(mu * (p1_f1[N_intro_plus1_f1] == t) + mu * (p2_f1[N_intro_plus1_f1] == t),
                    mu * (p1_f2[N_intro_plus1_f2] == t) + mu * (p2_f2[N_intro_plus1_f2] == t),
                    mu * (p1_f3[N_intro_plus1_f3] == t) + mu * (p2_f3[N_intro_plus1_f3] == t))
        portage_simu[from_56, t] = sapply(X = from_56, FUN = time_iteration_56toEnd)
      }
      
      portage_all_simu[,,simu] = portage_simu
    }
    cat("\n")
  }
  
  # Plot:
  
  indiv_farms=list(1:15, 16:30, 31:45)
  summa_prev = as.data.frame(matrix(NA, t_end*3, 11))
  colnames(summa_prev) = c("farm", "day", "val", "inf", "sup", "obs", "Pe", "Po", "T", "M", "ST")
  summa_prev$farm = c(rep(1,t_end), rep(2,t_end), rep(3,t_end))
  summa_prev$day = rep(1:t_end,3)
  for(t in 1:55){
    for (farm in 1:3){
      iter = which((summa_prev$farm == farm) & (summa_prev$day == t))
      prev_t = colSums(portage_all_simu[indiv_farms[[farm]],t,])/length(indiv_farms[[farm]])
      summa_prev[iter,"val"] = summary(100 * prev_t)[pred_val]
      summa_prev[iter,c("inf", "sup")] = 100 * quantile(prev_t, probs=c((100-pred_int)/200, 1-(100-pred_int)/200))
      if(t %in% sampling_days){
        summa_prev$obs[iter] = 100 * sum(portage_real[indiv_farms[[farm]],t+3] != 0) /length(indiv_farms[[farm]])
      }
      
      # To display AB exposure on top:
      
      for (class_i in c("Pe", "Po", "T", "M", "ST")){
        if(expo[1+15*(farm-1), t, class_i]){
          summa_prev[iter, class_i] = 1
        }
      }
    }
  }
  for(t in 56:t_end){
    for (farm in 1:3){
      iter = which((summa_prev$farm == farm) & (summa_prev$day == t))
      prev_t = colSums(portage_all_simu[intersect(from_56,indiv_farms[[farm]]),t,])/length(intersect(from_56,indiv_farms[[farm]]))
      summa_prev[iter,"val"] = summary(100 * prev_t)[pred_val]
      summa_prev[iter,c("inf", "sup")] = 100 * quantile(prev_t, probs=c((100-pred_int)/200, 1-(100-pred_int)/200))
      if(t %in% sampling_days){
        summa_prev$obs[iter] = 100 * sum(portage_real[intersect(from_56,indiv_farms[[farm]]),t+3] != 0) /length(intersect(from_56,indiv_farms[[farm]]))
      }
      
      # To display AB exposure on top:
      
      for (class_i in c("Pe", "Po", "T", "M", "ST")){
        if(expo[1+15*(farm-1), t, class_i]){
          summa_prev[iter, class_i] = 1
        }
      }
    }
  }
  
  library(ggplot2)
  
  if(plot_all_simu == F){
    p = ggplot(data = summa_prev, aes(x=day))
    p = p + geom_ribbon(aes(ymin=inf, ymax=sup), alpha=0.5, fill="#E69F00")
  }else{
    prev_all_simu = rbind(melt(colMeans(portage_all_simu[1:15,,], dims=1, na.rm=T)*100),
                          melt(colMeans(portage_all_simu[16:30,,], dims=1, na.rm=T)*100),
                          melt(colMeans(portage_all_simu[31:45,,], dims=1, na.rm=T)*100))
    colnames(prev_all_simu) = c("time", "simu", "val")
    prev_all_simu$farm = c(rep(1, nrow(prev_all_simu)/3), rep(2, nrow(prev_all_simu)/3), rep(3, nrow(prev_all_simu)/3))
    
    p = ggplot()
    p = p + ylim(0,115)
    p = p + geom_line(data = prev_all_simu, aes(x = time, y = val, group = simu), col = "grey")
    p = p + geom_line(data = summa_prev, aes(y=inf, x=day), col="indianred3")
    p = p + geom_line(data = summa_prev, aes(y=sup, x=day), col="indianred3")
  }
  
  p = p + geom_line(data = summa_prev, aes(y=val, x=day), col="black", size=1)
  if(mod %in% c(2, 4, 6, 12, 14, 16)){
    p = p + geom_vline(xintercept = 49, linetype = "dotted")
  }
  p = p + geom_point(data = summa_prev, aes(y = obs, x=day), size = 3, shape = 23,  fill = "red", color = "black")
  p = p + geom_point(data = summa_prev[! is.na(summa_prev$obs),], aes(y = -10, x=day), size = 3, shape = 4, color = "black")
  p = p + geom_rect(data = summa_prev, aes(xmin=day*Pe-0.5, xmax=day*Pe+0.5, ymin=105, ymax=107, fill="Penicillin"))
  p = p + geom_rect(data = summa_prev, aes(xmin=day*Po-0.5, xmax=day*Po+0.5, ymin=107, ymax=109, fill="Colistin"))
  p = p + geom_rect(data = summa_prev, aes(xmin=day*T-0.5, xmax=day*T+0.5, ymin=109, ymax=111, fill="Tetracycline"))
  p = p + geom_rect(data = summa_prev, aes(xmin=day*M-0.5, xmax=day*M+0.5, ymin=111, ymax=113, fill="Macrolide"))
  p = p + geom_rect(data = summa_prev, aes(xmin=day*ST-0.5, xmax=day*ST+0.5, ymin=113, ymax=115, fill="Sulfonamide-Trimethoprim"))
  p = p + scale_y_continuous(breaks = seq(0,100,25), limits = c(-12,115))
  p = p + scale_fill_manual(name = "Antibiotic class used:",
                            values = c("Penicillin" = "darkviolet", "Colistin" = "darkgreen", "Tetracycline" = "deepskyblue3", "Macrolide" = "darkorange", "Sulfonamide-Trimethoprim" = "brown3"))
  
  if(axis.title){
    p = p + xlab("Time (days)") + ylab("Prevalence of ESBL-EC in calves (%)")
    if(! save_fig){
      p = p + ggtitle("Observed and predicted carriage in farms")
    }
    p = p + theme_bw()
    p = p + theme(legend.position="bottom")
  }else{
    p = p + theme_void()
  }
  
  p = p + facet_wrap(~ farm, ncol=3, labeller = labeller(farm = c("1" = "Farm A", "2" = "Farm B", "3" = "Farm C")))
  p
  
  if(save_fig){
    ggsave(p, file = paste0("~/Figure simulations model", mod," all farms.pdf"), width=8, height=5)
  }
  
  return(list(portage_all_simu, p, summa_prev))
}

fig5 = simul_mod_allFarms(mod = 5, nsimu=1000, last_expo_class = last_expo_class, save_fig = T)

library(ggpubr)
p1 = simul_mod_perFarm(5, 1, nsimu=1000, for_save=T, last_expo_class=last_expo_class)[[2]]
p2 = simul_mod_perFarm(5, 2, nsimu=1000, for_save=T, last_expo_class=last_expo_class)[[2]]
p3 = simul_mod_perFarm(5, 3, nsimu=1000, for_save=T, last_expo_class=last_expo_class)[[2]]
com_leg = get_legend(simul_mod_perFarm(5, 1, nsimu=3, last_expo_class=last_expo_class)[[2]])

p_main = ggarrange(p1, p2, p3, com_leg, nrow = 2, ncol = 2, labels = c("Farm A", "Farm B", "Farm C"), label.x = 0.65, label.y = 0.98)
p_main
ggsave(p_main, file = "~/Figure individual fits.pdf", width=7, height=6.3)

