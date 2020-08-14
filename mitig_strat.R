# Simulations of mitigation strategies

# Needs data and functions from data_preparation.R and simul_mod.R

library(ggpubr)

# 1) Simulations of mitigation strategies

rm(last_expo_class)
val_col_init = c(0,3,6)
val_col_mid = c(0, 10)
val_arr = c(0,0.34,0.68)

htmap = expand.grid(col_init = val_col_init, col_mid = val_col_mid, arr = val_arr)
htmap = dplyr::filter(htmap, col_init == 6 | arr == 0.68)

cont_strat = function(colist_init, colist_mid, arriv){
  
  # Create new "last_expo_class" according to AB exposure defined for the simulation:
  
  #####
  expo[1:15,,] = 0
  expo[1:15, if(colist_init == 0){NULL}else{1:colist_init}, "Po"] = 1
  expo[1:15, if(colist_mid == 0){NULL}else{81:(80+colist_mid)}, "Po"] = 1
  
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
  #####
  
  result_simul = simul_mod_allFarms(mod = 5, nsimu = 1000, burnin = 0, deter_param = c("eta" = arriv, "N_intro_plus1_f1" = 1), path = "~/Model_", last_expo_class = last_expo_class)[[3]]
  
  result_simul = dplyr::filter(result_simul, farm == 1)
  
  list(mean_median_prev = mean(result_simul$val),
       mean_sup_prev = mean(result_simul$sup),
       mean_inf_prev = mean(result_simul$inf),
       final_median_prev = result_simul$val[which.max(result_simul$day)],
       final_sup_prev = result_simul$sup[which.max(result_simul$day)],
       final_inf_prev = result_simul$inf[which.max(result_simul$day)])
}

cl<-makeCluster(7)
clusterEvalQ(cl, library(rjags, ggplot2))
clusterExport(cl, c("to_mcmc.list", "expo","expo_tot","portage_real","farm_bel","box_bel","sampling_days","N_calves","t_end","from_56","simul_mod_allFarms"))
t = as.numeric(Sys.time())
values_table = clusterMap(cl, colist_init = htmap[,1], colist_mid = htmap[,2], arriv = htmap[,3], fun = cont_strat)
stopCluster(cl)
as.numeric(Sys.time()) - t

htmap = cbind(htmap, mean_median_prev = as.numeric(sapply(values_table, "[[", "mean_median_prev")), mean_sup_prev = as.numeric(sapply(values_table, "[[", "mean_sup_prev")), mean_inf_prev = as.numeric(sapply(values_table, "[[", "mean_inf_prev")), final_median_prev = as.numeric(sapply(values_table, "[[", "final_median_prev")), final_sup_prev = as.numeric(sapply(values_table, "[[", "final_sup_prev")), final_inf_prev = as.numeric(sapply(values_table, "[[", "final_inf_prev")))

######################################################################

# 2) Plot of results:

htmap$col_mid = as.factor(htmap$col_mid)
htmap$arr = factor(htmap$arr, labels = c("0" = "0%", "0.34" = "34%", "0.68" = "68% (baseline)"))
htmap$col_init = factor(htmap$col_init, labels = c("0" = "No exposure", "3" = "3 days", "6" = "6 days (baseline)"))

mitig_strat_plot = function(indic, fixed_var, plot_lines=F){
  if(indic == "mean_median_prev" & fixed_var == "col_init"){
    p = ggplot(data = dplyr::filter(htmap, col_init == names(which.max(table(htmap$col_init)))), aes(x = arr, y = mean_median_prev, ymin = mean_inf_prev, ymax = mean_sup_prev, col = col_mid, shape = col_mid))
    p = p + xlab("Prevalence on arrival") + ylab("Mean prevalence (%)")
  }
  if(indic == "mean_median_prev" & fixed_var == "arr"){
    p = ggplot(data = dplyr::filter(htmap, arr == names(which.max(table(htmap$arr)))), aes(x = col_init, y = mean_median_prev, ymin = mean_inf_prev, ymax = mean_sup_prev, col = col_mid, shape = col_mid))
    p = p + xlab("Initial exposure to selective antibiotics") + ylab("Mean prevalence (%)")
  }
  if(indic == "final_median_prev" & fixed_var == "col_init"){
    p = ggplot(data = dplyr::filter(htmap, col_init == names(which.max(table(htmap$col_init)))), aes(x = arr, y = final_median_prev, ymin = final_inf_prev, ymax = final_sup_prev, col = col_mid, shape = col_mid))
    p = p + xlab("Prevalence on arrival") + ylab("Prevalence at slaughter age (%)")
  }
  if(indic == "final_median_prev" & fixed_var == "arr"){
    p = ggplot(data = dplyr::filter(htmap, arr == names(which.max(table(htmap$arr)))), aes(x = col_init, y = final_median_prev, ymin = final_inf_prev, ymax = final_sup_prev, col = col_mid, shape = col_mid))
    p = p + xlab("Initial exposure to selective antibiotics") + ylab("Prevalence at slaughter age (%)")
  }
  
  p = p + geom_errorbar(position = position_dodge(0.4), size = 0.5, width = 0.01)
  p = p + geom_point(position = position_dodge(0.4), size = 3, stat="identity")
  
  if(indic == "final_median_prev" & fixed_var == "col_init"){
    p = p + annotate("point", x = names(which.max(table(htmap$arr))), y = 20.4, shape = 8, size = 3, colour = "black")
  }
  if(indic == "final_median_prev" & fixed_var == "arr"){
    p = p + annotate("point", x = names(which.max(table(htmap$col_init))), y = 20.4, shape = 8, size = 3, colour = "black")
  }
  
  p = p + scale_color_manual(values = c("0"="#1b9e77", "10"="#d95f02"), labels = c("0"="No mid-cycle exposure", "10"="10-day mid-cycle exposure"))
  p = p + scale_shape_manual(values = c("0"=16, "10"=15), labels = c("0"="No mid-cycle exposure", "10"="10-day mid-cycle exposure"))
  
  p = p + theme_bw()
  p = p + theme(legend.title = element_blank(), legend.direction = "vertical")
  p
}

p1 = mitig_strat_plot(indic = "mean_median_prev", fixed_var = "col_init")
p2 = mitig_strat_plot(indic = "mean_median_prev", fixed_var = "arr")
p3 = mitig_strat_plot(indic = "final_median_prev", fixed_var = "col_init")
p4 = mitig_strat_plot(indic = "final_median_prev", fixed_var = "arr")

p_all = ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, common.legend = T, legend = "top", labels = c("A", "B", "C", "D"), label.x = 0.9, label.y = 0.95)
p_all
ggsave(p_all, file = "~/Figure mitigation strategies.pdf", width=7, height=6)

