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

cont_strat = function(colist_init, colist_mid, arriv, nsim = 5000){
  
  # Create new "last_expo_class" according to AB exposure defined for the simulation:
  
  #####
  expo[,,] = 0
  expo[, if(colist_init == 0){NULL}else{1:colist_init}, "Po"] = 1
  expo[, if(colist_mid == 0){NULL}else{81:(80+colist_mid)}, "Po"] = 1
  
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
  
  res_sim = simul_mod_allFarms(mod = 5, nsimu = nsim, burnin = 0, deter_param = c("eta" = arriv, "N_intro_plus1_f1" = 1, "N_intro_plus1_f2" = 1, "N_intro_plus1_f3" = 1), path = "~/Model_", last_expo_class = last_expo_class)[[1]]
  
  indiv_farms = list(1:15, 16:30, 31:45)
  sum_tab = data.frame(farm = c(rep(1, nsim), rep(2, nsim), rep(3, nsim)),
                       simu = rep(1:nsim, 3),
                       mean_prev = NA,
                       final_prev = NA)
  
  for(f_i in 1:3){
    for(sim_i in 1:nsim){
      iter = which((sum_tab$farm == f_i) & (sum_tab$simu == sim_i))
      
      res_i = res_sim[indiv_farms[[f_i]],,sim_i]
      prev_i = 100*colMeans(res_i, na.rm=T)
      
      if(any(is.na(prev_i))){print("ERREUR: NAs in prevalence vector")}
      
      sum_tab$mean_prev[iter] = mean(prev_i)
      sum_tab$final_prev[iter] = prev_i[length(prev_i)]
    }
  }
  
  list(min_mean_prev_f1 = as.numeric(quantile(x=sum_tab$mean_prev[sum_tab$farm == 1], probs=0.025)),
       q1_mean_prev_f1 = as.numeric(quantile(x=sum_tab$mean_prev[sum_tab$farm == 1], probs=0.25)),
       median_mean_prev_f1 = median(x=sum_tab$mean_prev[sum_tab$farm == 1]),
       q3_mean_prev_f1 = as.numeric(quantile(x=sum_tab$mean_prev[sum_tab$farm == 1], probs=0.75)),
       max_mean_prev_f1 = as.numeric(quantile(x=sum_tab$mean_prev[sum_tab$farm == 1], probs=0.975)),
       
       min_mean_prev_f2 = as.numeric(quantile(x=sum_tab$mean_prev[sum_tab$farm == 2], probs=0.025)),
       q1_mean_prev_f2 = as.numeric(quantile(x=sum_tab$mean_prev[sum_tab$farm == 2], probs=0.25)),
       median_mean_prev_f2 = median(x=sum_tab$mean_prev[sum_tab$farm == 2]),
       q3_mean_prev_f2 = as.numeric(quantile(x=sum_tab$mean_prev[sum_tab$farm == 2], probs=0.75)),
       max_mean_prev_f2 = as.numeric(quantile(x=sum_tab$mean_prev[sum_tab$farm == 2], probs=0.975)),
       
       min_mean_prev_f3 = as.numeric(quantile(x=sum_tab$mean_prev[sum_tab$farm == 3], probs=0.025)),
       q1_mean_prev_f3 = as.numeric(quantile(x=sum_tab$mean_prev[sum_tab$farm == 3], probs=0.25)),
       median_mean_prev_f3 = median(x=sum_tab$mean_prev[sum_tab$farm == 3]),
       q3_mean_prev_f3 = as.numeric(quantile(x=sum_tab$mean_prev[sum_tab$farm == 3], probs=0.75)),
       max_mean_prev_f3 = as.numeric(quantile(x=sum_tab$mean_prev[sum_tab$farm == 3], probs=0.975)),
       
       min_final_prev_f1 = as.numeric(quantile(x=sum_tab$final_prev[sum_tab$farm == 1], probs=0.025)),
       q1_final_prev_f1 = as.numeric(quantile(x=sum_tab$final_prev[sum_tab$farm == 1], probs=0.25)),
       median_final_prev_f1 = median(x=sum_tab$final_prev[sum_tab$farm == 1]),
       q3_final_prev_f1 = as.numeric(quantile(x=sum_tab$final_prev[sum_tab$farm == 1], probs=0.75)),
       max_final_prev_f1 = as.numeric(quantile(x=sum_tab$final_prev[sum_tab$farm == 1], probs=0.975)),
       
       min_final_prev_f2 = as.numeric(quantile(x=sum_tab$final_prev[sum_tab$farm == 2], probs=0.025)),
       q1_final_prev_f2 = as.numeric(quantile(x=sum_tab$final_prev[sum_tab$farm == 2], probs=0.25)),
       median_final_prev_f2 = median(x=sum_tab$final_prev[sum_tab$farm == 2]),
       q3_final_prev_f2 = as.numeric(quantile(x=sum_tab$final_prev[sum_tab$farm == 2], probs=0.75)),
       max_final_prev_f2 = as.numeric(quantile(x=sum_tab$final_prev[sum_tab$farm == 2], probs=0.975)),
       
       min_final_prev_f3 = as.numeric(quantile(x=sum_tab$final_prev[sum_tab$farm == 3], probs=0.025)),
       q1_final_prev_f3 = as.numeric(quantile(x=sum_tab$final_prev[sum_tab$farm == 3], probs=0.25)),
       median_final_prev_f3 = median(x=sum_tab$final_prev[sum_tab$farm == 3]),
       q3_final_prev_f3 = as.numeric(quantile(x=sum_tab$final_prev[sum_tab$farm == 3], probs=0.75)),
       max_final_prev_f3 = as.numeric(quantile(x=sum_tab$final_prev[sum_tab$farm == 3], probs=0.975)))
}

cl<-makeCluster(7)
clusterEvalQ(cl, library(rjags, ggplot2))
clusterExport(cl, c("to_mcmc.list", "expo","expo_tot","portage_real","farm_bel","box_bel","sampling_days","N_calves","t_end","from_56","simul_mod_allFarms"))
t = as.numeric(Sys.time())
values_table = clusterMap(cl, colist_init = htmap[,1], colist_mid = htmap[,2], arriv = htmap[,3], fun = cont_strat)
stopCluster(cl)
as.numeric(Sys.time()) - t

htmap = cbind(htmap,
              min_mean_prev_f1 = as.numeric(sapply(values_table, "[[", "min_mean_prev_f1")),
              q1_mean_prev_f1 = as.numeric(sapply(values_table, "[[", "q1_mean_prev_f1")),
              median_mean_prev_f1 = as.numeric(sapply(values_table, "[[", "median_mean_prev_f1")),
              q3_mean_prev_f1 = as.numeric(sapply(values_table, "[[", "q3_mean_prev_f1")),
              max_mean_prev_f1 = as.numeric(sapply(values_table, "[[", "max_mean_prev_f1")),
              
              min_mean_prev_f2 = as.numeric(sapply(values_table, "[[", "min_mean_prev_f2")),
              q1_mean_prev_f2 = as.numeric(sapply(values_table, "[[", "q1_mean_prev_f2")),
              median_mean_prev_f2 = as.numeric(sapply(values_table, "[[", "median_mean_prev_f2")),
              q3_mean_prev_f2 = as.numeric(sapply(values_table, "[[", "q3_mean_prev_f2")),
              max_mean_prev_f2 = as.numeric(sapply(values_table, "[[", "max_mean_prev_f2")),
              
              min_mean_prev_f3 = as.numeric(sapply(values_table, "[[", "min_mean_prev_f3")),
              q1_mean_prev_f3 = as.numeric(sapply(values_table, "[[", "q1_mean_prev_f3")),
              median_mean_prev_f3 = as.numeric(sapply(values_table, "[[", "median_mean_prev_f3")),
              q3_mean_prev_f3 = as.numeric(sapply(values_table, "[[", "q3_mean_prev_f3")),
              max_mean_prev_f3 = as.numeric(sapply(values_table, "[[", "max_mean_prev_f3")),
              
              min_final_prev_f1 = as.numeric(sapply(values_table, "[[", "min_final_prev_f1")),
              q1_final_prev_f1 = as.numeric(sapply(values_table, "[[", "q1_final_prev_f1")),
              median_final_prev_f1 = as.numeric(sapply(values_table, "[[", "median_final_prev_f1")),
              q3_final_prev_f1 = as.numeric(sapply(values_table, "[[", "q3_final_prev_f1")),
              max_final_prev_f1 = as.numeric(sapply(values_table, "[[", "max_final_prev_f1")),
              
              min_final_prev_f2 = as.numeric(sapply(values_table, "[[", "min_final_prev_f2")),
              q1_final_prev_f2 = as.numeric(sapply(values_table, "[[", "q1_final_prev_f2")),
              median_final_prev_f2 = as.numeric(sapply(values_table, "[[", "median_final_prev_f2")),
              q3_final_prev_f2 = as.numeric(sapply(values_table, "[[", "q3_final_prev_f2")),
              max_final_prev_f2 = as.numeric(sapply(values_table, "[[", "max_final_prev_f2")),
              
              min_final_prev_f3 = as.numeric(sapply(values_table, "[[", "min_final_prev_f3")),
              q1_final_prev_f3 = as.numeric(sapply(values_table, "[[", "q1_final_prev_f3")),
              median_final_prev_f3 = as.numeric(sapply(values_table, "[[", "median_final_prev_f3")),
              q3_final_prev_f3 = as.numeric(sapply(values_table, "[[", "q3_final_prev_f3")),
              max_final_prev_f3 = as.numeric(sapply(values_table, "[[", "max_final_prev_f3")))

######################################################################

# 2) Plot of results:

htmap = dplyr::filter(htmap, (col_init == 6 & arr %in% c(0,0.34,0.68)) | (arr == 0.68 & col_init %in% c(0,3,6)))
htmap$col_mid = as.factor(htmap$col_mid)
htmap$arr = factor(htmap$arr, labels = c("0" = "0%", "0.34" = "34%", "0.68" = "68% (baseline)"))
htmap$col_init = factor(htmap$col_init, labels = c("0" = "No exposure", "3" = "3 days", "6" = "6 days (baseline)"))

mitig_strat_plot = function(indic, fixed_var, far=1, plot_lines=F){
  htmap_res = htmap[htmap$col_mid==10 & htmap$arr!="0%", c("col_init", "arr", paste0("median_mean_prev_f",far), paste0("median_final_prev_f",far))]
  print(htmap_res)
  
  if(indic == "mean_prev"){
    y_lab = "Mean prevalence over the cycle (%)"
  }else if(indic == "final_prev"){
    y_lab = "Prevalence at slaughter age (%)"
  }
  
  if(fixed_var == "col_init"){
    dat_plot = dplyr::filter(htmap, col_init == names(which.max(table(htmap$col_init))))
    p = ggplot(data = dat_plot, aes(x = arr))
    p = p + xlab("Prevalence on arrival")
  }
  
  if(fixed_var == "arr"){
    dat_plot = dplyr::filter(htmap, arr == names(which.max(table(htmap$arr))))
    p = ggplot(data = dat_plot, aes(x = col_init))
    p = p + xlab("Initial exposure to selective antibiotics")
  }
  
  p = p + geom_boxplot(stat = "identity", position = position_dodge(0.5), width=0.3,
                       aes(col = col_mid,
                           ymin = dat_plot[,paste0("min_",indic,"_f",far)],
                           lower = dat_plot[,paste0("q1_",indic,"_f",far)],
                           middle = dat_plot[,paste0("median_",indic,"_f",far)],
                           upper = dat_plot[,paste0("q3_",indic,"_f",far)],
                           ymax = dat_plot[,paste0("max_",indic,"_f",far)]))
  p = p + ylab(y_lab)
  p = p + ylim(0,100)
  
  if(indic == "final_prev" & fixed_var == "col_init" & far == 1){
    p = p + annotate("point", x = names(which.max(table(htmap$arr))), y = 20.4, shape = 8, size = 3.5, colour = "black")
  }
  if(indic == "final_prev" & fixed_var == "arr" & far == 1){
    p = p + annotate("point", x = names(which.max(table(htmap$col_init))), y = 20.4, shape = 8, size = 3.5, colour = "black")
  }
  
  p = p + scale_color_manual(values = c("0"="#1b9e77", "10"="#d95f02"), labels = c("0"="No mid-cycle exposure", "10"="10-day mid-cycle exposure"))
  
  p = p + theme_bw()
  p = p + theme(legend.title = element_blank(), legend.direction = "vertical")
  p
}

p1 = mitig_strat_plot(indic = "mean_prev", fixed_var = "col_init", far=1)
p2 = mitig_strat_plot(indic = "mean_prev", fixed_var = "arr", far=1)
p3 = mitig_strat_plot(indic = "final_prev", fixed_var = "col_init", far=1)
p4 = mitig_strat_plot(indic = "final_prev", fixed_var = "arr", far=1)

p_all = ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, common.legend = T, legend = "top", labels = c("A", "B", "C", "D"), label.x = 0.9, label.y = 0.95)
p_all
ggsave(p_all, file = "~/Figure mitigation strategies farm 1.pdf", width=7, height=6)

htmap[,c(1:3,grep("mean_prev_f1", colnames(htmap)))]
htmap[,c(1:3,grep("final_prev_f1", colnames(htmap)))]




# Sensitivity analysis: on farms B (2) and C (3)

p1 = mitig_strat_plot(indic = "mean_prev", fixed_var = "col_init", far=2)
p2 = mitig_strat_plot(indic = "mean_prev", fixed_var = "arr", far=2)
p3 = mitig_strat_plot(indic = "final_prev", fixed_var = "col_init", far=2)
p4 = mitig_strat_plot(indic = "final_prev", fixed_var = "arr", far=2)

p_all = ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, common.legend = T, legend = "top", labels = c("A", "B", "C", "D"), label.x = 0.9, label.y = 0.95)
p_all
ggsave(p_all, file = "~/Figure mitigation strategies farm 2.png", width=7, height=6)

htmap[,c(1:3,grep("mean_prev_f2", colnames(htmap)))]
htmap[,c(1:3,grep("final_prev_f2", colnames(htmap)))]

p1 = mitig_strat_plot(indic = "mean_prev", fixed_var = "col_init", far=3)
p2 = mitig_strat_plot(indic = "mean_prev", fixed_var = "arr", far=3)
p3 = mitig_strat_plot(indic = "final_prev", fixed_var = "col_init", far=3)
p4 = mitig_strat_plot(indic = "final_prev", fixed_var = "arr", far=3)

p_all = ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, common.legend = T, legend = "top", labels = c("A", "B", "C", "D"), label.x = 0.9, label.y = 0.95)
p_all
ggsave(p_all, file = "~/Figure mitigation strategies farm 3.png", width=7, height=6)

htmap[,c(1:3,grep("mean_prev_f3", colnames(htmap)))]
htmap[,c(1:3,grep("final_prev_f3", colnames(htmap)))]

