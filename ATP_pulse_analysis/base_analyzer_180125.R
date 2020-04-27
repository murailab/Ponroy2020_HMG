start_pulse = 100
end_pulse =299

#upper_filter = 1.15
#lower_filter = 0.85

activity_threshold =1.3


Ca_pulse <- function(files){
  
  file_name <- as.character(title_fun(files))
  
  save_prefix <- "/all_"
  start_file <- read_trim(files)
  
  #save_table_fun(input, save_prefix, file_name)
  plot_fun(start_file, save_prefix, file_name)
  
# filtered <- basefilter(start_file)
  
  # plot.ts(filtered, plot.type = "single") 
  
# filtered_save_prefix <- "/filtered_"
# save_plot_fun(filtered, filtered_save_prefix, file_name)
  
  shifted_file <- base_shift_fun(start_file)
  
  responders <- resp_filter_fun(shifted_file)
  
  
  
  responders_save_prefix <- "/responders_"
  if (NCOL(responders)>0){
  save_plot_fun(responders, responders_save_prefix, file_name)
  } else {
    save_table_fun(responders, responders_save_prefix, file_name)
  }
  
  AUC_list_all <- apply(shifted_file, 2, AUC_column)
  AUCall_save_prefix <- "/AUC_list_all"
  save_table_fun(AUC_list_all, AUCall_save_prefix, file_name)
  
  AUC_list <- apply(responders, 2, AUC_column)
  AUC_save_prefix <- "/AUC_list_"
  save_table_fun(AUC_list, AUC_save_prefix, file_name)

  #summarize data
  n_resp        <- NCOL(responders)
  n_ROIs        <- NCOL(shifted_file)
  perc_resp     <- n_resp/n_ROIs*100
  mean_amp_all  <- mean(apply(shifted_file[start_pulse:end_pulse, ], 2, max))
  sd_amp_all    <- sd(apply(shifted_file[start_pulse:end_pulse, ], 2, max))
  mean_AUC_all  <- mean(apply(AUC_list_all, 2, sum))
  sd_AUC_all    <- sd(apply(AUC_list_all, 2, sum))

  
 if (NCOL(responders)>1){

    mean_amp_resp <- mean(apply(responders[start_pulse:end_pulse, ], 2, max))
    sd_amp_resp   <- sd(apply(responders[start_pulse:end_pulse, ], 2, max))
    mean_AUC      <- mean(apply(AUC_list, 2, sum))
    sd_AUC        <- sd(apply(AUC_list, 2, sum))
  print("hurray")
  } else {
    mean_AUC      <- NA
    sd_AUC        <- NA
    mean_amp_resp <- NA
    sd_amp_resp   <- NA
    print("booo")
    print(file_name)
  }

  c(n_resp, n_ROIs, perc_resp, mean_amp_all, sd_amp_all, mean_amp_resp, sd_amp_resp, mean_AUC_all, sd_AUC_all, mean_AUC, sd_AUC)
}


path = getwd()

out_path = paste(path, "results", sep = "/")
file.names <- dir(path, pattern =".csv")
results <- sapply(file.names, Ca_pulse)

col_labels <- trim_names(file.names)
row_labels <- c("n_resp", "n_ROIs", "perc_resp", "mean_amp_all", "sd_amp_all", "mean_amp_resp", "sd_amp_resp","mean_AUC_all", "sd_AUC_all", "mean_AUC", "sd_AUC")

colnames(results) <- col_labels
rownames(results) <- row_labels

write.csv(results, paste(out_path, "results_summary.csv", sep = "/"))


library("MASS", lib.loc="/Library/Frameworks/R.framework/Versions/3.4/Resources/library")
write.matrix(results, paste(out_path, "results_matrix", sep = "/"))

