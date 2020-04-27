basefilter <-function(x){
  
  min <- apply(x[1:start_pulse, ], 2, min)
  max <- apply(x[1:start_pulse, ], 2, max)
  
  filter <- min > lower_filter & max < upper_filter
  
  #filtered_min <- 
  
  x[filter]
  
}

resp_filter_fun <- function(x){
  
  max <- apply(x[start_pulse:end_pulse, ], 2, max)
  filter_resp <- max > activity_threshold
  x[filter_resp]
}

read_trim <- function(x){
  
  input <- read.csv(x)
  input[, 2:NCOL(input)]
  
}


save_table_fun <-function(out_table, pre, name){
  save_path <- as.character(paste(out_path, pre, name, ".csv", sep=""))
  write.csv(out_table, save_path)
}


plot_fun <- function(out_table, pre, name){
  b  <- as.character(paste(out_path, pre, name, ".jpg", sep = ""))
  jpeg(b)
  plot.ts(out_table, plot.type = "single", main=name)
  dev.off()
}

save_plot_fun <- function(out_table, pre, name){
  save_table_fun(out_table, pre, name)
  plot_fun(out_table, pre, name)
}

title_fun <- function(X){
  title <-basename(X)
  as.character(substr(title, 1, 10))
}

trim_names <- function(x){
  n <- basename(x)
  substr(n, 1, 10)
}

AUC_column <- function(col){
  avg <- mean(col[1:start_pulse])
  start_vector <- col-avg
  
  a <- start_vector[start_pulse:end_pulse]
  b <- start_vector[(start_pulse+1):(end_pulse+1)]
  c <- a*b/2*2
  c
}


base_shift_fun <- function(x){
  min_x <- apply(x, 2, min)
  shift <- 1- min_x
  sweep(x, 2, shift, FUN = "+")
}