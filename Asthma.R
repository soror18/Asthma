# ==== Global GAM ====
library(mgcv)
library(gratia)
library(dplyr)

analyze_gam <- function(model, data, linear_vars, smooth_vars, output_prefix){

  delta_deviance <- function(model, data, var){
    fmla <- update(formula(model), paste(". ~ . -", var))
    mod_drop <- gam(fmla, data = data, method="REML", select=TRUE)
    dev_diff <- deviance(mod_drop) - deviance(model)
    return(dev_diff)
  }
  
  delta_results <- sapply(linear_vars,
                          function(v) delta_deviance(model, data, v))
  
  deriv_all <- derivatives(model, n = 500)
  
  smooth_summary <- summary(model)$s.table
  rownames(smooth_summary) <- gsub("s\\(|\\)", "", rownames(smooth_summary))
  
  tol <- 1e-8

  trend_list <- lapply(smooth_vars, function(var_name){
    
    d <- deriv_all %>%
      filter(.smooth == paste0("s(", var_name, ")")) %>%
      arrange(.data[[var_name]])
    
    x_vals <- d[[var_name]]
    mean_deriv <- mean(d$.derivative, na.rm=TRUE)
    
    trend <- ifelse(mean_deriv > 0, "positive", "negative")
    
    signs <- sign(d$.derivative)
    
    sign_change_idx <- which(diff(signs) != 0 &
                               (abs(d$.derivative[-1]) > tol |
                                  abs(d$.derivative[-length(d$.derivative)]) > tol))
    
    break_points <- if(length(sign_change_idx)==0)
      NA else x_vals[sign_change_idx + 1]
    
    if(var_name %in% rownames(smooth_summary)){
      f_val <- smooth_summary[var_name,"F"]
      p_val <- smooth_summary[var_name,"p-value"]
    } else {
      f_val <- NA
      p_val <- NA
    }
    
    data.frame(
      variable = var_name,
      delta_deviance = NA,
      mean_derivative = mean_deriv,
      trend = trend,
      break_points =
        ifelse(is.na(break_points), NA,
               paste(round(break_points,3), collapse=",")),
      F_value = f_val,
      p_value = p_val
    )
  })
  
  trend_table <- do.call(rbind, trend_list)
  
  lin_summary <- summary(model)$p.table
  matched_idx <- match(linear_vars, rownames(lin_summary))
  
  linear_table <- data.frame(
    variable = linear_vars,
    delta_deviance = delta_results,
    mean_derivative = NA,
    trend = NA,
    break_points = NA,
    F_value = lin_summary[matched_idx,"t value"],
    p_value = lin_summary[matched_idx,"Pr(>|t|)"]
  )

  final_table <- bind_rows(trend_table, linear_table)
  
  write.csv(final_table,
            paste0(output_prefix,"_final_table.csv"),
            row.names = FALSE)
  
  return(final_table)
}

data$ISO_A3 <- as.factor(data$ISO_A3)
data$Climate_type <- factor(data$Climate_type)
mod_full1 <- gam(
  Incidence_rate_0_19 ~ 
    s(GBD_PM_mean_log) + s(urban_index2) +
    Pop_0_19_ratio + s(Obesity_log) + s(Physicians_per_1000) +
    ppt_average_log + ws_average + Maximum_Temp_Variability_log +
    s(Land_TREES) + s(Land_GRASS_SHRUBS) +
    s(year) + factor(Climate_type) + s(ISO_A3, bs="re"),
  data = data, method="REML", select=TRUE
)
linear_vars1 <- c("Pop_0_19_ratio", "ws_average","Maximum_Temp_Variability_log","ppt_average_log") 
smooth_vars1 <- c("GBD_PM_mean_log","urban_index2","Obesity_log","Physicians_per_1000","Land_GRASS_SHRUBS","Land_TREES")

# ==== cluster ====
## ==== static ====
library(dplyr)
library(ggplot2)
library(tidyr)
library(cluster)
library(fmsb)
library(viridis)
library(e1071)  # Fuzzy C-Means

high_pm_threshold <- 3.2
q2_urban <- quantile(data$urban_index2, probs = 2/3, na.rm = TRUE)

cluster_input <- data %>%
  group_by(ISO_A3) %>%
  summarise(
    urb_mean = mean(urban_index2, na.rm = TRUE),
    urb_slope = coef(lm(urban_index2 ~ year))[2],
    pm_mean = mean(GBD_PM_mean_log, na.rm = TRUE),
    pm_slope = coef(lm(GBD_PM_mean_log ~ year))[2],
    high_pm_years = sum(GBD_PM_mean_log > high_pm_threshold),
    high_urban_years = sum(urban_index2 > q2_urban)
  ) %>% ungroup()

cluster_vars <- cluster_input[, -1]
scaled_vars <- scale(cluster_vars)

cluster_methods <- list(
  Kmeans = list(
    fun = function(x) kmeans(x, centers = 4, nstart = 25),
    assign = function(res) factor(res$cluster)
  ),
  PAM = list(
    fun = function(x) pam(x, k = 4),
    assign = function(res) factor(res$clustering)
  ),
  FCM = list(
    fun = function(x) cmeans(x, centers = 4, m = 2, iter.max = 100, verbose = FALSE),
    assign = function(res) factor(apply(res$membership, 1, which.max))
  )
)

  for(method_name in names(cluster_methods)) 
  cat("Processing", method_name, "...\n")

  res <- cluster_methods[[method_name]]$fun(scaled_vars)
  cluster_input[[paste0(method_name, "_cluster")]] <- cluster_methods[[method_name]]$assign(res)

  data <- data %>%
    left_join(cluster_input[, c("ISO_A3", paste0(method_name, "_cluster"))], by = "ISO_A3")
  
## ==== DTW;DTW+static ====
  library(dplyr)
  library(tidyr)
  library(cluster)
  library(fmsb)
  library(gridExtra)
  library(purrr)
  library(dtw)
  library(tibble)
  
  high_pm_threshold <- 3.2
  q2_urban <- quantile(data$urban_index2, probs = 2/3, na.rm = TRUE)

  data_seq <- data %>% arrange(ISO_A3, year) %>% select(ISO_A3, year, urban_index2, GBD_PM_mean_log)
  seq_list <- data_seq %>%
    group_by(ISO_A3) %>%
    summarise(
      urban_vec = list(urban_index2),
      pm_log_vec = list(GBD_PM_mean_log),
      highPM_vec = list(as.integer(GBD_PM_mean_log > high_pm_threshold)),
      highUrban_vec = list(as.integer(urban_index2 > q2_urban))
    ) %>%
    pmap(function(urban_vec, pm_log_vec, highPM_vec, highUrban_vec, ISO_A3) {
      list(
        urban_vec = urban_vec[[1]],
        pm_log_vec = pm_log_vec[[1]],
        highPM_vec = highPM_vec[[1]],
        highUrban_vec = highUrban_vec[[1]]
      )
    }) %>%
    set_names(unique(data_seq$ISO_A3))
  
  country_codes <- names(seq_list)
  n <- length(country_codes)

  dist_urban <- dist_pm <- dist_highPM <- dist_highUrban <- matrix(0, n, n)
  for (i in seq_len(n-1)) {
    for (j in (i+1):n) {
      d_u   <- dtw(seq_list[[i]]$urban_vec, seq_list[[j]]$urban_vec, distance.only = TRUE)$distance
      d_p   <- dtw(seq_list[[i]]$pm_log_vec, seq_list[[j]]$pm_log_vec, distance.only = TRUE)$distance
      d_hPM <- dtw(seq_list[[i]]$highPM_vec, seq_list[[j]]$highPM_vec, distance.only = TRUE)$distance
      d_hU  <- dtw(seq_list[[i]]$highUrban_vec, seq_list[[j]]$highUrban_vec, distance.only = TRUE)$distance
      dist_urban[i,j] <- dist_urban[j,i] <- d_u
      dist_pm[i,j]    <- dist_pm[j,i]    <- d_p
      dist_highPM[i,j]<- dist_highPM[j,i]<- d_hPM
      dist_highUrban[i,j]<- dist_highUrban[j,i]<- d_hU
    }
  }
  diag(dist_urban) <- diag(dist_pm) <- diag(dist_highPM) <- diag(dist_highUrban) <- 0
  
  D_dyn <- 0.3 * dist_urban + 0.3 * dist_pm + 0.2 * dist_highPM + 0.2 * dist_highUrban
  dist_dyn <- as.dist(D_dyn)
  
  pam_dyn <- pam(dist_dyn, k = 4)
  cluster_dyn_df <- tibble(ISO_A3 = country_codes,
                           DTW_only_cluster = factor(pam_dyn$clustering))

  data <- left_join(data, cluster_dyn_df, by = "ISO_A3")

  cluster_input <- data %>%
    group_by(ISO_A3) %>%
    summarise(
      urb_mean = mean(urban_index2, na.rm = TRUE),
      urb_slope = coef(lm(urban_index2 ~ year))[2],
      pm_mean = mean(GBD_PM_mean_log, na.rm = TRUE),
      pm_slope = coef(lm(GBD_PM_mean_log ~ year))[2],
      high_pm_years = sum(GBD_PM_mean_log > high_pm_threshold),
      high_urban_years = sum(urban_index2 > q2_urban)
    ) %>% ungroup()
  
  static_mat <- cluster_input %>% column_to_rownames("ISO_A3") %>% as.matrix()
  static_scaled <- scale(static_mat)
  dist_static <- as.matrix(dist(static_scaled))
  
  w_dyn <- 0.7
  w_stat <- 0.3
  D_total <- w_dyn * D_dyn + w_stat * dist_static
  dist_total <- as.dist(D_total)
  
  pam_total <- pam(dist_total, k = 4)
  cluster_total_df <- tibble(ISO_A3 = country_codes,
                             DTW_static_cluster = factor(pam_total$clustering))
  
  data <- left_join(data, cluster_total_df, by = "ISO_A3")
  
  # ==== D ====
  calc_coupling_df <- function(df, var1, var2, name1, name2){
    
    tmp <- df
    
    tmp$P <- (tmp[[var1]] - min(tmp[[var1]], na.rm=TRUE)) /
      (max(tmp[[var1]], na.rm=TRUE) - min(tmp[[var1]], na.rm=TRUE))
    
    tmp$A <- (max(tmp[[var2]], na.rm=TRUE) - tmp[[var2]]) /
      (max(tmp[[var2]], na.rm=TRUE) - min(tmp[[var2]], na.rm=TRUE))
    
    tmp$C <- 2 * sqrt(tmp$P * tmp$A) / (tmp$P + tmp$A)
    
    alpha <- 0.5
    beta <- 0.5
    tmp$T <- alpha * tmp$P + beta * tmp$A
    
    tmp$D <- sqrt(tmp$C * tmp$T)
    
    tmp$Var1 <- name1
    tmp$Var2 <- name2
    
    tmp <- tmp %>%
      dplyr::select(ISO_A3, Cluster, year, Var1, Var2, D)
    
    return(tmp)
  }
  
  
  D_PM_U <- calc_coupling_df(
    df,
    "urban_index2",
    "GBD_PM_mean_log",
    "Urbanization",
    "PM2.5"
  )
  
  D_PM_A <- calc_coupling_df(
    df,
    "GBD_PM_mean_log",
    "Incidence_rate_0_19",
    "PM2.5",
    "Asthma"
  )
  
  D_U_A <- calc_coupling_df(
    df,
    "urban_index2",
    "Incidence_rate_0_19",
    "Urbanization",
    "Asthma"
  )
  
  
  D_country_year <- dplyr::bind_rows(D_PM_U, D_PM_A, D_U_A)
  
  
  D_cluster_year <- D_country_year %>%
    dplyr::group_by(Cluster, year, Var1, Var2) %>%
    dplyr::summarise(
      D_mean = mean(D, na.rm = TRUE),
      .groups = "drop"
    )
  
  
  D_cluster_range <- D_cluster_year %>%
    dplyr::group_by(Cluster, Var1, Var2) %>%
    dplyr::summarise(
      D_min = min(D_mean, na.rm = TRUE),
      D_max = max(D_mean, na.rm = TRUE),
      .groups = "drop"
    )
  
  # ==== Cluster GAM ====
  library(mgcv)
  library(gratia)
  library(dplyr)

  analyze_gam <- function(model, data, linear_vars, smooth_vars, output_prefix){
    
    delta_deviance <- function(model, data, var){
      fmla <- update(formula(model), paste(". ~ . -", var))
      mod_drop <- gam(fmla, data = data, method = "REML", select = TRUE)
      dev_diff <- deviance(mod_drop) - deviance(model)
      return(dev_diff)
    }
    
    delta_results <- sapply(linear_vars,
                            function(v) delta_deviance(model, data, v))
    
    deriv_all <- derivatives(model)
    
    smooth_summary <- summary(model)$s.table
    rownames(smooth_summary) <- gsub("s\\(|\\)", "", rownames(smooth_summary))
    
    trend_list <- lapply(smooth_vars, function(var_name){
      
      d <- deriv_all %>%
        filter(.smooth == paste0("s(", var_name, ")")) %>%
        arrange(.data[[var_name]])
      
      x_vals <- d[[var_name]]
      
      mean_deriv <- mean(d$.derivative, na.rm = TRUE)
      
      trend <- ifelse(mean_deriv > 0, "positive", "negative")
      
      signs <- sign(d$.derivative)
      sign_change_idx <- which(diff(signs) != 0)
      
      break_points <- if(length(sign_change_idx) == 0)
        NA
      else
        x_vals[sign_change_idx + 1]
      
      segment_starts <- c(1, sign_change_idx + 1)
      segment_ends <- c(sign_change_idx, length(x_vals))
      
      local_trends <- sapply(seq_along(segment_starts), function(i){
        seg_mean <- mean(
          d$.derivative[segment_starts[i]:segment_ends[i]],
          na.rm = TRUE
        )
        if(is.na(seg_mean)) return(NA)
        if(seg_mean > 0) "positive" else "negative"
      })
      
      if(var_name %in% rownames(smooth_summary)){
        f_val <- smooth_summary[var_name, "F"]
        p_val <- smooth_summary[var_name, "p-value"]
      } else {
        f_val <- NA
        p_val <- NA
      }
      
      data.frame(
        variable = var_name,
        delta_deviance = NA,
        mean_derivative = mean_deriv,
        trend = trend,
        break_points =
          ifelse(is.na(break_points),
                 NA,
                 paste(round(break_points, 3), collapse = ",")),
        local_trends = paste(local_trends, collapse = ","),
        F_value = f_val,
        p_value = p_val,
        stringsAsFactors = FALSE
      )
    })
    
    trend_table <- do.call(rbind, trend_list)
    
    lin_summary <- summary(model)$p.table
    
    matched_idx <- match(linear_vars, rownames(lin_summary))
    
    linear_table <- data.frame(
      variable = linear_vars,
      delta_deviance = delta_results,
      mean_derivative = NA,
      trend = NA,
      break_points = NA,
      local_trends = NA,
      F_value = lin_summary[matched_idx, "t value"],
      p_value = lin_summary[matched_idx, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    
    final_table <- bind_rows(trend_table, linear_table) %>%
      arrange(match(variable, c(smooth_vars, linear_vars)))
    
    write.csv(
      final_table,
      paste0(output_prefix, "_final_table.csv"),
      row.names = FALSE
    )
    
    return(final_table)
  }
  
  # ==============================
  # cluster 1
  # ==============================
  data1$ISO_A3 <- as.factor(data1$ISO_A3)
  data1$Climate_type <- factor(data1$Climate_type)
  mod_full1 <- gam(
    Incidence_rate_0_19 ~ 
      s(GBD_PM_mean_log) + s(urban_index2) +
      Pop_0_19_ratio + Obesity_log + Physicians_per_1000 +
      ppt_average_log + ws_average + Maximum_Temp_Variability_log +
      s(Land_TREES) + Land_GRASS_SHRUBS +
      s(year) + factor(Climate_type) + s(ISO_A3, bs="re"),
    data = data1, method="REML", select=TRUE
  )
  linear_vars1 <- c("Pop_0_19_ratio","Obesity_log","Physicians_per_1000",
                    "ppt_average_log","ws_average","Maximum_Temp_Variability_log",
                    "Land_GRASS_SHRUBS","year")
  smooth_vars1 <- c("GBD_PM_mean_log","urban_index2","Land_TREES")
  
  # ==============================
  # cluster 2
  # ==============================
  data2$ISO_A3 <- as.factor(data2$ISO_A3)
  data2$Climate_type <- factor(data2$Climate_type)
  mod_full2 <- gam(
    Incidence_rate_0_19 ~ 
      s(GBD_PM_mean_log) + s(urban_index2) +
      Pop_0_19_ratio + Obesity_log + Physicians_per_1000 +
      ppt_average_log + ws_average + Maximum_Temp_Variability_log +
      s(Land_TREES) + Land_GRASS_SHRUBS +
      s(year) + factor(Climate_type) + s(ISO_A3, bs="re"),
    data = data2, method="REML", select=TRUE
  )
  linear_vars2 <- linear_vars1
  smooth_vars2 <- smooth_vars1

  # ==============================
  # cluster 3
  # ==============================
  data3$ISO_A3 <- as.factor(data3$ISO_A3)
  data3$Climate_type <- factor(data3$Climate_type)
  mod_full3 <- gam(
    Incidence_rate_0_19 ~ 
      s(GBD_PM_mean_log) + s(urban_index2)+
      Pop_0_19_ratio + Obesity_log + Physicians_per_1000 +
      ppt_average_log + ws_average + Maximum_Temp_Variability_log +
      Land_TREES + s(Land_GRASS_SHRUBS) +
      s(year) + factor(Climate_type) + s(ISO_A3, bs="re"),
    data = data3, method="REML", select=TRUE
  )
  linear_vars3 <- c("Pop_0_19_ratio","Obesity_log","Physicians_per_1000",
                    "ppt_average_log","ws_average","Maximum_Temp_Variability_log",
                    "Land_TREES")
  smooth_vars3 <- c("GBD_PM_mean_log","urban_index2","Land_GRASS_SHRUBS")
  
  # ==============================
  # cluster 4
  # ==============================
  data4$ISO_A3 <- as.factor(data4$ISO_A3)
  data4$Climate_type <- factor(data4$Climate_type)
  mod_full4 <- gam(
    Incidence_rate_0_19 ~ 
      s(GBD_PM_mean_log) + s(urban_index2) +
      Pop_0_19_ratio + Obesity_log + Physicians_per_1000 +
      ppt_average_log + ws_average + Maximum_Temp_Variability_log +
      Land_TREES + s(Land_GRASS_SHRUBS) +
      s(year) + factor(Climate_type) + s(ISO_A3, bs="re"),
    data = data4, method="REML", select=TRUE
  )
  smooth_vars4 <- c("GBD_PM_mean_log","urban_index2","Land_GRASS_SHRUBS")
  linear_vars4 <- c("Pop_0_19_ratio","Obesity_log","Physicians_per_1000",
                    "ppt_average_log","ws_average","Maximum_Temp_Variability_log",
                    "Land_TREES","year")
  
  # ==== Temporal stability ====
  core_vars <- data %>%
    dplyr::select(urban_index2, GBD_PM_mean_log, highPM, highUrban)
  
  cor_matrix <- cor(core_vars, use = "pairwise.complete.obs")
  
  colnames(cor_matrix) <- rownames(cor_matrix) <- c(
    "Urbanization Index",
    "PM2.5 (log)",
    "High PM",
    "High Urbanization"
  )
  
  data <- read.csv("E:/GBD/analysis/data0821.csv")
  
  q2_urban <- quantile(data$urban_index2, probs = 2/3, na.rm = TRUE)
  
  data <- data %>%
    dplyr::mutate(
      highPM = as.integer(GBD_PM_mean_log > 3.2),
      highUrban = as.integer(urban_index2 > q2_urban)
    )
  
  z_truncate <- function(x, limit = 3) {
    z <- scale(x)
    pmax(pmin(z, limit), -limit)
  }
  
  data <- data %>%
    dplyr::mutate(
      urban_z = z_truncate(urban_index2),
      pm_z = z_truncate(GBD_PM_mean_log),
      highPM_z = z_truncate(highPM),
      highUrban_z = z_truncate(highUrban)
    )
  
  cluster_centers <- data %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      urban_center_z = mean(urban_z, na.rm = TRUE),
      pm_center_z = mean(pm_z, na.rm = TRUE),
      highPM_center_z = mean(highPM_z, na.rm = TRUE),
      highUrban_center_z = mean(highUrban_z, na.rm = TRUE),
      .groups = "drop"
    )
  
  w_binary <- 0.3
  
  data <- data %>%
    dplyr::left_join(cluster_centers, by = "cluster") %>%
    dplyr::mutate(
      dist_self = sqrt(
        (urban_z - urban_center_z)^2 +
          (pm_z - pm_center_z)^2 +
          w_binary * ((highPM_z - highPM_center_z)^2 +
                        (highUrban_z - highUrban_center_z)^2)
      )
    )
  
  cluster_stats <- data %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      mu = mean(dist_self, na.rm = TRUE),
      sigma = sd(dist_self, na.rm = TRUE),
      .groups = "drop"
    )
  
  data <- data %>%
    dplyr::left_join(cluster_stats, by = "cluster") %>%
    dplyr::mutate(z_self = (dist_self - mu) / sigma)
  
  other_clusters <- cluster_centers %>%
    dplyr::rename(
      other_cluster = cluster,
      urban_c = urban_center_z,
      pm_c = pm_center_z,
      highPM_c = highPM_center_z,
      highUrban_c = highUrban_center_z
    )
  
  dist_df <- data %>%
    dplyr::select(
      ISO_A3, Country, year,
      urban_z, pm_z, highPM_z, highUrban_z, cluster
    ) %>%
    tidyr::crossing(other_clusters) %>%
    dplyr::filter(cluster != other_cluster) %>%
    dplyr::mutate(
      dist = sqrt(
        (urban_z - urban_c)^2 +
          (pm_z - pm_c)^2 +
          w_binary * ((highPM_z - highPM_c)^2 +
                        (highUrban_z - highUrban_c)^2)
      )
    )
  
  nearest_df <- dist_df %>%
    dplyr::group_by(ISO_A3, Country, year) %>%
    dplyr::slice_min(dist, n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      ISO_A3, Country, year,
      nearest_cluster = other_cluster,
      dist_min_other = dist
    )
  
  data <- data %>%
    dplyr::left_join(nearest_df, by = c("ISO_A3", "Country", "year")) %>%
    dplyr::mutate(
      deviate_flag = as.integer(z_self > 2 & dist_self > dist_min_other)
    )
  
  threshold <- 2
  
  data <- data %>%
    dplyr::mutate(
      urban_dev = urban_z - urban_center_z,
      pm_dev = pm_z - pm_center_z,
      highPM_dev = highPM_z - highPM_center_z,
      highUrban_dev = highUrban_z - highUrban_center_z,
      urban_flag = as.integer(abs(urban_dev) > threshold),
      pm_flag = as.integer(abs(pm_dev) > threshold),
      highPM_flag = as.integer(abs(highPM_dev) > threshold),
      highUrban_flag = as.integer(abs(highUrban_dev) > threshold)
    )
  
  deviate_years <- data %>%
    dplyr::filter(deviate_flag == 1) %>%
    dplyr::select(
      ISO_A3, Country, cluster, year,
      nearest_cluster, dist_self, dist_min_other, z_self,
      urban_flag, pm_flag, highPM_flag, highUrban_flag,
      urban_dev, pm_dev, highPM_dev, highUrban_dev
    )