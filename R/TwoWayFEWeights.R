library(dplyr)
library(fixest)
library(estimatr)

printf <- function(...)print(sprintf(...))
fn_ctrl_rename <- function(x) { paste("ctrl", x, sep="_") }
get_controls_rename <- function(controls) { unlist(lapply(controls, fn_ctrl_rename)) }
fn_treatment_rename <- function(x) { paste("OT", x, sep="_")}
get_treatments_rename <- function(treatments) { unlist(lapply(treatments, fn_treatment_rename)) }
fn_treatment_weight_rename <- function(x) { paste("weight_", x, sep = "") }
fn_random_weight_rename <- function(x) { paste("RW", x, sep="_")}
get_random_weight_rename <- function(ws) { unlist(lapply(ws, fn_random_weight_rename)) }

twowayfeweights_rename_var <- function(df, Y, G, T, D, D0, controls, treatments, random_weights) {
  controls_rename <- get_controls_rename(controls)
  treatments_rename <- get_treatments_rename(treatments)
  
  if (length(random_weights) > 0) {
    random_weight_rename <- get_random_weight_rename(random_weights)
    random_weight_df <- df[random_weights]
    colnames(random_weight_df) <- random_weight_rename
  }
  
  original_names = c(Y, G, T, D, controls, treatments)
  new_names = c("Y", "G", "T", "D", controls_rename, treatments_rename)
  
  if (!is.null(D0)) {
    original_names = c(original_names, D0)
    new_names = c(new_names, "D0")
  }
  
  df <- data.frame(df) %>% select_at(vars(original_names))
  colnames(df) <- new_names
  
  if (length(random_weights) > 0) {
    df <- cbind(df, random_weight_df)
  }
  
  df
}

twowayfeweights_normalize_var <- function(df, varname){
  var = sym(varname)
  sdf <- df %>%
    group_by(.data$G, .data$T) %>%
    summarise(tmp_mean_gt = mean(!!var), tmp_sd_gt = sd(!!var))
  
  tmp_sd_gt_sum = sum(sdf$tmp_sd_gt, na.rm=TRUE)
  if (tmp_sd_gt_sum > 0) {
    df <- df %>% 
      left_join(sdf, by=c("T", "G")) %>%
      mutate(!!var := .data$tmp_mean_gt) %>%
      select(-.data$tmp_mean_gt) %>%
      select(-.data$tmp_sd_gt)
  }
  
  list(retcode = (tmp_sd_gt_sum > 0), df = df)
}

twowayfeweights_transform <- function(df, controls, weights, treatments) {
  ret = twowayfeweights_normalize_var(df, "D")
  if (ret$retcode) {
    df <- ret$df
    printf("The treatment variable in the regression varies within some group * period cells.")
    printf("The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) apply to two-way fixed effects regressions")
    printf("with a group * period level treatment.")
    printf("The command will replace the treatment by its average value in each group * period.")
    printf("The results below apply to the two-way fixed effects regression with that treatment variable.")
  }
  
  for (control in controls) {
    ret = twowayfeweights_normalize_var(df, control)
    if (ret$retcode) {
      df <- ret$df
      printf("The control variable %s in the regression varies within some group * period cells.", control)
      printf("The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) apply to two-way fixed effects regressions")
      printf("with controls apply to group * period level controls.")
      printf("The command will replace replace control variable %s by its average value in each group * period.", control)
      printf("The results below apply to the regression with control variable %s averaged at the group * period level.", control)
    }
  }
  
  for (treatment in treatments) {
    ret = twowayfeweights_normalize_var(df, treatment)
    if (ret$retcode) {
      df <- ret$df
      printf("The other treatment variable %s in the regression varies within some group * period cells.", treatment)
      printf("The results in de Chaisemartin, C. and D'Haultfoeuille, X. (2020) apply to two-way fixed effects regressions")
      printf("with several treatments apply to group * period level controls.")
      printf("The command will replace replace other treatment variable %s by its average value in each group * period.", treatment)
      printf("The results below apply to the regression with other treatment variable %s averaged at the group * period level.", treatment)
    }
  }
  
  if (is.null(weights)) {
    df$weights <- 1
  } else {
    df$weights <- weights
  }
  
  df$Tfactor <- factor(df$T)
  TfactorLevels <- length(levels(df$Tfactor))
  df <- df %>% mutate(TFactorNum = as.numeric(factor(.data$Tfactor, labels = seq(1:TfactorLevels))))
  
  df
}

twowayfeweights_filter <- function(df, cmd_type, controls, treatments) {
  # Remove rows with NA values
  if (cmd_type != "fdTR") {
    df <- df %>%
      mutate(tag = rowSums(across(.cols = c("Y", "G", "T", "D", controls, treatments), .fns = is.na))) %>%
      filter(.data$tag == 0) %>%
      select(-.data$tag)
  } else {
    df <- df %>%
      mutate(tag1 = rowSums(across(.cols = c("D", "T", "Y"), .fns = is.na))) %>%
      mutate(tag2 = rowSums(across(.cols = c("D0"), .fns = is.na))) %>%
      filter(.data$tag1 == 0 | .data$tag2 == 0)
    
    if (length(controls) > 0) {
      df <- df %>%
        mutate(tag3 = rowSums(across(.cols = controls, .fns = is.na))) %>%
        filter(.data$tag1 == 1 | .data$tag3 == 0) %>%
        select(-.data$tag3)
    }
    df <- df %>% select(-.data$tag1, -.data$tag2)
  }
  
  df
}

twowayfeweights_calculate_fetr <- function(df, controls) {
  mean_D <- weighted.mean(df$D, df$weights, na.rm = TRUE)
  obs <- sum(df$weights)
  gdf <- df %>% group_by(.data$G, .data$T) %>% summarise(P_gt = sum(.data$weights))
  df <- df %>% 
    left_join(gdf, by=c("T", "G")) %>% 
    mutate(P_gt = .data$P_gt / obs) %>% 
    mutate(nat_weight = .data$P_gt * .data$D / mean_D)
  
  if (length(controls) == 0) {
    formula = "D ~ 1"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("D ~ ", formula, sep = "")
  }
  formula = paste(formula, " | G + Tfactor", sep = "")
  denom.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  df$eps_1 <- residuals(denom.lm)
  df$eps_1_E_D_gt <- df$eps_1 * df$D
  denom_W <- weighted.mean(df$eps_1_E_D_gt, df$weights, na.rm = TRUE)
  
  df <- df %>% 
    mutate(W = .data$eps_1 * mean_D / denom_W) %>% 
    mutate(weight_result = .data$W * .data$nat_weight) %>%
    select(-.data$eps_1, -.data$P_gt)
  
  if (length(controls) == 0) {
    formula = "Y ~ D"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("Y ~ D + ", formula, sep = "")
  }
  formula = paste(formula, " | G + Tfactor", sep = "")
  beta.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  beta <- as.numeric(coef(beta.lm)["D"])
  
  # * Keeping only one observation in each group * period cell
  # This should be done after this function
  # bys `group' `time': gen group_period_unit=(_n==1)	
  # 	drop if group_period_unit==0
  # 	drop group_period_unit
  df <- df %>%
    group_by(.data$G, .data$Tfactor) %>%
    filter(row_number(.data$D) == 1)
  
  list(df = df, beta = beta)
}

twowayfeweights_calculate_fdtr <- function(df, controls) {
  mean_D0 <- weighted.mean(df$D0, df$weights, na.rm = TRUE)
  obs <- sum(df$weights)
  gdf <- df %>% group_by(.data$G, .data$T) %>% summarise(P_gt = sum(.data$weights))
  df <- df %>% 
    left_join(gdf, by=c("G", "T")) %>% 
    mutate(P_gt = .data$P_gt / obs) %>% 
    mutate(nat_weight = .data$P_gt * .data$D0 / mean_D0)
  
  if (length(controls) == 0) {
    formula = "D ~ 1"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("D ~ ", formula, sep = "")
  }
  formula = paste(formula, " | Tfactor", sep = "")
  denom.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  df$eps_2 <- residuals(denom.lm)
  # df$eps_2 <- df$D - predict(denom.lm, df)
  
  df <- df %>% mutate(eps_2 = ifelse(is.na(.data$eps_2), 0, .data$eps_2))

  if (length(controls) == 0) {
    formula = "Y ~ D"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("Y ~ D + ", formula, sep = "")
  }
  formula = paste(formula, " | Tfactor", sep = "")
  beta.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  beta <- as.numeric(coef(beta.lm)["D"])
  
  df <- df %>% 
    arrange(.data$G, .data$TFactorNum) %>%
    group_by(.data$G) %>% 
    mutate(w_tilde_2 = ifelse(.data$TFactorNum + 1 == lead(.data$TFactorNum), .data$eps_2 - lead(.data$eps_2) * (lead(.data$P_gt) / .data$P_gt), NA)) %>%
    mutate(w_tilde_2 = ifelse(is.na(.data$w_tilde_2) | is.infinite(.data$w_tilde_2), .data$eps_2, .data$w_tilde_2)) %>%
    mutate(w_tilde_2_E_D_gt = .data$w_tilde_2 * .data$D0)
  
  denom_W <- weighted.mean(df$w_tilde_2_E_D_gt, df$P_gt, na.rm = TRUE)
  df <- df %>% 
    mutate(W = .data$w_tilde_2 * mean_D0 / denom_W) %>% 
    mutate(weight_result = .data$W * .data$nat_weight)
  df <- df %>%
    select(-.data$eps_2, -.data$P_gt, -.data$w_tilde_2, -.data$w_tilde_2_E_D_gt)
  
  list(df = df, beta = beta)
}

twowayfeweights_calculate_fes <- function(df, controls) {
  obs <- sum(df$weights)
  gdf <- df %>% group_by(.data$G, .data$T) %>% summarise(P_gt = sum(.data$weights))
  df <- df %>% 
    left_join(gdf, by=c("T", "G")) %>% 
    mutate(P_gt = .data$P_gt / obs)
  
  if (length(controls) == 0) {
    formula = "D ~ 1"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("D ~ ", formula, sep = "")
  }
  formula = paste(formula, " | G + Tfactor", sep = "")
  denom.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  df$eps_1 <- residuals(denom.lm)
  # df$eps_1 <- df$D - predict(denom.lm, df)
  
  df <- df %>% 
    mutate(eps_1_weight = .data$eps_1 * .data$weights) %>%
    arrange(.data$G, .data$Tfactor) %>%
    group_by(.data$G) %>%
    mutate(E_eps_1_g_ge_aux = rev(cumsum(rev(.data$eps_1_weight)))) %>%
    mutate(weights_aux = rev(cumsum(rev(.data$weights)))) %>%
    mutate(E_eps_1_g_ge = .data$E_eps_1_g_ge_aux / .data$weights_aux)
  
  if (length(controls) == 0) {
    formula = "Y ~ D"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("Y ~ D + ", formula, sep = "")
  }
  formula = paste(formula, " | G + Tfactor", sep = "")
  beta.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  beta <- as.numeric(coef(beta.lm)["D"])
  
  # * Keeping only one observation in each group * period cell
  #   bys `group' `time': gen group_period_unit=(_n==1)	
  # 	drop if group_period_unit==0
  # 	drop group_period_unit
  
  # df <- df %>% 
  #   group_by(.data$G, .data$Tfactor) %>%
  #   summarize(.data$P_gt, .data$nat_weight)
  
  df <- df %>% 
    arrange(.data$G, .data$Tfactor) %>%
    group_by(.data$G) %>% 
    mutate(delta_D = ifelse(.data$TFactorNum - 1 == lag(.data$TFactorNum), .data$D - lag(.data$D), NA)) %>%
    filter(!is.na(.data$delta_D)) %>%
    mutate(abs_delta_D = abs(.data$delta_D)) %>%
    mutate(s_gt = case_when(.data$delta_D > 0 ~ 1,
                            .data$delta_D < 0 ~ -1,
                            TRUE ~ 0)) %>%
    mutate(nat_weight = .data$P_gt * .data$abs_delta_D)
  
  P_S = sum(df$nat_weight, na.rm = TRUE)
  df <- df %>% 
    mutate(nat_weight = .data$nat_weight / P_S) %>%
    mutate(om_tilde_1 = .data$s_gt * .data$E_eps_1_g_ge / .data$P_gt)
  
  denom_W = weighted.mean(df$om_tilde_1, df$nat_weight, na.rm = TRUE)
  df <- df %>%
    mutate(W = .data$om_tilde_1 / denom_W) %>%
    mutate(weight_result = .data$W * .data$nat_weight) %>%
    select(-.data$eps_1, -.data$P_gt, -.data$om_tilde_1, -.data$E_eps_1_g_ge,
           -.data$E_eps_1_g_ge_aux, -.data$weights_aux, -.data$abs_delta_D, -.data$delta_D)
  
  list(df = df, beta = beta)
}

twowayfeweights_calculate_fds <- function(df, controls) {
  obs <- sum(df$weights)
  gdf <- df %>% group_by(.data$G, .data$T) %>% summarise(P_gt = sum(.data$weights))
  df <- df %>% 
    left_join(gdf, by=c("T", "G")) %>% 
    mutate(P_gt = .data$P_gt / obs)
  
  if (length(controls) == 0) {
    formula = "D ~ 1"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("D ~ ", formula, sep = "")
  }
  formula = paste(formula, " | Tfactor", sep = "")
  denom.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  df$eps_2 <- residuals(denom.lm)
  # df$eps_2 <- df$D - predict(denom.lm, df)
  
  if (length(controls) == 0) {
    formula = "Y ~ D"
  } else {
    formula = paste(controls, collapse = " + ")
    formula = paste("Y ~ D + ", formula, sep = "")
  }
  formula = paste(formula, " | Tfactor", sep = "")
  beta.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  beta <- as.numeric(coef(beta.lm)["D"])
  
  # * Keeping only one observation in each group * period cell
  #   bys `group' `time': gen group_period_unit=(_n==1)	
  # 	drop if group_period_unit==0
  # 	drop group_period_unit
  
  # df <- df %>% 
  #   group_by(.data$G, .data$Tfactor) %>%
  #   summarize(.data$P_gt, .data$nat_weight)
  
  df <- df %>%
    mutate(s_gt = case_when(.data$D > 0 ~ 1,
                            .data$D < 0 ~ -1,
                            TRUE ~ 0)) %>%
    mutate(abs_delta_D = abs(.data$D)) %>%
    mutate(nat_weight = .data$P_gt * .data$abs_delta_D)
  
  P_S = sum(df$nat_weight)
  df <- df %>% 
    mutate(nat_weight = .data$nat_weight / P_S) %>%
    mutate(W = .data$s_gt * .data$eps_2)
  denom_W = weighted.mean(df$W, df$nat_weight, na.rm = TRUE)
  df <- df %>% 
    mutate(W = .data$W / denom_W) %>% 
    mutate(weight_result = .data$W * .data$nat_weight) %>%
    select(-.data$eps_2, -.data$P_gt, -.data$abs_delta_D)
  
  list(df = df, beta = beta)
}

twowayfeweights_summarize_weights <- function(df, var_weight) {
  weight_plus <- df[[var_weight]][df[[var_weight]] > 0 & !is.na(df[[var_weight]])]
  nr_plus <- length(weight_plus)
  sum_plus <- sum(weight_plus, na.rm = TRUE)
  
  weight_minus <- df[[var_weight]][df[[var_weight]] < 0 & !is.na(df[[var_weight]])]
  nr_minus <- length(weight_minus)
  sum_minus <- sum(weight_minus, na.rm = TRUE)
  
  nr_weights <- nr_plus + nr_minus
  
  list(nr_plus = nr_plus, nr_minus = nr_minus, nr_weights = nr_weights, 
       sum_plus = sum_plus, sum_minus = sum_minus)
}

twowayfeweights_test_random_weights <- function(df, random_weights) {
  mat <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(mat) <- c("Coef", "SE", "t-stat", "Correlation")
  df_filtered <- df %>% filter(is.finite(.data$W))

  for (v in random_weights) {
    formula <- sprintf("%s ~ W", v)
    rw.lm = lm_robust(formula =as.formula(formula), data = df_filtered, weights = .data$nat_weight, clusters = .data$G, se_type = "stata")

    beta <- rw.lm$coefficients[["W"]]
    se <- rw.lm$std.error[["W"]]
    r2 <- rw.lm$r.squared
    mat[v, ] <- c(beta, se, beta/se, if (beta > 0) { sqrt(r2) } else { -sqrt(r2) })
    print(mat)
  }
  
  mat
}

twowayfeweights_result <- function(df, beta, random_weights) {
  ret <- twowayfeweights_summarize_weights(df, "weight_result")
  
  W_mean <- weighted.mean(df$W, df$nat_weight)
  W_sd <- sqrt(sum(df$nat_weight * (df$W - W_mean)^2, na.rm = TRUE))
  sensibility <- abs(beta) / W_sd
  
  df_result <- df %>% 
    select(.data$T, .data$G, .data$weight_result) %>% 
    rename(weight = .data$weight_result)
  
  ret$df_result = df_result
  ret$beta = beta
  ret$sensibility = sensibility
  
  if (length(random_weights) > 0) {
    ret$mat = twowayfeweights_test_random_weights(df, random_weights)
  }
  
  if (ret$sum_minus < 0) {
    df_sens <- df %>%
      filter(.data$weight_result != 0) %>%
      arrange(desc(.data$W)) %>% 
      mutate(P_k = 0, S_k = 0, T_k = 0)
    
    N = nrow(df_sens)
    df_sens$P_k[N] <- df_sens$nat_weight[N]
    df_sens$S_k[N] <- df_sens$weight_result[N]
    df_sens$T_k[N] <- df_sens$nat_weight[N] * df_sens$W[N]^2
    for (i in (N-1):1) {
      df_sens$P_k[i] <- df_sens$nat_weight[i] + df_sens$P_k[i+1]
      df_sens$S_k[i] <- df_sens$weight_result[i] + df_sens$S_k[i+1]
      df_sens$T_k[i] <- df_sens$nat_weight[i] * df_sens$W[i]^2 + df_sens$T_k[i+1]
    }
    df_sens <- df_sens %>% 
      mutate(sens_measure2 = abs(beta) / sqrt(.data$T_k + .data$S_k^2 / (1 - .data$P_k))) %>%
      mutate(indicator = as.numeric(.data$W < - .data$S_k / (1 - .data$P_k)))
    df_sens$indicator[1] = 0
    for (i in 2:N) {
      df_sens$indicator[i] = max(df_sens$indicator[i], df_sens$indicator[i-1])
    }
    total_indicator <- sum(df_sens$indicator)
    sensibility2 <- df_sens$sens_measure2[N - total_indicator + 1]
    ret$sensibility2 = sensibility2
  }
  ret
}

twowayfeweights_print_results <- function(cmd_type, r) {
  treat = case_when(cmd_type == "feTR" || cmd_type == "fdTR" ~ "ATT",
                    cmd_type == "feS" || cmd_type == "fdS" ~ "LATE", 
                    TRUE ~ "BLANK")
  assumption = case_when(cmd_type == "feTR" || cmd_type == "fdTR" ~ "Under the common trends assumption",
                         cmd_type == "feS" || cmd_type == "fdS" ~ "Under the common trends, treatment monotonicity, and if groups' treatment effect does not change over time", 
                         TRUE ~ "BLANK")
  
  printf("%s, ", assumption)
  printf("beta estimates a weighted sum of %d %ss.", r$nr_weights, treat)
  printf("%d %s receive a positive weight, and %d receive a negative weight.", r$nr_plus, treat, r$nr_minus) 
  printf("The sum of the positive weights is equal to %f.", r$sum_plus) 
  printf("The sum of the negative weights is equal to %f.", r$sum_minus)
  printf("beta is compatible with a DGP where the average of those %s is equal to 0,", treat)  
  printf("while their standard deviation is equal to %f.", r$sensibility)
  
  if (r$sum_minus < 0) {
    printf("beta is compatible with a DGP where those %s all are of a different sign than beta,", treat)  
    printf("while their standard deviation is equal to %f.", r$sensibility2)
  } else if (r$sum_minus == 0) {
    printf("All the weights are positive, so beta cannot be of a different sign than all those %s.", treat)
  }
  
  r$df_result
}


twowayfeweights_calculate_fetr_other_treatment <- function(df, controls, treatments) {
  mean_D <- weighted.mean(df$D, df$weights, na.rm = TRUE)
  obs <- sum(df$weights)
  gdf <- df %>% group_by(.data$G, .data$T) %>% summarise(P_gt = sum(.data$weights))
  df <- df %>% 
    left_join(gdf, by=c("T", "G")) %>% 
    mutate(P_gt = .data$P_gt / obs) %>% 
    mutate(nat_weight = .data$P_gt * .data$D / mean_D)
  
  vars = c(controls, treatments)
  formula = paste(vars, collapse = " + ")
  formula = paste("D ~ ", formula, sep = "")
  formula = paste(formula, " | G + Tfactor", sep = "")
  denom.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  df$eps_1 <- df$D - predict(denom.lm, df)
  df$eps_1_E_D_gt <- df$eps_1 * df$D
  denom_W <- mean(df$eps_1_E_D_gt, na.rm = TRUE)
  
  df <- df %>% 
    mutate(W = .data$eps_1 * mean_D / denom_W) %>% 
    mutate(weight_result = .data$W * .data$nat_weight)
  
  for (treatment in treatments) {
    varname = fn_treatment_weight_rename(treatment)
    df <- df %>% mutate(!!sym(varname) := .data$W * .data$P_gt * !!sym(treatment) / mean_D)
  }
  
  df <- df %>% select(-.data$eps_1, -.data$P_gt)
  
  formula = paste(vars, collapse = " + ")
  formula = paste("Y ~ D + ", formula, sep = "")
  formula = paste(formula, " | G + Tfactor", sep = "")
  beta.lm <- feols(as.formula(formula), data = df, weights = df$weights)
  beta <- as.numeric(coef(beta.lm)["D"])
  
  list(df = df, beta = beta)
}

twowayfeweights_result_other_treatment <- function(df, treatments, random_weights) {
  columns <- c("T", "G", "weight_result")
  ret <- twowayfeweights_summarize_weights(df, "weight_result")
  
  if (length(random_weights) > 0) {
    ret$mat = twowayfeweights_test_random_weights(df, random_weights)
  }
  
  for (treatment in treatments) {
    varname = fn_treatment_weight_rename(treatment)
    columns <- c(columns, varname)
    ret2 <- twowayfeweights_summarize_weights(df, varname)
    ret[[treatment]] <- ret2
  }
  df_result <- df %>% 
    select_at(vars(columns)) %>% 
    rename(weight = .data$weight_result)
  
  ret[["df_result"]] <- df_result
  ret
}

twowayfeweights_print_result_other_treatment <- function(r, treatments) {
  printf("Under the common trends assumption, beta estimates the sum of several terms.")
  printf("The first term is a weighted sum of %d ATTs of the treatment.", r$nr_weights)
  printf("%d ATTs receive a positive weight, and %d receive a negative weight.", r$nr_plus, r$nr_minus)
  printf("The sum of the positive weights is equal to %f.", r$sum_plus)
  printf("The sum of the negative weights is equal to %f.", r$sum_minus)
  
  for (treatment in treatments) {
    r2 = r[[treatment]]
    printf("The next term is a weighted sum of %d ATTs of treatment %s included in the other_treatments option.", r2$nr_weights, treatment)
    printf("%d ATTs receive a positive weight, and %d receive a negative weight.", r2$nr_plus, r2$nr_minus)
    printf("The sum of the positive weights is equal to %f.", r2$sum_plus)
    printf("The sum of the negative weights is equal to %f.", r2$sum_minus)
  }
  
  r$df_result
}


#' twowayfeweights
#'
#' @param df the data frame for input
#' @param Y the name of Y variable
#' @param G the name of group variable
#' @param T the name of time variable
#' @param D the name of treatment variable
#' @param cmd_type the type of command, including fetr, fdtr, fes, fds
#' @param D0 the name of the mean of the treatment in group g and at period t
#' @param controls the list of names of control variables, empty if not specified
#' @param weights a column of data that replaces the default weight
#' @param other_treatments the list of other treatment variables to include in the regression other than D
#' @param test_random_weights weights when this option is specified, the command estimates the correlation
#'                            between each variable in varlist and the weights
#' @covariance
#' @average_effect
#' @param parallel parallelly perform bootstrap
#'
#' @return twowayfeweights returns data frame that contains the following columns
#'         T: time variable
#'         G: group variable
#'         weight: the result of the weight
#'
#' @export

twowayfeweights <- function(df, Y, G, T, D, cmd_type, D0 = NULL, controls = c(), weights = NULL, other_treatments = c(), test_random_weights = c()) {
  if (!is.null(D0) && cmd_type != "fdTR") {
    printf("Type fdTR requires D0 defined")
    return(c())
  }
  
  if (length(other_treatments) > 0 && cmd_type != "feTR") {
    printf("When the other_treatments option is specified, you need to specify the type(feTR) option.")
    return(c())
  }
  
  controls_rename <- get_controls_rename(controls)
  treatments_rename <- get_treatments_rename(other_treatments)
  random_weight_rename <- get_random_weight_rename(test_random_weights)
  
  df_renamed <- twowayfeweights_rename_var(df, Y, G, T, D, D0, controls, other_treatments, test_random_weights)
  df_transformed <- twowayfeweights_transform(df_renamed, controls_rename, weights, treatments_rename)
  df_filtered <- twowayfeweights_filter(df_transformed, cmd_type, controls_rename, treatments_rename)
  
  if (length(other_treatments) == 0) {
    res <- if (cmd_type == "feTR") {
      twowayfeweights_calculate_fetr(df_filtered, controls_rename)
    } else if (cmd_type == "fdTR") {
      twowayfeweights_calculate_fdtr(df_filtered, controls_rename)
    } else if (cmd_type == "feS") {
      twowayfeweights_calculate_fes(df_filtered, controls_rename)
    } else if (cmd_type == "fdS") {
      twowayfeweights_calculate_fds(df_filtered, controls_rename)
    }
    
    res <- twowayfeweights_result(res$df, res$beta, random_weight_rename)
    df_result <- twowayfeweights_print_results(cmd_type, res)
    df_result
  } else {
    res <- twowayfeweights_calculate_fetr_other_treatment(df_filtered, controls_rename, treatments_rename)
    res <- twowayfeweights_result_other_treatment(res$df, treatments_rename, random_weight_rename)
    df_result <- twowayfeweights_print_result_other_treatment(res, treatments_rename)
    df_result
  }
}