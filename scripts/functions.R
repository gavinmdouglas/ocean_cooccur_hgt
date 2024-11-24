
split_multi_category_rows <- function(in_df, category_col, delimiter = ",", num_cores = 1) {

  multi_category_row_i <- grep(delimiter, in_df[, category_col])

  in_df_nonmulti <- in_df[-multi_category_row_i, , drop = FALSE]
  in_df_multi <- in_df[multi_category_row_i, , drop = FALSE]

  in_df_multi_merged_raw <- parallel::mclapply(1:nrow(in_df_multi),
                                               function(i) {
                                                 row_subset <- in_df_multi[i, , drop = FALSE]
                                                 category_split <- base::strsplit(x = row_subset[, category_col], split = ",")[[1]]
                                                 split_row <- row_subset[rep(1, length(category_split)), , drop = FALSE]
                                                 split_row[, category_col] <- category_split
                                                 return(split_row)
                                               },
                                               mc.cores = num_cores)

  in_df_multi_merged <- do.call(rbind, in_df_multi_merged_raw)

  return(rbind(in_df_nonmulti, in_df_multi_merged))
}


run_glmms_per_grouping <- function(in_tab, subset_col_i, dependent_var_string,
                                   backup1_dependent_var_string,
                                   backup2_dependent_var_string,
                                   convert_to_integer_no_yes=FALSE,
                                   num_cores=1) {

  subset_colnames <- colnames(in_tab)[subset_col_i]

  raw_coef_out <- list()
  glmm_calls <- list()

  for (subset_colname in subset_colnames) {
    message('Working on: ', subset_colname)
    if (convert_to_integer_no_yes) {
      in_tab[, subset_colname] <- as.integer(factor(in_tab[, subset_colname], levels=c('No', 'Yes'))) - 1
    }

    tmp_glmm_out <- tryCatch({

      in_formula <- as.formula(paste(subset_colname, " ~ ", dependent_var_string, sep = ""))

      glmmTMB(formula = in_formula,
              data = in_tab,
              family = "binomial",
              control = glmmTMBControl(optimizer = nlminb,
                                       parallel = num_cores,
                                       profile = TRUE,
                                       optCtrl = list(iter.max = 2000,
                                                      eval.max = 2000)))

    }, condition = function(cond) {

      diff_approach_out <- tryCatch({

        in_formula <- as.formula(paste(subset_colname, " ~ ", backup1_dependent_var_string, sep = ""))

        glmmTMB(formula = in_formula,
                data = in_tab,
                family = "binomial",
                control = glmmTMBControl(optimizer = nlminb,
                                         parallel = num_cores,
                                         profile = TRUE,
                                         optCtrl = list(iter.max = 2000,
                                                        eval.max = 2000)))
      }, condition = function(cond) {

        in_formula <- as.formula(paste(subset_colname, " ~ ", backup2_dependent_var_string, sep = ""))

        glmmTMB(formula = in_formula,
                data = in_tab,
                family = "binomial",
                control = glmmTMBControl(optimizer = nlminb,
                                         parallel = num_cores,
                                         profile = TRUE,
                                         optCtrl = list(iter.max = 2000,
                                                        eval.max = 2000)))
      })

      diff_approach_out
    })

    x_summary <- summary(tmp_glmm_out)
    raw_coef <- data.frame(x_summary$coefficients$cond)
    raw_coef$BH <- p.adjust(raw_coef$Pr...z.., 'BH')
    orig_col <- colnames(raw_coef)
    raw_coef$level <- subset_colname
    raw_coef$variable <- rownames(raw_coef)
    raw_coef <- raw_coef[, c('level', 'variable', orig_col)]
    rownames(raw_coef) <- NULL
    raw_coef_out[[subset_colname]] <- raw_coef
    glmm_calls[[subset_colname]] <- tmp_glmm_out$call
  }

  glmm_coefs <- do.call(rbind, raw_coef_out)
  rownames(glmm_coefs) <- NULL

  return(list(glmm_coefs, glmm_calls))
}

