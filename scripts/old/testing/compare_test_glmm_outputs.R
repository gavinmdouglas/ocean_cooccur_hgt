rm(list = ls(all.names = TRUE))

hgt_by_cooccur <- readRDS('~/tmp/glmm_test/hgt_by_cooccur.rds')

hgt_by_cooccur_w_dup_table_and_random <- readRDS('~/tmp/glmm_test/hgt_by_cooccur_w_dup_table_and_random.rds')

hgt_by_cooccur_w_dup_table <- readRDS('~/tmp/glmm_test/hgt_by_cooccur_w_dup_table.rds')

hgt_by_cooccur$logLik
hgt_by_cooccur$coefficients

hgt_by_cooccur_w_dup_table$logLik
hgt_by_cooccur_w_dup_table$coefficients

hgt_by_cooccur_w_dup_table_and_random$logLik
hgt_by_cooccur_w_dup_table_and_random$coefficients
