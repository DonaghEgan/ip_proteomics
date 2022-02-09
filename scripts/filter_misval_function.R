filter_missval_own <- function(protein.df, exp_design.df, thr = 0) {
  
  # Make assay values binary (1 = valid value)
  bin_data <- protein.df
  idx <- is.na(protein.df)
  bin_data[!idx] <- 1
  bin_data[idx] <- 0
  
  # Filter se on the maximum allowed number of
  # missing values per condition (defined by thr)
  keep <- bin_data %>%
    data.frame() %>%
    rownames_to_column() %>%
    gather(ID, value, -rowname) %>%
    left_join(., data.frame(exp_design.df), by = c("ID" = "label")) %>%
    group_by(rowname, condition) %>%
    summarize(miss_val = n() - sum(value)) %>%
    filter(miss_val <= thr) %>%
    spread(condition, miss_val)
  protein_fltrd <- protein.df[keep$rowname, ]
  return(protein_fltrd)
}

test <- filter_missval_own(protein_data, experimental_design) 


