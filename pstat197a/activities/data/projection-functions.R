projection_fn <- function(.dtm, .prop){
  # coerce feature matrix to sparse
  dtm_mx <- .dtm %>%
    as.matrix() %>%
    as('sparseMatrix')
  
  # compute svd
  svd_out <- sparsesvd(dtm_mx)
  
  # select number of projections
  var_df <- tibble(var = svd_out$d^2) %>%
    mutate(pc = row_number(),
           cumulative = cumsum(var)/sum(var))
  
  n_pc <- which.min(var_df$cumulative < .prop)
  
  # extract loadings
  loadings <- svd_out$v[, 1:n_pc] %>% as.matrix()
  
  # extract scores
  scores <- (dtm_mx %*% svd_out$v[, 1:n_pc]) %>% as.matrix()
  
  # adjust names
  colnames(loadings) <- colnames(scores) <- paste('pc', 1:n_pc, sep = '')
  
  # output
  out <- list(n_pc = n_pc,
              var = var_df,
              projection = loadings,
              data = as_tibble(scores))
  
  return(out)
}

reproject_fn <- function(.dtm, .projection_fn_out){
  as_tibble(as.matrix(.dtm) %*% .projection_fn_out$projection)
}
