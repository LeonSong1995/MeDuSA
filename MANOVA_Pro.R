#############################################################################################################
#' @title MANOVA_Pro: Combines multiple analysis of variance (MANOVA) with polynomial regression to detect differences in cell-state abundance among groups.

#' @param MeDuSA_obj The MeDuSA object. 
#' @param degree A numeric variable used to specify the polynomial degrees.
#' @param condition  A character vector containing the biological condition for each bulk sample.

#' @export

MANOVA_Pro <- function(MeDuSA_obj,degree = NULL,condition = NULL){

  ### check the input object
  if (!class(MeDuSA_obj) %in% "MeDuSA_Object"){
    stop("Please input the MeDuSA Object.")
  }
  
  if (is.null(condition)){
    stop("Please provide the biological condition for each sample.")
  }
  
  if (is.null(degree)){
    stop("Please input the polynomial degrees.")
  }

  ### load the estimated results
  abundance_estimate = MeDuSA_obj@Estimation$cell_state_abundance
  if(ncol(abundance_estimate) != length(condition)){
    stop("The number of bulk samples does not match the number of biological conditions.")
  }

  state = MeDuSA_obj@Estimation$TimeBin

  ### conduct the polynomial regression analysis 
  coef = t(lm(abundance_estimate~poly(state,k=degree))$coefficients)
  colnames(coef) = paste0('coef',seq(1,ncol(coef)))

  ### conduct the MANOVA
  out = c(summary(manova(coef~condition))$stats[1,])

  return(out)
}
