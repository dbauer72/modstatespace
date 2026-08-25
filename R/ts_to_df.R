

#' Convert a [`tsibble`] object to a [`matrix`]
#'
#' The matrix returned contains the variables in columns, while the time column makes up 
#' the rownames. 
#'
#' @param ts_obj the [`tsibble::tsibble`] object
#' @param var_name string indicating the column of ts_obj containing variable names
#' 
#' @return matrix; containing the data in matrix form. 
#'
#' @export
#'
tsibble_to_mat <- function(ts_obj,var_name = "y"){
  
  time_col <- tsibble::index(ts_obj)
  time_stamps <- as.data.frame(ts_obj[,time_col])
  var_col <- tsibble::key(ts_obj)
  vdf <- as.data.frame(ts_obj |> dplyr::distinct( vars))
  uvars <- unique(as.data.frame(ts_obj)[,as.character(var_col)])
  s <- length(uvars)
  T <- length(unique(time_stamps[[1]]))
  
  dat_mat <- matrix(NA,T,s)
  for (i in 1:s){
    df <- as.data.frame( ts_obj |> dplyr::filter(vars == uvars[i]))
    dat_mat[,i] <- df[,var_name]
  }

  # set the row and column names 
  rownames(dat_mat) <- as.character(as.data.frame(ts_obj[1:T,time_col])[[1]])
  colnames(dat_mat) <- uvars
  return(dat_mat)
}
