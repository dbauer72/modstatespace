

#' Plot impulse responses and covariance sequences from a collection of systems
#'
#' This function plots impulse responses nd covariance sequences for a number of models 
#' selecting the input/output channels to include.
#' The data is assumed to be stored in a tibble containing the columns lag,in_chan,out_chan,syst.
#'
#' @param x     tibble containing the data
#' @param inp   vector of integers containing the numbers of the input channel
#' @param outp  vector of integers containing the numbers of the output channel
#' @param main  title string
#' @param ylab  y-axis label string 
#'
#' @return ggplot2 object 
#'
plot_IR_tibble = function(x, inp=1, outp=1, main='',ylab = '') {

  if (inherits(x,"tbl_df") == FALSE){
    stop('plot_IR_tibble requires a tibble as input')
  }
  
  # filter inputs and outputs
  prod_arrow <- function(a=1,b=1){paste(a,'->',b)}
  x_filt <- x |> dplyr::filter(in_chan %in% inp,out_chan %in% outp) |> dplyr::mutate(in_out = prod_arrow(in_chan,out_chan))
  
  # generate plot 
  p <- ggplot2::ggplot(data = x_filt, mapping = ggplot2::aes(x = lag, y = values, group = Systems, color = Systems)) +
    ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::facet_grid(in_chan ~ out_chan) +
    ggplot2::labs( x = "Lag", y = ylab, title = main)

  # add faceting 
  
  return(p)
}


utils::globalVariables(c("in_chan","out_chan","lag","values","vars","Systems"))