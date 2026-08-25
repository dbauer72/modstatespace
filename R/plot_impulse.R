
#'
#' plot_IR plots the Impulse Response for a number of systems 
#'
#' @param syst_vec vector of systems or StSp object, if only one system
#' @param channel_in integer: SSM input element of the impulse response to plot. Default is `1`.
#' @param channel_out integer: SSM output element of the impulse response to plot. Default is `1`.
#' @param L integer, number of lags to include. 
#'
#' @return the `ggplot2` object of the plot
#'
#' @export
#'
plot_IR <- function(syst_vec, channel_in=1, channel_out=1, L =20 ){
  
  # put all systems into a vector. 
  if (is.null(names(syst_vec)) == FALSE){
    syst_v <- list(1)
    syst_v[[1]] <- syst_vec
  } else {
    syst_v <- syst_vec
  }
  
  # dimensions 
  M <- length(channel_in)
  N <- length(channel_out)
  S <- length(syst_v)

  # initialise data vectors 
  in_chan <- vector(mode= "numeric",length= M*N*(L+1)*S)
  out_chan <- vector(mode= "numeric",length= M*N*(L+1)*S)
  in_syst <- vector(mode= "numeric",length= M*N*(L+1)*S)
  lag <- vector(mode= "numeric",length= M*N*(L+1)*S)
  values <- vector(mode= "numeric",length= M*N*(L+1)*S)

  ind_in <- as.vector(outer(channel_in,rep(1,N)))
  ind_out <- as.vector(outer(rep(1,M),channel_out))
  
  labs_syst <- list()
  for (s in 1:S){
    labs_syst[[s]] <- paste('Syst',s)
  }
  # cycle over systems 
  for (s in 1:S){
    syst <- syst_v[[s]]
    # calculate IRs 
    IR <- impulse(syst,M=L,plots=0)
   
    # fill into tibble 
    for (l in 1:(L+1)){
      ir_l <- IR[[1]][[l]][channel_out,channel_in]
      ind <- (s-1)*M*N*(L+1)+(l-1)*M*N+c(1:(M*N))
      lag[ind] <- as.integer(l-1)
      in_syst[ind] <- labs_syst[[s]]
      in_chan[ind] <- as.integer(ind_in)
      out_chan[ind] <- as.integer(ind_out)
      values[ind] <- as.vector(ir_l)
    }
  }

  # write into data.frame
  data <- tibble::tibble(lag = lag, Systems = in_syst, in_chan = in_chan, out_chan = out_chan, values = values)
  
  p <- plot_IR_tibble(data, inp=channel_in, outp=channel_out, main = 'Impulse Response', ylab = "IR[i,j]")

  return(p) 
  
} 


  

