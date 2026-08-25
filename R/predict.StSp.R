#' @title prediction function 

#' calculates the predictions 
#' @name predict.StSp
#' @param syst StSp system object used for predictions. 
#' @param y time series data
#' @param h prediction horizon
#' @param T_est_val time point for separating estimation and validation data 
#' 
#' @return predictions
#'
#' @export 
#' @method predict StSp
predict.StSp <- function (syst, y, h=1,T_est_val = 1) {
  
  # convert tsibble to data matrix 
  y_mat <- t(tsibble_to_mat(y)) 
  
  # find point to switch from estimation to validation data. 
  time <- colnames(y_mat)
  time_est <- time[1:T_est_val]
  time_val <- time[T_est_val:dim(y_mat)[2]]
  
  # use Kalman filter to estimate state at latest time point 
  kf_out <- kf(y= y_mat[,time_est], syst= syst, pe=TRUE, init=0)
  xhat <- kf_out$xhat[,length(time_est)]
  Sigt <- kf_out$sigmat[[length(kf_out$sigmat)]]
  Kt <- kf_out$Kt[[length(kf_out$Kt)]]
  Pt <- kf_out$Pt
  At <- syst$A
  Ct <- syst$C

  Abar <- At - Kt %*% Ct
  yhat <- y_mat[,time_val]*0
  resh <- yhat
  res1 <- resh
  
  # advance prediction 
  xhatp1 <- At %*% xhat + Kt %*% kf_out$res[,length(time_est)]
  xhatph <- xhatp1

  # calculate prediction error according to system (NOTE: estimation uncertainty not taken into account)
  Sig_th <- vector("list",length(time_val))
  Sig_th[[1]] <- Sigt + Ct %*% Pt %*% t(Ct) 
  
  Pt <- Abar %*% Pt %*% t(Abar)
  Qt <- Kt %*% Sigt %*% t(Kt)
  CAth <- Ct
  
  yhat[,1] <- Ct %*% xhatp1
  resh[,1]<- y_mat[,time_val[1]] - yhat[,1]
  res1[,1] <- resh[,1]
  if (h>1){
    for (jh in 1:(h-1)){
      Sig_th[[jh+1]] <- Sig_th[[jh]] + CAth %*% Qt %*% t(CAth) + CAth %*% Pt %*% t(CAth) 
      CAth <- CAth %*% At
      Pt <- Abar %*% Pt %*% t(Abar)
      #xhatp1 <- Ft %*% xhatp1 + Kt %*% res1[,jh]
      yhat[,jh+1] <- CAth %*% xhatp1
      resh[,jh+1] <- y_mat[,time_val[jh+1]] - yhat[,jh+1]
      #res1[,jh+1] <- y_mat[,time_val[jh+1]] - 
    }
  }
  # roll over validation sample 
  
  for (th in (h+1):length(time_val)){
    xhatp1 <- At %*% xhatp1 + Kt %*% res1[,th-h]
    yhat[,th] <- CAth %*% xhatp1
    resh[,th] <- y_mat[,time_val[th]]-yhat[,th]
    res1[,th-h+1] <- y_mat[,time_val[th-h+1]]- Ct %*% xhatp1
    Sig_th[[th]] <- Sig_th[[h]]
  }
  
  # provide return values 
  
  return(list(yhat = yhat,res1=res1, resh=resh,Sig_th=Sig_th))
}


#' 
#'  plots the predictions 
#'  
#' @title plot predictions  
#' 
#' @name plot_predict
#' @param pr list containing the predictions as 'resh' and prediction variances as 'Sig_th'.
#' @param window pair of integers specifying the window to plot. 
#' 
#' @return ggplot2 object containing the plot 
#'
#' @export 
#'
plot_predict <- function (pr, window = NULL) {

  # calculate standard errors 
  yhat <- pr$yhat
  res <- pr$resh
  
  Th <- dim(yhat)[2]
  if (is.null(window) == FALSE){
    Th <- diff(window)+1
    time_indices <- seq(from=window[1],to=window[2],by=1)
  } else {
    time_indices <- 1:Th
  }
  
  s <- dim(yhat)[1]
  stds <- matrix(0,s,Th)
  for (j in 1:Th){
    vcov <- pr$Sig_th[[time_indices[j]]]
    stds[,j] <- sqrt(diag(vcov))
  }
  
  # calculate elements of data frame 
  resids <- matrix(0,Th*s,1)
  yhats <- matrix(0,Th*s,1)
  ys <- matrix(0,Th*s,1)
  time <- ys
  comp <- time
  stand_errors <- ys  
      
  # function for converting strings to dates (or integers, or strings)
  convert_time_vec <- function(x) {
    date_converted <- tryCatch({as.Date(x)}, 
                               error = function(e) {
                                  rep(as.Date(NA), length(x))
                              })
    
    # Check if date conversion succeeded entirely without NAs
    if (!any(is.na(date_converted)) && length(x) > 0) {
      # If all elements successfully converted to dates without NA
      return(list(x = date_converted, isdate = TRUE, isnum = FALSE))
    }
    
    int_converted <- suppressWarnings(as.integer(x))
    
    if (!any(is.na(int_converted)) && length(x) > 0) {
      # If all elements successfully converted to integers without NA
      return(list(x = int_converted, isdate = FALSE, isnum = TRUE))
    }
    
    return(list(x = x, isdate = FALSE, isnum = FALSE))
  }
  
  for (js in 1:s){
    # cycle over dimensions 
    ind <- Th*(js-1)+c(1:Th)
    yhats[ind] <- yhat[js,time_indices]
    resids[ind] <- res[js,time_indices]
    conv <-  convert_time_vec(colnames(yhat[,time_indices]))
    if ((conv$isnum == FALSE) && (conv$isdate == FALSE)){
      time[ind] <- ind
    } else {
      time[ind] <- conv$x
    }
    comp[ind] <- rownames(yhat)[js]
    stand_errors[ind] <- stds[js,]
  }
  
  yupper <- yhats + 1.96*stand_errors
  ylower = yhats - 1.96*stand_errors
  ys <- yhats+ resids
  # fill in the data 
  pr_df <- tibble::tibble(data.frame(ys = ys, yhats = yhats, resids = resids, comp= comp, time=time, stds = stand_errors, ylow = ylower, yup = yupper))
  
  # plot data frame using ggplot2 functions 
  if (conv$isdate == TRUE){
    p <- ggplot2::ggplot() +
      ggplot2::geom_line(data = pr_df, mapping = ggplot2::aes(x = as.Date(time), y = ys), color = 'black') +
      ggplot2::geom_line(data = pr_df, mapping = ggplot2::aes(x = as.Date(time), y = yhats), color = 'red') +
      ggplot2::geom_ribbon(data = pr_df, mapping = ggplot2::aes(x = as.Date(time), ymin=ylower, ymax=yupper), alpha=0.2) + 
      ggplot2::facet_grid( ~ comp) + ggplot2::labs( x = "Time", y = "series", title = "Prediction plot") 
  } else {
    p <- ggplot2::ggplot() +
      ggplot2::geom_line(data = pr_df, mapping = ggplot2::aes(x = time, y = ys), color = 'black') +
      ggplot2::geom_line(data = pr_df, mapping = ggplot2::aes(x = time, y = yhats), color = 'red') +
      ggplot2::geom_ribbon(data = pr_df, mapping = ggplot2::aes(x = time, ymin=ylower, ymax=yupper), alpha=0.2) + 
      ggplot2::facet_grid( ~ comp) + ggplot2::labs( x = "Time", y = "series", title = "Prediction plot") 
  }
  plot(p)
  
  # return output 
  return(p)
}

