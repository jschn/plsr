#'@importFrom stats biplot sd var
#'@importFrom graphics abline barplot hist par plot segments

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Be aware that plsr 0.0.1 contains experimental and partly untested code.\nUse cautiously.")
}

#'Biplot for plsr Objects
#'
#'Produces a biplot from a plsr object
#'
#' @param x The plsr object
#' @param side The side for which the biplot should be generated. Can be "X" (default) to generate
#'    a biplot of the loadings of X onto the latent space or "Y" for the loadings of Y.
#' @param LVs Vector of length two which specifies the latent variables to be plotted against each other.
#'    For example, the default LVs=c(1,2) will plot latent variable 1 against latent variable 2.
#' @param ... optional arguments to be passed to biplot.default.
#' @examples
#' \dontrun{
#' biplot(x)
#' }
#' @export
biplot.plsr = function(x,side="X",LVs=c(1,2),...){
  if (side=="X"){
    V = x$decomposition$V
    LX = x$decomposition$LX
    biplot(LX[,LVs],V[,LVs],xlab = colnames(V[,LVs])[1],ylab = colnames(V[,LVs])[2],...)
  }
  if (side=="Y"){
    U = x$decomposition$U
    LY = x$decomposition$LY
    biplot(LY[,LVs],U[,LVs],xlab = colnames(U[,LVs])[1],ylab = colnames(U[,LVs])[2],...)
  }
}

bootstrap_saliences <- function(data,indices,X_ncol, V) {
  data = data[indices,]

  #data contains both X and Y, so need to separate them again
  X_boot = data[,1:X_ncol]
  Y_ind1 = (X_ncol+1)
  Y_ind2  = ncol(data)
  Y_boot = data[,Y_ind1:Y_ind2]

  R_boot = t(Y_boot)%*%X_boot
  svd_sol=svd(R_boot)
  U_boot = svd_sol$u
  V_boot = svd_sol$v
  D_boot = svd_sol$d

  #following McIntosh and Lobaugh 2004 (https://doi.org/10.1016/j.neuroimage.2004.07.020)
  #calculate procrustes rotation
  svd_boot=svd(t(V)%*%V_boot)
  N = svd_boot$u
  P = svd_boot$v
  Q = N%*%t(P)

  #apply rotation
  V_boot_rot = V_boot %*% diag(D_boot) %*% Q
  U_boot_rot = U_boot %*% diag(D_boot) %*% Q

  #order of linearized matrices: D (just diagonal),U,V
  out = c(sqrt(colSums(V_boot_rot**2)),c(U_boot_rot),c(V_boot_rot))
  return(out)
}

#' Calculates the precision of the p-value estimated by permutation testing
#'
#' Following Ojala&Garriga (2010): "Permutation tests for studying classifier performance"
#'
#' @param p The p value
#' @param k Number of permutation iterations
#' @return The precision given \code{p} and \code{k}.
#' @examples
#' permutation_precision(0.05,1000)
#' @export
permutation_precision = function(p,k){
  return(sqrt((p*(1-p)/k)))
}


#' Calculate variance explained for plsr object
#'
#' Calculates explained variance per component of the original data sets X and Y.
#'
#' @param plsr_obj A plsr object
#' @return A list containing the elements ExpVarX and ExpVarY, which contain the explained variances for X and Y respectively
#' @examples
#' \dontrun{
#' explained_variance(plsr_object)
#' }
explained_variance=function(plsr_obj){
  Y= plsr_obj$orig_data$Y
  X = plsr_obj$orig_data$X
  LX= plsr_obj$decomposition$LX
  LY= plsr_obj$decomposition$LY
  U= plsr_obj$decomposition$U
  V= plsr_obj$decomposition$V
  vars_y = matrix(NA, ncol=ncol(LY), nrow=1)
  vars_x = matrix(NA, ncol=ncol(LX), nrow=1)
  Var_X = apply(X,2,var)
  Var_Y = apply(Y,2,var)
  for (c in 1:ncol(LX)){
    X_hat=LX[,1:c]%*%t(V[,1:c]) # original X calculated from component c
    X_res= X-X_hat #residual of X-X_hat
    Var_X_res = apply(X_res,2,var) # variance per column of residual

    Y_hat=LY[,1:c]%*%t(U[,1:c])
    Y_res = Y-Y_hat
    Var_Y_res = apply(Y_res,2,var)

    var_exp_x = (sum(Var_X) - sum(Var_X_res))/sum(Var_X)
    var_exp_y = (sum(Var_Y) - sum(Var_Y_res))/sum(Var_Y)

    vars_x[1,c] = var_exp_x
    vars_y[1,c] = var_exp_y
  }
  return(list(ExpVarX=vars_x,ExpVarY=vars_y))
}

#' Constructor for plsr objects
#'
#' @param decomp List of singular value decomposition results
#' @param perm List of permutation testing results
#' @param bootstrp List of bootstrapping results
#' @param sclng List of scaling paramters applied to original data
#' @param org_dat List of original data
#' @param cl Call of pls function
#' @examples
#' \dontrun{
#' plsr_obj=new_plsr(d,p,b,s,o,c)
#' }
#' @export
new_plsr=function(decomp=list(),perm=list(), bootstrp=list(), sclng=list(),org_dat=list(),cl=list()){
  #TODO: stopifnots here?
  structure(list(decomposition=decomp, permutation=perm, bootstrapping=bootstrp, scaling=sclng,orig_data=org_dat,call=cl),class="plsr")
}

#' Print loadings of plsr object
#'
#' This will print the loading matrices V and U that project from original data spaces X and Y to latent space.
#' @param x A plsr object
#' @param mat Which matrix to print (U or V), if NULL (default) will print both
#' @examples
#' \dontrun{
#' loadings(x)
#' loadings(x,"U")
#' loadings(x,"V")
#' }
#' @export
loadings.plsr=function(x,mat=NULL){
  if (mat=='U' | is.null(mat)){
    cat("Loading matrix for Y-side (U)\n")
    print(x$decomposition$U)
  }
  if (mat=='V' | is.null(mat)){
    cat("Loading matrix for X-side (V)\n")
    print(x$decomposition$V)
  }
}

#'  Plot function for plsr objects
#'
#'  Plots information about a plsr object. The following plots will be generated:
#' \itemize{
#' \item barplot of p-values of latent variables estimated via permuatation testing
#' \item Histograms of the distributions of latent variables derived via permutation testing
#' \item A plot showing the effect of the first latent variable on the original data spaces
#' \item Several plots to visualize bootstrapping results
#' }
#'
#' @param x The plsr object.
#' @param ... Further arguments.
#' @examples
#' \dontrun{
#' plot(x)
#' }
#'
#' @export
plot.plsr = function(x, ...){
  plot_perm_results(x)
  invisible(readline(prompt="Press [enter] for next plot"))
  plot_perm_distr(x)
  invisible(readline(prompt="Press [enter] for next plot"))
  plot_latent_variables(x)
  invisible(readline(prompt="Press [enter] for next plot"))
  plot_boot_results(x)
  invisible(readline(prompt="Press [enter] for next plot"))
}

plot_boot_results = function(plsr_obj, sig_threshold=1.96){
  warning("Bootstrapping functionality is still untested!")
  #saliences
  U = plsr_obj$decomposition$U
  V = plsr_obj$decomposition$V

  #standard errors
  U_SE = plsr_obj$bootstrapping$U_SE
  V_SE = plsr_obj$bootstrapping$V_SE

  #salience elements in standard error units
  U_in_SE = U/U_SE
  V_in_SE = V/V_SE

  molten_U_in_SE = reshape2::melt(U_in_SE)
  molten_U_in_SE$Var2=factor(molten_U_in_SE$Var2,levels = rev(levels(factor(molten_U_in_SE$Var2))))
  molten_U_in_SE_sig = molten_U_in_SE
  molten_U_in_SE_sig$value = molten_U_in_SE_sig$value > sig_threshold


  molten_V_in_SE = reshape2::melt(V_in_SE)
  molten_V_in_SE$Var2=factor(molten_V_in_SE$Var2,levels=rev(levels(factor(molten_V_in_SE$Var2))))
  molten_V_in_SE_sig = molten_V_in_SE
  molten_V_in_SE_sig$value = molten_V_in_SE_sig$value > sig_threshold

  print(ggplot2::ggplot(data=molten_U_in_SE,ggplot2::aes_string(x="Var1",y="Var2",fill="value"))+ggplot2::geom_tile(color="black")+ggplot2::scale_fill_continuous()+
    ggplot2::ggtitle("Standardized U"))
  invisible(readline(prompt="Press [enter] for next plot"))
  print(ggplot2::ggplot(data=molten_V_in_SE,ggplot2::aes_string(x="Var1",y="Var2",fill="value"))+ggplot2::geom_tile(color="black")+ggplot2::scale_fill_continuous()+
    ggplot2::ggtitle("Standardized V"))
  invisible(readline(prompt="Press [enter] for next plot"))
  print(ggplot2::ggplot(data=molten_U_in_SE_sig,ggplot2::aes_string(x="Var1",y="Var2",fill="value"))+ggplot2::geom_tile(color="black")+ggplot2::scale_fill_discrete()+
    ggplot2::ggtitle("Significant U Elements"))
  invisible(readline(prompt="Press [enter] for next plot"))
  print(ggplot2::ggplot(data=molten_V_in_SE_sig,ggplot2::aes_string(x="Var1",y="Var2",fill="value"))+ggplot2::geom_tile(color="black")+ggplot2::scale_fill_discrete()+
    ggplot2::ggtitle("Significant V Elements"))
}


#' Plot explained variance of plsr object
#'
#' Calculates and plots the variance explained in the original data X and Y by each additional latent variable.
#'
#' @param plsr_obj The plsr object.
#'
#' @examples
#' \dontrun{plot_explained_variance(plsr_obj)}
#'
plot_explained_variance=function(plsr_obj){
  exp_list = explained_variance(plsr_obj)
  par(mfrow = c(1,2))
  barplot(exp_list$ExpVarX, main = "Explained Variance of X per number of LVs",
          names.arg = 1:length(exp_list$ExpVarX), xlab = "Number of LVs", ylab = "Explained Variance")
  barplot(exp_list$ExpVarY, main = "Explained Variance of Y per number of LVs",
          names.arg = 1:length(exp_list$ExpVarX), xlab = "Number of LVs", ylab = "Explained Variance")

  par(mfrow = c(1,1))
}

#' Plots latent variables
#'
#' This function will plot the effects of increasing and decreasing one or several latent variables by the specified standard deviation.
#'
#' @param plsr_obj A plsr object
#' @param lv_num An integer or list of integer specifying which latent variables to plot.
#' @param sd Range in standard deviations from +[sd] to -[sd].
#' @param FUN A vector containing two functions, which will be used for plotting. Default is c(barplot,barplot).
#' @param args1 Arguments for the plotting function in \code{FUN[1]}
#' @param args2 Arguments for the plotting function in \code{FUN[2]}
#' @examples
#' \dontrun{
#' plot_latent_variables(plsr_obj)
#' plot_latent_variables(plsr_obj,lv=1:2, sd=2, FUN=c(customplot,barplot))
#' }
#'
#' @export
plot_latent_variables = function(plsr_obj, lv_num=1, sd=3, FUN=c(barplot,barplot), args1=NULL, args2=NULL){
  U = plsr_obj$decomposition$U
  V = plsr_obj$decomposition$V
  D = plsr_obj$decomposition$D
  scaling = plsr_obj$scaling

  steps = c(-sd,0,sd)
  num_steps = length(steps)
  old_setting = par()$mfrow
  par(mfrow=c(2,num_steps))
  F1 = FUN[[1]]
  F2 = FUN[[2]]

  for (i in 1:num_steps){
    plot_title = paste(steps[i],"SDs")
    Y_side=rowSums(cbind((U%*%sqrt(D)*steps[i])[,lv_num]))*scaling$Y_scale+scaling$Y_mean
    args = c(list(Y_side,main=plot_title),args1)
    do.call(F1,args)
  }
  for (i in 1:num_steps){
    X_side = rowSums(cbind((V*steps[i])[,lv_num]))*scaling$X_scale+scaling$X_mean
    args = c(list(X_side),args2)
    do.call(F2,args)
  }
  par(mfrow=old_setting)
}

#' Plot permuation results for plsr object
#'
#' Plots the p-values for the latent variables estimated through permutation testing.
#'
#' @param plsr_obj A plsr_obj.
#' @param ... Additional arguments passed to \code{barplot}.
#' @param alpha The significance threshold used. Will be indicated in the plot by a horizontal line.
#'    If NULL (default), the alpha value of the plsr object will be used.
#' @param main The title of the plot.
#' @param lwd The line width of the line indicating alpha.
#' @param col The color of the line indicating alpha.
#'
#' @examples
#' \dontrun{
#' plot_perm_results(plsr_obj)
#' }
#' @export
plot_perm_results=function(plsr_obj,...,alpha = NULL,main = "Permutation Testing Results",lwd=2,col="red"){
  #TODO: maybe always limit to 0 to 1. When all p values are low or high, you cannot see any difference and also not the significance line
  lv_names = paste("LV", 1:nrow(plsr_obj$decomposition$D))
  barplot(plsr_obj$permutation$p_values,names.arg = lv_names,main=main,ylab="p-value",...)
  if (is.null(alpha)) abline(h=plsr_obj$permutation$alpha,lwd=lwd, col=col) else abline(h=alpha,lwd=lwd, col=col)
}

#' Plots null distributions constructed via permutation testing
#'
#' Plots histograms of the null distribution for values of singular values of latent variables constructed via permutation testing.
#' @param plsr_obj A plsr object.
#' @param ... Further parameters to be passed to \code{hist}.
#' @param lwd Line width of vertical line indicating the estimated value of the singular value.
#' @param bar_col Color of the bars in the histograms.
#' @param line_col Color of the vertial line indicating the estimated value of the singular value.
#' @examples
#' \dontrun{plot_perm_distr(plsr_obj)}
#' @export
plot_perm_distr = function(plsr_obj,..., lwd=2,bar_col="grey", line_col="red"){
  n_comp = ncol(plsr_obj$decomposition$D)
  plot_dim=ceiling(sqrt(n_comp))
  old_setting = par()$mfrow
  par(mfrow=c(plot_dim,plot_dim))

  for (i in 1:n_comp){
    #ensure that the real singular value can always be drawn in the plot by setting the xlim parameter accordingly
    xlim_range = range(plsr_obj$permutation$D_perm[,i])
    if (xlim_range[2]<plsr_obj$decomposition$D[i,i]){
      xlim_range[2]= plsr_obj$decomposition$D[i,i]+1
    }
    hist_title = paste("Null Distribution of Singular Value", i)
    xlab_name = paste("Values of Singular Value", i)
    hist(plsr_obj$permutation$D_perm[,i],...,main = hist_title, xlim = xlim_range, xlab = xlab_name, col = bar_col)
    abline(v=plsr_obj$decomposition$D[i,i],lwd=lwd, col=line_col)
  }
  par(mfrow=old_setting)
}

#' Print plsr object
#'
#' Prints information about a plsr object.
#'
#' @param x A plsr object.
#' @param ... Further arguments.
#'
#' @export
print.plsr=function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("X:",ncol(x$orig_data$X), "variables and", nrow(x$orig_data$X), "observations")
  cat("\n")
  cat("Y:",ncol(x$orig_data$Y), "variables and", nrow(x$orig_data$Y), "observations")
  cat("\n")
}

# TODO: check math of backward prediction
# forward(backward(x)) = x should hold
#' Predict from a plsr object
#'
#' This function can be used to make predictions from one original data space to the other.
#' Prediction direction can be forward, meaning X to Y direction and backward, meaning Y to X prediction.
#'
#' @param object A plsr object.
#' @param new_data The data from which you want to predict.
#' @param direction The direction of prediction. Default is "forward" meaning X to Y. Every other argument will result in backward prediction.
#' @param ... Additional arguments.
#'
#' @examples
#' \dontrun{predict(plsr_obj,rnorm(10),"forward")}
#' @export
predict.plsr=function(object,new_data,direction="forward",...){
  if (class(new_data)=="data.frame"){
    new_data=as.matrix(new_data)
  }

  if ((class(new_data)!="matrix")){
    #if new_data is neither matrix nor data.frame then assume it's a vector and turn it into a 1 by length(vector) matrix
    new_data = matrix(new_data, nrow=1)
  }
  U=object$decomposition$U #matrix for projection of tracking data into latent space
  D=object$decomposition$D #covariance of latent components on diagonal
  V=object$decomposition$V #matrix for projection of ratings into latent space
  scaling = object$scaling

  if (direction=="forward"){#X to Y
    x_centered = sweep(new_data,2,scaling$X_mean,"-")
    x_scaled = sweep(x_centered,2,scaling$X_scale,"/")
    x_v_space = x_scaled%*%V #vector in latent space

    pred = U%*%sqrt(D)%*%t(x_v_space)
    pred_scaled_back = pred*scaling$Y_scale+scaling$Y_mean

    return(pred_scaled_back)
  }
  else{
    warning("Backward prediction is untested and might give incorrect results!")
    y_centered = sweep(new_data,2,scaling$Y_mean,"-")
    y_scaled = sweep(y_centered,2,scaling$Y_scale,"/")
    y_u_space = y_scaled%*%U #vector in latent space, TODO: this should be transposed (maybe?)

    pred = y_u_space%*%sqrt(D)%*%t(V)
    pred_scaled_back = pred*scaling$X_scale+scaling$X_mean

    return(pred_scaled_back)
  }

}

#TODO: options to discard results from permutation and bootstrap steps to keep plsr objects small

#' Run partial least squares analysis
#'
#' This is the main function of the plsr package. It will calculate a partial least squares solution
#' for the provided data and perform permutation testing and bootstrapping on the resulting latent variables.
#' Results will be saved as a plsr object.
#'
#' @param X A matrix of m observations on n_x variables.
#' @param Y A matrix of m observations on n_y dimensions.
#' @param n_perm Number of permutation iterations. Default is 100.
#' @param n_boot Number of bootstrap iterations. Default is 100.
#' @param scale Scaling of X and Y (Boolean).
#' @param verbose Provides additional output.
#' @param alpha The significance level for permutation testing.
#' @return A plsr Object.
#' @examples
#' X=matrix(rnorm(300), ncol = 3)
#' Y=matrix(rnorm(1000), ncol = 10)
#' pls(X,Y)
#' pls(X,Y, n_perm = 10, n_boot = 10)
#' @export
pls = function(X,Y,n_perm=100,n_boot=100, scale=T, verbose=F, alpha=0.05){
  #TODO: should this also work if X or Y only have one dimensions? IF so, need to make it work
  if (n_perm<10) n_perm=10
  if (n_boot<10) n_boot=10

  #standardize matrices
  if (scale){
    X = scale(X)
    Y = scale(Y)
  }
  #SVD
  R = t(Y)%*%X
  svd_sol=svd(R)
  U = svd_sol$u
  V = svd_sol$v
  D = svd_sol$d

  #name rows and columns of V
  rownames(V)=colnames(X)
  colnames(V) = paste("LV",1:ncol(V),sep="")

  #Projection of original values into singular vector space: scores
  LX = X%*%V
  LY = Y%*%U

  D_perm_vals = matrix(NA,nrow = n_perm, ncol = length(D)) #will store singular values obtained during permutation


  permutation_progress <- utils::txtProgressBar(min = 0, max = n_perm, style = 2)
  cat("Permuting...\n\n")
  #permutation loop
  for(p in 1:n_perm){

    utils::setTxtProgressBar(permutation_progress,p)

    #let's permute X
    X_perm = X[sample(nrow(X)),]
    #SVD
    R_perm = t(Y)%*%X_perm
    svd_sol=svd(R_perm)
    U_perm = svd_sol$u
    V_perm = svd_sol$v
    D_perm = svd_sol$d

    #following McIntosh and Lobaugh 2004 (https://doi.org/10.1016/j.neuroimage.2004.07.020)
    #calculate procrustes rotation to align solutions
    svd_perm=svd(t(V)%*%V_perm)
    N = svd_perm$u
    P = svd_perm$v
    Q = N%*%t(P)

    #apply rotation
    V_perm_rot = V_perm %*% diag(D_perm) %*% Q
    U_perm_rot = U_perm %*% diag(D_perm) %*% Q

    #now need to calculate the singular values from rotated matrix again and store in D_perm_vals
    #"can be obtained by calculating the columnwise square root of the sums-of-squares"
    D_perm_vals[p,]=sqrt(colSums(V_perm_rot**2))
  }

  close(permutation_progress)
  cat("Done!\n")

  #calculate p values here: percentage of singular value distribution that is bigger than our obtained svs
  p_vals = c()
  for (v in 1:length(D)){
    #use conservative estimate, avoid p=0
    new_pval=(sum(D_perm_vals[,v]>=D[v])+1)/(n_perm+1)
    p_vals=c(p_vals,new_pval)
  }
  #bootstrapping with boot package
  cat("Bootstrapping...\n")
  warning("Bootstrapping functionality is still untested! No guarantee for correctness! Use with extra care.")
  boot_results = boot::boot(data = cbind(X,Y),statistic = bootstrap_saliences, R=n_boot,X_ncol=ncol(X) , V=V)
  cat("Done!\n")
  #TODO: make this pretty
  #separating results for the three matrices again
  D_boot_vals = boot_results$t[,1:length(D)]
  U_boot_vals = boot_results$t[,(length(D)+1):((nrow(U)*ncol(U))+length(D))]
  V_boot_vals = boot_results$t[,(((nrow(U)*ncol(U))+length(D))+1):(ncol(boot_results$t))]

  sv_se = c()
  for (v in 1:ncol(D_boot_vals)){
    new_sv_se=(sqrt((n_boot-1)/n_boot) * sd(D_boot_vals[,v]))#using uncorrected sd here because we have access to the population
    sv_se=c(sv_se,new_sv_se)
  }

  u_se = rep(NA, ncol(U)*nrow(U))
  for (v in 1:ncol(U_boot_vals)){
    u_se[v]=(sqrt((n_boot-1)/n_boot) * sd(U_boot_vals[,v])) #using uncorrected sd here because we have access to the population
  }

  v_se = rep(NA, ncol(V)*nrow(V))
  for (v in 1:ncol(V_boot_vals)){
    v_se[v]=(sqrt((n_boot-1)/n_boot) * sd(V_boot_vals[,v])) #using uncorrected sd here because we have access to the population
  }

  #reshape linearized V and U standard error matrices again
  V_SE = matrix(v_se,nrow=nrow(V))
  U_SE = matrix(u_se,nrow=nrow(U))

  #organize output
  precision = permutation_precision(p_vals,n_perm)
  decomp = list(U=U,V=V,D=diag(D),LX=LX,LY=LY)
  perm = list(n_perm=n_perm,D_perm= D_perm_vals, p_values=p_vals, alpha=alpha, precision=precision)
  bootstrp = list(n_boot=n_boot,sv_se=sv_se, D_boot= D_boot_vals, U_boot= U_boot_vals,V_boot= V_boot_vals, U_SE = U_SE,V_SE=V_SE)

  X_mean = attr(X,"scaled:center")
  X_scale = attr(X,"scaled:scale")
  Y_mean = attr(Y, "scaled:center")
  Y_scale = attr(Y, "scaled:scale")

  sclng = list(X_mean=X_mean,X_scale=X_scale,Y_mean=Y_mean,Y_scale=Y_scale)
  org_dat = list(X=X,Y=Y)
  output = new_plsr(decomp, perm, bootstrp, sclng,org_dat, match.call())

  if (sum(p_vals<alpha)!=sum((p_vals+precision)<alpha)) warning("Some p-values are not stable under this precision.\n Try to increase the number of permutation steps (n_perm).")
  if (sum(p_vals<alpha)!=sum((p_vals-precision)<alpha)) warning("Some p-values are not stable under this precision.\n Try to increase the number of permutation steps (n_perm).")

  return(output)
}

#' Summary of plsr object
#'
#' @param object A plsr object.
#' @param ... Further arguments.
#' @examples
#' \dontrun{
#' summary(plsr_object)
#' }
#' @export
summary.plsr=function(object, ...){
  cat("Permutation iterations:", object$permutation$n_perm,"\n")
  cat("P_values:\n")
  print(object$permutation$p_values)
  cat("\n")
  cat("Bootstrap iterations:", object$bootstrapping$n_boot, "\n")
  cat("\n")
  if(dim(object$decomposition$U)[1]<11 && dim(object$decomposition$U)[2]<11){
    cat("Loading matrix for Y-side (U)\n")
    print(object$decomposition$U)
  }
  else{
    cat("Loading matrix for Y-side (U) too big...ommited\n\n")
  }
  if(dim(object$decomposition$V)[1]<11 && dim(object$decomposition$V)[2]<11){
    cat("Loading matrix for X-side (V)\n\n")
    print(object$decomposition$V)
  }
  else{
    cat("Loading matrix for X-side (V) too big...ommited\n")
  }
}
