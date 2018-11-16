biplot.plsr = function(plsr_obj,direction = "forward",...){
  V = plsr_obj$decomposition$V
  LX = plsr_obj$decomposition$LX
  biplot(LX,V,...)
}

bootstrap_saliences <- function(data,indices,X_ncol, V) {
  #TODO: find out if I can pass things by reference: would be good for X,Y,U and V
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

#' Calculates explained variance per component of the original data sets X and Y
#'
#' @param plsr_obj A plsr object
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


new_plsr=function(decomp=list(),perm=list(), bootstrp=list(), sclng=list(),org_dat=list(),cl=list()){
  #TODO: stopifnots here?

  structure(list(decomposition=decomp, permutation=perm, bootstrapping=bootstrp, scaling=sclng,orig_data=org_dat,call=cl),class="plsr")
}

loadings.plsr=function(plsr_obj){
  print("this should display loadings")
}

plot.plsr = function(plsr_obj){
  plot_perm_results(plsr_obj)
  invisible(readline(prompt="Press [enter] for next plot"))
  plot_latent_variables(plsr_obj)
  invisible(readline(prompt="Press [enter] for next plot"))
  plot_boot_results(plsr_obj)
  invisible(readline(prompt="Press [enter] for next plot"))
}

plot_boot_results = function(plsr_obj, sig_threshold=1.96){
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

  print(ggplot2::ggplot(data=molten_U_in_SE,ggplot2::aes(x=Var1,y=Var2,fill=value))+ggplot2::geom_tile(color="black")+ggplot2::scale_fill_continuous()+
    ggplot2::ggtitle("Standardized U"))
  invisible(readline(prompt="Press [enter] for next plot"))
  print(ggplot2::ggplot(data=molten_V_in_SE,ggplot2::aes(x=Var1,y=Var2,fill=value))+ggplot2::geom_tile(color="black")+ggplot2::scale_fill_continuous()+
    ggplot2::ggtitle("Standardized V"))
  invisible(readline(prompt="Press [enter] for next plot"))
  print(ggplot2::ggplot(data=molten_U_in_SE_sig,ggplot2::aes(x=Var1,y=Var2,fill=value))+ggplot2::geom_tile(color="black")+ggplot2::scale_fill_discrete()+
    ggplot2::ggtitle("Significant U Elements"))
  invisible(readline(prompt="Press [enter] for next plot"))
  print(ggplot2::ggplot(data=molten_V_in_SE_sig,ggplot2::aes(x=Var1,y=Var2,fill=value))+ggplot2::geom_tile(color="black")+ggplot2::scale_fill_discrete()+
    ggplot2::ggtitle("Significant V Elements"))
}

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
#' @param plsr_obj A plsr object
#' @param lv_num An integer or list of integer specifying which latent variables to plot
#' @param sd Range in standard deviations from +[sd] to -[sd]
#' @param frame Which timestep to plot
plot_latent_variables = function(plsr_obj,lv_num=1,sd=3,frame=1){
  U = plsr_obj$decomposition$U
  V = plsr_obj$decomposition$V
  D = plsr_obj$decomposition$D
  scaling = plsr_obj$scaling

  #steps = -sd:sd
  steps = c(-sd,0,sd)
  num_steps = length(steps)
  old_setting = par()$mfrow
  par(mfrow=c(2,num_steps))

  for (i in 1:num_steps){
    plot_title = paste(steps[i],"SDs")
    plot_frame(rowSums(cbind((U%*%sqrt(D)*steps[i])[,lv_num]))*scaling$Y_scale+scaling$Y_mean,lim=F,title=plot_title,single_frame = frame)
    #plot_frame(rowSums(cbind((U*steps[i])[,lv_num]))*scaling$Y_scale+scaling$Y_mean,lim=F,title=plot_title,single_frame = frame)
  }
  for (i in 1:num_steps){
    #barplot(rowSums(cbind((V%*%sqrt(D)*steps[i])[,lv_num]))*scaling$X_scale+scaling$X_mean,cex.axis = 1.5,font=2)
    barplot(rowSums(cbind((V*steps[i])[,lv_num]))*scaling$X_scale+scaling$X_mean,cex.axis = 1.5,font=2)
  }

  par(mfrow=old_setting)


}

plot_perm_results=function(plsr_obj,sig_level = 0.05){

  #TODO: maybe always limit to 0 to 1. When all p values are low or high, you cannot see andy difference and also not the significance line
  lv_names = paste("lv", 1:nrow(plsr_obj$decomposition$D))
  barplot(plsr_obj$permutation$p_values,names.arg = lv_names, main = "Permutation Testing Results")
  abline(h=sig_level, col="red")
}

print.plsr=function(plsr_obj){
  cat("\n")
  cat("plsr object\n") #TODO: dimensions of orig data
  cat("\n")
  cat("Call:\n")
  print(plsr_obj$call)
  cat("\n")
}

#TODO: stimmt noch nicht ganz
#forward(backward(x)) = x sollte gelten. Tut es aber nicht
predict.plsr=function(plsr_obj,new_data,direction="forward"){
  if (class(new_data)=="data.frame"){
    new_data=as.matrix(new_data)
  }

  if ((class(new_data)!="matrix")){
    #if new_data is neither matrix nor data.frame then assume it's a vector and turn it into a 1 by length(vector) matrix
    new_data = matrix(new_data, nrow=1)
  }
  U=plsr_obj$decomposition$U #matrix for projection of tracking data into latent space
  D=plsr_obj$decomposition$D #covariance of latent components on diagonal
  V=plsr_obj$decomposition$V #matrix for projection of ratings into latent space
  scaling = plsr_obj$scaling


  if (direction=="forward"){#X to Y
    x_centered = sweep(new_data,2,scaling$X_mean,"-")
    x_scaled = sweep(x_centered,2,scaling$X_scale,"/")
    x_v_space = x_scaled%*%V #vector in latent space

    pred = U%*%sqrt(D)%*%t(x_v_space)
    pred_scaled_back = pred*scaling$Y_scale+scaling$Y_mean

    return(pred_scaled_back)
  }
  else{
    y_centered = sweep(new_data,2,scaling$Y_mean,"-")
    y_scaled = sweep(y_centered,2,scaling$Y_scale,"/")
    y_u_space = y_scaled%*%U #vector in latent space

    pred = y_u_space%*%sqrt(D)%*%t(V)
    pred_scaled_back = pred*scaling$X_scale+scaling$X_mean

    return(pred_scaled_back)
  }

}



#TODO: options to discard results from permutation and bootstrap steps
#TODO: need to provide examples for doc

#' Run partial least square analysis
#'
#' @param X A matrix of m observations on n dimensions
#' @param Y A matrix of m observations on n dimensions
#' @param n_perm Number of permutation iterations
#' @param n_boot Number of bootstrap iterations
#' @param scale Scaling of X and Y (Boolean)
#' @param verbose Provides additional output
#' @return A PLS Object
pls = function(X,Y,n_perm=10,n_boot=10, scale=T, verbose=F){
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


  permutation_progress <- txtProgressBar(min = 0, max = n_perm, style = 2)
  cat("Permuting...\n\n")
  #permutation loop
  for(p in 1:n_perm){

    setTxtProgressBar(permutation_progress,p)

    #let's permute X
    X_perm = X[sample(nrow(X)),]
    #SVD
    R_perm = t(Y)%*%X_perm
    svd_sol=svd(R_perm)
    U_perm = svd_sol$u
    V_perm = svd_sol$v
    D_perm = svd_sol$d

    #calculate procrustes rotation
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
  #TODO: answer this
  #now do I compare only the first LV to the first column of D_perm and the second LV to the second column..or what?

  #calculate p values here: percentage of singular value distribution that is bigger than our obtained svs
  p_vals = c()
  for (v in 1:length(D)){
    #use conservative estimate, avoid p=0
    new_pval=(sum(D_perm_vals[,v]>=D[v])+1)/(n_perm+1)
    p_vals=c(p_vals,new_pval)
  }
  #bootstrapping with boot package
  cat("Bootstrapping...\n")
  boot_results = boot::boot(data = cbind(X,Y),statistic = bootstrap_saliences, R=n_boot,X_ncol=ncol(X) , V=V)
  cat("Done!\n")
  #TODO: make this pretty
  #separating results for the three matrices again
  D_boot_vals = boot_results$t[,1:length(D)]
  U_boot_vals = boot_results$t[,(length(D)+1):((nrow(U)*ncol(U))+length(D))]
  V_boot_vals = boot_results$t[,(((nrow(U)*ncol(U))+length(D))+1):(ncol(boot_results$t))]


  sv_se = c()
  for (v in 1:ncol(D_boot_vals)){
    new_sv_se=(sqrt((n_boot-1)/n_boot) * sd(D_boot_vals[,v]))#using uncorrected sd here because we have access to the population, correct?
    sv_se=c(sv_se,new_sv_se)
  }

  u_se = rep(NA, ncol(U)*nrow(U))
  for (v in 1:ncol(U_boot_vals)){
    u_se[v]=(sqrt((n_boot-1)/n_boot) * sd(U_boot_vals[,v])) #using uncorrected sd here because we have access to the population, correct?
  }

  v_se = rep(NA, ncol(V)*nrow(V))
  for (v in 1:ncol(V_boot_vals)){
    v_se[v]=(sqrt((n_boot-1)/n_boot) * sd(V_boot_vals[,v])) #using uncorrected sd here because we have access to the population, correct?
  }

  #reshape linearized V and U standard error matrices again
  V_SE = matrix(v_se,nrow=nrow(V))
  U_SE = matrix(u_se,nrow=nrow(U))

  #organize output
  decomp = list(U=U,V=V,D=diag(D),LX=LX,LY=LY)
  perm = list(D_perm= D_perm_vals, p_values=p_vals)
  bootstrp = list(sv_se=sv_se, D_boot= D_boot_vals, U_boot= U_boot_vals,V_boot= V_boot_vals, U_SE = U_SE,V_SE=V_SE)

  X_mean = attr(X,"scaled:center")
  X_scale = attr(X,"scaled:scale")
  Y_mean = attr(Y, "scaled:center")
  Y_scale = attr(Y, "scaled:scale")

  sclng = list(X_mean=X_mean,X_scale=X_scale,Y_mean=Y_mean,Y_scale=Y_scale)
  org_dat = list(X=X,Y=Y)
  output = new_plsr(decomp, perm, bootstrp, sclng,org_dat, match.call())

  return(output)
}

summary.plsr=function(plsr_obj){
  cat("\n")
  cat("plsr object\n\n") #TODO: dimensions of orig data
  cat("Permutation iterations:", plsr_obj$call$n_perm,"\n")#TODO: save this as actual value
  cat("P_values:\n")
  print(plsr_obj$permutation$p_values)
  cat("\n")
  cat("Bootstrap iterations:", plsr_obj$call$n_boot, "\n")
  cat("\n")
  if(dim(plsr_obj$decomposition$U)[1]<11 && dim(plsr_obj$decomposition$U)[2]<11){
    cat("Loading matrix for Y-side (U)\n")
    print(plsr_obj$decomposition$U)
  }
  else{
    cat("Loading matrix for Y-side (U) too big...ommited\n\n")
  }
  if(dim(plsr_obj$decomposition$V)[1]<11 && dim(plsr_obj$decomposition$V)[2]<11){
    cat("Loading matrix for X-side (V)\n\n")
    print(plsr_obj$decomposition$V)
  }
  else{
    cat("Loading matrix for X-side (V) too big...ommited\n")
  }
}
