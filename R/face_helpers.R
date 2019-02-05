concat_into_vector=function(frame_matrix){
  #will concat all rows of matrix into vector
  return(as.vector(t(frame_matrix)))
}

plot_frame=function(frame_vec,lim=T, single_frame=NULL,...){

  connections<- NULL;rm(connections); utils::data(connections); #to get rid of R CMD check complaints about using connections


  if (is.numeric(single_frame)){
    a = (single_frame-1)*136+1
    b = single_frame*136
    frame_vec=frame_vec[a:b]
  }

  frame_vec=concat_into_vector(frame_vec)
  x <- frame_vec[seq(1,length(frame_vec),2)]
  y <- frame_vec[seq(2,length(frame_vec),2)]*-1

  if(lim){
    plot(
      x,
      y,
      ylim=c(-0.2,0.13),xlim=c(-0.2,0.2),
      xlab="X", ylab="Y", ...
    )
  }
  else{
    plot(
      x,
      y,
      xlab="X", ylab="Y", ...
    )
  }

  segments(
    x[connections[,1]],
    y[connections[,1]],
    x[connections[,2]],
    y[connections[,2]]
  )
}

plot_latent_variables_face=function(plsr_obj,frame,...){
  args1=list(single_frame=frame, lim=T)
  args2=list(cex.axis = 1.5,font=2)
  plot_latent_variables(plsr_obj,FUN=c(plot_frame,barplot),args1=args1,args2=args2,...)
}
