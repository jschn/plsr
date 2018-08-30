concat_into_vector=function(frame_matrix){
  #will concat all rows of matrix into vector
  return(as.vector(t(frame_matrix)))
}

plot_frame=function(frame_vec,title="",lim=T, single_frame=NULL){

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
      xlab="X", ylab="Y", main=title
    )
  }
  else{
    plot(
      x,
      y,
      xlab="X", ylab="Y", main=title
    )
  }

  segments(
    x[connections[,1]],
    y[connections[,1]],
    x[connections[,2]],
    y[connections[,2]]
  )
}
