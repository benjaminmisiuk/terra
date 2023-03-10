k.aspect <- function(x, k=3, unit="degrees", type="aspect", method="queen"){ 
  
  require(raster)    
  
  #j is the number of cells on either side of the focal cell; l is used to generate the focal matrix
  j <- (k/2)-0.5
  l <- j-1
  
  if(method=="queen"){
    
    #create matrix weights for x-component
    xl.end <- matrix(c(1, rep(0, times=k-1)), ncol=k, nrow=1)
    xr.end <- matrix(c(rep(0, times=k-1), 1), ncol=k, nrow=1)
    
    x.mids <- matrix(0, ncol=k, nrow=l)
    
    xl.mid <- matrix(c(2, rep(0, times=k-1)), ncol=k, nrow=1)
    xr.mid <- matrix(c(rep(0, times=k-1), 2), ncol=k, nrow=1)
    
    xl.mat <- rbind(xl.end, x.mids, xl.mid, x.mids, xl.end)
    xr.mat <- rbind(xr.end, x.mids, xr.mid, x.mids, xr.end)
    
    #create matrix weights for y-component
    yt.end <- matrix(c(1, rep(0, times=k-1)), ncol=1, nrow=k)
    yb.end <- matrix(c(rep(0, times=k-1), 1), ncol=1, nrow=k)
    
    y.mids <- matrix(0, ncol=l, nrow=k)
    
    yt.mid <- matrix(c(2, rep(0, times=k-1)), ncol=1, nrow=k)
    yb.mid <- matrix(c(rep(0, times=k-1), 2), ncol=1, nrow=k)
    
    yt.mat <- cbind(yt.end, y.mids, yt.mid, y.mids, yt.end)
    yb.mat <- cbind(yb.end, y.mids, yb.mid, y.mids, yb.end)
    
    #use focal statistics for e, w, n, s components of the k-neighbourhood
    dz.dx.l <- focal(x, xl.mat, fun=sum)
    dz.dx.r <- focal(x, xr.mat, fun=sum)
    
    dz.dy.t <- focal(x, yt.mat, fun=sum)
    dz.dy.b <- focal(x, yb.mat, fun=sum)
    
    #calculate dz/dx and dz/dy using the components. 8*j is the weighted run, or distance between ends: 4*j*2, or (4 values in each row)*(length of the side)*(2 sides)
    dz.dx <- (dz.dx.r-dz.dx.l)/(8*j*res(x)[1])
    dz.dy <- (dz.dy.b-dz.dy.t)/(8*j*res(x)[2])
  }
  
  if(method=="rook"){
    
    #create matrix weights for x-component
    x.ends <- matrix(0, ncol=k, nrow=j)
    
    xl.mid <- matrix(c(1, rep(0, times=k-1)), ncol=k, nrow=1)
    xr.mid <- matrix(c(rep(0, times=k-1), 1), ncol=k, nrow=1)
    
    xl.mat <- rbind(x.ends, xl.mid, x.ends)
    xr.mat <- rbind(x.ends, xr.mid, x.ends)
    
    #create matrix weights for y-component
    y.ends <- matrix(0, ncol=j, nrow=k)
    
    yt.mid <- matrix(c(1, rep(0, times=k-1)), ncol=1, nrow=k)
    yb.mid <- matrix(c(rep(0, times=k-1), 1), ncol=1, nrow=k)
    
    yt.mat <- cbind(y.ends, yt.mid, y.ends)
    yb.mat <- cbind(y.ends, yb.mid, y.ends)
    
    #use focal statistics for e, w, n, s components of the k-neighbourhood
    dz.dx.l <- focal(x, xl.mat, fun=sum)
    dz.dx.r <- focal(x, xr.mat, fun=sum)
    
    dz.dy.t <- focal(x, yt.mat, fun=sum)
    dz.dy.b <- focal(x, yb.mat, fun=sum)
    
    #calculate dz/dx and dz/dy using the components. 2*j is the run: (2 sides)*(length of each side)
    dz.dx <- (dz.dx.r-dz.dx.l)/(2*j*res(x)[1])
    dz.dy <- (dz.dy.b-dz.dy.t)/(2*j*res(x)[2])
  }
  
  aspect.k <- (180/pi)*atan2(dz.dy, -1*dz.dx)
  
  f <- function(x) ifelse(x<0, 90-x, ifelse(x>90, 360-x+90, 90-x))
  aspect.k <- (calc(aspect.k, f))*(pi/180)
  
  #for 0-360
  if(type=="aspect"){
    
    if(unit=="radians"){
      return(aspect.k)
    }
    
    if(unit=="degrees"){
      return(aspect.k*(180/pi))
    }
  }
  
  #for eastness and northness
  if(type=="components"){
    return(list(sin(aspect.k), cos(aspect.k)))
  }
}
