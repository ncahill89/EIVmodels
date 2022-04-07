igp_data<-function(dat)
{
  ############# Set up the grid for the GP ###################
    xgrid = seq(min(dat$x_st), max(dat$x_st), length.out = 50)
    Ngrid = length(xgrid)

  ###Change data to lower zero for integration
  minx = min(dat$x_st)
  x = dat$x_st-minx
  xstar = xgrid - minx

  Dist <- fields::rdist(xstar) ###Distance matrix required for the model
  D <- cbind(x,dat$y) ###Combine the x,y data for the model

  ########Initialize quadrature for the integration########
  N <- nrow(dat)
  L = 30    ## this sets the precision of the integration quadrature (higher is better but more computationally expensive)
  index=1:L
  cosfunc=cos(((2*index-1)*pi)/(2*L))

  quad1=array(dim=c(nrow=N,ncol=Ngrid,L))
  quad2=array(dim=c(nrow=N,ncol=Ngrid,L))

  for(j in 1:Ngrid)
  {   for(k in 1:N)
  {
    quad1[k,j,]=abs((x[k]*cosfunc/2)+(x[k]/2)-xstar[j])^1.99
    quad2[k,j,]=((x[k]/2)*(pi/L))*(sqrt(1-cosfunc^2))
  }
  }


  return(list(xstar = xstar,
              N = N,
              Ngrid = Ngrid,
              Dist = Dist,
              quad1 = quad1,
              quad2 = quad2,
              cosfunc = cosfunc,
              ppi = pi,
              L = L))
}
