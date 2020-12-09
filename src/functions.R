# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

get.elbow.points.indices <- function(x, y, threshold) {
  d1 <- diff(y) / diff(x) # first derivative
  d2 <- diff(d1) / diff(x[-1]) # second derivative
  indices <- which(abs(d2) > threshold)  
  return(indices)
}



find_modes<- function(x) {
  modes <- NULL
  for ( i in 2:(length(x)-1) ){
    if ( (x[i] > x[i-1]) & (x[i] > x[i+1]) ) {
      modes <- c(modes,i)
    }
  }
  if ( length(modes) == 0 ) {
    modes = 'This is a monotonic distribution'
  }
  return(modes)
}


# ROC optimal performance
opt.cut = function(perf, pred){
  cut.ind = mapply(FUN=function(x, y, p){
    d = (x - 0)^2 + (y-1)^2
    ind = which(d == min(d))
    c(sensitivity = y[[ind]], specificity = 1-x[[ind]], 
      cutoff = p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}

# root mean squared error
rmse <- function(sm) 
  sqrt(mean(sm$residuals^2))

turn_plot <- function(call_param, audio){
  #layout(rbind(c(1,2), c(1,2)))
  ggplot(audio,aes(x=abs(turn),y=call_param,colour=bat))+geom_point(alpha = 0.1)+
    geom_smooth(method="lm",alpha=0.3)+scale_y_continuous(limits = c(0,quantile(audio$call_param, .95)))
}

speed_plot <- function(call_param, audio){
  ggplot(audio,aes(x=speed,y=call_param,colour=bat))+geom_point(alpha = 0.1)+
    geom_smooth(method="lm",alpha=0.3)+scale_y_continuous(limits = c(0,quantile(audio$call_param, .95)))
}


range01 <- function(x){(x-min(x))/(max(x)-min(x))}

## better diff time function
dtime <- function(t, ...)
  as.numeric(difftime(t[-1], t[-length(t)], ...))


## basic euclidean distance
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

## degrees to radians
deg2rad <- function(deg) {(deg * pi) / (180)}

## test is step lengths and turning angles are from a normal distribution
jarque_bera <- function(m){
  #pseudo residuals
  with(momentuHMM::pseudoRes(m),{
    list(
      step=tseries::jarque.bera.test(stepRes[complete.cases(stepRes)]),
      angleRes=tseries::jarque.bera.test(angleRes[complete.cases(angleRes)])
    )})
}


## Plot confidence intervals around coefficients for covariates in HMM analysis
plot_CI <- function(m) {
  
  coef_ci <- function(df_ci,header) {
    plotrix::plotCI(1:nrow(df_ci),df_ci[,1], ui=df_ci[,3], li=df_ci[,2], xaxt='n',
                    xlab="Covariate",ylab="95% Confidence Interval",lwd=2)
    
    axis(1,1:nrow(df_ci),rownames(df_ci))
    title(header)
    grid()
    abline(h=0,col=2,lty=2)
  }
  pp <- mapply(function(beta,lower,upper,rows)
  {
    df <- data.frame(beta,lower,upper)
    rownames(df) <- rows
    df
  },
  as.data.frame(m$mle$beta),
  as.data.frame(moveHMM::CI(m)$beta$lower),
  as.data.frame(moveHMM::CI(m)$beta$upper),
  MoreArgs = list(rows=rownames(m$mle$beta)),
  SIMPLIFY = F)
  names(pp) <- colnames(m$mle$beta)
  n <- ncol(m$mle$beta)
  op <- par(mfcol=c(2,n/2),mar=c(3,4,1,2))
  options(warn=-1)
  invisible(mapply(coef_ci,pp,names(pp)))
  options(warn=0)
  par(op)
}

# extract environmental data for move data.frame
extract_env <- function(moves, date, layer, env_rasts){
  ids <- c('erdMWchla3day', # Chlorophyll-a, Aqua MODIS, NPP, 0.0125?, West US, EXPERIMENTAL, 2002-present (3 Day Composite)
           'erdMBchla3day', # Chlorophyll-a, Aqua MODIS, NPP, 0.025 degrees, Pacific Ocean, 2006-present, EXPERIMENTAL (3 Day Composite)
           'erdMWsstd3day', # SST, Aqua MODIS, NPP, 0.0125?, West US, Day time (11 microns), 2002-present (3 Day Composite)
           'bathy', 'slope')
  
  idx <- which(ids == layer)
  dates <- unique(date(moves$t))
  didx <- which(dates == date)
  var_day <- paste(ids[idx], dates[didx], sep = "_")
  var_day <- gsub("-", ".", var_day)
  
  if(idx == 4){
    var_day <- 'bathy'
  }
  if(idx == 5){
    var_day <- 'slope'
  }
  rast <- env_rasts[[var_day]]
  midx <- try(which(date(moves$t) == dates[didx]))
  
  XY <- data.frame(ID = 1, X = moves$X[midx], Y = moves$Y[midx])
  coordinates(XY) <- c("X", "Y")
  proj4string(XY) <- CRS("+proj=utm +zone=12 ellps=WGS84")
  
  # if(idx == 1){ moves$chla[midx] <- raster::extract(rast, XY)}
  # if(idx == 2){ moves$NPP[midx] <- raster::extract(rast, XY)}
  # if(idx == 3){ moves$SST[midx] <- raster::extract(rast, XY)}
  # if(idx == 4){ moves$bathy[midx] <- raster::extract(rast, XY)}
  return(raster::extract(rast, XY))
}

# extract env data for sp object
extract_env_ssf <- function(moves, date, layer, env_rasts){
  
  ids <- c('erdMWchla3day', # Chlorophyll-a, Aqua MODIS, NPP, 0.0125?, West US, EXPERIMENTAL, 2002-present (3 Day Composite)
           'erdMBchla3day', # Chlorophyll-a, Aqua MODIS, NPP, 0.025 degrees, Pacific Ocean, 2006-present, EXPERIMENTAL (3 Day Composite)
           'erdMWsstd3day', # SST, Aqua MODIS, NPP, 0.0125?, West US, Day time (11 microns), 2002-present (3 Day Composite)
           'bathy', 'slope')
  
  idx <- which(ids == layer)
  dates <- unique(date(moves$time))
  didx <- which(dates == date)
  var_day <- paste(ids[idx], dates[didx], sep = "_")
  var_day <- gsub("-", ".", var_day)
  
  if(idx == 4){
    var_day <- 'bathy'
  }
  if(idx == 5){
    var_day <- 'slope'
  }
  rast <- env_rasts[[var_day]]
  midx <- which(date(moves$time) == date(date))
  XY <- data.frame(ID = 1, X = coordinates(moves)[midx,1], Y = coordinates(moves)[midx,2])
  coordinates(XY) <- c("X", "Y")
  proj4string(XY) <- CRS("+proj=utm +zone=12 ellps=WGS84")
  
  # if(idx == 1){ moves$chla[midx] <- raster::extract(rast, XY)}
  # if(idx == 2){ moves$NPP[midx] <- raster::extract(rast, XY)}
  # if(idx == 3){ moves$SST[midx] <- raster::extract(rast, XY)}
  # if(idx == 4){ moves$bathy[midx] <- raster::extract(rast, XY)}
  return(raster::extract(rast, XY))
}

# prep data for SSF
ssf_ready <- function(move_sf){
  #time ordering for caterpillar
  #habitat used
  used <- rep(1, nrow(move_sf))
  stratum <- 1:nrow(move_sf)
  used_sf <- sf::st_sf(data.frame(move_sf,used, stratum))
  # used_sf <- subset(used_sf)
  
  #null
  step_sample <- 10
  used <- rep(0, (nrow(used_sf)-2) * step_sample)
  #caterpillar only works with projected
  cater_sf <- caterpillar(used_sf,step_sample,used_sf$time, do_plot = T)
  null_sf <- sf::st_sf(data.frame(used , cater_sf))
  
  all_sf <- rbind(null_sf,used_sf)
  return(all_sf)
}

#' caterpillar
#'
#' @param Z a complex vector of locations ordered by time
#' @param do_plot whether to plot the result
#'
#' @return a matrix with a dimension of (n-2)*(n-2) of pseudo absences
#' (n-2)*n.pseudos
#' @export
#'
#' @examples
#' #simulate
#' X <- cumsum(arima.sim(n=3000, model=list(ar=.7)))
#' Y <- cumsum(arima.sim(n=3000, model=list(ar=.7)))
#' Z <- X + 1i*Y
#'
#' op <- par(mfrow=c(1,2))
#' znull<-caterpillar_elie(Z[1:10],T)
#' znull<-caterpillar(Z[1:10],T)
#' par(op)
#'
#' print(system.time(caterpillar_elie(Z)))
#' print(system.time(caterpillar(Z)))

caterpillar <- function(used_sf,n.pseudos=10, time = NA_real_, do_plot=F) {
  #if move data aren't regularly sampled must divide by difftime
  
  if(!inherits(used_sf,"sf"))
    stop("used_sf must be of type sf")
  
  geom <- sf::st_geometry(used_sf)
  xy <- do.call(rbind,geom)
  
  Z <- complex(re = xy[,1],  imag=xy[,2])
  n <- length(Z)
  S <- Mod(diff(Z))
  Phi <- Arg(diff(Z))
  Theta <- diff(Phi)
  RelSteps <- complex(mod = S[-1], arg=Theta)
  
  Z_inner <- Z[-c(1,n)] #2:9
  Azimuths <- complex(mod = 1, arg=Phi[-(n-1)] ) # [-9]
  
  #Add the rotated steps to the last step
  z_null <-apply(outer(Azimuths,RelSteps),2,'+',Z_inner)
  #the first and n-1 sample of pseudo-absences
  z_null <- subset(z_null,select=sample(2:ncol(z_null), n.pseudos))
  
  ####
  #ELIE CODE
  #Add the rotated steps to the last step
  #z_null <- matrix(0, ncol=n-2, nrow=n-2)
  #for(i in 1:length(Z1))
  # z_null[i,] <- Z1[i] + RelSteps * Azimuths[i]
  ####
  
  #Make the fuzzy catterpillar plot
  if(do_plot)
  {
    plot(Z, type="o", col=1:10, pch=19, asp=1)
    invisible(
      lapply(1:nrow(z_null), function(i)
        segments(rep(Re(Z_inner[i]), n-2),
                 rep(Im(Z_inner[i]), n-2),
                 Re(z_null[i,]), Im(z_null[i,]), col=i+1)
      )
    )
  }
  # +1 because, stratum starts from 1 in used_df
  stratum <- rep(1:nrow(z_null), ncol(z_null)) + 1
  if (all(!is.na(time)))
    time <- rep(time[-c(1,length(time))],ncol(z_null))
  
  #as.vector is by col. use t to make it by row
  z_null <- as.vector(z_null)
  cater <- data.frame(stratum=stratum, x=Re(z_null), y=Im(z_null),time=time )
  
  return(sf::st_as_sf(cater,coords=c('x','y'), crs=st_crs(used_sf) ))
}


scantrack <- function(track, ...){
  layout(rbind(c(1,2), c(1,3)))
  with(track, plot(x,y,type = "o", main = Individual[1], col = date(time)))
  lines(Islands)
  lines(Coast)
  with(track, plot(time,x))
  with(track, plot(time,y))
  layout(1)
}