#used by pairsplus
panel.hist <- function(x,right=FALSE,diagCol=5,linefun=mean, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x,plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y,col=diagCol,   ...)
  abline(v=linefun(x),col=2)
}

#used by pairsplus
panel.cor <- function(x, y, digits=2, prefix="", corcex=0.5,method="pearson", ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- (cor(x, y,method=method,use="pairwise.complete.obs"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(corcex)) corcex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = corcex * sqrt(abs(r))+0.5)
}

#used by pairsplus
panel.scatter <- function(x,y,fitcurve='linear',fitcol=2,crossx=0,crossy=0,...) {
  points(x,y,...)

  if(fitcurve=='linear')
  {
    abline(lsfit(x,y),col=fitcol)

  } else if(fitcurve=='crosshairs') {

    abline(h=crossy,col=fitcol)
    abline(v=crossx,col=fitcol)

  } else if(fitcurve=='spline') {
    abline(smooth.spline(x,y),col=fitcol,...)
  }
}

#' pairs plot
#' @description fancy pairs plot of a dataframe. uses `pairs` internally.
#' @param x a numeric dataframe
#' @param diag.panel
#' @param diagCol
#' @param fitcurve
#' @param ...
#'
#' @return
#' @export
#'
#' @examples pairsPlus(cars)
pairsPlus<-function(x,diag.panel=panel.hist,diagCol=4,fitcurve='linear',...)
{pairs(x,diag.panel=diag.panel,upper.panel=panel.cor,lower.panel=panel.scatter,...)}
