## Author Xiao Xiang
## Created: 10-12-2018
## Latest Update: 10-12-2018
## Updates------- 10-12---
# betaPrior function

### Binomial 
##Beta
# Set beta prior

betaPrior <- function(method = "mean", location = NULL, scale = NULL){
  # methods: mean, mean-sd, mode
  if (method == "mode" & scale <= 2){
    stop("Incorrect parameter for mode, please ensure scale > 2")
  }
  params <- switch(method,
                   "mean" = c(location*scale, (1-location)*scale),
                   "mean-sd" = c(location*(location*(1-location)/scale**2-1), (1-location)*(location*(1-location)/scale**2-1)),
                   "mode" = c(location*(scale-2)+1, (1-location)*(scale-2)+1))
  if (any(params < 0)){
    stop("Incorrect parameters, please check!")
  }
  return(params)
}

### diagnositic
## diagMCMC for coda sample from the book doing bayesian data analysis
openGraph = function( width=7 , height=7 , mag=1.0 , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    tryInfo = try( X11( width=width*mag , height=height*mag , type="cairo" , 
                        ... ) )
    if ( class(tryInfo)=="try-error" ) {
      lineInput = readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
      graphics.off() 
      X11( width=width*mag , height=height*mag , type="cairo" , ... )
    }
  } else { # Windows OS
    tryInfo = try( windows( width=width*mag , height=height*mag , ... ) )
    if ( class(tryInfo)=="try-error" ) {
      lineInput = readline("WARNING: Previous graphics windows will be closed because of too many open windows.\nTO CONTINUE, PRESS <ENTER> IN R CONSOLE.\n")
      graphics.off() 
      windows( width=width*mag , height=height*mag , ... )    
    }
  }
}

saveGraph = function( file="saveGraphOutput" , type="pdf" , ... ) {
  if ( .Platform$OS.type != "windows" ) { # Mac OS, Linux
    if ( any( type == c("png","jpeg","jpg","tiff","bmp")) ) {
      sptype = type
      if ( type == "jpg" ) { sptype = "jpeg" }
      savePlot( file=paste0(file,".",type) , type=sptype , ... )     
    }
    if ( type == "pdf" ) {
      dev.copy2pdf(file=paste0(file,".",type) , ... )
    }
    if ( type == "eps" ) {
      dev.copy2eps(file=paste0(file,".",type) , ... )
    }
  } else { # Windows OS
    file=paste0(file,".",type) 
    savePlot( file=file , type=type , ... )
  }
}


# Functions for computing limits of HDI's:

HDIofMCMC = function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

HDIofICDF = function( ICDFname , credMass=0.95 , tol=1e-8 , ... ) {
  # Arguments:
  #   ICDFname is R's name for the inverse cumulative density function
  #     of the distribution.
  #   credMass is the desired mass of the HDI region.
  #   tol is passed to R's optimize function.
  # Return value:
  #   Highest density iterval (HDI) limits in a vector.
  # Example of use: For determining HDI of a beta(30,12) distribution, type
  #   HDIofICDF( qbeta , shape1 = 30 , shape2 = 12 )
  #   Notice that the parameters of the ICDFname must be explicitly named;
  #   e.g., HDIofICDF( qbeta , 30 , 12 ) does not work.
  # Adapted and corrected from Greg Snow's TeachingDemos package.
  incredMass =  1.0 - credMass
  intervalWidth = function( lowTailPr , ICDFname , credMass , ... ) {
    ICDFname( credMass + lowTailPr , ... ) - ICDFname( lowTailPr , ... )
  }
  optInfo = optimize( intervalWidth , c( 0 , incredMass ) , ICDFname=ICDFname ,
                      credMass=credMass , tol=tol , ... )
  HDIlowTailPr = optInfo$minimum
  return( c( ICDFname( HDIlowTailPr , ... ) ,
             ICDFname( credMass + HDIlowTailPr , ... ) ) )
}

HDIofGrid = function( probMassVec , credMass=0.95 ) {
  # Arguments:
  #   probMassVec is a vector of probability masses at each grid point.
  #   credMass is the desired mass of the HDI region.
  # Return value:
  #   A list with components:
  #   indices is a vector of indices that are in the HDI
  #   mass is the total mass of the included indices
  #   height is the smallest component probability mass in the HDI
  # Example of use: For determining HDI of a beta(30,12) distribution
  #   approximated on a grid:
  #   > probDensityVec = dbeta( seq(0,1,length=201) , 30 , 12 )
  #   > probMassVec = probDensityVec / sum( probDensityVec )
  #   > HDIinfo = HDIofGrid( probMassVec )
  #   > show( HDIinfo )
  sortedProbMass = sort( probMassVec , decreasing=TRUE )
  HDIheightIdx = min( which( cumsum( sortedProbMass ) >= credMass ) )
  HDIheight = sortedProbMass[ HDIheightIdx ]
  HDImass = sum( probMassVec[ probMassVec >= HDIheight ] )
  return( list( indices = which( probMassVec >= HDIheight ) ,
                mass = HDImass , height = HDIheight ) )
}

DbdaAcfPlot = function( codaObject , parName=varnames(codaObject)[1] , plColors=NULL ) {
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain = length(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat = NULL
  yMat = NULL
  for ( cIdx in 1:nChain ) {
    acfInfo = acf(codaObject[,c(parName)][[cIdx]],plot=FALSE) 
    xMat = cbind(xMat,acfInfo$lag)
    yMat = cbind(yMat,acfInfo$acf)
  }
  matplot( xMat , yMat , type="o" , pch=20 , col=plColors , ylim=c(0,1) ,
           main="" , xlab="Lag" , ylab="Autocorrelation" )
  abline(h=0,lty="dashed")
  EffChnLngth = effectiveSize(codaObject[,c(parName)])
  text( x=max(xMat) , y=max(yMat) , adj=c(1.0,1.0) , cex=1.25 ,
        labels=paste("ESS =",round(EffChnLngth,1)) )
}

DbdaDensPlot = function( codaObject , parName=varnames(codaObject)[1] , plColors=NULL ) {
  if ( all( parName != varnames(codaObject) ) ) { 
    stop("parName must be a column name of coda object")
  }
  nChain = length(codaObject) # or nchain(codaObject)
  if ( is.null(plColors) ) plColors=1:nChain
  xMat = NULL
  yMat = NULL
  hdiLims = NULL
  for ( cIdx in 1:nChain ) {
    densInfo = density(codaObject[,c(parName)][[cIdx]]) 
    xMat = cbind(xMat,densInfo$x)
    yMat = cbind(yMat,densInfo$y)
    hdiLims = cbind(hdiLims,HDIofMCMC(codaObject[,c(parName)][[cIdx]]))
  }
  matplot( xMat , yMat , type="l" , col=plColors , 
           main="" , xlab="Param. Value" , ylab="Density" )
  abline(h=0)
  points( hdiLims[1,] , rep(0,nChain) , col=plColors , pch="|" )
  points( hdiLims[2,] , rep(0,nChain) , col=plColors , pch="|" )
  text( mean(hdiLims) , 0 , "95% HDI" , adj=c(0.5,-0.2) )
  EffChnLngth = effectiveSize(codaObject[,c(parName)])
  MCSE = sd(as.matrix(codaObject[,c(parName)]))/sqrt(EffChnLngth) 
  text( max(xMat) , max(yMat) , adj=c(1.0,1.0) , cex=1.25 ,
        paste("MCSE =\n",signif(MCSE,3)) )
}

diagMCMC = function( codaObject , parName=varnames(codaObject)[1] ,
                     saveName=NULL , saveType="jpg" ) {
  DBDAplColors = c("skyblue","black","royalblue","steelblue")
  openGraph(height=5,width=7)
  par( mar=0.5+c(3,4,1,0) , oma=0.1+c(0,0,2,0) , mgp=c(2.25,0.7,0) , 
       cex.lab=1.5 )
  layout(matrix(1:4,nrow=2))
  # traceplot and gelman.plot are from CODA package:
  require(coda)
  coda::traceplot( codaObject[,c(parName)] , main="" , ylab="Param. Value" ,
                   col=DBDAplColors ) 
  tryVal = try(
    coda::gelman.plot( codaObject[,c(parName)] , main="" , auto.layout=FALSE , 
                       col=DBDAplColors )
  )  
  # if it runs, gelman.plot returns a list with finite shrink values:
  if ( class(tryVal)=="try-error" ) {
    plot.new() 
    print(paste0("Warning: coda::gelman.plot fails for ",parName))
  } else { 
    if ( class(tryVal)=="list" & !is.finite(tryVal$shrink[1]) ) {
      plot.new() 
      print(paste0("Warning: coda::gelman.plot fails for ",parName))
    }
  }
  DbdaAcfPlot(codaObject,parName,plColors=DBDAplColors)
  DbdaDensPlot(codaObject,parName,plColors=DBDAplColors)
  mtext( text=parName , outer=TRUE , adj=c(0.5,0.5) , cex=2.0 )
  if ( !is.null(saveName) ) {
    saveGraph( file=paste0(saveName,"Diag",parName), type=saveType)
  }
}