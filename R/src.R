#' @title Simulating Unbiased and Frequency-Biased Cultural Transmission 
#' @description
#' @param N Number of individuals.
#' @param timesteps Number of generations (timesteps).
#' @param mu Innovation rate.
#' @param wamUp Number of generations executed to reach equilibrium conditions (discarded from analysis)
#' @param top Number indicating the most common variants considered for the turn-over rate analysis (e.g. if top is equal to 5 the turnover rate is computed as the average number of new variants introduced in the top 1, 2, 3, ,4 and 5 ranked variants).
#' @param bias Frequency bias parameter. If negative the transmission is conformist, if positive the transmission is anti-conformist, and if equal to 0 the transmission is unbiased (see details below).
#' @param raw Logical variable indicating whether the matrix of variants frequencies are returned or not. Default is FALSE.
#' @details
#' @return A list containing: 1) a matrix of variants frequencies, with rows equal to  generations and columns to specific variants (\code{rawMatrix}, available on when the argument \code{raw} is set to \code{TRUE}.; 2) a data.frame with the observed turnover rates (column: \code{obs.z}), the expected rates according to  Bentley et al 2007 (column: \code{exp.z.Bentley}) and Evans and Giometto 2012 (column: \code{exp.z.EvansGiometto}); 3) the fitted value for the parameter $b$.

#' @references
#' @examples
#' @import stats
#' @import utils


transmission<-function(N=500,timesteps=501,mu=0.01,warmUp=300,top=NA,bias=0,raw=FALSE)
{

	##Initialise Agents (each with a different mental template)
	agents <- 1:N
	traitCounter <- N + 1 #counter for inifinite allele innovation 
	rawList <- vector("list",length= timesteps-warmUp)
	if (length(bias)==1){bias=rep(bias,timesteps)} #time-varyng frequency bias; default is none

	for (t in 1:timesteps)
	{ 	

		samplePool = agents
		if (length(unique(samplePool))>1)
		{
			sampleTraits = as.numeric(names(table(samplePool)))
			sampleProp = table(samplePool)/sum(table(samplePool))
			sampleTraitsProb = sampleProp^(1-bias[t]) / sum(sampleProp^(1-bias[t]))
			# Cultural Transmission:
			agents = sample(sampleTraits,size=N,replace=TRUE,prob=sampleTraitsProb)
		}
		# Innovation
		index = which(runif(N)<=mu)
		if (length(index)>0)
		{

			newTraits<-traitCounter:c(traitCounter+length(index)-1)
			agents[index]=newTraits
			traitCounter=max(newTraits)+1
		}


		#Store 
		if (t>warmUp) {rawList[[t-warmUp]] = agents}	    
	}

	#Transform 
	cases <- unique(unlist(rawList))
	rawMatrix <- t(sapply(rawList,instances,cases=cases))

	#Compute Diversity Indices
	tF <- 1 - apply(rawMatrix,1,function(x){sum(c(x/sum(x))^2)})
	tF.exp <- 1 - (1/(2*N*mu + 1))

	#TurnoverRate
	if (is.na(top)) {top = min(apply(rawMatrix,1,function(x){sum(x>0)}))} 
	zMatrix=turnover(rawMatrix,top=top)
	z = apply(zMatrix,2,mean) #average turn-over per top list of size y
	z.frame=data.frame(y=1:top,obs.z=z)
	
	aN = 1.38 * (mu^0.55) * N^0.13 

	z.frame$exp.z.Bentley=z.frame$y*sqrt(mu)
	z.frame$exp.z.N= aN*(z.frame$y)^0.86	/ 2 #see page 2 on Giometto and Evans on defintion of z
	z.frame$zfitN= z.frame$obs.z/aN

	#Empirical estimate of b
        modN=lm(log(zfitN+0.000001)~log(y),data=z.frame)
        bN=as.numeric(coefficients(modN)[2])

	if (!raw){rawMatrix=NULL}

	return(list(rawMatrix=rawMatrix,
		    z.frame=z.frame,
		    tF=tF,
		    tF.exp = tF.exp,
		    bN=bN))
}



#Utility Function for counting cases
instances <- function(x,cases)
{
	x <- c(x,cases)
	return(table(x) - 1)
}

#Compute turnover rate from frequency matrix
turnover <- function(mat,top)
{
	z<-matrix(NA,nrow=c(nrow(mat)-1),ncol=top)
	for (i in 1:c(nrow(mat)-1))
	{
		t1 = names(sort(mat[i,],decreasing=TRUE))
		t2 = names(sort(mat[i+1,],decreasing=TRUE))
		z[i,]=sapply(1:top,function(x,t1,t2){return(x-length(intersect(t1[1:x],t2[1:x])))},t1=t1,t2=t2)
	}
	return(z)
}


reScale = function(x,lo,hi)
{
return(((hi-lo)/(max(x)-min(x)))*(x-max(x))+hi)
}





heteroPopTransmission<-function(N=500,timesteps=301,mu=0.01,bsd=0)
{

	##Initialise Agents (each with a different mental template)
	tFseries = numeric(length=timesteps)
	agents <- 1:N
	traitCounter <- N + 1 #counter for inifinite allele innovation 

	biasedSampling = function(b,s,p)
	{
		p = p^(1-b) / sum(p^(1-b))
		return(sample(s,size=1,prob=p))
	}

	for (t in 1:timesteps)
	{ 	
		samplePool = agents
		if (length(unique(samplePool))>1)
		{
			sampleTraits = as.numeric(names(table(samplePool)))
			sampleProp = table(samplePool)/sum(table(samplePool))
			bias = rnorm(N,mean=0,sd=bsd)
			agents = sapply(bias,biasedSampling,s=sampleTraits,p=sampleProp) 
		}
		# Innovation
		index = which(runif(N)<=mu)
		if (length(index)>0)
		{

			newTraits<-traitCounter:c(traitCounter+length(index)-1)
			agents[index]=newTraits
			traitCounter=max(newTraits)+1
		}

		tFseries[t]=1 - sum((table(agents)/N)^2)

	}

	#Compute Diversity Indices
	return(tFseries)
}
