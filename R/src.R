#' @title Simulating Unbiased and Frequency-Biased Cultural Transmission 
#' @description Function for simulating unbiased/frequency-biased cultural transmission.
#' @param N Number of individuals.
#' @param timesteps Number of generations (timesteps).
#' @param mu Innovation rate.
#' @param wamUp Number of generations executed to reach equilibrium conditions (discarded from analysis)
#' @param top Number indicating the most common variants considered for the turn-over rate analysis (e.g. if \code{top} is equal to 5 the turnover rate is computed as the average number of new variants introduced in the top 1, 2, 3, ,4 and 5 ranked variants at each timestep).
#' @param bias Frequency bias parameter. Either a vector with length equal to the number of timesteps (to allow for a time-varying frequency bias), or single value (for a fixed bias). 
#' @param raw Logical variable indicating whether the matrix of variants frequencies are returned or not. Default is FALSE.
#' @details The function simulates unbiased (\code{bias}=0), confomist (\code{bias}<0), and anti-conformist (\code{bias}>0) cultural transmission with a user-defined population size and mutation rates, returning Simpson's diversity index of variants frequency, observed and expected turnover rates, and estimates of the exponent \eqn{x} following the procedure described in Acerbi and Bentley 2014 (notice that Acerbi and Bentley refer to this exponent as \eqn{b}). Frequency bias is modelled using the following formula:

#' \deqn{\pi_i=\frac{\left(\frac{m_i}{N}\right)^{1-\beta}}{\sum\limits_{j=1}^k \left(\frac{m_j}{N}\right)^{1-\beta}}(1-\mu)}

#' with \eqn{pi_i} equal to the probability of copying the variant \eqn{i}, \eqn{m_i} to the number of individuals possessing the variant \eqn{i}, \eqn{N} to the population size, \eqn{mu} to the mutation rate, and \eqn{\beta} to the frequency bias parameter (i.e. \code{bias}). 

#' @return A list containing: 1) a matrix of variants frequencies, with rows equal to  generations and columns to specific variants (\code{rawMatrix}, available on when the argument \code{raw} is set to \code{TRUE}.; 2) a data.frame with the observed turnover rates (column: \code{obs.z}), the expected rates according to  Bentley et al 2007 (column: \code{exp.z.Bentley}) and Evans and Giometto 2012 (column: \code{exp.z.EvansGiometto}); 3) the fitted value for the parameter \code{x}; 4) a vector with a time-series of observed Simpson's diversity index (\code{obs.div}); and 5) the expected Simpson's diversity under neutrality, equal to \eqn{ 1 - (1/(2*N*\mu + 1))} (\code{exp.div}).

#' @references 
#' Acerbi, A., Bentley, A.R., 2014. Biases in cultural transmission shape the turnover of populat traits. \emph{Evolution and Human Behavior}, 35, 228–236.\cr
#' Bentley, R.A., Lipo, C.P., Herzog, H.A., Hahn, M.W., 2007. Regular rates of popular culture change reflect random copying. \emph{Evolution and Human Behavior}, 28, 151–158. \cr
#' Crema, E.R., Kandler, A., Shennan, S., 2016. Revealing patterns of cultural transmission from frequency data: equilibrium and non-equilibrium assumptions. \emph{Scientific Reports} 6, 39122. https://doi.org/10.1038/srep39122 \cr
#' Evans, T.S., Giometto, A., 2011. Turnover Rate of Popularity Charts in Neutral Models. arXiv: [physics:soc-ph].
#' 
#' @import stats
#' @import utils
#' @export

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
	obs.div <- 1 - apply(rawMatrix,1,function(x){sum(c(x/sum(x))^2)})
	exp.div <- 1 - (1/(2*N*mu + 1))

	#TurnoverRate
	if (is.na(top)) {top = min(apply(rawMatrix,1,function(x){sum(x>0)}))} 
	zMatrix=turnover(rawMatrix,top=top)
	z = apply(zMatrix,2,mean) #average turn-over per top list of size y
	z.frame=data.frame(y=1:top,obs.z=z)
	
	aN = 1.38 * (mu^0.55) * N^0.13 

	z.frame$exp.z.Bentley=z.frame$y*sqrt(mu)
	z.frame$exp.z.EvansGiometto= aN*(z.frame$y)^0.86	/ 2 #see page 2 on Giometto and Evans on defintion of z
	z.frame$zfitN= z.frame$obs.z/aN

	#Empirical estimate of b
        modN=lm(log(zfitN+0.000001)~log(y),data=z.frame)
        b=as.numeric(coefficients(modN)[2])

	if (!raw){rawMatrix=NULL}
	
	return(list(rawMatrix=rawMatrix,
		    z.frame = z.frame[,1:4],
		    x = b,
		    obs.div = obs.div,
		    exp.div = exp.div))
}

#' @import utils
#' @keywords internal
#' @export

#Utility Function for counting cases
instances <- function(x,cases)
{
	x <- c(x,cases)
	return(table(x) - 1)
}


#' @import utils
#' @keywords internal
#' @export

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




#' @import utils
#' @keywords internal
#' @export

reScale = function(x,lo,hi)
{
return(((hi-lo)/(max(x)-min(x)))*(x-max(x))+hi)
}




#' @title Simulating Frequency-Biased Cultural Transmission with an Heterogenous Population.
#' @description Simulates unbiased/frequency-biased transmission with the bias parameter of each individual randomly drawn from a normal distribution with user-defined mean and standard deviation. 
#' @param N Number of individuals.
#' @param timesteps Number of generations (timesteps).
#' @param mu Innovation rate.
#' @param bmean Mean of the frequency bias parameter b.
#' @param bsd Standard deviation of the frequency bias parameter b.
#' @details The function simulates cultural transmission within a population of heterogenous learners, each with a different levels of frequency bias at each generation. More formally the probability of an individual $x$ copying a variant $i$ is given by the following equation.

#' \deqn{\pi_i=\frac{\left(\frac{m_i}{N}\right)^{1-\beta_x}}{\sum\limits_{j=1}^k \left(\frac{m_j}{N}\right)^{1-\beta_x}}(1-\mu)}

#' where \eqn{\beta_x} is randomly drawn from a normal distribution with mean \code{bmean} and standard deviation \code{bsd}. 
#'
#' @return A vector containing the Simpson's diversity index at each timestep/generation.
#' 
#' @import stats
#' @import utils
#' @export

heteroPopTransmission<-function(N=500,timesteps=301,mu=0.01,bmean=0,bsd=0)
{

	##Initialise Agents (each with a different mental template)
	diversity.ts = numeric(length=timesteps)
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

		diversity.ts[t]=1 - sum((table(agents)/N)^2)

	}
	return(diversity.ts)
}
