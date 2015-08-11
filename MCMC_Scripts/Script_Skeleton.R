###Read in data
{
  if (Sys.info()[8] == "dpritch")
    setwd("/home/dpritch/Documents/Projects/Dunson Day Specific/Software")
  else
    setwd("/Users/Sam/Documents/Sam/School/Graduate/Dave/DSP_Software/")
}
load("Data/PracticeDat.RData")
source("MCMC_Scripts/Script_SkeletonHelper.R")

###Create data objects
K <- 5 #fertile window
fwLen <- length(X) / length(Y) # this is equal to K for this dataset
ncyc <- length(Y) #number of cycles
N <- ncyc*K #number observation days
n <- length(unique(id)) #number of individuals
H <- dim(U)[2] #number of covariates


idIdx <- lapply(X=unique(id), FUN=function(x) which(id == x))
  # index to expand from baseline to cycles
dayExpanIdx <- rep(x=seq(from=1, to=n), times=table(id))
phiPropDist <- "normal"


###Matrix objects
#eyen <- diag(n)
#onen <- rep(1, n)

###Initial values
phi <- 1
g <- rep(1, H) #cant use gamma, since we use the gamma function
beta <- log(g)
xi <- rep(1, n)

###Hyperparameters
a <- rep(3, H)
b <- rep(1, H)
p <- rep(0.5, H)
A <- matrix(nrow=H, ncol=2)
A[,1] <- rep(0, H)
A[,2] <- rep(Inf, H)
hypPhi <- list(c1=1, c2=1)

###Metropolis objects
#phi_propvar <- 1
#acceptance_phi <- 0
delta <- 0.2
numAcceptPhi <- 0
  

###Distribution functions

##Sample from multinomial
rmult <- function(ns, P) {
  Ys <- matrix(0, nrow(P), ncol(P)) 
	if (sum(ns)>0)
		Ys[ns>0,1] <- rbinom(sum(ns>0), ns[ns>0], P[ns>0,1])
	ns <- ns-Ys[ ,1]
	if (sum(ns)>0) Ys[ns>0,2] <- rbinom(sum(ns>0), ns[ns>0], P[ns>0,2]/(1-P[ns>0,1]))
		ns<- ns-Ys[,2]
	for(j in 3:(ncol(P)-1)) {
		if(sum(ns)>0)
			Ys[ns>0,j] <- rbinom(sum(ns>0), ns[ns>0], P[ns>0,j]/(1-P[ns>0,1:(j-1)]%*%rep(1, j-1)))
		ns <- ns-Ys[ ,j]
	}
	Ys[ ,ncol(P)] <- ns
	Ys
}

##Sample from one-inflated gamma truncated by region A
rI1gammaT <- function(n, p, a, b, A){
  z <- rbinom(n, 1, p)
  z+(1-z)*qgamma(runif(n, pgamma(A[1],a,b), pgamma(A[2],a,b)), a, b)
}

##Log gamma constant
lC <- function(a, b) {
	b*log(a)-lgamma(a)
}

###MCMC Objects
nsims <- 10000
output <- "/outputfolder/"
write(1:H, file=paste(output,"GAMMA.csv",sep=""), sep=",", ncolumns=H)
write(1:n, file=paste(output,"XI.csv",sep=""), sep=",", ncolumns=n)
write(1, file=paste(output,"PHI.csv",sep=""), sep=",", ncolumns=1)

###Time MCMC Sampler
begin <- Sys.time()

###Begin Sampler
for (s in 1:nsims) {
	
	##Sample latent variable Z
	Z <- numeric(length=N)
	for (i in 1:ncyc) {
		
		#Subset data to individual i, cycle j
		idcyc <- (K*(i-1)+1):(i*K)
		u <- U[idcyc, ]
		x <- X[idcyc]
		y <- Y[i]
		ID <- unique(id[idcyc])

		#Sample Z_ij
		if (y!=0) {
			
			#Sample z_ij (scalar Z_ij) from zero truncated poisson (https://stat.ethz.ch/pipermail/r-help/2005-May/070678.html)
			mu <- xi[ID]*exp(u%*%beta)
			mean_pois <- as.numeric(x%*%mu)
			unif <- runif(1)
			tt <- -log(1-unif*(1-exp(mean_pois)))
			T1 <- mean_pois-tt
			z_ij <- rpois(1, T1)+1
			
			#Sample Z_ij (vector Z_ij)
			ns <- z_ij
			P <- matrix(x*mu/mean_pois, nrow=1, ncol=K)
			Z[idcyc] <- rmult(ns, P)
				
		#End sample Z_ij
		}
	
	##End latent Z loop
	}

	##Gamma_h full conditional
	for (h in 1:H) {
		
		#Set hyperparameters
		ah <- a[h]
		bh <- b[h]
		ph <- p[h]
		Ah <- A[h, ]
		
		#Calculate posterior hyperparameters
		
			#ah_tilde
			uh <- U[ ,h]		
			ah_t <- ah+sum(uh*Z)	
			
			#bh_tilde
			temp <- g^U
			prod_gamma <- apply(temp[ ,-h], 1, prod)
			xi_i <- unlist( mapply(rep, xi, as.numeric(table(id))) )
			bh_t <- bh+sum(xi_i[X==1]*prod_gamma[X==1]) #sum over only X==1
			
			#ph_tilde
			cons1 <- exp( lC(ah,bh) - lC(ah_t,bh_t) + (bh_t-bh) )
			cons2 <- ( pgamma(Ah[2], ah_t, bh_t) - pgamma(Ah[1], ah_t, bh_t) ) / ( pgamma(Ah[2], ah, bh) - pgamma(Ah[1], ah, bh) )
			den<- cons1 * cons2 * ( 1 - ph ) + ph
			ph_t <- ph / den
		
		#Sample from one-inflated gamma truncated to region Ah
		g[h] <- rI1gammaT(1, ph_t, ah_t, bh_t, Ah)		

	##End gamma sampler
	}
	
	##Update beta
	beta <- log(g)
	
	##Write gamma to file
	write(g,file=paste(output,"GAMMA.csv",sep=""),sep=",",ncolumns=H,append=TRUE)
	
	##Xi_i Full Conditional
  betaCoef <- beta
  uProdBeta <- U %*% betaCoef
  W <- Z
  
	xi <- sampXi(W, uProdBeta, phi, idIdx)
  xiDay <- xi[dayExpanIdx]
	write(xi,file=paste(output,"XI.csv",sep=""),sep=",",ncolumns=n,append=TRUE)

		
	##Metroplois Step For Phi
	phiProp <- sampPhiProp(phi, phiPropDist, delta)
  phiLogR <- getPhiLogR(xi=xi, phiCurr=phi, phiProp=phiProp, hypPhi=hypPhi)
  
	if (log(runif(1)) < phiLogR) {
	  phi <- phiProp
	  #acceptance_phi <- acceptance_phi+1 #keep track of acceptance rate
    numAcceptPhi <- numAcceptPhi + 1
	}
	write(phi, file=paste(output,"PHI.csv",sep=""), sep=",", ncolumns=1, append=TRUE)
	
	##Verbose (we can personalize the verbose that we output...)
	if (s!=nsims) print(paste("Completed Percentage: ",round((s/nsims)*100,digits=0),"%",sep=""))
	print(paste("Gamma: ",round(gamma[1],digits=3),", ",round(gamma[2],digits=3),sep=""))
	print(paste("Phi: ",round(phi,digits=3),sep=""))
	print(paste("Deviance: ",round(dev,digits=3),sep=""))
	print(paste("Acceptance Rate: Phi (",100*round(mean(acceptance_phi)/s,digits=2),"%)",sep=""))
	print("######################################################################################################")
	###Time MCMC Sampler
	if (s==nsims) {
		after<-Sys.time()
		time<-after-begin
		print(paste("Run Time: ",round(time,digits=2)," ",attr(time,"unit"),sep=""))
		print("######################################################################################################")

	}
	
###End MCMC Sampler	
}


