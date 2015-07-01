###Load Libraries (may need a few extra functions...)

###Read in Data (format data into the Y, U and X matrices)
{
  if (Sys.info()[8] == "dpritch")
    setwd("/home/dpritch/Documents/Projects/Dunson Day Specific/Software")
  else
    setwd("/Users/Sam/Documents/Sam/School/Graduate/Dave/DSP_Software/")
}
load("Data/PracticeDat.RData")

###Create Data Objects (define objects like n, n_i)
K<-5 #fertile window
ncyc<-length(Y) #number of cycles
N<-dim(U)[1] #number observation days
n<-length(unique(id)) #number of individuals
H<-dim(U)[2] #number of covariates

###Matrix Objects (define objects )
eyen<-diag(n)
onen<-rep(1,n)

###Initial Values
phi<-1
g<-rep(1,H)
b<-log(g)
xi<-rep(1,n)

###Hyperparameters
a<-rep(3,H)
b<-rep(1,H)
p<-rep(0.5,H)
A<-matrix(nrow=H,ncol=2)
A[,1]<-rep(0,H)
A[,2]<-rep(Inf,H)

###Metropolis Objects
phi_propvar<-1
acceptance_phi<-0

###Distribution Functions

## Sample from multinomial
rmult<- function(ns, P){
  Ys<- matrix(0,nrow(P),ncol(P))
  if(sum(ns)>0)
    Ys[ns>0,1]<- rbinom(sum(ns>0),ns[ns>0],P[ns>0,1])
  ns<- ns-Ys[,1]
  if(sum(ns)>0) Ys[ns>0,2]<- rbinom(sum(ns>0),ns[ns>0],P[ns>0,2]/(1-P[ns>0,1]))
  ns<- ns-Ys[,2]
  for(j in 3:(ncol(P)-1)){
    if(sum(ns)>0)
      Ys[ns>0,j]<- rbinom(sum(ns>0),ns[ns>0], 
                          P[ns>0,j]/(1-P[ns>0,1:(j-1)]%*%rep(1,j-1)))
    ns<- ns-Ys[,j]
  }
  Ys[,ncol(P)]<- ns
  Ys
}

#Sample from one-inflated gamma truncated by region A
rI1gammaT <- function(n,p,a,b,A){
  z <- rbinom(n,1,p)
  z+(1-z)*qgamma(runif(n,pgamma(A[1],a,b),pgamma(A[2],a,b)),a,b)
}

##Gamma Constant
C <- function(a,b) {out<-b^a/gamma(a); return(out)}

###MCMC Objects
nsims<-10000
output<-"/outputfolder/"
write(1:H,file=paste(output,"GAMMA.csv",sep=""),sep=",",ncolumns=H)
write(1:n,file=paste(output,"XI.csv",sep=""),sep=",",ncolumns=n)
write(1,file=paste(output,"PHI.csv",sep=""),sep=",",ncolumns=1)

###Time MCMC Sampler
begin<-Sys.time()

###Begin Sampler
for (s in 1:nsims) {
	
	##Sample Latent Variable Z
	Z <- numeric(length=N)
	for (i in 1:ncyc) {
		
		#Define data for individual i, cycle j
		idcyc <- (K*(i-1)+1):(i*K)
		u <- U[idcyc,]
		x <- X[idcyc]
		y <- Y[i]
		ID <- unique(id[idcyc])

		#Sample Z_ij
		if (y!=0) {
			
			#Sample z_ij (scalar Z_ij) from zero truncated poisson (https://stat.ethz.ch/pipermail/r-help/2005-May/070678.html)
			mu <- xi[ID]*exp(u%*%b)
			mean_pois <- as.numeric(x%*%mu)
			unif <- runif(1)
			tt <- -log(1-unif*(1-exp(mean_pois)))
			T1 <- mean_pois-tt
			z_ij <- rpois(1,T1)+1
			
			#Sample Z_ij (vector Z_ij)
			ns <- z_ij
			P <- matrix(x*mu/mean_pois,nrow=1,ncol=K)
			Z[idcyc]<-rmult(ns,P)
				
		#End Sample Z_ij
		}
	
	##End latent Z loop
	}

	##Gamma_h full conditional
	for (h in 1:H) {
		
		#Set hyperparameters
		ah <- a[h]
		bh <- b[h]
		ph <- p[h]
		Ah <- A[h,]
		
		#Calculate posterior hyperparameters
		
			#ah_tilde
			uh <- U[,h]		
			ah_t <- ah+sum(uh*Z)	
			
			#bh_tilde
			temp <- g^U
			temp <- temp[,-h]
			prod_gamma <- apply(temp,1,prod)
			xi_i <- unlist(mapply(rep,xi,as.numeric(table(id))))
			bh_t <- bh+sum(xi_i*prod_gamma)
			
			#ph_tilde
			cons1 <- C(ah,bh)/C(ah_t,bh_t)
			cons2 <- (pgamma(Ah[2],ah_t,bh_t)-pgamma(Ah[1],ah_t,bh_t))/(pgamma(Ah[2],ah,bh)-pgamma(Ah[1],ah,bh))
			cons <- cons1*cons2
			num <- (ph*exp(-(bh_t-bh)))
			den <- num+(1-ph)*cons
			ph_t <- num/den
		
		#Sample from one-inflated gamma truncated to region Ah
		g[h] <- rI1gammaT(1,ph_t,ah_t,bh_t,Ah)		
	}
	
	##Update beta
	b <- log(g)
	
	##Write gamma to file
	write(g,file=paste(output,"GAMMA.csv",sep=""),sep=",",ncolumns=H,append=TRUE)
	
	##Xi_i Full Conditional
	for (i in 1:n) {
		xi[i]<-rgamma()
	}
	write(xi,file=paste(output,"XI.csv",sep=""),sep=",",ncolumns=n,append=TRUE)
		
	##Metroplois Step For Phi
			
		#Sample from proposal distribution
		phi_prop<-rnorm(1,phi,phi_propvar)

		#We will have to either transform phi or only accept if it is positive (we can discuss this more later...)
		
		#Compute log acceptance ratio		
		log.r<-

		##Update phi
		if (log(runif(1))<log.r) {
			phi<-phi_prop
			acceptance_phi<-acceptance_phi+1 #keep track of acceptance rate
		}
	write(phi,file=paste(output,"PHI.csv",sep=""),sep=",",ncolumns=1,append=TRUE)
	
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

