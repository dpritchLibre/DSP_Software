###Load Libraries (may need a few extra functions...)

###Read in Data (format data into the Y, U and X matrices)

###Create Data Objects (define objects like n, n_i)

###Matrix Objects (define objects )
eyen<-diag(n)
onen<-rep(1,n)

###Initial Values
phi<-1
gamma<-rep(1,H)
xi<-rep(1,n)

###Hyperparameters


###Metropolis Objects
phi_propvar<-1
acceptance_phi<-0

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

	##Gamma_h full conditional (we will probably want to loop over this)
	for (h in 1:H) {
		gamma[h]<-		
	}
	write(gamma,file=paste(output,"GAMMA.csv",sep=""),sep=",",ncolumns=H,append=TRUE)
	
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

