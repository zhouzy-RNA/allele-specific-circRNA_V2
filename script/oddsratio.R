args<-commandArgs(T)
input <-paste0(args[1],".linear_circular_RNA.out.3.need.ref_alt") 
output <-paste0(args[1],".linear_circular_RNA.out.3.need.ref_alt.oddsratio") 

oddsratioWald.proc <- function(n00, n01, n10, n11, alpha = 0.05){	 
	myOR <- (n00 * n11)/(n01 * n10)
	siglog <- sqrt((1/n00) + (1/n01) + (1/n10) + (1/n11))
	logOR<-log(myOR)
	loglo<-exp(logOR-1.96*siglog)
	loghi<-exp(logOR+1.96*siglog)
	oframe<-data.frame(LowerCI = loglo, OR = myOR, UpperCI = loghi, alpha = alpha)
	oframe
}

my_data = read.table(file=input, header = FALSE, colClasses = c('character', 'character', 'character', 'character', 'character', 'numeric', 'numeric', 'numeric', 'numeric'))

new_data = my_data

rowcount<-dim(new_data)[1]

new_data$oddsratio    <- rep(0, rowcount) 
new_data$ci_low      <- rep(0, rowcount) 
new_data$ci_high     <- rep(0, rowcount) 


for(XX in 1:rowcount){ 
	AA <- mRNA_ref     <- new_data[XX,6]
	BB <- mRNA_alt     <- new_data[XX,7]
	CC <- circRNA_ref  <- new_data[XX,8]
	DD <- circRNA_alt  <- new_data[XX,9]
	
	if(AA+BB==0 || AA+CC==0 || BB+DD==0 || CC+DD==0)
	{ next	}
	
	if(AA==0 || BB==0 || CC==0 || DD==0)
	{ AA=AA+0.5;BB=BB+0.5;CC=CC+0.5;DD=DD+0.5 }

	if(AA>0 && BB>0 && CC>0 && DD>0) {
	
		newdf<- oddsratioWald.proc(AA,BB,CC,DD)
		if(newdf$OR<1){ 

			newdf <- oddsratioWald.proc(BB,AA,DD,CC)
		} 		
	
		new_data[XX,10]  <- newdf$OR   
		new_data[XX,11]  <-	newdf$LowerCI   
		new_data[XX,12]  <-	newdf$UpperCI  
	}
}
write.table(new_data, file = output, row.names = F, quote = F, sep="\t")
