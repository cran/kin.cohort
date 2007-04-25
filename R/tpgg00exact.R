`tpgg00exact` <-
function(f,siz){
#
# f:vector with allele frequencies
# siz: vector with respective number of genotypes in proband
#
	if (length(f)==1) {
      if (siz[1]==2)
 	      tpm<-cp1exact(f)    #prob of carrier case
      else if (siz[1]==3)
 	      tpm<-cp32exact(f)
      else
         stop("only 2 or 3 genotypes allowed")

   }	else {
      cp<-list(NULL,NULL)
      for (l in 1:2){

         if (siz[l]==2)
     	      cp[[l]] <- cp1exact(f[l])                 #  2 genotypes -> 2
         else if (siz[l]==3)
    	      cp[[l]] <- cp32exact(f[l])                #  3 genotypes -> 2
         else
            stop("only 2 or 3 genotypes allowed")
      }

      ngeno.pro<-siz[1]*siz[2]
   	tpm = array(data=0, dim=c(ngeno.pro, 4, 2))
      gen1<- as.numeric(gl(siz[1],siz[2])) # 11 22 33
      gen2<- rep(1:siz[2],siz[1])          # 12 12 12
      ngeno.rel<- 4
      car1<- as.numeric(gl(2,2))           # 11 22
      car2<- rep(1:2,2)                    # 12 12
      for (k in 1:2) {
         for (i in 1:ngeno.pro) {
              for (j in 1:ngeno.rel) {
                   tpm[i,j,k] <- cp[[1]][gen1[i],car1[j] ,k]*cp[[2]][gen2[i],car2[j] ,k]

      } } }

   }
   tpm
}

