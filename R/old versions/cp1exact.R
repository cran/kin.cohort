`cp1exact` <-
function(p){
# function for computing mendelian conditional probabilities
# p = allele frequency
# tp = 2 by 2 by 2 array of conditional probabilities
#      tp(i,j,k) gives the conditional probability of genotype j
#       given the proband has genotype i for k-th type of relationship
#       i = 1(non-carrier), 2(carrier)
#       j = 1(non-carriers), 2(carriers)
#       k = 1(parent or offspring), 2(sibling)

q = 1-p
tp = array(data=0, dim=c(2, 2, 2))

# parent or offspring
tp[1,1,1] = q
tp[1,2,1] = p  # 1-tp[1,1,1] 

tp[2,1,1] = p*q^2/(p^2+2*p*q)
tp[2,2,1] = p*(1+p*q)/(p^2+2*p*q) # 1-tp[2,1,1] 

# sibling
tp[1,1,2] = (1/4)*(1+q)^2
tp[1,2,2] = (1/4)*p^2+(1/2)*p*(1+q) # 1-tp[1,1,2] 

tp[2,1,2] = ((1/4)*p^2*q^2 + (1/2)*p*q^2*(1+q))/(p^2+2*p*q)
tp[2,2,2] = ((1/4)*p^2*(1+p)^2 + p^2*q*(1+p)+p*q*(1+p*q))/(p^2+2*p*q) # 1-tp[2,1,2] 

tp
}

