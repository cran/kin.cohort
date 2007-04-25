`cp32exact` <-
function(p){
# function for computing mendelian conditional probabilities
#
# colapses 3 genotypes into noncarrier/carrier
#
# p = allele frequency
# tp = 2 by 2 by 2 array of conditional probabilities
#      tp(i,j,k) gives the conditional probability of genotype j
#       given the proband has genotype i for k-th type of relationship
#       i = 1(non-carrier), 2(carrier)
#       j = 1(non-carriers), 2(carriers)
#       k = 1(parent or offspring), 2(sibling)
q <- 1-p
tp <- array(data=0, dim=c(3, 2, 2))

# parent or offspring
tp[1,1,1] = q;       #  aa|aa
tp[1,2,1] = p;       #  Aa|aa + (0 x AA|aa )

tp[2,1,1] = q/2      #  aa|Aa                          #  p*q^2/(2*p*q);
tp[2,2,1] = 1-q/2    #  (p/2 x AA|Aa) + (1/2 x Aa|Aa)  #  p*q*(1+p)/(2*p*q);

tp[3,1,1] = 0;       #  aa|AA
tp[3,2,1] = 1;       #  (q x Aa|AA) + (p x AA|AA)

# sibling
tp[1,1,2] = (1/4)*(1+q)^2;               # aa|aa
tp[1,2,2] = (1/4)*p^2 + (1/2)*p*(1+q);   # AA|aa + Aa|aa

tp[2,1,2] = (1/4)*q*(1+q);                 ### (1/4)*q*(1+q) + 1/16  aa|Aa
tp[2,2,2] = (1/4)*p*(1+p) + (1/2)*(1+p*q); ### (1/4)*p*(1+p) + 1/16  AA|Aa   (3+4*p*q)/8  Aa|Aa

tp[3,1,2] = (1/4)*q^2                           # aa|AA          # 1-tp[3,2,2];
tp[3,2,2] = (p + (1/4)*q^2) + (p*q + (1/2)*q^2);# AA|AA + Aa|AA  # (1/4)*(1+p)^2+(1/2)*q*(1+p);

tp
}

