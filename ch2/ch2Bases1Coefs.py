from sympy import *

#this script computes the coefficients of t^{2m+1}

alpha0=0
alpha1=1
alphaAll=[alpha0,alpha1]

for k in range(1,7):
    tmp=0
    for j in range(0,k+1):
        tmp+=alphaAll[j]*alphaAll[k-j]
    tmp*=-Rational(1,k+1)
    alphaAll.append(tmp)

print(alphaAll)