from sympy import *

#this script generates the coefficient mu_{0}^{m,n}

h=symbols("h",cls=Symbol,real=True)

m=1
n=1

mu0=0
for j in range(0,m-n+1):
    mu0+=binomial(n-1+j,j)*(1+h)**j

mu0*=(-h)**n

pprint(factor(mu0.expand()))