from sympy import *

k=3

lmdStr="lambda0"
for i in range(1,k):
    lmdStr+=",lambda"+str(i)

lmdAll=symbols(lmdStr,cls=Symbol)

eps=symbols("epsilon",cls=Symbol)

uStr="u0"
for i in range(1,k):
    uStr+=",u"+str(i)


uAll=symbols(uStr,cls=Symbol)


RSumPart=0
for m in range(0,k):
    RSumPart+=lmdAll[m]*uAll[k-1-m]

dbSum=0
for m in range(0,k):
    tmp=0
    for j in range(0,m+1):
        tmp+=uAll[j]*uAll[m-j]
    dbSum+=uAll[k-1-m]*tmp
dbSum*=eps


RSumPart+=dbSum

pprint(RSumPart.expand())