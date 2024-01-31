from sympy import *
from sympy.simplify.fu import TR8
import re
import numpy as np


n=1
x,eps,h,C1=symbols("x,epsilon,h,C1",cls=Symbol)

u0=sqrt(2)*sin(n*pi*x)
# u1=h*eps*sqrt(2)/16*sin(3*n*pi*x)/(n**2*pi**2)

lmd0=n**2*pi**2-3*eps/2

lmdUnknown=symbols("lambda",cls=Symbol)
uAll=[u0]
lmdAll=[lmdUnknown]

def R(k,uAll,lmdAll):
    """

    :param k:
    :param uAll:
    :param lmdAll:
    :return: R
    """
    ret=diff(uAll[k-1],(x,2))
    # pprint(ret)

    # print("k="+str(k))
    for m in range(0,k):
        # pprint("lmdAll[m]="+str(lmdAll[m]))
        # pprint("uAll[k-1-m]="+str(uAll[k-1-m]))
        ret+=lmdAll[m]*uAll[k-1-m]


    dbSum = 0
    for m in range(0, k):
        tmp = 0
        for j in range(0, m + 1):
            tmp += uAll[j] * uAll[m - j]
        dbSum += uAll[k - 1 - m] * tmp
    dbSum *= eps

    ret+=dbSum

    ret1=TR8(ret.expand())
    return expand(ret1)

RAll=[0]

def extractCoefs(R):
    """

    :param R: RHS of deformation equation
    :return:
    """

    expandedR=expand(R)
    # pprint(expandedR)
    expandedRStr=str(expandedR)
    occurrences=re.findall(r"sin\(\d*\*?pi\*x\)",expandedRStr)
    occurrences=list(set(occurrences))
    # print(occurrences)
    numStr=[re.findall(r"(\d+)\*pi\*x",elem) for elem in occurrences]
    nums=[int(item[0]) for item in numStr if len(item)>0]
    maxCoef=np.max(nums)
    retPairs=[]
    k=0
    while (2*k+1)*n<=maxCoef:
        coefTmp=expandedR.coeff(sin((2*k+1)*n*pi*x))
        if coefTmp!=0:
            retPairs.append([k,coefTmp])
        k+=1

    # pprint(retPairs)
    return retPairs

def solveRHS(retPairs):
    """

    :param retPairs: [k, coef of sin[(2k+1)n*pi*x] ]
    :return:
    """
    pair0=retPairs[0]
    # pprint(pair0)
    lmdNext=solve(pair0[1],lmdUnknown)[0]
    # pprint(retPairs)
    # pprint(lmdNext)
    LInvRCoefs=[]
    for k in range(1,len(retPairs)):
        cTmp=retPairs[k][1]
        newCTmp=-cTmp/((4*k**2+4*k)*n**2*pi**2)*h
        LInvRCoefs.append([k,newCTmp])
    return lmdNext,LInvRCoefs

# R1=R(1,uAll, lmdAll)
# coefPairs=extractCoefs(R1)
#
# lmdNext,LInvRCoefs=solveRHS(coefPairs)
# # pprint(lmdNext)
# pprint(coefPairs)




#init u1
k=1
R1=R(k,uAll,lmdAll)
RAll.append(R1)

RSum1=sum(RAll)
coefPairs1=extractCoefs(RSum1)

lmd0,LInvRCoefs1=solveRHS(coefPairs1)
u1=0
for item in LInvRCoefs1:
    k,coefTmp=item
    u1+=coefTmp*sin((2*k+1)*n*pi*x)

u1+=C1*sin(n*pi*x)

uAll.append(u1)

lmdAll.insert(-1,lmd0)
# pprint(lmd0)



# # #interation
M=1
for j in range(2,M+1):
    RNew = R(j, uAll, lmdAll)
    RAll.append(RNew)
    RSum = sum(RAll)
    coefPairs = extractCoefs(RSum)
    lmdNext, LInvRCoefs = solveRHS(coefPairs)
    uNext = 0
    for item in LInvRCoefs:
        k, coefTmp = item
        # pprint(k)
        # pprint(coefTmp)
        uNext += coefTmp * sin((2 * k + 1) * n * pi * x)
    uAll.append(uNext)

    lmdAll.insert(-1, lmdNext)


# lmdTruncated=sum(lmdAll)
pprint(uAll[0])
uSumTrucated=sum(uAll)

coefPairsTruncated=extractCoefs(uSumTrucated)
pprint(uAll[0]+uAll[1])
# eqn=0
# for item in coefPairsTruncated:
#     _,coefTmp=item
#     eqn+=coefTmp**2/2
#
# eqn=eqn.expand()
# # pprint(eqn)
#
#
# hVal=-1
# epsVal=-5
#
# val=eqn.subs([(eps,epsVal),(h,hVal)])
#
# p=Poly(val,C1)
#
# C1Roots=np.roots( p.all_coeffs())
#
# print(C1Roots)