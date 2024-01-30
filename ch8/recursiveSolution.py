from sympy import *
from sympy.simplify.fu import TR8
import re
import numpy as np


n=1
x,eps,h,C1=symbols("x,epsilon,h,C1",cls=Symbol)

u0=sqrt(2)*sin(n*pi*x)
u1=C1*sin(n*pi*x)+h*eps*sqrt(2)/16*sin(3*n*pi*x)/(n**2*pi**2)

lmd0=n**2*pi**2-3*eps/2

lmdUnknown=symbols("lambda",cls=Symbol)
uAll=[u0,u1]
lmdAll=[lmd0,lmdUnknown]

def R(k,uAll,lmdAll):
    """

    :param k:
    :param uAll:
    :param lmdAll:
    :return: R
    """
    ret=diff(uAll[k-1],(x,2))

    for m in range(0,k):
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
    return ret1

RAll=[0,TR8(R(1,[u0],[lmd0]).expand())]

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
    lmdNext=solve(pair0[1],lmdUnknown)
    # pprint(retPairs)
    pprint(lmdNext)
    LInvRCoefs=[]
    for k in range(1,len(retPairs)):
        cTmp=retPairs[k][1]
        newCTmp=-cTmp/((4*k**2+4*k)*n**2*pi**2)*h
        LInvRCoefs.append([k,newCTmp])
    return lmdNext,LInvRCoefs

R1=RAll[1]
coefs=extractCoefs(R1)
lmdNext,LInvRCoefs=solveRHS(coefs)

pprint(lmdNext)
pprint(LInvRCoefs)
