from sympy import *
from sympy.simplify.fu import TR8
import re
import numpy as np
from datetime import datetime

n=1
x,eps,h,C1=symbols("x,epsilon,h,C1",cls=Symbol)

u0=sqrt(2)*sin(n*pi*x)
# u1=h*eps*sqrt(2)/16*sin(3*n*pi*x)/(n**2*pi**2)

# lmd0=n**2*pi**2-3*eps/2

lmdUnknown=symbols("lambda",cls=Symbol)
uAll=[u0]
lmdAll=[lmdUnknown]

def R(k,uAll,lmdAll):
    """

    :param k:
    :param uAll:
    :param lmdAll:
    :return: R, using u0,u1,...,u_{k-1}
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

def extractCoefs(RHS):
    """

    :param RHS: RHS of deformation equation, without h
    :return:
    """

    expandedR=expand(RHS)
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

    :param retPairs: [k, coef of sin[(2k+1)n*pi*x] ],result from extractCoefs(RHS)
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
        newCTmp=-cTmp/((4*k**2+4*k)*n**2*pi**2)
        LInvRCoefs.append([k,newCTmp])
    return lmdNext,LInvRCoefs

# R1=R(1,uAll, lmdAll)
# coefPairs=extractCoefs(R1)
#
# lmdNext,LInvRCoefs=solveRHS(coefPairs)
# # pprint(lmdNext)
# pprint(coefPairs)




#init u1
# k=1
# R1=R(k,uAll,lmdAll)
# RAll.append(R1)
#
# RSum1=sum(RAll)
# coefPairs1=extractCoefs(RSum1*h)
#
# lmd0,LInvRCoefs1=solveRHS(coefPairs1)
# R1New=expand(RAll[1].subs([(lmdUnknown,lmd0)]))
# RAll[1]=R1New
#
#
#
# u1=0
# for item in LInvRCoefs1:
#     k,coefTmp=item
#     u1+=coefTmp*sin((2*k+1)*n*pi*x)

# u1+=C1*sin(n*pi*x)

# uAll.append(u1)
#
# lmdAll.insert(-1,lmd0)



#
# k=2
# R2=R(k,uAll,lmdAll)
# RAll.append(R2)
#
# RSum2=sum(RAll)
# tmp=expand(h*RSum2)
# # pprint(tmp.coeff(sin(pi*x)))
# coefPairs2=extractCoefs(RSum2*h)
# lmd1,LInvRCoefs2=solveRHS(coefPairs2)
# R2New=RAll[2].subs([(lmdUnknown,lmd1)])
# RAll[2]=R2New
# u2=0
# for item in LInvRCoefs2:
#     k,coefTmp=item
#     u2+=coefTmp*sin((2*k+1)*n*pi*x)
# pprint((u0+u1+u2).coeff(sin(5*pi*x)))


# # #iteration
M=5
for j in range(1,M+1):
    tIterStart=datetime.now()
    print("iteration "+str(j))
    RNew = R(j, uAll, lmdAll)
    RAll.append(RNew)
    RSum = sum(RAll)
    coefPairs = extractCoefs(RSum*h)
    lmdNext, LInvRCoefs = solveRHS(coefPairs)
    RNewReplaced=RAll[j].subs([(lmdUnknown,lmdNext)])
    RAll[j]=RNewReplaced
    uNext = 0
    for item in LInvRCoefs:
        k, coefTmp = item
        # pprint(k)
        # pprint(coefTmp)
        uNext += coefTmp * sin((2 * k + 1) * n * pi * x)
    uAll.append(uNext)

    lmdAll.insert(-1, lmdNext)
    tIterEnd=datetime.now()
    print("iteration "+str(j)+" time: ",tIterEnd-tIterStart)


# lmdTruncated=sum(lmdAll[:-1])
# pprint(uAll[2])
uSumTrucated=sum(uAll)+C1*sin(n*pi*x)

# #
coefPairsTruncated=extractCoefs(uSumTrucated)
# pprint(coefPairsTruncated)
# # pprint(uAll[0]+uAll[1])
eqn=0
for item in coefPairsTruncated:
    _,coefTmp=item
    eqn+=coefTmp**2/2
# #
eqn-=1
eqn=eqn.expand()
# pprint(eqn)
# #
# #
hVal=-1
epsVal=-10
# #
val=eqn.subs([(eps,epsVal),(h,hVal)])
# pprint(val)
# #
p=Poly(val,C1)
# #
C1Roots=np.roots( p.all_coeffs())
# # #
C10=C1Roots[0]

lmdSumTRuncated=sum(lmdAll[:-1])

lmdSumTRuncatedVal=lmdSumTRuncated.subs([(eps,epsVal),(h,hVal),(C1,C10)]).evalf()

pprint(lmdSumTRuncatedVal/(n**2*np.pi**2))