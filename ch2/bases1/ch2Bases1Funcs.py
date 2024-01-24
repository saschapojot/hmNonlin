import numpy as np
from sympy import *
from datetime import datetime
import matplotlib.pyplot as plt



#this script computes the Vm terms
h=symbols("h",cls=Symbol,real=True)

tau=symbols("tau",cls=Symbol)

t=symbols("t",cls=Symbol)

V0Func=tau

V0Val=V0Func.subs(tau,t)

def chi(m):
    if m==0 or m==1:
        return 0
    if m >=2:
        return 1

VFuncsAll=[V0Func]#variable is tau
RFuncsAll=[0]#variable is tau
VValsAll=[V0Val]
tSymbolStart=datetime.now()
mEnd=4
for m in range(1,mEnd+1):
    RmFunc=diff(VFuncsAll[m-1],tau)
    for j in range(0,m):
        RmFunc+=VFuncsAll[j]*VFuncsAll[m-1-j]
    RmFunc-=(1-chi(m))
    RmFunc=apart(RmFunc,t)
    RFuncsAll.append(RmFunc)
    Vm=chi(m)*VValsAll[m-1]+h*integrate(RmFunc,(tau,0,t))
    VValsAll.append(Vm)
    VFuncsAll.append(Vm.subs(t,tau))
tSymbolEnd=datetime.now()
print("symbolic computation time: ",tSymbolEnd-tSymbolStart)
def combineSameOrders(m):
    """

    :param m: order of Vm
    :return:
    """
    Vm=VValsAll[m]
    maxOrder=2*m+1
    outExpression=0
    for j in range(1,maxOrder+1):
        coef=expand(Vm).coeff(t**j)
        coef=factor(coef)
        # pprint("j="+str(j)+", coef="+str(coef))
        outExpression+=coef*t**j
    return outExpression


VTruncated=sum(VValsAll)


V=lambdify([t,h],VTruncated,"numpy")

hValsAll=[-1,-1/2,-1/10]
colors=["black","lightcoral","forestgreen","darkslategrey","navy"]

tValsAll=np.linspace(0,10,100)
plt.figure()

for i in range(0,len(hValsAll)):
    est=[V(tVal,hValsAll[i]) for tVal in tValsAll]
    # print(est)
    plt.plot(tValsAll,est,c=colors[i])

plt.scatter(tValsAll,np.tanh(tValsAll),color="purple",s=6)

plt.savefig("ch2Fig2.1.png")