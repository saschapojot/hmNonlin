import numpy as np
from sympy import *
from datetime import datetime
import matplotlib.pyplot as plt
# from findiff import FinDiff
#this script computes the Vm terms for bases 1/(1+t)^{m}

h=symbols("h",cls=Symbol,real=True)

tau=symbols("tau",cls=Symbol)

t=symbols("t",cls=Symbol)

V0Func=1-1/(1+tau)
V0Val=V0Func.subs(tau,t)
H=1/(1+tau)
def chi(m):
    if m==0 or m==1:
        return 0
    if m >=2:
        return 1
VFuncsAll=[V0Func]#variable is tau
RFuncsAll=[0]#variable is tau
VValsAll=[V0Val]
tSymbolStart=datetime.now()
mEnd=40
for m in range(1,mEnd+1):
    tStart=datetime.now()
    RmFunc=diff(VFuncsAll[m-1],tau)
    for j in range(0,m):
        RmFunc+=VFuncsAll[j]*VFuncsAll[m-1-j]
    RmFunc-=(1-chi(m))
    RmFunc=apart(RmFunc,tau)
    RFuncsAll.append(RmFunc)
    Vm=chi(m)*VValsAll[m-1]+h*integrate(RmFunc*H,(tau,0,t))/(1+t)
    VValsAll.append(Vm)
    VFuncsAll.append(Vm.subs(t,tau))
    tEnd=datetime.now()
    print("computing order "+str(m)+": ",tEnd-tStart)
tSymbolEnd=datetime.now()

print("symbolic computation time: ",tSymbolEnd-tSymbolStart)

#use apart() to do partial  fraction expansion

# pprint(apart(VValsAll[1],t))#
VTruncated=sum(VValsAll)
Vddt=diff(VTruncated,(t,2))
Vdddt=diff(VTruncated,(t,3))
# print(VTruncated)
Vddt0=Vddt.subs(t,0)
Vdddt0=Vdddt.subs(t,0)

V2=lambdify(h,Vddt0,"numpy")
V3=lambdify(h,Vdddt0,"numpy")


hValsAll=np.linspace(-1.9,-0.1,20)
colors=["black","lightcoral","forestgreen","darkslategrey","navy","skyblue"]

V02=[V2(h) for h in hValsAll]
V03=[V3(h) for h in hValsAll]

plt.figure()
plt.plot(hValsAll,V02,color="blue",label="$V^{(2)}(0)$")
plt.plot(hValsAll,V03,color="red",label="$V^{(3)}(0)$")
plt.legend(loc="best")
plt.title("h curve, order "+str(m))
plt.savefig("ch2bases2hcurve.png")






# plt.figure()
#
# plt.plot(hValsAll,Vdd0,color="red",label="$V^{''}(0)$")
# plt.plot(hValsAll,Vddd0,color="blue",label="$V^{'''}(0)$")
#
# plt.xlabel("$\hbar$")
# plt.legend(loc="best")
# plt.savefig("bases2hcurve.png")
#
#
#
#
#
# plt.figure()
# plt.plot(hValsAll,Vdd0,color="red",label="$V^{''}$")
# #
