from sympy import *

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

for m in range(1,4):
    RmFunc=diff(VFuncsAll[m-1],tau)
    for j in range(0,m):
        RmFunc+=VFuncsAll[j]*VFuncsAll[m-1-j]
    RmFunc-=(1-chi(m))
    RFuncsAll.append(RmFunc)
    Vm=chi(m)*VValsAll[m-1]+h*integrate(RmFunc,(tau,0,t))
    VValsAll.append(Vm)
    VFuncsAll.append(Vm.subs(t,tau))

pprint(expand(VValsAll[3]))

