# -*- coding: utf-8 -*-

def PlanEq(U,V,degU,degV,K,u,v):
    eqprec = (degU+1)*(degV+1)
    A = Matrix(K,eqprec,eqprec)
    Ut = U.add_bigoh(eqprec)
    Vt = V.add_bigoh(eqprec)
    #print("pow...")
    powU=[U^0 for i in range(degV+1)]
    for i in range(1,degV+1):
        powU[i]=powU[i-1]*U
    powV=[V^0 for i in range(degU+1)]
    for i in range(1,degU+1):
        powV[i]=powV[i-1]*V
    #print("prods...")
    for i in range(degV+1):
        for j in range(degU+1):
            prod = powU[i]*powV[j]
            for k in range(eqprec):
                A[k,i*(degU+1)+j] = K(prod[k])
    #print("Ker...")      
    noyau = A.right_kernel().basis()
    #if not len(noyau)==1:
        #print("PlanEq : dimension du noyau "+str(len(noyau))+"\n")
    eq = sum([sum([noyau[0][i*(degU+1)+j]*u**i * v**j for j in range(degU+1)]) for i in range(degV+1)])
    eq = eq//gcd(eq,eq.derivative(u))
    eq = eq//gcd(eq,eq.derivative(v)) # remove square factors and factors that only depend on one variable
    #print("eq epluchee:",eq)
    return eq

def NewtonEq(eqU,U,V,t,hiqprec):
    DeqU=eqU.derivative()
    #print("Newton")
    prec=eqU(V).valuation()
    while prec<hiqprec:
        #print("precision",prec)
        prec*=2
        if prec>hiqprec:
            prec=hiqprec
        #print("esperee",prec)
        V=V.polynomial()+O(t^prec)
        test=eqU.map_coefficients(lambda b: b+O(t^prec))(V)
        testD=DeqU.map_coefficients(lambda b: b+O(t^prec))(V)
        V-=test/testD
        #print("val test",test.valuation())
        #print("val testD",testD.valuation())
        prec=2*(test.valuation()-testD.valuation())
        #print("obtenu",prec)
    return V

def PowerSeriesNthRoot(f,n,a0,x):
    if n==1:
        return f
    g = a0+O(x**(f.prec()))
    prec = 1
    while prec<f.prec():
        ft = f+O(x**(2*prec))
        gt = g+O(x**(2*prec))
        dg = (g-ft/g**(n-1))/n
        g -= dg.polynomial()
        prec *= 2
    return g

def liftFp2Z(x,p):
    n = x.lift()
    if 2*n > p:
        n -= p
    return n

def SerRedModp(f,Fp,Fpt):
    return Fpt([Fp(a) for a in f.list()]).add_bigoh(f.prec())

def NewtonGamma0_0(f,G,hiqprec,p): # f cuspform / Z of level G, compute eqn between f/dj and j mod p, use it to get more terms of f mod p by Newton, and finally lift to Z (/!\ take p large enough !)
    Fp=GF(p)
    Fpt.<t>=Fp[[]]
    Fpty.<y>=Fpt[]
    Fpuv.<u,v>=Fp[]
    Qq.<q>=QQ[[]]
    indG=G.index()
    g0=G.genus()
    nptes=len(G.cusps())
    degU=indG
    degV=2*g0-2+indG+nptes
    #print("E4 E6")
    E4=1+240*Fpt([sigma(n,3) for n in range(1,hiqprec)]).shift(1).add_bigoh(hiqprec+1)
    E6=1-504*Fpt([sigma(n,5) for n in range(1,hiqprec)]).shift(1).add_bigoh(hiqprec+1)
    U=(1-E6^2/E4^3)/1728
    J=1/U
    dJshift=J.derivative().shift(1)
    #print("qexp")
    fp=SerRedModp(f.qexp((degU+1)*(degV+1)),Fp,Fpt)
    V=(fp/dJshift).power_series()
    #print("eq")
    eq=PlanEq(U,V,degU,degV,Fp,u,v)
    #print("eq(U,.)")
    eqU=[0 for j in range(degU+1)]
    Ui=1
    for i in range(degV+1):
        for j in range(degU+1):
            eqU[j] += eq[i,j]*Ui
        Ui*=U
    eqU=Fpty(eqU)
    #print("Newton")
    V=NewtonEq(eqU,U,V,t,hiqprec)
    #print("fini, prec V:",V.prec())
    fp=(V*dJshift).power_series()
    #print("prec fp",fp.prec())
    #print("QQ")
    return Qq([liftFp2Z(fp[i],p) for i in range(fp.prec())]).add_bigoh(fp.prec())


def NewtonGamma0_1(f,f0,G,hiqprec,p): # same as above, but uses f0 (a cuspform / Z of level G) instead of dj
    Fp=GF(p)
    Fpt.<t>=Fp[[]]
    Fpty.<y>=Fpt[]
    Fpuv.<u,v>=Fp[]
    Qq.<q>=QQ[[]]
    indG=G.index()
    g0=G.genus()
    degU=indG
    degV=2*g0-2
    #print("E4 E6")
    E4=1+240*Fpt([sigma(n,3) for n in range(1,hiqprec)]).shift(1).add_bigoh(hiqprec+1)
    E6=1-504*Fpt([sigma(n,5) for n in range(1,hiqprec)]).shift(1).add_bigoh(hiqprec+1)
    U=(1-E6^2/E4^3)/1728
    #print("qexp")
    fp=SerRedModp(f.qexp((degU+1)*(degV+1)+1),Fp,Fpt)
    f0p=SerRedModp(f0,Fp,Fpt)
    V=fp/f0p
    #print("eq")
    eq=PlanEq(U,V,degU,degV,Fp,u,v)
    #print("eq(U,.)")
    eqU=[0 for j in range(degU+1)]
    Ui=1
    for i in range(degV+1):
        for j in range(degU+1):
            eqU[j] += eq[i,j]*Ui
        Ui*=U
    eqU=Fpty(eqU)
    #print("Newton")
    V=NewtonEq(eqU,U,V,t,hiqprec)
    #print("fini, prec V:",V.prec())
    fp=V*f0p
    #print("prec fp",fp.prec())
    #print("QQ")
    return Qq([liftFp2Z(fp[i],p) for i in range(fp.prec())]).add_bigoh(fp.prec())


def NewtonGamma0_2(f,f0,f1,G,hiqprec,p): # same as above, but uses f1/f0 instead of j
    Fp=GF(p)
    Fpt.<t>=Fp[[]]
    Fpty.<y>=Fpt[]
    Fpuv.<u,v>=Fp[]
    Qq.<q>=QQ[[]]
    g0=G.genus()
    degU=2*g0-2
    degV=2*g0-2
    #print("Red mod p")
    f0p=SerRedModp(f0,Fp,Fpt)
    f1p=SerRedModp(f1,Fp,Fpt)
    #print("U")
    U=f1p/f0p
    #print("qexp")
    fp=Fpt(f.qexp((degU+1)*(degV+1)+1))
    V=fp/f0p
    #print("eq")
    eq=PlanEq(U,V,degU,degV,Fp,u,v)
    #print("eq(U,.)")
    eqU=[0 for j in range(degU+1)]
    Ui=1
    for i in range(degV+1):
        for j in range(degU+1):
            eqU[j] += eq[i,j]*Ui
        Ui*=U
    eqU=Fpty(eqU)
    #print("Newton")
    V=NewtonEq(eqU,U,V,t,hiqprec)
    #print("fini, prec V:",V.prec())
    fp=V*f0p
    #print("prec fp",fp.prec())
    #print("QQ")
    return Qq([liftFp2Z(fp[i],p) for i in range(fp.prec())]).add_bigoh(fp.prec())


# Returns C such that if a in Q(zetaN) has all |cmplx embs|<B, then all |coeffs of a| (as a polynomial in zetaN) are < B*C.
def CycloBound(N):
    phiN=euler_phi(N)
    S=Integers(N).list_of_elements_of_multiplicative_group()
    V=Matrix(CC,phiN,phiN)
    for i in range(phiN):
        for j in range(phiN):
            V[j,i]=exp(2*I*RR.pi()*((i*S[j])%N)/N)
    #print(V.det())
    V=V^-1
    return max([sum([abs(V[i,j]) for j in range(phiN)]) for i in range(phiN)])


def mu_r(Fp,r):
    g=Fp.multiplicative_generator()
    g=g^((Fp.order()-1)/r)
    return [g^a for a in Integers(r).list_of_elements_of_multiplicative_group()]

def SerRedModP(f,Fp,z,Fpt):
    return Fpt([sum([Fp(a[i])*z^i for i in range(len(a.list()))]) for a in f.list()]).add_bigoh(f.prec())

def liftFP2Zzeta(X,L,zeta,p):
    A=[liftFp2Z(a,p) for a in sum([X[i]*L[i] for i in range(len(L))]).list()]
    return sum([A[i]*zeta^i for i in range(len(A))])


@parallel
def NewtonGamma1_2_sub(fq,f0p,U,r,hiqprec,z,degU,degV,t,u,v):
    Fpt=f0p.parent()
    Fp=Fpt.base_ring()
    Fpty.<y>=Fpt[]
    fp=SerRedModP(fq,Fp,z,Fpt)
    V=fp/f0p
    valV=V.valuation()
    a0=V[valV]
    V=V^r
    #print("eq")
    eq=PlanEq(U,V,degU,degV,Fp,u,v)
    #print("eq(U,.)")
    eqU=[0 for j in range(degU+1)]
    Ui=1
    for i in range(degV+1):
        for j in range(degU+1):
            eqU[j] += eq[i,j]*Ui
        Ui*=U
    eqU=Fpty(eqU)
    V=NewtonEq(eqU,U,V,t,hiqprec)
    #print("root")
    fp=PowerSeriesNthRoot(V.shift(-r*valV),r,a0,t).shift(valV)*f0p
    return fp

def NewtonGamma1_2(f,f0,f1,r,G,G0,hiqprec,p): # f of level G and nebentypus killed by r and coeffs in Q(zeta_r), f0,f1 of level G0 >= G and triv nebentypus, p must be 1 mod r.
    Fp=GF(p)
    Fpx.<x>=Fp[]
    Phi=Fpx(cyclotomic_polynomial(r))
    mu=[z[0] for z in Phi.roots()]
    mu = sorted(mu,key=lambda x:x.lift()) # Arrange in some order
    lagrange=[Phi//(x-z) for z in mu]
    lagrange=[lagrange[i]/lagrange[i](mu[i]) for i in range(len(mu))] # lagrange interpolents, in the same order
    Fpt.<t>=Fp[[]]
    Fpuv.<u,v>=Fp[]
    K=CyclotomicField(r)
    Kq.<q>=K[[]]
    zeta=K.gen()
    degU=2*G0.genus()-2
    degV=2*G.genus()-2
    #print("Red mod P")
    f0p=SerRedModp(f0,Fp,Fpt)
    f1p=SerRedModp(f1,Fp,Fpt)
    #print("U")
    U=f1p/f0p
    #print("qexp")
    fq=Kq(f.qexp((degU+1)*(degV+1)+1))
    todo=[(fq,f0p,U,r,hiqprec,z,degU,degV,t,u,v) for z in mu_r(Fp,r)]
    #print("Eq & Newton en parallele")
    res=NewtonGamma1_2_sub(todo)
    listfp=[]
    for X in res:
        listfp.append((X[1],X[0][0][5])) # list of pairs (fp,z), to be Chinese-lifted
    #print("Retour")
    listprec=[fp[0].prec() for fp in listfp]
    #print("precs",listprec)
    prec=min(listprec)
    listfp = sorted(listfp,key=lambda x:x[1].lift()) # sort in the same order as mu_r and lagrange
    #print("Relevement chinois")
    #print("Out")
    return Kq([liftFP2Zzeta([fp[0][i] for fp in listfp],lagrange,zeta,p) for i in range(prec)]).add_bigoh(prec)

