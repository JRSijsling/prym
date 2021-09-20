# -*- coding: utf-8 -*-

def PrimaryDecomp(n):
    return [p^e for (p,e) in factor(n)]

def IntDecomp2(n,Q): # returns nQ,nQ' s.t. nQ*nQ'=n and nQ=gcd(n,Q^oo)
    nQ=1
    nR=1
    for (p,e) in factor(n):
        if p.divides(Q):
            nQ*=p^e
        else:
            nR*=p^e
    return (nQ,nR)

def CharDecomp(eps):
    d=dict()
    for e in eps.decomposition():
        d[e.modulus()]=e
    return d

def CharDecomp2(eps,Q): # Let N=eps.modulus(). If Q||N, let R=N//Q; returns decomp of eps as epsQ * epsR
    N=eps.modulus()
    R=N//Q
    if N%Q or gcd(Q,R)>1:
        raise ValueError("must have Q||N")
    dec=eps.decomposition()
    R=eps.modulus()//Q
    epsQ=prod([e.extend(Q) for e in dec if e.modulus().divides(Q)])
    if epsQ==1: # Make sure we get a char, even if empty product above
        epsQ=DirichletGroup(1)[0]
    epsR=prod([e.extend(R) for e in dec if e.modulus().divides(R)])
    if epsR==1: # Make sure we get a char, even if empty product above
        epsR=DirichletGroup(1)[0]
    return (epsQ,epsR)

def CharFlip(eps,Q): # Return epsR / epsQ
    N=eps.modulus()
    R=N//Q
    if N%Q or gcd(Q,R)>1:
        raise ValueError("must have Q||N")
    G = eps.parent()
    res = G.one()
    for e in eps.decomposition():
        if Q%e.modulus():
            res = res*G(e)
        else:
            res = res/G(e)
    return res

def GaussSum(eps,field):
    # disagreement with SAGE
    if eps.modulus()==1:
        return 1
    if field.is_exact():
        return eps.gauss_sum()
    return eps.gauss_sum_numerical(prec=field.prec())
