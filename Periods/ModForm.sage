# -*- coding: utf-8 -*-

load("PrimaryDecomp.sage")

class ModForm: # Virtual class, represents constant * B(d) * modform
    
    def __init__(self):
        #print("Base Constructor")
        self.Bd=1
    
    def a(self,n):
        d = self.Bd
        if n%d == 0:
            return self.a_before_Bd(n//d)
        return self.field(0)

    def qexp(self,n):
        res = zero_vector(self.field,n)
        d = self.Bd
        for i in range((n-1)//d+1):
            res[d*i] = self.a_before_Bd(i)
        return res
    
    def B(self,d):
        res = copy(self)
        res.Bd *= d
        res.N *= d
        res.eps = res.eps.extend(res.N)
        return res
    
    def Diamond(self,d):
        res = copy(self)
        res.constant *= self.field(self.eps(d))
        return res
    
    def W(self,Q):
        d = self.Bd
        N = self.newN # Level at which f is new
        if gcd(N,d) > 1:
            raise NotImplementedError("I know how to compute f|B(d)|W_Q  only for d corpime to the level of f.")
        Qd = gcd(Q,d)
        QN = gcd(Q,N)
        Q1 = Q//(Qd*QN)
        if gcd(Q,Q1*d*N//Q) > 1:
            raise ValueError("Forbidden value of Q in W_Q (Q and N/Q must be coprime)")
        g = self.W_initial_level(QN).B(Q1*d//Qd)
        field = self.field
        eps = self.eps.primitive_character().extend(N)
        epsQN,_ = CharDecomp2(eps,QN)
        g.constant *= field(Q1/Qd)^(self.k/2) * ((field(epsQN(d)) * field(eps(Qd))).conjugate())
        return g
