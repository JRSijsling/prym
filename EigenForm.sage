# -*- coding: utf-8 -*-

load("ModForm.sage")

class EigenForm(ModForm):
    def __init__(self,coeffs,N=1,k=2,eps=1,field=CC):
        #ModForm.__init__(self)
        self.Bd = 1
        self.coeffs = coeffs
        self.prec=len(coeffs)
        if coeffs[1] != 1:
            for i in range(self.prec):
                self.coeffs[i] /= coeffs[1]
        self.constant = coeffs[1]
        self.prec=len(coeffs)
        self.N=N
        self.newN=N
        self.k=k
        self.eps=eps
        if eps==1:
            self.eps=DirichletGroup(N)[0]
        self.field=field
    
    def a_before_Bd(self,n):
        return self.constant * self.coeffs[n]
    
    def Hecke_list_vap(self,hecke_list): # eigenval of self // Hecke operator represented by hecke_list, where [... , [a,n] , ...] represents ... + a*T(n) + ...
        v=0
        for x in hecke_gen_list:
            if gcd(self.Bd,x[1])>1:
                raise NotImplementedError("f|Bd|Tn not implemented yet when gcd(d,n)>1.")
            v+= x[0]*self.coeffs[x[1]]
        return v

    def W_eigenval_primary(self,q): # W_q pseudo-eigenval of f (as opposed to f|Bd). /!\ Assumes q||N is prime
        N = self.N // self.Bd
        epsq = CharDecomp(self.eps.primitive_character().extend(N))[q]
        if epsq.is_trivial():
            return -q^(1-self.k/2)*((self.coeffs[q]/self.coeffs[1]).conjugate())
        return GaussSum(epsq,self.field) * ((self.coeffs[q]/self.coeffs[1]).conjugate()) / (q^(self.k/2))

    def W_eigenval(self,Q): # W_Q pseudo-eigenval of f (as opposed to f|Bd). /!\ Assumes Q||N is squarefree
        N = self.N // self.Bd
        dec = CharDecomp(self.eps.primitive_character().extend(N))
        return prod([self.field(dec[q](Q/q))*self.W_eigenval_primary(q) for q in PrimaryDecomp(Q)])

    def W_initial_level(self,Q): # W_Q operator action of f (as opposed to f|Bd). /!\ Assumes Q||newN is squarefree
        N = self.newN
        if not N%Q==0:
            raise ValueError("Forbidden Q in W_Q(Newform) : Q must divide N.")
        R = N // Q
        if gcd(Q,R) > 1:
            raise ValueError("Forbidden Q in W_Q(Newform) : Q and N/Q must be coprime.")
        eps = self.eps.primitive_character().extend(N)
        field = self.field
        if Q==1:
            g = EigenForm(self.coeffs, N, self.k, eps, field)
            g.constant = self.constant
            return g
        l = self.W_eigenval(Q)
        if Q==N: # Speed optimisation, not necessary
            g = EigenForm([z.conjugate() for z in self.coeffs], N, self.k, eps.bar(), field)
        else:
            epsQ,epsR = CharDecomp2(eps,Q)
            epsQbar = epsQ.bar()
            prec = len(self.coeffs)
            res = zero_vector(field,prec)
            for n in range(1,prec):
                nQ,nR = IntDecomp2(n,Q)
                res[n] = field(epsR(nQ)) * field(epsQbar(nR)) * (self.coeffs[nQ].conjugate()) * self.coeffs[nR]
            g = EigenForm(res, N, self.k, epsQbar.extend(N)*epsR.extend(N), field)
        g.constant = self.constant * l
        return g
    ## source : T. Asai, On the Fourier coefficients of automorphic forms at various cusps and some applications to Rankin’s convolution, thm2
    
    def Fricke_twist_eigenval(self,chi): # pseudo-eigenval of self twisted by chi wrt. Fricke (of level N*l^2, l=chi.modulus())
        if self.Bd != 1:
            raise NotImplementedError("Not implemented for f|Bd.")
        if not chi.is_primitive():
            raise NotImplementedError("Only implemented for primitive chi.")
        l = chi.modulus()
        return self.W_eigenval(self.newN) * self.field(self.eps(l)) * self.field(chi(-self.newN)) * GaussSum(chi,self.field) / GaussSum(chi.bar(),self.field)
    ## source : Twists of Newforms and Pseudo-Eigenvalues of W-Operators, Atkin + Li, p. 228
 
    def Integral(self,chi): # 2 i pi * the integral of self along sum_{x mod l} chibar(-x) {oo,x/l}, where chi.modulus()=l. /!\ Assumes weight==2.
        if self.k != 2:
            raise NotImplementedError("Only implemented for weight 2 (got weight "+str(self.k)+")")
        d = self.Bd
        l = chi.modulus()
        field = self.field
        sqrtN = field(sqrt(self.newN))
        lFchi = self.Fricke_twist_eigenval(chi)
        q = exp(-2*field.pi()/(l*sqrtN))
        S = 0
        for n in range(len(self.coeffs)-1,0,-1):
            a = self.coeffs[n]*self.field(chi(n))
            S += (a-lFchi*(a.conjugate()))/n
            S *= q
        return field(chi(d)/d) * S * l / GaussSum(chi,field)

    @parallel
    def _par_Integral(self,chi,tag):
        return self.Integral(chi)
    def Integrals(self,list_chi):
        n = len(list_chi)
        ints = zero_vector(self.field,n)
        todo = [(list_chi[i],i) for i in range(n)]
        res = self._par_Integral(todo)
        for X in res:
            ints[X[0][0][1]] = X[1]
        return ints

    
    def UpMatBdSpace(self,p,D):
        # Given the list D of divisors of some N2 and a prime p,
        # returns matrix of the U(p) operator of level N*N2 on the space with basis the self|B(d), d|N2, with the same order as D
        N2 = max(D)
        if self.Bd != 1:
            raise ValueError("UpMatBdSpace should not be called on F|B(d). Something fishy is happening.")
        if gcd(N2,self.N)>1:
            raise ValueError("UpMatBdSpace must be called with D = divisors of some N2 coprime to N.")
        nD = len(D)
        field = self.field
        U = Matrix(field,nD,nD)
        for j in range(nD):
            d = D[j]
            # Compute (*) = self | B(d) | U(p) (U at level NN')
            if d%p==0:
                # Case p|d : (*) = self | B(d/p)
                U[D.index(d//p),j] = 1
            else:
                # Case p !| d: (*) = a_p self | B(d) - p^(k-1) eps(p) self | B(dp)
                U[j,j] = self.coeffs[p]
                if self.N % p:
                    U[D.index(d*p),j] = -field(p^(self.k-1)*self.eps(p)) # This is 0 if p|N.
        return U

    def TnMatBdSpace(self,n,D):
        # Given the list D of divisors of some N2 coprime to N = new level of self and some integer n,
        # returns matrix of the T(n) operator of level N*N2 on the space with basis the self|B(d), d|N2, with the same order as D
        if self.Bd != 1:
            raise ValueError("TnMatBdSpace should not be called on F|B(d). Something fishy is happening.")
        N = self.N
        N2 = max(D)
        if gcd(N2,self.N)>1:
            raise ValueError("TnMatBdSpace must be called with D = divisors of some N2 coprime to N.")
        n2,n1 = IntDecomp2(n,N*N2) # n = n1*n2 coprime, n2=gcd(n,N2^oo)
        # So at level N*N2, T(n) = T(n1)*U(n2), and T(n1) commutes with B(d) for d|N2.
        fa2 = factor(n2)
        nD = len(D)
        Tn1 = self.coeffs[n1] * identity_matrix(self.field,nD,nD) # Matrix of T(n1)
        # n2 = prod(p_i^e_i) -> U(n2) = prod(U(p_i)^e_i).
        return Tn1 * prod([(self.UpMatBdSpace(p,D))^e for (p,e) in factor(n2)])
        

#@parallel
#def CalcNewformAndIntegrate(forms,prec,Pj,eps,list_chi,field):
#    d=len(Pj)
#    print eps
#    coeffs=[sum([Pj[i]*forms[i][m] for i in range(d)]) for m in range(prec)]
#    F=EigenForm(coeffs,eps.modulus(),2,eps,field)
#    ints=[]
#    for chi in list_chi:
#        ints.append(F.Integral(chi))
#    return [F.coeffs,F.eps.element(),ints]
