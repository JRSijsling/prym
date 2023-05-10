load("LinAlg.sage")
load("HeckeGen.sage")
load("qexp.sage")
load("EigenForm.sage")

@parallel
def CalcEigenFormC(tag,qexps,coords,hiqprec,field): # tag used to recover nebentypus (weird pickling bug : cannot do parallel with Dirichlet characters ???)
    d = len(coords)
    A = zero_vector(field,hiqprec)
    for i in range(1,hiqprec):
        for j in range(d):
            A[i] += coords[j]*field(qexps[j][i])
    return A

@parallel
def AllIntegrals(j,Newforms,list_chi):
    return [Newforms[j].Integral(chi) for chi in list_chi]

def GammaHPeriods(N,H,prec):
    if not N.is_squarefree():
        raise NotImplementedError("Non-squarefree level")
    zeroplus = RR(10)**-prec
    field = ComplexField(prec*log(10)/log(2))
    G = GammaH(N,H)
    g = G.genus()
    print("Genus",g)
    Eps = G.characters_mod_H() # Dir chars trivial on H
    ExpH = lcm([e.order() for e in Eps])
    KH = CyclotomicField(ExpH)
    Eps = [e.change_ring(KH) for e in Eps]
    S = [CuspForms(e) for e in Eps]
    hecke_gen_list = Find_Hecke_Gen(S,KH) # Generator of Hecke alg
    print("Gen Hecke",hecke_gen_list)
    
    print("Finding basis of homology consisting of twisted winding elements")
    p = 2*N + 1
    while not is_prime(p):
        p += 2*N
    M = ModularSymbols(G,2)
    B = M.boundary_map().matrix().transpose() # Matrix of boundary operator
    matS,list_proj_S = Ker_with_echlist(B,2*g,0) # Basis of cupsidal subspace,
                                        # and pivot list to help projection
    mat_proj_S = matrix(QQ,2*g,M.dimension()) # Matrix to help write down a cuspidal elt of M in terms of basis given by matS
    # If v in S, we can reduce it along the permuted-echelon basis of S using the pivots of list_proj_S
    # The matrix of this reduction process is :
    for i in range(2*g):
        mat_proj_S[i,list_proj_S[i]]=1
    # Tp-(p+1) maps M to S and is invertible on S
    TM = (p+1-M.hecke_matrix(p).transpose()) # matrix of p+1-T(p) on M
    TS = mat_proj_S* TM * matS # Induced operator on S
    matDrinfeld = TS**(-1) * mat_proj_S * TM # Explicit Manin-Drinfeld matrix
    fini = false
    list_chi = []
    pmax = 0
    K = QQ
    c = 1 # Along the loop below, K=Q(zeta_c), where c is minimal for all the characters chi to assume values in K
    M = ModularSymbols(G,2,base_ring=QQ)   
    Tgen_symb_M = Hecke_Gen_Matrix(M,hecke_gen_list)
    Tgen_symb = mat_proj_S * Tgen_symb_M * matS # Matrix of Hecke generator on S
    mat_winds = matrix(QQ,2*g,0)
    m = 0 # Modulus of the characters considered
    while not fini:
        m +=1
        if gcd(N,m) > 1:
            continue
        print("Now considering characters mod",m)
        # Enlarge K ?
        em = euler_phi(m) # chars mod m assume values in Q(zeta_em)
        if not c % em == 0: # Is Q(zeta_em) included in K=Q(zeta_c) ?
            c = lcm(c,em)
            K = CyclotomicField(c)
            print("Enlarging K to "+str(K))
            mat_winds = mat_winds.base_extend(K)
            matDrinfeld = matDrinfeld.base_extend(K)
            list_chi = [chi.base_extend(K) for chi in list_chi]
        for chi in DirichletGroup(m,K):
            if fini:
                break
            if not chi.is_primitive(): # Only want primitive chars
                continue
            s = Matrix(K,1,M.dimension()) # Coordinates of twisted winding element on basis of M
            for x in range(m):
                s += chi.bar()(-x) * (Matrix(M.modular_symbol([infinity,x/m]).list()).base_extend(K))
            s = matDrinfeld.base_extend(K) * (s.transpose()) # Now express it on the basis of S
            if s.is_zero():
                continue
            list_chi.append(chi)
            Ts = matrix(K,2*g,g) # Coordinates of T^i s on basis of s for i < 2g, where T is gen of Hecke
            for j in range(g):
                for i in range(2*g):
                    Ts[i,j] = s[i,0]
                s = Tgen_symb.base_extend(K) * s
            mat_winds = mat_winds.augment(Ts)
            max_rank,mat_S2wind,S2wind_imax = ExtractBasis(mat_winds) # If max_rank, then we have spanned all S, and mat_S2wind expresses the basis of S on the T^k s_chi.
            if max_rank:
                fini=true

    print("q-expansions to high accuracy")
    # q-expansion accuracy needed to get the periods with accuracy < zeroplus
    mmax = max([chi.modulus() for chi in list_chi])
    hiqprec=ceil(1.2 * log(10.) * prec * mmax * sqrt(N) / (2*RR.pi()) )
    print("q-adic accuracy",hiqprec)
    hinewbound = max([sigma(n,0)*ceil(sqrt(n)) for n in range(1,hiqprec)]) # Bound on the coeffs of a newform in this range
    Shi = [] # list of triples (basis,P,eps), where basis is the q-echelon basis of S_2(eps)^new, and P expresses newforms in this space on this basis
    M0 = [M for M in divisors(N) if M > 10 and CuspForms(M).new_submodule().dimension() >= 2]
    if len(M0)==0:
        raise NotImplementedError("Less than 2 newforms on Gamma0(M) for all M|N !")
    M0 = min(M0) # Most of the work done at this level, the other levels will use the first 2 forms from this level to find small modular equations
    S0 = CuspForms(M0).new_submodule()
    d = S0.dimension() # At least 2
    print("First, processing the",d,"forms of level Gamma0("+str(M0)+")")
    P = HeckeEigen(S0,field)
    Q = P**-1 # expresses the basis of S0 in terms of newforms of S0
    eps = S0.character() # Trivial character
    S0 = S0.basis()
    G0 = Gamma0(M0)
    S0hi = []
    f0 = S0[0]
    f0den = lcm([a.denominator() for a in f0.qexp(G0.sturm_bound()+1).coefficients()]) # denom of f
    B = ceil(hinewbound*f0den*sum([abs(Q[i,0]) for i in range(d)])) # Bound on coeffs of fden*f
    p = next_prime(2*B+2)
    f0 = NewtonGamma0_0(f0den*f0,G0,hiqprec,p)/f0den
    S0hi.append(f0)
    f1 = S0[1]
    f1den = lcm([a.denominator() for a in f1.qexp(G0.sturm_bound()+1).coefficients()]) # denom of f
    B = ceil(hinewbound*f1den*sum([abs(Q[i,1]) for i in range(d)])) # Bound on coefs of fden*f
    p = next_prime(2*B+2)
    f1 = NewtonGamma0_1(f1den*f1,f0,G0,hiqprec,p)/f1den
    S0hi.append(f1)
    for j in range(2,d):
        f = S0[j]
        fden = lcm([a.denominator() for a in f.qexp(G0.sturm_bound()+1).coefficients()]) # denom of f
        B = ceil(hinewbound*fden*sum([abs(Q[i,j]) for i in range(d)])) # Bound on coeffs of fden*f
        p = next_prime(2*B+2)
        S0hi.append(NewtonGamma0_2(fden*f,f0,f1,G0,hiqprec,p)/fden)
    Shi.append((S0hi,P,eps))
    
    for M in divisors(N):
        if M < 11:
            continue
        print("Now processing newforms of level Gamma_1("+str(M)+")")
        EpsM = [eps.primitive_character().extend(M) for eps in Eps if M % eps.conductor() == 0] # Characters that descend mod M
        EpsM = [eps.change_ring(CyclotomicField(eps.order())) for eps in EpsM] # Individually minimise value field
        Stodo = [CuspForms(eps).new_submodule() for eps in EpsM]
        Stodo = [s for s in Stodo if s.dimension()] # Drop empty spaces
        #  q-exps of forms on Gamma1(M). Do only one of each conjugate pair of nebentypus.
        MM0 = lcm(M,M0) # common level for the current forms and the forms f0 and f1
        G0 = Gamma0(MM0) # f0 and f1 are also modular of this level
        while len(Stodo) > 0:
            S = Stodo.pop() # Take a space
            eps = S.character()
            r = eps.order()
            if r==1 and M==M0:
                continue # Level M0, trivial character has already been done
            if r > 2:
                # eps != epsbar, remove epsbar space from todo list
                for i in range(len(Stodo)):
                    if Stodo[i].character() == eps.bar():
                        Stodo.pop(i)
                        break
            C = CycloBound(r)
            P = HeckeEigen(S,field)
            Q = P**-1
            d = S.dimension()
            S = S.basis()
            Sepshi = []
            GH = GammaH(MM0,[x for x in range(MM0) if gcd(x,MM0)==1 and x%M in eps.kernel()]) # forms in S are also modular of that level
            for j in range(d):
                f = S[j]
                fden = lcm([a.denominator() for a in f.qexp(GH.sturm_bound()+1).coefficients()]) # denom de f
                B = ceil(C*hinewbound*fden*sum([abs(Q[i,j]) for i in range(d)])) # Borne sur les coeffs de fden*f
                p = next_prime(2*B+2)
                while ExpH>1 and not (p%ExpH)==1:
                    p = next_prime(p+1)
                Sepshi.append(NewtonGamma1_2(fden*f,f0,f1,r,GH,G0,hiqprec,p)/fden)
            Shi.append((Sepshi,P,eps))
    print("Got exact q-expansions, embedding them into C")
    todo=[]
    hiqprec=min([min([f.prec() for f in s[0]]) for s in Shi])
    nebens = [s[2] for s in Shi]
    for i in range(len(Shi)):
        qexps,P,eps = Shi[i]
        for j in range(len(qexps)):
            todo.append((i,qexps,P.column(j),hiqprec,field))
    res = CalcEigenFormC(todo)
    Newforms = []
    for X in res:
        coeffs = X[1]
        tag = X[0][0][0]
        eps = nebens[tag]
        M = eps.modulus()
        Newforms.append(EigenForm(coeffs,M,2,eps,field))
        if eps.order()>2:
            Newforms.append(EigenForm([z.conjugate() for z in coeffs],M,2,eps.bar(),field))
    nNew = len(Newforms)
    print("Computing integrals")
    todo = [(j,Newforms,list_chi) for j in range(nNew)]
    Ints = ["NO DATA" for j in range(nNew)]
    res = AllIntegrals(todo)
    for X in res:
        Ints[X[0][0][0]] = X[1]
    print("Computing periods")
    matPeriodes = Matrix(field,g,2*g) # Rows = eigenforms, columns = homology loops
    i = 0
    # Must apply B(d) operators too. Order of the eigenforms:
    # newform_0, .. , newform_0 | B(N/M0), newform_1, .. , newform_1 | B(N/M1), ..
    # where Mi = level of newform_i, and the divisors of N/Mi are arranged in increasing order
    for n in range(nNew):
        F = Newforms[n]
        M = F.N
        D = divisors(N//M)
        D.sort()
        nD = len(D)
        nchi = len(list_chi)
        T = sum([t[0]*F.TnMatBdSpace(t[1],D) for t in hecke_gen_list]) # Matrix of Hecke gen on the F|B(d), d|N/M
        Tpow = T.powers(g)
        #ints = F.Integrals(list_chi) # Compute integrals along modular symobls s_chi in parallel
        for d in range(nD):
            # Periods of F|B(D[d])
            for j in range(2*g):
                per = 0
                for k in range(nchi):
                    for m in range(g):
                        for d2 in range(nD):
                            per += field(mat_S2wind[k*g+m,j]) * Tpow[m][d2,d] * field(list_chi[k](D[d2])/D[d2]) * Ints[n][k]
                matPeriodes[i,j] = per
            i += 1
    return matPeriodes
