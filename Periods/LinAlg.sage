# -*- coding: utf-8 -*-

def IsApprox0(z,zeroplus):
    return bool(z.abs() <= zeroplus)

def IsApproxZeroVector(v,n,zeroplus):
    for i in range(n):
        if not IsApprox0(v[i],zeroplus):
            return false
    return true
    
def MatrixMax(A):
    max = A[0,0].abs()
    for i in range(A.nrows()):
        for j in range(A.ncols()):
            if (A[i,j].abs() > max):
                max = A[i,j].abs()
    return max
    
def arr(x):
    return floor(x+1/2)


def QR(A,zeroplus): # returns Q unitary, R upper triangular, and P permutation matrix s.t. Q*A*P=R
    m = A.nrows()
    n = A.ncols()
    s = min(m,n)
    field = A.base_ring()
    R = copy(A)
    P = identity_matrix(field,n)
    Q = identity_matrix(field,m)
    for i in range(s):
        # Recherche de la colonne de plus grande norme
        jmax = i
        max = field(0)
        for j in range(i,n):
            norm2 = sum([R[k,j].abs()**2 for k in range(i,m)])
            if norm2 > max:
                jmax = j
                max = norm2
        # Echange des colonnes i et jmax
        if i != jmax:
            for k in range(m):
                x = R[k,i]
                R[k,i] = R[k,jmax]
                R[k,jmax] = x
            for k in range(n):
                x = P[k,i]
                P[k,i] = P[k,jmax]
                P[k,jmax] = x
        # Y a-t-il vraiment du travail à faire sur cette colonne ?
        boulot = false
        for k in range(i+1,m):
            if R[k,i].abs()>zeroplus:
                boulot = true
                break
        if not boulot:
            continue
        # Calcul du vecteur u de Householder
        u = [R[k,i] for k in range(i,m)]
        z = field(1)
        if u[0] != field(0):
            z = u[0]/u[0].abs()
        norm2 = sum([ui.abs()**2 for ui in u])
        u[0] -= z*sqrt(norm2)
        norm2 = sum([ui.abs()**2 for ui in u])
        # Application de Householder
        for j in range(i,n):
            lambd = field(2)*sum([u[k-i].conjugate()*R[k,j] for k in range(i,m)])/norm2
            for k in range(i,m):
                R[k,j] -= lambd*u[k-i]
        for j in range(m):
            lambd = field(2)*sum([u[k-i].conjugate()*Q[k,j] for k in range(i,m)])/norm2
            for k in range(i,m):
                Q[k,j] -= lambd*u[k-i]
    return Q,R,P


def Ker(A,dim,zeroplus): # Return K s.t. A*K=0
    # If dim too large, the rightmost cols of K won't be in the Ker
    Q,R,P = QR(A.transpose(),zeroplus) # Q*tA*P = R -> tP*A*tQ=L -> A*tQ = P*L -> prendre les dernières colonnes de Q 
    n = A.ncols()
    K = matrix(A.base_ring(),n,dim)
    for j in range(dim):
        for i in range(n):
            K[i,j] = Q[n-1-j,i]
    return K

def EqnMatrix(A,rang,zeroplus): # Return K s.t. K*A=0
    Q,R,P = QR(A,zeroplus) # Q*A*P = R, R est mince donc a des 0s en bas -> Q*A=R*tP a des 0 en bas -> prendre les dernières lignes de Q
    m = A.nrows()
    K = matrix(A.base_ring(),m-rang,m)
    for j in range(m):
        for i in range(m-rang):
            K[i,j] = Q[i+rang,j]
    return K

def ImgCoords(A,rang,zeroplus): # return U s.t. the columns of A*U form a basis of Im(A).
    Q,R,P = QR(A.transpose(),zeroplus) # A*tQ=P*L -> les premières colonnes de tQ sont les combinaisions des colonnes de A qui forment une base de Im(A)
    n = A.ncols()
    U = matrix(A.base_ring(),n,rang)
    for j in range(rang):
        for i in range(n):
            U[i,j] = Q[j,i]
    return U
    

def Ker_with_echlist(A,dim,zeroplus): # (Not sure) Returns K=kernel of dimension dim in ech form with rows permuted, and list of len=dim telling where the pivots(=1) are in the cols of K
    m=A.nrows()
    n=A.ncols()
    B=copy(A)
    P=identity_matrix(A.base_ring(),n)
    P1loc = list(range(n))
    i = 0
    j = 0
    while(i < m and j < n):
        # Recherche du pivot le plus grand sur Li
        max = B[i,j].abs()
        jmax = j
        for k in range(j+1,n):
            if(B[i,k].abs() > max):
                max = B[i,k].abs()
                jmax = k
        # La ligne était-elle nulle ?
        if (IsApprox0(max,zeroplus)):
            i += 1
            continue
        # Echange des colonnes j et jmax
        if (jmax != j):
            for k in range(m):
                z = B[k,j]
                B[k,j] = B[k,jmax]
                B[k,jmax] = z
            for k in range(n):
                z = P[k,j]
                P[k,j] = P[k,jmax]
                P[k,jmax] = z
            k = P1loc[j]
            P1loc[j] = P1loc[jmax]
            P1loc[jmax] = k
        # Transvection : on tue Li
        inv = 1/B[i,j]
        for k in range(j+1,n):
            lambd = B[i,k] * inv
            # Ck <- Ck - lambd * Cj
            for l in range(i,m):
                B[l,k] -= lambd * B[l,j]
            for l in range(n):
                P[l,k] -= lambd * P[l,j]
        i += 1
        j+= 1
    # dim ker = n-j si pas d'erreur d'arrondi
    #if (dim != n-j):
     #   print "Dans Ker : dimension de noyau inattendue"
     #  print n-j
     #  print "Dimension attendue :"
     #  print dim
    # Extraction des résultats
    K = matrix(A.base_ring(),n,dim)
    for i in range(n):
        for j in range(dim):
            K[i,j] = P[i,j+n-dim]
    return K,[P1loc[i] for i in range(n-dim,n)]
 
def HasNonZeroKer(M,zeroplus):
    K = Ker(M,1)
    v = M*K
    for i in range(M.nrows()):
        if(not IsApprox0(v[i,0]),zeroplus):
            #print v[i,0]
            return false
    return true

def prettymatrix(M,zeroplus):
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            if IsApprox0(M[i,j].real(),zeroplus):
                M[i,j] = M[i,j].imag() * I
            if IsApprox0(M[i,j].imag(),zeroplus):
                M[i,j] = M[i,j].real()
                


def ExtractBasis(A): # Attention : optimisé pour les domaines exacts
    # (Not sure) returns b,Q,max where b=true iff A is surjective, Q is a matrix s.t. A*Q = Id (i.e. a section), and imax isi minimal such that the first imax cols of A span all of the image space.
    m = A.ncols()
    n = A.nrows()
    B = copy(A)
    P = identity_matrix(A.base_ring(),m)
    imax = 0
    if n > m:
        return false,P,imax
    for i in range(n):
        if B[i,i] == 0: # Il faut trouver un élément non nul sur cette ligne i, à droite
            for ix in range(i+1,m):
                if not (B[i,ix] == 0):
                    break
            if B[i,ix] == 0: # Pas trouvé, A n'est pas de rang n, on renvoie une erreur
                    return false,P,0
            if ix > imax:
                imax = ix
            # Echange des colonnes i et ix
            for j in range(i,n):
                x = B[j,i]
                B[j,i] = B[j,ix]
                B[j,ix] = x
            for j in range(m):
                x = P[j,i]
                P[j,i] = P[j,ix]
                P[j,ix] = x
        else:
            if i > imax:
                imax = i
        # Maintenant, B[i,i] n'est pas nul, on divise, transvecte pour tuer toute la ligne i
        Bii_inv = 1/B[i,i]
        B[i,i] = 1
        for k in range(i+1,n):
            B[k,i] = B[k,i] * Bii_inv
        for k in range(m):
            P[k,i] = P[k,i] * Bii_inv
        for j in range(m):
            if not (j == i):
                lambd = B[i,j]
                B[i,j] = 0
                for k in range(i+1,n):
                    B[k,j] = B[k,j] - lambd*B[k,i]
                for k in range(m):
                    P[k,j] = P[k,j] - lambd*P[k,i]
    Q = matrix(A.base_ring(),m,n)
    for i in range(m):
        for j in range(n):
            Q[i,j] = P[i,j]
    return true,Q,imax

                
def ImgExtract(A): # Renvoie une liste d'indices de colonnes de A engendrant Im A
#                    /!\  Optimisé pour domaines exacts
    m = A.nrows()
    n = A.ncols()
    B = copy(A)
    cols = range(n)
    i = 0
    j = 0
    while i < m:
        if B[i,j] == 0:
            # Recherche d'une tete de col non nulle
            k = j+1
            for k in range(j+1,n):
                if not B[i,k] == 0:
                    break
            # Cette ligne était-elle nulle ?
            if B[i,k] == 0:
                i = i+1
                continue
            # Non, on échange colonnes j et k
            for l in range(i,m):
                buff = B[l,j]
                B[l,j] = B[l,k]
                B[l,k] = buff
            buff = cols[j]
            cols[j] = cols[k]
            cols[k] = buff
        # Puis on transvecte
        inv = 1/B[i,j]
        for k in range(j+1,n):
            lambd = B[i,k]*inv
            # Ck <- Ck - lambd * Cj
            for l in range(i,m):
                B[l,k] = B[l,k] - lambd * B[l,j]
        i = i+1
        j = j+1
    return [cols[k] for k in range(j)]





def VectAppend(Mat,rg,nlignes,piv,v,tol): 
# Mat de forme
#
# x 0 0 0
# x 0 0 0
# x x 0 0
# x x x 0
# x x x 0
#
# Les hauts des colonnes sont donnés par le vecteur piv, rg = nb de colonnes non nulles
# le param nlignes est superflu ?
# 
# "Ajoute" v à Mat en le réduisant et en décalant au besion les colonnes pour conserver une forme échelonnée
# Renvoie bool, newMat, newrg, newpiv
# Où bool = true ssi. le rg a augmenté en ajoutant v

    w = copy(v)
    top = 0
    for j in range(rg+1):
        for i in range(top,piv[j]):
            if w[i].abs() > tol:
                # Décalage des colonnes de j à rg-1
                for k in range(rg-j):
                    jk = rg-1-k
                    for ik in range(piv[jk],nlignes):
                        Mat[ik,jk+1] = Mat[ik,jk]
                    piv[jk+1] = piv[jk]
                piv[j] = i
                for ik in range(i,nlignes):
                    Mat[ik,j] = w[ik]
                return true,Mat,rg+1,piv
        if j<rg:
            top = piv[j]
            lambd = w[top]/Mat[top,j]
            for i in range(top,nlignes):
                w[i] -= lambd * Mat[i,j]
    return false,Mat,rg,piv 
