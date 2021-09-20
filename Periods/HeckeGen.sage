# -*- coding: utf-8 -*-

def VAP_distinctes(matlist,common_base_field):
    size = 0
    polys = []
    for M in matlist:
        charP = M.charpoly()
        if charP.gcd(charP.derivative()).degree() > 0:
            return false # Une matrice a une valeur propre multiple
        polys.append(charP)
        size = size + 1
    for i in range(size):
        for j in range(i):
            if (polys[i].base_extend(common_base_field)).gcd(polys[j].base_extend(common_base_field)).degree() > 0:
                return false # Deux matrices ont une valeur propre commune
    return true


def Find_Hecke_Gen(list_S_by_nt,common_base_field,avoid_N=1):
    # Trouve un générateur sur Q de l'algèbre de Hecke formé de T(n) pour n premier à avoid_N
    list_hecke_mats = []
    list_n = []
    new_n = 1
    nb_trials = 1
    while true:
        # Increase n
        new_n += 1
        while gcd(new_n,avoid_N) > 1:
            new_n += 1
        # Calcul des matrices du nouvel opérateur de Hecke
        list_hecke_mats_n = []
        for S in list_S_by_nt:
            list_hecke_mats_n.append(S.hecke_matrix(new_n))
        # On commence par tester si le nouvel opérateur est générateur à lui tout seul
        if VAP_distinctes(list_hecke_mats_n,common_base_field):
            return [[1,new_n]]
        # Ce n'était pas le cas; on essaie donc des combinaisons linéaires
        if new_n > 2:
            for trials in range(nb_trials):
                list_mats = copy(list_hecke_mats_n)
                coeffs = []
                for matlist in list_hecke_mats:
                    coeff = ZZ.random_element()
                    i = 0
                    for S in list_S_by_nt:
                        list_mats[i] = list_mats[i] + coeff * matlist[i]
                        i = i+1
                    coeffs.append(coeff)
                if VAP_distinctes(list_mats,common_base_field):
                    gene = []
                    i = 0
                    for list in list_hecke_mats:
                        gene.append([coeffs[i],i+2])
                        i = i+1
                    gene.append([1,new_n])
                    return gene
        list_hecke_mats.append(list_hecke_mats_n)
        nb_trials = nb_trials * 2

def Hecke_Gen(space,list):
    T = space.hecke_algebra().zero()
    for x in list:
        T = T + x[0] * space.hecke_operator(x[1])
    return T
    
def Hecke_Gen_Matrix(space,list):
    T = matrix(space.base_ring(),space.dimension(),space.dimension())
    for x in list:
        T = T + x[0] * space.hecke_matrix(x[1])
    return T.transpose()


def Hecke_Eigen_Mat(space,list,field=CC,zeroplus=1e-53): # Returns complex matrix expressing base forms of space in terms of normalised eigenforms
    d=space.dimension()
    T=Hecke_Gen_Matrix(space,list)
    f = T.charpoly()
    if not f.is_squarefree():
        raise ValueError("charpoly is not squarefree !")
    f = f.base_extend(field)
    vap = [z[0] for z in f.roots()]
    P = Matrix(field,d,d)
    for j in range(d):
        v = Ker(T.base_extend(field)-vap[j],1,zeroplus)
        for i in range(d):
            P[i,j] = v[i,0]/v[0,0]
    return P

def Hecke_newf_vap(coeffs,hecke_gen_list): # eigenval of newf // hecke_gen
    v=0
    for x in hecke_gen_list:
        v+= x[0]*coeffs[x[1]]
    return v

def HeckeEigen(S,field): # Let S be a semi-simple Hecke module, returns matrix expressing eigenforms of S on basis of S
    gen = Find_Hecke_Gen([S],S.base_ring())
    return Hecke_Eigen_Mat(S,gen,field)
