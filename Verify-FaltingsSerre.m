/* Faltings-Serre verification of output */
_<x> := PolynomialRing(Rationals());
K<nu> := NumberField(x^2-x-15);
r := nu;
ZK := Integers(K);
ClS, mClS := SUnitGroup(2*ZK);
ds := [mClS(ClS![a1,a2,a3]) : a1,a2,a3 in [0,1]];
Ls := [K] cat [AbsoluteField(ext<K | Polynomial([-d,0,1])>) : d in ds[2..#ds]];
ZLs := [Integers(L) : L in Ls];
[#RayClassGroup(64*ZL,[1..#RealPlaces(NumberField(ZL))]) mod 3 : ZL in ZLs];
L := Ls[5];
ZL := Integers(L);
R, m := RayClassGroup(2*ZL);
_, m3 := quo<R | 3*R>;
M := NumberField(AbelianExtension(m3^-1*m));
Mabs := AbsoluteField(M);
M := RelativeField(K,Mabs);
M;
M0abs := OptimizedRepresentation(AbsoluteField(Subfields(M)[2][1]));
RelativeField(K,M0abs);

ZMabs := Integers(Mabs);
ClS, mClS := SUnitGroup(2*Integers(Mabs));
ds := [mClS(ClS.i) : i in [1..9]];
h := function(c); if c eq 1 then return GF(2)!0; else return GF(2)!1; end if; end function;
hds := Matrix([[h(LegendreSymbol(ZMabs!d,pp)) : d in ds] : pp in PrimesUpTo(200,Mabs) | Norm(\
pp) mod 2 eq 1]);
Rank(hds) eq Ncols(hds);
[Factorization(Norm(pp)) : pp in PrimesUpTo(200,Mabs)];

