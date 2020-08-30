/* Another quick way to certify the existence of a divisor and hence verify the endomorphisms */
SetVerbose("EndoFind", 3);
SetVerbose("EndoCheck", 3);
K<nu> := BaseNumberFieldExtra(Polynomial([-15,-1,1]), 200);
_<x,y> := PolynomialRing(K,2);                            
F := (-195*nu + 859)*x^4 + (-64*nu + 282)*(-x^2)*y^2 + (148*nu - 652)*(-x^2)*y + (-74*nu + 326)*(-x^2) + (-5*nu + 24)*y^4 + (-16*nu + 68)*y^3 + (96*nu - 422)*y^2 + (-148*nu + 652)*y - 47*nu + 207;
X := PlaneCurve(F);
Pi := PeriodMatrix(X);
GeoEndoRep := GeometricEndomorphismRepresentationCC(Pi);
T := GeoEndoRep[3][1];
_, seq := NumberFieldExtra(Eltseq(T), K);
TK := Matrix(K, 3, 3, seq);
print TK;
P0 := X![nu+3,1,0];
DivisorFromMatrixAmbientSplit(X, P0, X, P0, TK);

