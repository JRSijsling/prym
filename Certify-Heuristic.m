/* A heuristic check of the endomorphism algebra */
SetVerbose("EndoFind", 2);

prec := 300;

/* Intermediate calculation for diagonalizing involution on a plane quartic
F := RationalsExtra(prec);
S<x,y,z> := PolynomialRing(F, 3);
R := PolynomialRing(F, 2);
h := hom< R -> S | [y,z] >;

D := [-3..3];
mons2 := [ h(mon) : mon in MonomialsOfWeightedDegree(R, 2) ];
mons4 := [ h(mon) : mon in MonomialsOfWeightedDegree(R, 4) ];

repeat
    f := x^4 + x^2*&+[ Random(D)*mon : mon in mons2 ] + &+[ Random(D)*mon : mon in mons4 ];
    I := DixmierOhnoInvariants(f);
until I[#I] ne 0 and Evaluate(f, [0,0,1]) eq 0;
print f;

repeat
    U := Matrix(F, [[ Random(D) : i in [1..3] ] : j in [1..3]]);
until Determinant(U) ne 0;
f := f^U;
U0 := DiagonalMatrix(F, [-1,1,1]);
print U^(-1)*U0*U;

X := PlaneCurve(f);
*/

K<nu> := BaseNumberFieldExtra(Polynomial([-15,-1,1]), 200);
_<x,y> := PolynomialRing(K,2);                            
f := (-195*nu + 859)*x^4 + (-64*nu + 282)*(-x^2)*y^2 + (148*nu - 652)*(-x^2)*y + (-74*nu + 326)*(-x^2) + (-5*nu + 24)*y^4 + (-16*nu + 68)*y^3 + (96*nu - 422)*y^2 + (-148*nu + 652)*y - 47*nu + 207;
X := PlaneCurve(f);


print "";
print "Curve:";
print X;

print "";
print "Check endomorphisms:";
Xriem := RiemannSurface(f, RealPlaces(K)[1] : Precision := 200);
P := BigPeriodMatrix(Xriem);
EndoRep := GeometricEndomorphismRepresentation(P, K);

