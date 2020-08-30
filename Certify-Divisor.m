/* Certifying the existence of a divisor and hence verify the endomorphisms */
SetVerbose("EndoCheck", 3);

F := Rationals();
R<x,y,z> := PolynomialRing(F, 3);
f := x^4 - 2*x^2*y^2 - 3*x^2*z^2 + y^4 - 2*y^3*z + y*z^3;
f := x^4 + 212/43*x^3*y + 48/43*x^3*z + 360/43*x^2*y^2 + 252/43*x^2*y*z - 324/43*x^2*z^2 + 224/43*x*y^3 + 612/43*x*y^2*z - 1512/43*x*y*z^2 + 108/43*x*z^3 + 24/43*y^4 + 12*y^3*z - 1512/43*y^2*z^2 + 108/43*y*z^3;
X := PlaneCurve(f);
P0 := X ! [0, 0, 1];
P0 := X ! [2, -2, -1/3];

T := DiagonalMatrix(F, [-1,1,1]);
T := Matrix(F, [
[ 7/3,    2,    2],
[-4/3,   -1,   -2],
[-8/9, -4/3, -1/3]
]);

R<t> := PolynomialRing(Rationals());
F<nu> := NumberField(t^2 - t - 15);
R<x,y,z> := PolynomialRing(F, 3);
f := (-195*nu + 859)*x^4 + (-64*nu + 282)*(-x^2)*y^2 + (148*nu - 652)*(-x^2)*y*z + (-74*nu + 326)*(-x^2)*z^2 + (-5*nu + 24)*y^4 + (-16*nu + 68)*y^3*z + (96*nu - 422)*y^2*z^2 + (-148*nu + 652)*y*z^3 + (- 47*nu + 207)*z^4;
X := PlaneCurve(f);
P0 := X![nu+3,1,0];

fu := 1 - nu;
T := DiagonalMatrix(F, [-1,1,1]);
T := Matrix(F, [
[            1,             0,             0],
[            0,           3/2, 1/2*(-fu - 3)],
[            0, 1/2*(-fu + 4),          -3/2]
]);
print T;

print "Curve:";
print X;
print "Calculating divisor:";
time test, D := DivisorFromMatrixAmbientSplit(X, P0, X, P0, T : LowerBound := 1);
eqs := DefiningEquations(D);
R<y2,y1,x2,x1> := Parent(eqs[1]);
print "Divisor:";
print D;

print "Calculating Cantor representation...";
time test, fs := CantorFromMatrixAmbientSplit(X, P0, X, P0, T : LowerBound := 1);
R<x,y> := Parent(fs[1]);
print "Cantor representation:";
print fs;

Y := X;
fs := [ X`KU ! f : f in fs ];
ceqs := Y`cantor_eqs;

print "";
print "Check 0:";
print [ Evaluate(ceq, fs) : ceq in ceqs ];
