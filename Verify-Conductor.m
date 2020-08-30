/* Conductor of various curves occuring in the argument */
_<t> := PolynomialRing(Rationals());
F<nu> := NumberField(t^2 - t - 15);
R<x,y,z> := PolynomialRing(F, 3);
f := (-195*nu + 859)*x^2 + (-64*nu + 282)*(-x)*y^2 + (148*nu - 652)*(-x)*y*z + (-74*nu + 326)*(-x)*z^2 + (-5*nu + 24)*y^4 + (-16*nu + 68)*y^3*z + (96*nu - 422)*y^2*z^2 + (-148*nu + 652)*y*z^3 + (- 47*nu + 207)*z^4;
f := (10*nu + 33)*x^4 + (-12*nu - 44)*x^3*z + (-14*nu - 48)*x^2*y^2 + (2*nu + 26)*x^2*z^2 + (8*nu + 28)*x*y^2*z + (4*nu - 16)*x*z^3 + (5*nu + 17)*y^4 + (-4*nu - 14)*y^2*z^2 + (nu - 6)*z^4;
f := (5*nu + 17)*x^4 + (-14*nu - 48)*x^2*y^2 + (8*nu + 28)*x^2*y*z + (-4*nu - 14)*x^2*z^2 + (10*nu + 33)*y^4 + (-12*nu - 44)*y^3*z + (2*nu + 26)*y^2*z^2 + (4*nu - 16)*y*z^3 + (nu - 6)*z^4;

R<t> := PolynomialRing(F);
a := (-195*nu + 859);
b := (-64*nu + 282)*(-1)*t^2 + (148*nu - 652)*(-1)*t + (-74*nu + 326)*(-1);
c := (-5*nu + 24)*t^4 + (-16*nu + 68)*t^3 + (96*nu - 422)*t^2 + (-148*nu + 652)*t + (- 47*nu + 207);

a := (5*nu + 17);
b := (-14*nu - 48)*t^2 + (8*nu + 28)*t + (-4*nu - 14);
c := (10*nu + 33)*t^4 + (-12*nu - 44)*t^3 + (2*nu + 26)*t^2 + (4*nu - 16)*t + (nu - 6);

p := b^2 - 4*a*c;
X := HyperellipticCurve(p);
E := EllipticCurve(X);

print Conductor(E);
print Conductor(EllipticCurve(t^3 + nu*t^2 + (1 - nu)*t + 1));
p := (10*nu + 34)*(t^3 + nu*t^2 + (1 - nu)*t + 1);
X := HyperellipticCurve(p);
E := EllipticCurve(X);
print Conductor(E);

