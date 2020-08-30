/* Find places of good reduction */
_<nu,x,y,z> := PolynomialRing(Integers(), 4);
r := nu;
F := (-195*r + 859)*x^4 + (-64*r + 282)*x^2*y^2 + (-5*r + 24)*y^4 + (148*r - 
    652)*x^2*y*z + (-16*r + 68)*y^3*z + (-74*r + 326)*x^2*z^2 + (96*r - 
    422)*y^2*z^2 + (-148*r + 652)*y*z^3 + (-47*r + 207)*z^4;
// F := x^4 + (-2*nu + 6)*x^2*y^2 + (10*nu + 34)*x^2*y*z + (-74*nu - 252)*x^2*z^2 + (-5*nu + 24)*y^4 +
//    (-18*nu - 32)*y^3*z + (299*nu + 1018)*y^2*z^2 + (-2760*nu - 9398)*y*z^3 + (9429*nu + 
//    32107)*z^4;
J := ideal<Parent(nu) | nu^2-nu-15, F, Derivative(F,x), Derivative(F,y), Derivative(F,z)>;
J;
Saturation(J, ideal<Parent(x) | [x,y,z]>);
