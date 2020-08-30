/* Further decrease of size by twisting curve: also determines the right twist */
S := CuspForms(KroneckerCharacter(61));
f := Newforms(S)[1][1];
Mf := BaseField(f);
_<t> := PolynomialRing(ComplexField());
R<x> := PolynomialRing(Rationals());
K<nu> := NumberField(x^2-x-15);
r := nu;
ZK := Integers(K);
_<x,y,z> := PolynomialRing(K,3);
// F := (-5*r + 24)*x^4 + (-64*r + 282)*x^2*y^2 + (-195*r + 859)*y^4 + (-16*r +
//     68)*x^3*z + (148*r - 652)*x*y^2*z + (96*r - 422)*x^2*z^2 + (-74*r +
//     326)*y^2*z^2 + (-148*r + 652)*x*z^3 + (-47*r + 207)*z^4;
F := (-195*r + 859)*x^4 + (-64*r + 282)*x^2*y^2 + (-5*r + 24)*y^4 + (148*r - 
    652)*x^2*y*z + (-16*r + 68)*y^3*z + (-74*r + 326)*x^2*z^2 + (96*r - 
    422)*y^2*z^2 + (-148*r + 652)*y*z^3 + (-47*r + 207)*z^4;
// F := x^4 + (-2*nu + 6)*x^2*y^2 + (4*nu + 14)*x^2*y*z + (-12*nu - 42)*x^2*z^2 + (-5*nu + 24)*y^4 + 
//     (22*nu - 142)*y^3*z + (50*nu + 163)*y^2*z^2 + (-176*nu - 658)*y*z^3 + (224*nu + 997)*z^4;
X := Curve(ProjectiveSpace(Parent(x)), F);
tpps := [];
badpps := [ 1*ZK ];
for pp in PrimesUpTo(200,K) do
  p := Norm(pp);
  if not IsPrime(p) then continue; end if;
  Xpp := Curve(Reduction(X, pp));
  if IsIrreducible(Xpp) and IsNonsingular(Xpp) then
    print pp;
    Lpp := Numerator(ZetaFunction(Xpp));
    print Factorization(Lpp);
    LpCC := &*[Polynomial([1,-r[1],KroneckerCharacter(61)(p)*p]) : r in
Roots(CharacteristicPolynomial(Coefficient(f,p)),ComplexField())];
    Lp := Polynomial(Integers(),[Round(c) : c in Coefficients(LpCC)]);
    print Factorization(Lp);
    if not IsDivisibleBy(Lpp,Lp) then
      if not IsDivisibleBy(Evaluate(Lpp,-Parent(Lpp).1),Lp) then
        print "     =====> What?!?  Is this a prime of bad reduction?";
        Append(~badpps, pp);
      else
        Append(~tpps, pp);
      end if;
    end if;
  else
    Append(~badpps, pp);
  end if;
end for;
ClS, mClS := SUnitGroup(2*&*badpps : GRH := true);
ds := [mClS(ClS![a1,a2,a3]) : a1,a2,a3 in [0,1]];
[< {LegendreSymbol(ZK!d,pp) : pp in tpps}, K!d > : d in ds];
Norm(-5*nu + 22);

delta := -5*nu+22;
F := (-195*r + 859)/delta^2*x^4 + (-64*r + 282)/delta*x^2*y^2 + (-5*r + 24)*y^4 + (148*r - 
    652)/delta*x^2*y*z + (-16*r + 68)*y^3*z + (-74*r + 326)/delta*x^2*z^2 + (96*r - 
    422)*y^2*z^2 + (-148*r + 652)*y*z^3 + (-47*r + 207)*z^4;
F *:= (5*nu+17);

_<nu,x,y,z> := PolynomialRing(Integers(),4);
F := (5*nu + 17)*x^4 + (-14*nu - 48)*x^2*y^2 + (8*nu + 28)*x^2*y*z + (-4*nu - 
    14)*x^2*z^2 + (10*nu + 33)*y^4 + (-12*nu - 44)*y^3*z + (2*nu + 26)*y^2*z^2 +
    (4*nu - 16)*y*z^3 + (nu - 6)*z^4;
J := ideal<Parent(nu) | nu^2-nu-15, F, Derivative(F,2), Derivative(F,3),
Derivative(F,4)>;
Saturation(J, ideal<Parent(nu) | [x,y,z]>);

