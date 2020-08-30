/* Final simplification attempt */
Attach("binary_quartics.m");
SetSeed(1);

/* Define base number field K of discriminant 61 */
R<x> := PolynomialRing(Rationals());
f := x^2 - x - 15;
K<r> := NumberField(f);
ZZK := Integers(K);
U, phiU := UnitGroup(ZZK);

/* Base elliptic curve over K */
R<t> := PolynomialRing(K);
g0 := t^3 + (-13*r - 44)*t + 4*(13*r + 44);
j0 := BinaryQuarticInvariants(g0);
print j0;

/*
D := [-2..2];
while true do
    g := t^3 + &+[ (K ! [ Random(D) : i in [1..2] ])*t^i : i in [0..2] ];
    if Discriminant(g) ne 0 then
        j := BinaryQuarticInvariants(g);
        if j eq j0 then
            print g;
            Disc := Discriminant(g);
            print Factorization(ideal< ZZK | Disc >);
        end if;
    end if;
end while;
*/

f := t^3 + (-13*r - 44)*t + 4*(13*r + 44);
j := BinaryQuarticInvariants(f);
print j;
print Factorization(Conductor(EllipticCurve(f)));

g := t^3 + r*t^2 + (-r + 1)*t + 1;
j := BinaryQuarticInvariants(g);
print j;
print Factorization(Conductor(EllipticCurve(g)));

/* Reduce size of defining equation of elliptic curve */
print IsGL2Equivalent(f, g, 4);
T := Matrix(K, 2,2, [ 1, 1/3*r, 0, 1/9*(2*r - 6) ]);

num := T[1,1]*t + T[1,2];
den := T[2,1]*t + T[2,2];
subst := num/den;
print subst;

g0 := Numerator(Evaluate(f, subst)*den^4);
lg0 := LeadingCoefficient(g0);
Fac := Factorization(ideal< ZZK | lg0/2 >);
test, gen := IsPrincipal(Fac[1][1]);
Fac := Factorization(ideal< ZZK | gen^2*lg0/2 >);
test, rt := IsSquare(gen^2*lg0/2);

gen2 := 2;
gen31 := r - 4;
gen32 := r + 3;
subst := 1/2*(r + 2)*t + 1/2*(r + 5);
g0 := Evaluate(f, subst);
test, c := IsPrincipal(ideal< ZZK | Coefficients(g0) >);
c := LeadingCoefficient(g0);
print c;
Fac := Factorization(ideal< ZZK | c >);
print Fac;

u := gen2^3*gen31^(-6)*c;
print u @@ phiU;

R<x,y,w> := PolynomialRing(K, 3);
f1 := y^2 - (x^3 + (-13*r - 44)*x + 4*(13*r + 44));
f2 := w^2 - (6*y + (-r - 6)*x^2 + (15*r + 51)*x - 45*r - 153);
xnew := 1/2*(r + 2)*x + 1/2*(r + 5);
ynew := gen2^(-2)*gen31^3*phiU(-1*U.2)*y;
wnew := w/(phiU(0*U.2)*2*gen31^(-2)*gen32^(-1));

f1new := Evaluate(f1, [ xnew, ynew, wnew ]);
test, gcd := IsPrincipal(ideal< ZZK | Coefficients(f1new) >);
f1new /:= gcd;
f1new /:= MonomialCoefficient(f1new, y^2);
print f1new;

f2new := Evaluate(f2, [ xnew, ynew, wnew ]);
test, gcd := IsPrincipal(ideal< ZZK | Coefficients(f2new) >);
f2new /:= gcd;
//f2new /:= MonomialCoefficient(f2new, w^2);
print f2new;

print Factorization(ideal< ZZK | MonomialCoefficient(f2new, y) >);

f1new := (-10*r - 34)*x^3 + (-44*r - 150)*x^2 + (34*r + 116)*x + y^2 - 10*r - 34;
f2new := (54*r + 189)*x^2 + (-30*r - 114)*x + (6*r - 42)*y + w^2 + 15*r + 57;
f2new := (2029*r + 6909)*x^2 + (-1170*r - 3984)*x + (-74*r - 252)*y + w^2 + 585*r + 1992;
f2new := (7*r + 24)*x^2 + (-4*r - 14)*x - 2*y + (-2*r + 9)*w^2 + 2*r + 7;

R<x,w> := PolynomialRing(K, 2);
y := (-1/(6*r - 42))*((54*r + 189)*x^2 + (-30*r - 114)*x + w^2 + 15*r + 57);
y := (-1/(-74*r - 252))*((2029*r + 6909)*x^2 + (-1170*r - 3984)*x + w^2 + 585*r + 1992);
y := (1/2)*((7*r + 24)*x^2 + (-4*r - 14)*x + w^2 + 2*r + 7);

F := (-10*r - 34)*x^3 + (-44*r - 150)*x^2 + (34*r + 116)*x + y^2 - 10*r - 34;

X := Curve(AffineSpace(R), F);
X := ProjectiveClosure(X);
F := DefiningPolynomial(X);
test, gcd := IsPrincipal(ideal< ZZK | Coefficients(F) >);
F /:= gcd;
_<x,y,z> := Parent(F);
print F;

exit;

I, W := DixmierOhnoInvariants(F);
print Factorization(ideal< ZZK | I[#I] >);
print WPSNormalize(W, I);

F := (10*r + 33)*x^4 + (-12*r - 44)*x^3*z + (-14*r - 48)*x^2*y^2 + (2*r + 26)*x^2*z^2 + (8*r + 28)*x*y^2*z + (4*r - 16)*x*z^3 + (5*r + 17)*y^4 + (-4*r - 14)*y^2*z^2 + (r - 6)*z^4;
test, gcd := IsPrincipal(ideal< ZZK | Coefficients(F) >);
print gcd;

I, W := DixmierOhnoInvariants(F);
print Factorization(ideal< ZZK | I[#I] >);
print WPSNormalize(W, I);
