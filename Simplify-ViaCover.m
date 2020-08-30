/* Simplification as in the article as a degree-2 cover of a simplified
 * elliptic curve */
prec := 700;
SetDefaultRealFieldPrecision(prec);

R<x> := PolynomialRing(Rationals());
f := x^2 - x - 15;
K<r> := BaseNumberFieldExtra(f, prec);
ZZK := Integers(K);
S<x,y,z> := PolynomialRing(K, 3);

G := (452307960*r + 1540165069)*x^2 + (-14730336518385462242*r - 50158634757806812984)*x*y^2 + (747661483315157764468*r + 2797238357213231310296)*x*y*z + (198608956544567351576412098230*r - 874897416261413161610312252752)*x*z^2 + (89122297587520274092967110785*r + 402226432605545967984627611047)*y^4 + (-207862540760215810112186254961508006508*r + 915660327845926417512727220157458898884)*y^3*z + (20839159933760437997281998731149856623306601029542*r - 91799101026289912313510476667809593826711397843838)*y^2*z^2 + (-1021468147879718335000162531352567743523620796694515003473900*r + 4499694709403129492382709444842357688433107047259479799411796)*y*z^3 + (-8899911420057458297429757435907534450845121391846075027327752395706239*r + 39205220852079577121308484470146580541045867832945120613168590393242095)*z^4;
NG := #Sprint(G);
print NG;

/*
for mon in Monomials(G) do
    print mon;
    print Factorization(ideal< ZZK | MonomialCoefficient(G, mon) >);
end for;
*/

/*
F := Evaluate(G, [ x^2, y, z ]);
print F;
*/

T<t> := PolynomialRing(K);
h := hom< S -> T | [ 1,t,1] >;
a := h(&+[ MonomialCoefficient(G, mon)*mon : mon in [ x^2 ] ]);
b := h(&+[ MonomialCoefficient(G, mon)*mon : mon in [ x*y^2, x*y*z, x*z^2 ] ]);
c := h(&+[ MonomialCoefficient(G, mon)*mon : mon in [ y^4, y^3*z, y^2*z^2, y*z^3, z^4 ] ]);
f := b^2 - 4*a*c;
print f;

/*
for i := 0 to 4 do
    print i;
    print Factorization(ideal< ZZK | Coefficient(f, i) >);
end for;
*/

j := BinaryQuarticInvariants(f);
print j;
g := t^3 + (-13*r - 44)*t + 4*(13*r + 44);
j := BinaryQuarticInvariants(g);
print j;

print IsGL2Equivalent(f, g, 4);
T := Matrix(K, 2,2, [ 1, 1/518325113739780541910511913*(-189941303676068175727952941*r - 974554372499935413043934834), 1/518325113739780541910511913*(-642276438367034082451628*r - 2187031456608220067294393), 1/518325113739780541910511913*(7982940752511962968405349*r + 27182909815319262180034390) ]);

g0 := f^T;
test, rt := IsSquare(LeadingCoefficient(g0));
print test;

num := T[1,1]*t + T[1,2];
den := T[2,1]*t + T[2,2];
subst := num/den;
g0 := Numerator(Evaluate(f, subst)*den^4);
test, rt := IsSquare(LeadingCoefficient(g0));
print test;
c := rt;
assert g eq (1/c)^2*Numerator(Evaluate(f, subst)*den^4);


g := t^3 + (-13*r - 44)*t + 4*(13*r + 44);
E := EllipticCurve(g);
K<x,y> := FunctionField(E);
num := x + 1/518325113739780541910511913*(-189941303676068175727952941*r - 974554372499935413043934834);
den := 1/518325113739780541910511913*(-642276438367034082451628*r - 2187031456608220067294393)*x + 1/518325113739780541910511913*(7982940752511962968405349*r + 27182909815319262180034390);
xnew := num/den;
ynew := c*y/den^2;
bnew := (-14730336518385462242*r - 50158634757806812984)*xnew^2 + (747661483315157764468*r + 2797238357213231310296)*xnew + 198608956544567351576412098230*r - 874897416261413161610312252752;
a := 452307960*r + 1540165069;
f := (ynew - bnew)/(2*a);

D := Support(Divisor(f))[1];
test, g := IsPrincipal(D - 4*Divisor(E ! 0));
print E;
print g;

R<x> := PolynomialRing(Rationals());
f := x^2 - x - 15;
K<r> := NumberField(f);

R<x,y,w> := PolynomialRing(K, 3);
f1 := y^2 - (x^3 + (-13*r - 44)*x + 4*(13*r + 44));
f2 := w^2 - (6*y + (-r - 6)*x^2 + (15*r + 51)*x - 45*r - 153);
X := Curve(AffineSpace(R), [ f1, f2 ]);

phi := CanonicalMap(X);
F := DefiningPolynomial(Image(phi));
I, W := DixmierOhnoInvariants(F);
print WPSNormalize(W, I);


R<x,w> := PolynomialRing(K, 2);
y := (1/6)*(w^2 - ((-r - 6)*x^2 + (15*r + 51)*x - 45*r - 153));
f := y^2 - (x^3 + (-13*r - 44)*x + 4*(13*r + 44));
X := Curve(AffineSpace(R), f);
X := ProjectiveClosure(X);
_<x,y,z> := CoordinateRing(Ambient(X));
F := DefiningPolynomial(X);
F := Evaluate(F, [y,x,z]);
F /:= MonomialCoefficient(F, x^4);
print F;

