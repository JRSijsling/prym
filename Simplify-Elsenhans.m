/* Trying to decrase size of defining equation by generalizing Elsenhans' minimal model algortithms */

function RandomFraction(B)

 num := Random([-B..B]);
 repeat
  den := Random([-B..B]);
 until den ne 0;
 return num/den;
end function;


function ReducePolynomial(f, v, R0)

 if IsZero(f) then
  return R0 ! 0;
 end if;
 return R0 ! &+[ Evaluate(MonomialCoefficient(f, mon), v) * Monomial(R0, Exponents(mon)) : mon in Monomials(f) ];
end function;


function LiftElement(c0, v, K)

 FF := ResidueClassField(v);
 //c0 := FF ! c0;

 /* Try given generator */
 gen := K.1;
 if Valuation(gen, v) ge 0 then
  gen0 := Evaluate(gen, v);
 else
  gen := K ! 0; gen0 := FF ! 0;
 end if;

 if not sub< FF | gen0 > eq FF then
  print "Finding other generator";
  repeat
   done := false;
   gen := K ! [ RandomFraction(10) : i in [1..Degree(K)] ];
   if Valuation(gen, v) ge 0 then
    gen0 := Evaluate(gen, v);
    if sub< FF | gen0 > eq FF then
     done := true;
    end if;
   end if;
  until done;
 end if;

 w := Matrix([ Eltseq(c0) ]);
 A := Matrix([ Eltseq(gen0^i) : i in [0..(Degree(FF) - 1)] ]);
 v := Eltseq(Solution(A, w));
 c := &+[ (Rationals() ! Integers() ! v[i+1])*gen^i : i in [0..(Degree(FF) - 1)] ];
 return K ! c;
end function;


function LiftPolynomial(f0, v, R)

 K := BaseRing(R);
 return R ! &+[ LiftElement(MonomialCoefficient(f0, mon0), v, K) * Monomial(R, Exponents(mon0)) : mon0 in Monomials(f0) ];
end function;


function TryReducibleReduction(f, pi)

 R := Parent(f);
 K := BaseRing(R); ZZK := Integers(K); v := Place(ideal< ZZK | pi >);
 FF := ResidueClassField(v); R0 := PolynomialRing(FF, Rank(R));

 res_l := []; gcd_l := []; subs_l := [];
 f0 := ReducePolynomial(f,v,R0);
 fac0 := Factorization(f0);
 lf0 := [i[1] : i in fac0 | (TotalDegree(i[1]) eq 1)];

/* linear factors -- test weight vector [0,0,1] */
 for act0 in lf0 do
  j := 1;
  while MonomialCoefficient(act0,R0.j) ne 1 do j := j+1; end while;
  subs := [R.i : i in [1..3]];
  subs[j] := -LiftPolynomial(act0, v, R) + (pi + 1)*R.j;
  res := Evaluate(f,subs);

  test, gcd := IsPrincipal(ideal< ZZK | Coefficients(res) >);
  v0 := Valuation(gcd, v);
  assert v0 ge 1;

  Append(~res_l, res / gcd);
  Append(~gcd_l, gcd);
  Append(~subs_l, subs);
 end for;
 return res_l, subs_l, gcd_l;
end function;


function NormalizePoint(pt0, pi, R)

 K := BaseRing(R); ZZK := Integers(K); v := Place(ideal< ZZK | pi >);
 FF := ResidueClassField(v); R0 := PolynomialRing(FF, Rank(R));

 j := 1; while pt0[j] eq 0 do j := j + 1; end while;
 rep := [pt0[i] / pt0[j] : i in [1..3]];
 subs := [LiftElement(rep[i], v, K) * R.1 : i in [1..3] ]; /* Move rep to [1,0,0] */
 if j eq 1 then
  subs[2] := subs[2] + pi*R.2; subs[3] := subs[3] + pi*R.3;
 end if;
 if j eq 2 then
  subs[1] := subs[1] + pi*R.2; subs[3] := subs[3] + pi*R.3;
 end if;
 if j eq 3 then
  subs[1] := subs[1] + pi*R.2; subs[2] := subs[2] + pi*R.3;
 end if;
 return subs;
end function;


function LocalMinimizationStepPlaneQuartic(f,pi);
 print "Entering minimization step";

 R := Parent(f);
 K := BaseRing(R); ZZK := Integers(K); v := Place(ideal< ZZK | pi >);
 FF := ResidueClassField(v); R0 := PolynomialRing(FF, Rank(R));

 res_l, subs_l, gcd_l := TryReducibleReduction(f,pi);
 ind := [j : j in [1..#gcd_l] | Valuation(gcd_l[j],v) ge 2 ];
 if #ind gt 0 then
  print "Weight vector (0,0,1)";
  return res_l[ind[1]], subs_l[ind[1]];
 end if;

/* Reducible case does not work -- try weight vector (0,1,3): */
 f0 := ReducePolynomial(f,v,R0);
 fac0 := Factorization(f0);
 //print fac0;

 lf0 := [i[1] : i in fac0 | (TotalDegree(i[1]) eq 1)];
/* Search for critical points on lines -- try weight vector [0,1,3] */
 for act0 in lf0 do
/* Parametrize the line */
  T := PolynomialRing(K,2);
  T0 := PolynomialRing(FF,2);
  co := Matrix([[MonomialCoefficient(act0,R0.i)] : i in [1..3]]);
  ker := Kernel(co);
  assert Dimension(ker) eq 2;

  bm := BasisMatrix(ker);
  para := [ &+[LiftElement(bm[j,i], v, K) * T.j : j in [1..2]] : i in [1..3]];
  eqns := [Evaluate(f,para) / pi] cat [Evaluate(Derivative(f,i),para) : i in [1..3]];
  eqns0 := [ReducePolynomial(eqn,v,T0) : eqn in eqns];
  S0 := Scheme(ProjectiveSpace(T0), eqns0);

/* Scheme of critical points */
  cpt0 := Points(S0);
  for i := 1 to #cpt0 do
   pt0 := [cpt0[i][j] : j in [1..2]];
   pt0 := [Evaluate(ReducePolynomial(para[k],v,T0), pt0) : k in [1..3]];
   subs := NormalizePoint(pt0, pi, R);
   ttc := Evaluate(f,subs);

   test, gcd := IsPrincipal(ideal< ZZK | Coefficients(ttc) >);
   v0 := Valuation(gcd, v);
   assert v0 ge 2;

   if v0 ge 3 then
    print "Weight vector (0,1,1)";
    return ttc / gcd, subs;
   end if;
   res_l, subs_l, gcd_l := TryReducibleReduction(ttc / gcd,pi);
   for j := 1 to #res_l do
    if Valuation(gcd_l[j],v) ge 3 then
     print "Weight vector (0,1,2)";
     return res_l[j], [Evaluate(subs[k],subs_l[j]) : k in [1..3]];
    end if;
    res_l2, subs_l2, gcd_l2 := TryReducibleReduction(res_l[j],pi);
    for k := 1 to #res_l2 do
     if Valuation(gcd_l[j] * gcd_l2[k],v) ge 4 then
      subs := [Evaluate(subs[m],subs_l[j]) : m in [1..3]];
      subs := [Evaluate(subs[m],subs_l2[k]) : m in [1..3]];
      print "Weight vector (0,1,3)";
      return res_l2[j], subs;
     end if;
    end for;
   end for;
  end for;
 end for;

/* Complete weight vector [0,1,1] -- critical point is no longer on a rational line
   The double of a smooth conic does not lead to a critical point.  */
 if fac0[1][2] eq 2 and TotalDegree(fac0[1][1]) eq 2 then
/* Singular point of double conic */
  eqns0 := [fac0[1][1]] cat [Derivative(fac0[1][1],i) : i in [1..3]];
 else
  nl0 := &*([fac0[i][1] : i in [1..#fac0] | TotalDegree(fac0[i][1]) gt 1] cat [1]);
/* Singular point, but not on a line -- Points on lines are alreaddy treated. */
  eqns0 := [nl0] cat [Derivative(fac0[1][1],i) : i in [1..3]];
 end if;

/* Scheme of singular points */
 S0 := Scheme(ProjectiveSpace(R0), eqns0);
 cpt0 := Points(S0);
 //print cpt0;

 for i := 1 to #cpt0 do
  subs := NormalizePoint([cpt0[i][k] : k in [1..3]], pi, R);
  //print subs;
  res := Evaluate(f,subs);

  test, gcd := IsPrincipal(ideal< ZZK | Coefficients(res) >);
  v0 := Valuation(gcd, v);
  assert v0 ge 1;

  if v0 ge 3 then
   print "Weight vector (0,1,1)";
   return res / gcd, subs;
  end if;
 end for;

 print "Stopping";
 return f, [Parent(f).i : i in [1..3]];
end function;


function MinimizePlaneQuartic(f,pi)

 R := Parent(f);
 K := BaseRing(R); ZZK := Integers(K); v := Place(ideal< ZZK | pi >);
 FF := ResidueClassField(v); R0 := PolynomialRing(FF, Rank(R));
 res := f;

 test, gcd := IsPrincipal(ideal< ZZK | Coefficients(res) >);
 res := res / gcd;
 print #Sprint(res);
 subs := [R.i : i in [1..3]];

 repeat
  res,subs_n := LocalMinimizationStepPlaneQuartic(res,pi);
  test, gcd := IsPrincipal(ideal< ZZK | Coefficients(res) >);
  res := res / gcd;
  print #Sprint(res);
  subs := [Evaluate(a,subs_n) : a in subs];
 until subs_n eq [Parent(subs_n[1]).i : i in [1..3]];

 return R ! res, subs;
end function;


R<x> := PolynomialRing(Rationals());
f := x^2 - x - 15;
K<r> := NumberField(f);
ZZK := Integers(K);
U, phiU := UnitGroup(ZZK);

R<x,y,z> := PolynomialRing(K, 3);
F := x^4 + (-2*r + 6)*x^2*y^2 + (4*r + 14)*x^2*y*z + (-12*r - 42)*x^2*z^2 + (-5*r + 24)*y^4 + (22*r - 142)*y^3*z + (50*r + 163)*y^2*z^2 + (-176*r - 658)*y*z^3 + (224*r + 997)*z^4;

u := phiU(U.1 + U.2);
F := u^2*x^4 + (-2*r + 6)*u*x^2*y^2 + (4*r + 14)*u*x^2*y*z + (-12*r - 42)*u*x^2*z^2 + (-5*r + 24)*y^4 + (22*r - 142)*y^3*z + (50*r + 163)*y^2*z^2 + (-176*r - 658)*y*z^3 + (224*r + 997)*z^4;
F := 4*u^2*x^4 + (-2*r + 6)*2*u*x^2*y^2 + (4*r + 14)*2*u*x^2*y*z + (-12*r - 42)*2*u*x^2*z^2 + (-5*r + 24)*y^4 + (22*r - 142)*y^3*z + (50*r + 163)*y^2*z^2 + (-176*r - 658)*y*z^3 + (224*r + 997)*z^4;
I, W := DixmierOhnoInvariants(F);

Fac := Factorization(ideal< ZZK | I[#I] >);
test, gen := IsPrincipal(Fac[1][1]);
F := MinimizePlaneQuartic(F, gen);
I, W := DixmierOhnoInvariants(F);
Fac := Factorization(ideal< ZZK | I[#I] >);

F := MinimizePlaneQuartic(F, 2);
I, W := DixmierOhnoInvariants(F);
Fac := Factorization(ideal< ZZK | I[#I] >);
print Fac;

print F;
test, d := IsPrincipal(ideal< ZZK | Coefficients(F) >);
F /:= d;

