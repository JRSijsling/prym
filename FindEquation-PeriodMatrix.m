/* Transform reconstructed curve to a curve yielding the exact period matrix
 * supplied by Pari. This simplifies the curve, while on the other hand making
 * the base field larger */
SetSeed(1);

//prec := 1000;
prec := 900;
SetDefaultRealFieldPrecision(prec);

/* Base field and polynomial */
R<x> := PolynomialRing(Rationals());
f := x^2 - x - 15;
K<r> := BaseNumberFieldExtra(f, prec);
ZZK := Integers(K);
S<x,y,z> := PolynomialRing(K, 3);
R := PolynomialRing(K, 2);
h := hom< S -> R | [ R.1, R.2, 1 ] >;

F := (27307965693753899931844585236559742132400430059592529296875*r - 120295011428448504225733593228325497613945662503643310546875)*x^4 + (15417121920249437332293413294278870057376067736083984375000*r - 67915943086870853749604617208563573060155694168945312500000)*x^3*y + (78155549129617409866649767950780305755219439092604744648437500*r - 344284388997432521172802616153482553574957360522973648429687500)*x^3*z + (100105738249792236317588696268951821445075880279541015625000*r - 441069816545920639692977413563534497155933003143310546875000)*x^2*y^2 + (735112762927886451396520398634505757328277201868280775390625000*r - 3238172376160861232188446236989422030972340656324589470703125000)*x^2*y*z + (-102406986297877921599885210309422036317052427062492050152495468750*r + 451087137720113362748411271932725154057029475094494337658177593750)*x^2*z^2 + (143317362045937921558244335092859076882143398284912109375000*r - 633418727084374262902242501534473697373541107177734375000000)*x*y^3 + (1291675339755613084789726398017218770449977816358732666015625000*r - 5686673889245649903305198470167012633953489889797505615234375000)*x*y^2*z + (-2547878921616522466688765929910856861590052498746221081267359375000*r + 11221851574248983592677030057213718115814523478782917672749828125000)*x*y*z^2 + (960331960293989693226695177815015337351612194107715679474135968330500*r - 4230047329592862592568258939678860633785747664420122687788717511427500)*x*z^3 + (578913976305813570655728246683201809823951005935668945312500*r - 2587637032498458172893071975310078671316399812698364257812500)*y^4 + (2022372202657451125365319358163520552262148412540435791015625000*r - 8807481250464730037170334161971249120345760289513977050781250000)*y^3*z + (4403060098956258692422966349774886816876884373995333903666015625000*r - 19501261277960566038981496469990693924435530155081352203808593750000)*y^2*z^2 + (340033390583783926373864349706435529674679220409111954567230953125000*r - 1448163723213624150552410834442029712732482479704491483960470168250000)*y*z^3 + (2885421569863669681858331608379178163318264181842855111741619345682549763*r - 12719889818818547544412693016029624997890364334618820880961485211306780239)*z^4;
F /:= MonomialCoefficient(F, x^4);
sigma := InfinitePlaces(K)[1];
print Evaluate(K.1, sigma);

/* Calculate period matrix (now stored in Matrices.m)
CC<I> := K`CC;
S := RiemannSurface(h(F), sigma : Precision := Precision(CC));
Q := ChangeRing(S`BigPeriodMatrix, CC);
//PrintFileMagma("Matrices.m", Q);
*/

FCC := EmbedPolynomialExtra(F);
CC<I> := K`CC;

load "Matrices.m";
P := ChangeRing(P, CC);
Q := ChangeRing(Q, CC);

/* Find isomorphisms between the given period matrix P and the new one Q */
SetVerbose("EndoFind", 3);
isos := SymplecticIsomorphismsCC(P, Q);
print isos;

G, _, Gphi := AutomorphismGroup(K);
g := &+[ Gphi(G.1)(MonomialCoefficient(F, mon))*mon : mon in Monomials(F) ];

/* Get period matrix corresponding to P */
T := Matrix(CC, 3,3, Eltseq(isos[1][1]));
gCC := EmbedPolynomialExtra(g);
fCC := TransformForm(gCC, T);

CC := BaseRing(Parent(fCC));
coeffs := [ c : c in Coefficients(fCC) | Abs(c) gt CC`epscomp ];
min, ind := Minimum([ Abs(c) : c in coeffs ]);
fCC /:= coeffs[1];

monsCC := Monomials(fCC);
coeffsCC := [ MonomialCoefficient(fCC, monCC) : monCC in monsCC ];
exps := [ Exponents(monCC) : monCC in monsCC ];

/* Find equation corresponding to P */
L, coeffs, hKL := NumberFieldExtra(coeffsCC, K : UpperBound := 24);
