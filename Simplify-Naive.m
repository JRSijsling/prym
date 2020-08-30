/* Trying to decrase size of defining equation naively */
SetSeed(1);

//prec := 1000;
prec := 700;
SetDefaultRealFieldPrecision(prec);

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

/*
CC<I> := K`CC;
S := RiemannSurface(h(F), sigma : Precision := Precision(CC));
Q := ChangeRing(S`BigPeriodMatrix, CC);
//PrintFileMagma("Matrices.m", Q);

FCC := EmbedPolynomialExtra(F);
CC<I> := K`CC;
*/

/*
load "Matrices.m";
P := ChangeRing(P, CC);
Q := ChangeRing(Q, CC);
*/

/*
SetVerbose("EndoFind", 3);
isos := SymplecticIsomorphismsCC(Q, Q);
print isos;

_, Aseq := NumberFieldExtra(Eltseq(isos[2][1]), K);
print Aseq;

// Blech
G, _, Gphi := AutomorphismGroup(K);
A := Matrix(K, 3,3, [ Gphi(G.1)(a) : a in Aseq ] );

D, T := JordanForm(A);
G := F^(T^(-1));
*/

G := (-65445264537800236469976444735926939637815778378926359092933451464954692978146462727760*r - 29955400541188572029626701763847982934806324465134178052023791247997245556470116262672)*x^2 + (-1180146280931482937269385929317311915764124090599471414822423750980853814924986050000*r + 1820879899157188105247286217011530937403461058145631487080618701200576074121336875000)*x*y^2 + (42418111333452084075436945686953758338059978425990053654469598120954728439213493750000*r + 130123821967149510828307385748124225173798995456806076428218977566213742117149031250000)*x*y*z + (-1993962875342883789386781856933986932527283525187281265943717862191427269322745703125000*r - 6783108638526739992762830300281308985367415726160455629476304121628001745782683593750000)*x*z^2 + (-11809291691427425895597494546797108808848704775233407453991103833916763947939453125*r + 44416110298005943743283340277443479418116011132144713119700488278698232847558593750)*y^4 + (232849372471017419416378273196273257492063643155515780277894880341839819785156250000*r + 896566282831639277235011706212700067133695539760953113179617273281576939721679687500)*y^3*z + (-22400785391982212558206824843889583446026166208132783727678383624356783770446777343750*r - 77052172522031745417354703918249036516001358478231269971492716506523902415466308593750)*y^2*z^2 + (789953890477247259943957753335133817150795692454149044023730317348331028518676757812500*r + 2692721658569653499723163511856467259047166473700413735567254943243497244262695312500000)*y*z^3 + (-7158853249435648167052411859935368891722653544163405500514790505770116810798645019531250*r - 24374951337057966740816558946386331424785346359424668192439929575711716902256011962890625)*z^4;

Fac := Factorization(ideal< ZZK | MonomialCoefficient(G, x^2) >);
Fac := Reverse(Fac);
test, gen := IsPrincipal(&*[ Fac[i][1] : i in [1..6] ]);
G := Evaluate(G, [ 1063^2*x / gen, y, z ]);
test, d := IsPrincipal(ideal< ZZK | Coefficients(G) >);
G /:= d;

test, gen := IsPrincipal(ideal< ZZK | [ 5, r ] >);
G := Evaluate(G, [ x, gen^3*y, z ]);
G := Evaluate(G, [ gen^12*x, y, z ]);
test, d := IsPrincipal(ideal< ZZK | Coefficients(G) >);
G /:= d;

gen := 2;
G := Evaluate(G, [ gen^(-2)*x, y, z ]);
test, d := IsPrincipal(ideal< ZZK | Coefficients(G) >);
G /:= d;

test, gen := IsPrincipal(ideal< ZZK | [ 3, r ] >);
G := Evaluate(G, [ gen^(-12)*x, y, gen^12*z ]);
test, d := IsPrincipal(ideal< ZZK | Coefficients(G) >);
G /:= d;

test, gen := IsPrincipal(ideal< ZZK | [ 3, 5*r + 1 ] >);
G := Evaluate(G, [ gen^(-21)*x, y, gen^21*z ]);
test, d := IsPrincipal(ideal< ZZK | Coefficients(G) >);
G /:= d;

G /:= MonomialCoefficient(G, x^2);
test, d := IsPrincipal(ideal< ZZK | Coefficients(G) >);
G /:= d;

test, gen := IsPrincipal(ideal< ZZK | [ 5, 4*r + 1 ] >);
G := Evaluate(G, [ gen^11*x, gen^3*y, z ]);
test, d := IsPrincipal(ideal< ZZK | Coefficients(G) >);
G /:= d;

test, gen := IsPrincipal(ideal< ZZK | [ 47, 12 + r ] >);
G := Evaluate(G, [ gen^3*x, gen*y, z ]);
test, d := IsPrincipal(ideal< ZZK | Coefficients(G) >);
G /:= d;

G := (452307960*r + 1540165069)*x^2 + (-14730336518385462242*r - 50158634757806812984)*x*y^2 + (747661483315157764468*r + 2797238357213231310296)*x*y*z + (198608956544567351576412098230*r - 874897416261413161610312252752)*x*z^2 + (89122297587520274092967110785*r + 402226432605545967984627611047)*y^4 + (-207862540760215810112186254961508006508*r + 915660327845926417512727220157458898884)*y^3*z + (20839159933760437997281998731149856623306601029542*r - 91799101026289912313510476667809593826711397843838)*y^2*z^2 + (-1021468147879718335000162531352567743523620796694515003473900*r + 4499694709403129492382709444842357688433107047259479799411796)*y*z^3 + (-8899911420057458297429757435907534450845121391846075027327752395706239*r + 39205220852079577121308484470146580541045867832945120613168590393242095)*z^4;

/*
for mon in Monomials(G) do
    print mon;
    print Factorization(ideal< ZZK | MonomialCoefficient(G, mon) >);
end for;
*/

NG := #Sprint(G);
print NG;

F := Evaluate(G, [ x^2, y, z ]);
print F;

/*
// For Stoll--Cremona reduction after Bouyer--Streng (which fails):
T<t> := PolynomialRing(K);
h := hom< S -> T | [ 1,t,1] >;
a := h(&+[ MonomialCoefficient(G, mon)*mon : mon in [ x^2 ] ]);
b := h(&+[ MonomialCoefficient(G, mon)*mon : mon in [ x*y^2, x*y*z, x*z^2 ] ]);
c := h(&+[ MonomialCoefficient(G, mon)*mon : mon in [ y^4, y^3*z, y^2*z^2, y*z^3, z^4 ] ]);
f := (b^2 - 4*a*c)/(4*a^2);
print f;
*/


