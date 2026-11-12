(* Final analytic finite parts from the two source notebooks, 2026-11-12.
   Positive soft branch and original single-topology normalization are retained.
   Gauge: nu = Sqrt[(mA/H)^2-1/4] > 0; top: nu = mt/H > 0.
   No regulator is used. See README.md for provenance and validation. *)
BeginPackage["SignalPlots`Coefficients`"];
GaugeCoefficients::usage = "GaugeCoefficients[nu,p:40] gives complex CBub, CTri, A1, A2, A3, B, D as an Association.";
TopCoefficients::usage = "TopCoefficients[nu,p:50] gives complex A1, A2, A3, B, D as an Association.";
GaugeCoefficientExpressions::usage = "GaugeCoefficientExpressions[nu] gives the exact analytic expressions, including parameter derivatives.";
TopCoefficientExpressions::usage = "TopCoefficientExpressions[nu] gives the exact analytic expressions.";
AngularCoefficient::usage = "AngularCoefficient[coefficients,x,y,phi3] sums the five complex angular terms; x,y are cosines in [-1,1].";
Begin["`Private`"];

angularKeys = {"A1", "A2", "A3", "B", "D"};
gaugeKeys = Join[{"CBub", "CTri"}, angularKeys];
f0[a_, b_] := Gamma[a+b-3/2] Gamma[3/2-a] Gamma[3/2-b]/
  (Gamma[3-a-b] Gamma[a] Gamma[b]);

(* F2/F3: lower-parameter derivatives of regularized 3F2(1).
   F4: the genuinely mixed lower-parameter derivative. Native Wolfram
   analytic continuation at z=1 is preserved, including in numerical calls. *)
gaugeLeadingMasters[v_] := {
  HypergeometricPFQ[{1,1/2+I v,2 I v},{1+2 I v,1+2 I v},1],
  Derivative[{0,0,0},{1,0},0][HypergeometricPFQRegularized][
    {1+2 I v,1/2-I v,1},{1,2},1],
  Derivative[{0,0,0},{1,0},0][HypergeometricPFQRegularized][
    {2+2 I v,1/2-I v,1},{1,3},1]};
gaugeMasters[v_] := Append[gaugeLeadingMasters[v],
  Derivative[{0,0,0},{1,1},0][HypergeometricPFQRegularized][
    {2 I v,1/2-I v,1},{1,1},1]];

(* Compact exact finite P_L/P_T from the final explicit-coefficients section.
   The longitudinal finite contact term is included in delta. *)
gaugeHard[v_, masters_List] := Module[{u=I v,w=Exp[Pi v],g=Gamma[2 I v],
    h=Gamma[1/2+I v],psi=PolyGamma[0,2 I v],psi1=PolyGamma[1,1/2+I v],t,delta},
  t = (g/h^2) ((1-w^-2) (psi^2-psi1-4 u psi/(1+2 u)-2/(1+2 u))
      -4 I Pi/(1+w^2) (psi-2 u/(1+2 u))
      +Pi^2 (1+4 (w^2-1)/(1+w^2)^2))
    -I ((2+u) w-u/w)/(1+2 u)
    +I (-w masters[[1]]/(4 u^2)+2 u masters[[2]]/w+masters[[4]]/w);
  delta = 4 g/((1-2 u) (3+2 u) h^2)
      ((1-w^-2) ((1+4 u)/(1+2 u)+2 u psi)-4 I Pi u/(1+w^2))
    +2 I/((1-2 u) (3+2 u)) ((8+10 u+3 u^2) w-(6+8 u+3 u^2)/w)
    -2 I u (1+2 u) masters[[3]]/((1-2 u) w);
  g {(1+2 u)/(1-2 u) t+delta,t}
];

(* Cauchy evaluation of the exact F4 derivative, as in the source notebook.
   For F(z)=3F2reg(a,b,1;1+z,1+z;1), F4=[z^2]F(z)-d_c^2 2F1reg(a,b;c;1)|1.
   Radius is a differentiation contour, NOT epsilon in a regulated integral.
   The number of contour samples grows with the requested precision. *)
mixedMaster[v_,p_Integer,n_Integer:0] := Module[
    {vv=N[v,p+20],a,b,rad=1/128,nn,vals,gauss,pure},
  a=2 I vv; b=1/2-I vv;
  nn=If[n==0,8 Ceiling[(p+8)/(8 Log[10,64])],n];
  vals=Table[With[{z=N[rad Exp[2 Pi I j/nn],p+20]},
    HypergeometricPFQ[{a,b,1},{1+z,1+z},1]/Gamma[1+z]^2 Exp[-4 Pi I j/nn]],
    {j,0,nn-1}];
  gauss=Gamma[1-a-b]/(Gamma[1-a] Gamma[1-b]);
  pure=gauss ((PolyGamma[0,1-a-b]-PolyGamma[0,1-a]-PolyGamma[0,1-b])^2
    +PolyGamma[1,1-a-b]-PolyGamma[1,1-a]-PolyGamma[1,1-b]);
  N[Total[vals]/(nn rad^2)-pure,p]
];

gaugeAssemble[v_,{pl_,pt_}] := Module[{p1,r,loop,k,box},
  p1=(2+I v)/(2 v) Sinh[Pi v] Gamma[2+2 I v];
  r=3 (3 I-2 v)^2/(I-2 v)^2;
  loop={
    (-3 I+2 v) (-4-5 I v+2 v^2)/(2 (-2 I+v) (-I+2 v)^2),
    -v (-3 I+2 v)/(2 (-2 I+v) (-I+2 v)),
    3 v (-3 I+2 v)/(2 (-2 I+v) (-I+2 v)),
    -v (-I+v) (-3 I+2 v)/((-2 I+v) (-I+2 v)^2),
    (-I+v) (-3 I+2 v)^2/(2 (-2 I+v) (-I+2 v)^2)};
  k=Gamma[-I v]^4 f0[-I v,-I v]/(4 Pi)^(7/2);
  box=(pl-pt)^2 loop+pt (2 pl+pt) r/3 {1,0,0,0,0};
  AssociationThread[gaugeKeys,k Join[{p1^2 r,p1 (pl+2 pt) r/3},box]]
];
GaugeCoefficientExpressions[v_] := gaugeAssemble[v,gaugeHard[v,gaugeMasters[v]]];

(* One ordinary 3F2(1) remains in the exact top finite part. *)
topMaster[v_] := HypergeometricPFQ[{1,1+I v,1+2 I v},{2+2 I v,2+2 I v},1];
topHard[v_,ff_] := Module[{ee,t,r1,r2},
  ee=Pi^2 Gamma[1+2 I v]/(Cosh[Pi v] Gamma[1+I v]^2);
  t=Cosh[Pi v] Gamma[1+2 I v];
  r1=(5+2 I v) (4+9 I v-4 v^2)/(2 (1+I v) (2+I v) (1+2 I v));
  r2=I v (5+2 I v)/(2 (1+I v) (2+I v));
  t {ee-2 ff/(1+2 I v)^2-r1,-(1+2 I v) ee+2 ff/(1+2 I v)+r2}
];
topAssemble[v_,{p1_,p2_}] := Module[{rho=1/2-I v,bb},
  bb=Gamma[rho]^4 f0[rho,rho]/
    (4 (4 Pi)^(7/2) (rho-3) (rho-2) (4 rho-7) (4 rho-5));
  AssociationThread[angularKeys,bb {
    (rho-1) (2 rho-3) (4 (rho-3) (3 p1^2+2 p1 p2)+(2 rho-5) p2^2),
    -p2^2 (rho-1) (2 rho-3)^2,
    3 p2^2 (rho-1) (2 rho-3)^2,
    -2 p2^2 (rho-1) (rho-2) (2 rho-3),
    p2^2 (rho-2) (2 rho-5) (2 rho-3)}]
];
TopCoefficientExpressions[v_] := topAssemble[v,topHard[v,topMaster[v]]];

validInputQ[v_,p_] := IntegerQ[p] && p>=16 && TrueQ[v>0] &&
  TrueQ[Im[v]==0] && FreeQ[v,_Real?MachineNumberQ] && Precision[v]>=p;
inputFailure[] := Failure["Precision",<|"MessageTemplate"->
  "Use real nu>0, integer precision>=16, and exact or genuinely high-precision input (e.g. 3/5, not 0.6)."|>];
finiteAssociationQ[a_] := AssociationQ[a] && VectorQ[Values[a],NumericQ] &&
  FreeQ[a,Indeterminate|DirectedInfinity[_]];
GaugeCoefficients[v_?NumericQ,p_Integer:40] := Module[{vv,mv,result},
  If[!validInputQ[v,p],Return[inputFailure[]]];
  vv=N[v,p+40];
  mv=Append[N[gaugeLeadingMasters[vv],p+20],mixedMaster[v,p+20]];
  result=N[gaugeAssemble[vv,gaugeHard[vv,mv]],p];
  If[finiteAssociationQ[result],result,Failure["GaugeEvaluation",<|"Nu"->v|>]]
];
TopCoefficients[v_?NumericQ,p_Integer:50] := Module[{vv,result},
  If[!validInputQ[v,p],Return[inputFailure[]]];
  vv=N[v,p+30];
  result=N[topAssemble[vv,topHard[vv,N[topMaster[vv],p+20]]],p];
  If[finiteAssociationQ[result],result,Failure["TopEvaluation",<|"Nu"->v|>]]
];
AngularCoefficient[c_Association,x_,y_,phi_] := Lookup[c,angularKeys] . {
  1,x^2+y^2,x^2 y^2,4 x y Sqrt[1-x^2] Sqrt[1-y^2] Cos[phi],
  (1-x^2) (1-y^2) Cos[2 phi]};
End[];
EndPackage[];
