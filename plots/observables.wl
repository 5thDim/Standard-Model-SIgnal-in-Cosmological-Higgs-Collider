(* Curvature clock observables in SMS in CHC, Section 5.
   Lu et al. (2019), arXiv:1907.07390, Fig. 4 benchmark and Section 2 Pzeta.
   The user-selected Higgs background is h0/H = 1.
   Nu and external momentum directions remain independent function inputs. *)
BeginPackage["SignalPlots`Observables`", {"SignalPlots`Coefficients`"}];
Lu2019Parameters::usage = "Lu2019Parameters[] returns exact benchmark parameters, with h0OverH=1.";
AssembleClockCoefficients::usage = "AssembleClockCoefficients[species,rawCoefficients,parameters:Lu2019Parameters[]] assembles the five angular coefficients of C and its normalized complex clock amplitude. Species is W, Z, or t; raw coefficients come from GaugeCoefficients or TopCoefficients.";
ClockCoefficients::usage = "ClockCoefficients[species,nu,precision:60,parameters:Lu2019Parameters[]] evaluates raw coefficients once and assembles the clock coefficients. Reuse this result when scanning angles.";
ClockObservables::usage = "ClockObservables[assembled,theta1,theta3,phi3] returns C, ComplexAmplitude, tNL, and Phase. Angles are in radians. ComplexAmplitude=tNL Exp[I Phase]; top Phase includes the overall minus sign. At zero amplitude Phase is Missing[ZeroAmplitude].";
AssembleClockMassScan::usage = "AssembleClockMassScan[species,rawNuScan,{theta1,theta3,phi3},parameters:Lu2019Parameters[]] converts stored raw coefficients to complex curvature clock amplitudes versus m/H, without evaluating special functions.";
Lu2019MassBenchmark::usage = "Lu2019MassBenchmark[parameters:Lu2019Parameters[]] gives Lu (2019), Eqs. (11),(45), loop-dominated W/Z masses and the background top mass, in Hubble units.";
AssembleClockShape::usage = "AssembleClockShape[speciesAssemblies,{theta1,theta3,phi3}] stores the W/Z/top complex amplitudes and frequencies at fixed masses and directions.";
ClockCollapsedSignal::usage = "ClockCollapsedSignal[shape,r1,r2] evaluates the real leading collapsed trispectrum and T/E, with r1=ks/k12, r2=ks/k34 and E=(r1 r2/(r1+r2))^3.";
Begin["`Private`"];

angularKeys = {"A1", "A2", "A3", "B", "D"};
parameterKeys = {"g", "cW", "yt", "lambda", "Rh", "Pzeta", "h0OverH"};
Lu2019Parameters[] := <|
  "g" -> 3/5, "cW" -> 22/25, "yt" -> 1, "lambda" -> 1/100,
  "Rh" -> 7/50, "Pzeta" -> 2/10^9, "h0OverH" -> 1|>;

validParametersQ[p_Association] := ContainsAll[Keys[p], parameterKeys] &&
  AllTrue[Lookup[p, parameterKeys], NumericQ[#] && TrueQ[Im[#] == 0] &] &&
  TrueQ[p["g"] > 0 && 0 < p["cW"] <= 1 && p["yt"] > 0 &&
    p["lambda"] >= 0 && 0 <= p["Rh"] <= 1 &&
    p["Pzeta"] > 0 && p["h0OverH"] >= 0];

AssembleClockCoefficients[species_String, raw_Association,
    parameters_Association:Lu2019Parameters[]] := Module[
    {required, missing, g, cw, yt, rh, pzeta, h, q, basis, factor, sign},
  If[!MemberQ[{"W", "Z", "t"}, species],
    Return[Failure["Species", <|"Species" -> species|>]]];
  If[!validParametersQ[parameters],
    Return[Failure["Parameters", <|"RequiredKeys" -> parameterKeys|>]]];
  required = If[species == "t", angularKeys, Join[{"CBub", "CTri"}, angularKeys]];
  missing = Complement[required, Keys[raw]];
  If[missing =!= {}, Return[Failure["CoefficientKeys", <|"MissingKeys" -> missing|>]]];
  {g, cw, yt, rh, pzeta, h} = Lookup[parameters,
    {"g", "cW", "yt", "Rh", "Pzeta", "h0OverH"}];
  If[species == "t",
    basis = KeyTake[raw, angularKeys];
    factor = 3 yt^4 rh^4/(32 Pi^2 pzeta);
    sign = -1,
    q = If[species == "W", g^2 h^2, g^2 h^2/cw^2];
    basis = AssociationThread[angularKeys, q^2 Lookup[raw, angularKeys]];
    AssociateTo[basis, "A1" -> raw["CBub"] - 2 q raw["CTri"] + q^2 raw["A1"]];
    factor = If[species == "W", g^4 rh^4/(128 Pi^2 pzeta),
      g^4 rh^4/(256 Pi^2 cw^4 pzeta)];
    sign = 1
  ];
  <|"Species" -> species, "Parameters" -> parameters,
    "RawCoefficients" -> raw,
    "AngularCoefficients" -> basis,
    "ComplexAmplitudeCoefficients" -> Map[sign factor # &, basis],
    "AmplitudePrefactor" -> factor, "PhaseSign" -> sign|>
];

ClockCoefficients[species_String, nu_?NumericQ, precision_Integer:60,
    parameters_Association:Lu2019Parameters[]] := Module[{raw, assembled},
  If[!MemberQ[{"W", "Z", "t"}, species],
    Return[Failure["Species", <|"Species" -> species|>]]];
  If[!validParametersQ[parameters],
    Return[Failure["Parameters", <|"RequiredKeys" -> parameterKeys|>]]];
  raw = If[species == "t", TopCoefficients[nu, precision], GaugeCoefficients[nu, precision]];
  If[FailureQ[raw], Return[raw]];
  assembled = AssembleClockCoefficients[species, raw, parameters];
  If[FailureQ[assembled], Return[assembled]];
  Join[assembled, <|"Nu" -> nu,
    "MassOverH" -> If[species == "t", nu, Sqrt[nu^2 + 1/4]],
    "Precision" -> precision|>]
];

ClockObservables[assembled_Association, theta1_, theta3_, phi3_] := Module[{c, amplitude},
  If[!ContainsAll[Keys[assembled],
      {"Species", "AngularCoefficients", "AmplitudePrefactor", "PhaseSign"}],
    Return[Failure["ClockAssembly", <|"MessageTemplate" ->
      "Use the result of AssembleClockCoefficients or ClockCoefficients."|>]]];
  c = AngularCoefficient[assembled["AngularCoefficients"], Cos[theta1], Cos[theta3], phi3];
  amplitude = assembled["PhaseSign"] assembled["AmplitudePrefactor"] c;
  <|"Species" -> assembled["Species"], "C" -> c,
    "ComplexAmplitude" -> amplitude, "tNL" -> Abs[amplitude],
    "Phase" -> If[TrueQ[amplitude == 0], Missing["ZeroAmplitude"], Arg[amplitude]]|>
];

AssembleClockMassScan[species_String, scan_Association, angles_List,
    parameters_Association:Lu2019Parameters[]] := Module[{model, rows, assembled, value, nu},
  model=If[species=="t","Top","Gauge"];
  If[!MemberQ[{"W","Z","t"},species] || Lookup[scan,"Kind",None]=!="NuScan" ||
      Lookup[scan,"Model",None]=!=model || !KeyExistsQ[scan,"Rows"],
    Return[Failure["ClockMassScan",<|"Species"->species|>]]];
  If[Length[angles]=!=3 || !AllTrue[angles,NumericQ[#] && TrueQ[Im[#]==0]&],
    Return[Failure["ClockAngles",<||>]]];
  If[!validParametersQ[parameters],Return[Failure["Parameters",<||>]]];
  rows=Table[
    nu=row["Nu"];
    assembled=AssembleClockCoefficients[species,row["Coefficients"],parameters];
    If[FailureQ[assembled],Return[assembled]];
    value=ClockObservables[assembled,Sequence@@angles];
    <|"Nu"->nu,"MassOverH"->If[species=="t",nu,Sqrt[nu^2+1/4]],
      "C"->value["C"],"ComplexAmplitude"->value["ComplexAmplitude"]|>,
    {row,scan["Rows"]}];
  If[!VectorQ[Lookup[rows,"ComplexAmplitude"],NumericQ] ||
      !FreeQ[rows,Indeterminate|DirectedInfinity[_]],
    Return[Failure["ClockMassEvaluation",<|"Species"->species|>]]];
  <|"SchemaVersion"->1,"Kind"->"ClockMassScan","Species"->species,
    "Parameters"->parameters,"Angles"->angles,"Precision"->scan["Precision"],
    "CoefficientSourceHash"->scan["CoefficientSourceHash"],"Rows"->rows|>
];

Lu2019MassBenchmark[parameters_Association:Lu2019Parameters[]] := Module[{mh2,mw2,masses,nus},
  If[!validParametersQ[parameters] || !TrueQ[parameters["lambda"]>0],
    Return[Failure["MassParameters",<||>]]];
  mh2=Sqrt[6 parameters["lambda"]/Pi^3];
  mw2=3 parameters["g"]^2/(8 Pi^2 mh2);
  masses=<|"W"->Sqrt[mw2],"Z"->Sqrt[mw2]/parameters["cW"],
    "t"->parameters["yt"] parameters["h0OverH"]/Sqrt[2]|>;
  If[!TrueQ[masses["W"]>1/2 && masses["Z"]>1/2 && masses["t"]>0],
    Return[Failure["NonOscillatoryMasses",<|"MassOverH"->masses|>]]];
  nus=<|"W"->Sqrt[masses["W"]^2-1/4],"Z"->Sqrt[masses["Z"]^2-1/4],"t"->masses["t"]|>;
  <|"Parameters"->parameters,"HiggsMassSquaredOverH2"->mh2,"MassOverH"->masses,
    "Nu"->nus,"VectorMassPrescription"->"Lu2019LoopDominated"|>
];
AssembleClockShape[assemblies_Association,angles_List] := Module[{species={"W","Z","t"},p,oscillators},
  If[!ContainsAll[Keys[assemblies],species] || Length[angles]=!=3 ||
      !AllTrue[angles,NumericQ[#] && TrueQ[Im[#]==0]&],Return[Failure["ClockShape",<||>]]];
  p=assemblies["W"]["Parameters"];
  If[!AllTrue[species,assemblies[#]["Parameters"]===p && TrueQ[assemblies[#]["Nu"]>0]&],
    Return[Failure["ClockShapeAssemblies",<||>]]];
  oscillators=AssociationMap[Function[s,Join[
    KeyTake[assemblies[s],{"Nu","MassOverH"}],ClockObservables[assemblies[s],Sequence@@angles]]],species];
  <|"Kind"->"ClockCollapsedShape","Parameters"->p,"Angles"->angles,
    "Oscillators"->oscillators,"ClockAssemblies"->assemblies|>
];
ClockCollapsedSignal[shape_Association,r1_?NumericQ,r2_?NumericQ] := Module[{q,e,parts,o},
  If[Lookup[shape,"Kind",None]=!="ClockCollapsedShape" || !TrueQ[0<r1<=1 && 0<r2<=1],
    Return[Failure["ClockShapeRatios",<||>]]];
  q=r1 r2; e=(q/(r1+r2))^3; o=shape["Oscillators"];
  parts=AssociationMap[Function[s,If[s=="t",q,1] Re[
    o[s]["ComplexAmplitude"] Exp[2 I o[s]["Nu"] Log[q/4]]]],{"W","Z","t"}];
  <|"Envelope"->e,"SpeciesReduced"->parts,"TOverE"->Total[Values[parts]],
    "T"->e Total[Values[parts]]|>
];
End[];
EndPackage[];
