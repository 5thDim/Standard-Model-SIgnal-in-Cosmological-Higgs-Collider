BeginPackage["SignalPlots`Data`",{"SignalPlots`Coefficients`"}];
DefaultGaugeGrid::usage = "DefaultGaugeGrid[] is the exact original 100-point nu grid.";
DefaultTopGrid::usage = "DefaultTopGrid[] is the exact original 171-point nu grid.";
GenerateNuScan::usage = "GenerateNuScan[model,grid,p] retains every complex coefficient at every nu.";
GenerateAngularData::usage = "GenerateAngularData[model,nu,phi3,cosGrid,p] computes coefficients ONCE and stores {x,y,complex C} rows, y outermost.";
SavePlotData::usage = "SavePlotData[path,data] writes portable high-precision WXF data.";
LoadPlotData::usage = "LoadPlotData[path] reloads WXF data without recomputing special functions.";
NuPlotData::usage = "NuPlotData[scan,measure] projects complex rows to LogMagnitude or principal Phase; phase curves are split at jumps > Pi.";
AngularPlotData::usage = "AngularPlotData[data,measure,unit:1] projects complex rows to Magnitude or principal Phase. Magnitude is divided only by the explicit display unit.";
CachedNuScan::usage = "CachedNuScan[path,model,grid,p,recompute:False] reuses matching WXF data or computes and saves it.";
CachedAngularData::usage = "CachedAngularData[path,model,nu,phi3,grid,p,recompute:False] reuses matching WXF data or computes and saves it.";
GenerateClockAngularData::usage = "GenerateClockAngularData[assembled,nu,phi3,thetaGrid,p:60] stores normalized complex clock amplitudes on a theta1/theta3 grid. The assembled coefficients are reused without special-function evaluation.";
CachedClockAngularData::usage = "CachedClockAngularData[path,assembled,nu,phi3,thetaGrid,p,recompute:False] caches complex curvature clock data, including physical parameters and assembly/source signatures.";
ClockMassPlotData::usage = "ClockMassPlotData[clockMassScan,\"Magnitude\"|\"Phase\"] projects stored complex amplitudes versus m/H. Zero amplitudes and principal-phase jumps are not joined.";
GenerateClockShapeData::usage = "GenerateClockShapeData[shape,logGrid,p:60] evaluates signed T, T/E and species contributions on x=-log10(r1), y=-log10(r2). Matrices have y outermost.";
CachedClockShapeData::usage = "CachedClockShapeData[path,shape,logGrid,p,recompute:False] caches both signed shape maps with physics, grid and source metadata.";
GenerateReducedClockShapeData::usage = "GenerateReducedClockShapeData[shape,uniformLogGrid,p:60] computes only T/E, using its dependence on x+y to evaluate each distinct sum once. It does not generate T or E grids.";
CachedReducedClockShapeData::usage = "CachedReducedClockShapeData[path,shape,uniformLogGrid,p,recompute:False] caches a T/E-only shape map separately from full T and T/E data.";
Begin["`Private`"];
coefficientSignature=FileHash[FileNameJoin[{DirectoryName[$InputFileName],"coefficients.wl"}],"SHA256"];
clockDataSignature=FileHash[$InputFileName,"SHA256"];
observableSignature=FileHash[FileNameJoin[{DirectoryName[$InputFileName],"observables.wl"}],"SHA256"];
DefaultGaugeGrid[] := Range[1/100,1,1/100];
DefaultTopGrid[] := Join[Range[1/10,1,1/100],Range[21/20,5,1/20]];
coefficientFunction["Gauge"] := GaugeCoefficients;
coefficientFunction["Top"] := TopCoefficients;
coefficientFunction[_] := $Failed;
GenerateNuScan[model_String,grid_List,p_Integer:40] := Module[{f,rows},
  f=coefficientFunction[model];
  If[f===$Failed,Return[Failure["Model",<|"Model"->model|>]]];
  rows=Table[<|"Nu"->v,"Coefficients"->f[v,p]|>,{v,grid}];
  If[!AllTrue[rows,AssociationQ[#["Coefficients"]]&],
    Return[Failure["ScanEvaluation",<|"Model"->model|>]]];
  <|"SchemaVersion"->1,"CoefficientSourceHash"->coefficientSignature,"Kind"->"NuScan","Model"->model,"Precision"->p,
    "NuGrid"->grid,"CoefficientNames"->Keys[First[rows]["Coefficients"]],"Rows"->rows|>
];
GenerateAngularData[model_String,v_,phi_,grid_List,p_Integer:60] := Module[
    {f,c,ng,ph,rows},
  f=coefficientFunction[model];
  If[f===$Failed,Return[Failure["Model",<|"Model"->model|>]]];
  If[!AllTrue[Join[grid,{phi}],FreeQ[#,_Real?MachineNumberQ] && Precision[#]>=p&] ||
      !AllTrue[grid,TrueQ[-1<=#<=1]&],
    Return[Failure["AngularGrid",<|"MessageTemplate"->"Use exact or high-precision phi3 and cosine grid in [-1,1]."|>]]];
  c=f[v,p]; If[FailureQ[c],Return[c]];
  ng=N[grid,p]; ph=N[phi,p];
  rows=Flatten[Table[{x,y,N[AngularCoefficient[c,x,y,ph],p]},
    {y,ng},{x,ng}],1];
  <|"SchemaVersion"->1,"CoefficientSourceHash"->coefficientSignature,"Kind"->"AngularGrid","Model"->model,"Nu"->v,
    "Phi3"->phi,"Precision"->p,"CosGrid"->grid,"GridDimensions"->{Length[grid],Length[grid]},
    "Coefficients"->c,"Rows"->rows|>
];
SavePlotData[path_String,data_Association] := (
  If[!DirectoryQ[DirectoryName[path]],CreateDirectory[DirectoryName[path],CreateIntermediateDirectories->True]];
  Export[path,data,"WXF"]);
LoadPlotData[path_String] := Import[path,"WXF"];
cached[path_,metadata_,make_,recompute_] := Module[{data},
  If[!TrueQ[recompute] && FileExistsQ[path],
    data=LoadPlotData[path];
    If[AssociationQ[data] && Lookup[data,"CoefficientSourceHash",None]===coefficientSignature &&
      AllTrue[Keys[metadata],Lookup[data,#,None]===metadata[#]&],Return[data]]];
  data=make[]; If[AssociationQ[data],SavePlotData[path,data]]; data
];
CachedNuScan[path_String,model_String,grid_List,p_Integer,recompute_:False] :=
  cached[path,<|"Kind"->"NuScan","Model"->model,"NuGrid"->grid,"Precision"->p|>,
    Function[GenerateNuScan[model,grid,p]],recompute];
CachedAngularData[path_String,model_String,v_,phi_,grid_List,p_Integer,recompute_:False] :=
  cached[path,<|"Kind"->"AngularGrid","Model"->model,"Nu"->v,"Phi3"->phi,
      "CosGrid"->grid,"Precision"->p|>,
    Function[GenerateAngularData[model,v,phi,grid,p]],recompute];

clockMetadata[assembly_,v_,phi_,grid_,p_] := <|
  "Kind"->"ClockAngularGrid","Species"->assembly["Species"],"Nu"->v,"Phi3"->phi,
  "Precision"->p,"ThetaGrid"->grid,"Parameters"->assembly["Parameters"],
  "AssemblyHash"->Hash[assembly["ComplexAmplitudeCoefficients"],"SHA256"],
  "ObservableSourceHash"->observableSignature,"DataSourceHash"->clockDataSignature|>;
GenerateClockAngularData[assembly_Association,v_,phi_,grid_List,p_Integer:60] := Module[
    {metadata,theta,u,w,basis,values,n,rows},
  If[!ContainsAll[Keys[assembly],{"Species","Parameters","ComplexAmplitudeCoefficients"}],
    Return[Failure["ClockAssembly",<||>]]];
  If[KeyExistsQ[assembly,"Nu"] && !TrueQ[assembly["Nu"]==v],
    Return[Failure["ClockNuMismatch",<|"AssemblyNu"->assembly["Nu"],"RequestedNu"->v|>]]];
  If[!IntegerQ[p] || p<16 || Length[grid]<2 || !AllTrue[Differences[grid],TrueQ[#>0]&] ||
      !AllTrue[Join[grid,{phi,v}],FreeQ[#,_Real?MachineNumberQ] && Precision[#]>=p&] ||
      !TrueQ[v>0] || !AllTrue[grid,TrueQ[0<=#<=Pi]&],
    Return[Failure["ClockAngularGrid",<|"MessageTemplate"->
      "Use nu>0 and an increasing exact/high-precision theta grid in [0,Pi]."|>]]];
  metadata=clockMetadata[assembly,v,phi,grid,p];
  theta=N[grid,p]; n=Length[grid];
  u=N[Cos[grid]^2,p]; w=N[Sin[grid]^2,p];
  basis=Lookup[assembly["ComplexAmplitudeCoefficients"],{"A1","A2","A3","B","D"}];
  values=ConstantArray[basis[[1]],{n,n}]+basis[[2]] Outer[Plus,u,u]+
    basis[[3]] Outer[Times,u,u]+
    4 basis[[4]] N[Cos[phi],p] Outer[Times,N[Cos[grid] Sin[grid],p],N[Cos[grid] Sin[grid],p]]+
    basis[[5]] N[Cos[2 phi],p] Outer[Times,w,w];
  rows=Flatten[MapIndexed[{theta[[#2[[2]]]],theta[[#2[[1]]]],#1}&,values,{2}],1];
  If[!FreeQ[rows,Indeterminate|DirectedInfinity[_]],
    Return[Failure["ClockAngularEvaluation",<|"Species"->assembly["Species"]|>]]];
  Join[<|"SchemaVersion"->1,"CoefficientSourceHash"->coefficientSignature|>,metadata,
    <|"GridDimensions"->{n,n},"ClockAssembly"->assembly,"Rows"->rows|>]
];
CachedClockAngularData[path_String,assembly_Association,v_,phi_,grid_List,p_Integer,recompute_:False] :=
  cached[path,clockMetadata[assembly,v,phi,grid,p],
    Function[GenerateClockAngularData[assembly,v,phi,grid,p]],recompute];
NuPlotData[scan_Association,measure_String] := Module[{names,series,project},
  names=scan["CoefficientNames"];
  project=Switch[measure,"LogMagnitude",Log[Abs[#]]&,"Phase",Arg,
    _,Return[Failure["Measure",<|"Measure"->measure|>]]];
  series=Table[({#["Nu"],project[#["Coefficients"][key]]}& /@ scan["Rows"]),{key,names}];
  <|"Model"->scan["Model"],"Measure"->measure,"CoefficientNames"->names,
    "Series"->series,"Segments"->If[measure=="Phase",
      (Split[#,Abs[Last[#1]-Last[#2]]<=Pi&]& /@ series),List /@ series]|>
];
AngularPlotData[data_Association,measure_String,unit_:1] := Module[{project},
  project=Switch[measure,"Magnitude",Abs[#]/unit&,"Phase",Arg,
    _,Return[Failure["Measure",<|"Measure"->measure|>]]];
  {#[[1]],#[[2]],project[#[[3]]]}& /@ data["Rows"]
];
ClockMassPlotData[data_Association,measure_String] := Module[{project,blocks,segments},
  If[Lookup[data,"Kind",None]=!="ClockMassScan",
    Return[Failure["ClockMassData",<||>]]];
  project=Switch[measure,"Magnitude",Abs,"Phase",Arg,
    _,Return[Failure["Measure",<|"Measure"->measure|>]]];
  blocks=Select[SplitBy[data["Rows"],TrueQ[#["ComplexAmplitude"]==0]&],
    !TrueQ[First[#]["ComplexAmplitude"]==0]&];
  segments=({#["MassOverH"],project[#["ComplexAmplitude"]]}& /@ #)& /@ blocks;
  If[measure=="Phase",segments=Join@@(Split[#,Abs[Last[#1]-Last[#2]]<=Pi&]& /@ segments)];
  <|"Species"->data["Species"],"Measure"->measure,"Segments"->segments|>
];
shapeMetadata[shape_,grid_,p_] := <|"Kind"->"ClockShapeGrid","Precision"->p,
  "LogGrid"->grid,"Parameters"->shape["Parameters"],"Angles"->shape["Angles"],
  "ShapeHash"->Hash[shape,"SHA256"],"ObservableSourceHash"->observableSignature,
  "DataSourceHash"->clockDataSignature|>;
GenerateClockShapeData[shape_Association,grid_List,p_Integer:60] := Module[
    {r,q,e,logq,o,parts,reduced,raw},
  If[Lookup[shape,"Kind",None]=!="ClockCollapsedShape" || Length[grid]<2 || p<16 ||
      !AllTrue[Differences[grid],TrueQ[#>0]&] ||
      !AllTrue[grid,TrueQ[#>=0] && FreeQ[#,_Real?MachineNumberQ] && Precision[#]>=p&],
    Return[Failure["ClockShapeGrid",<||>]]];
  r=N[10^-grid,p]; q=Outer[Times,r,r]; e=(q/Outer[Plus,r,r])^3;
  logq=-N[Log[10],p] Outer[Plus,N[grid,p],N[grid,p]]-N[Log[4],p];
  o=shape["Oscillators"];
  parts=AssociationMap[Function[s,N[If[s=="t",q,1] Re[
    o[s]["ComplexAmplitude"] Exp[2 I N[o[s]["Nu"],p] logq]],p]],{"W","Z","t"}];
  reduced=Total[Values[parts]]; raw=e reduced;
  If[!FreeQ[{e,parts,raw,reduced},Indeterminate|DirectedInfinity[_]],
    Return[Failure["ClockShapeEvaluation",<||>]]];
  Join[<|"SchemaVersion"->1,"CoefficientSourceHash"->coefficientSignature|>,shapeMetadata[shape,grid,p],
    <|"GridDimensions"->{Length[grid],Length[grid]},"Shape"->shape,
      "Envelope"->e,"SpeciesReduced"->parts,"TOverE"->reduced,"T"->raw|>]
];
CachedClockShapeData[path_String,shape_Association,grid_List,p_Integer,recompute_:False] :=
  cached[path,shapeMetadata[shape,grid,p],Function[GenerateClockShapeData[shape,grid,p]],recompute];
GenerateReducedClockShapeData[shape_Association,grid_List,p_Integer:60] := Module[
    {n,sums,logq,o,series,total,matrix},
  If[Lookup[shape,"Kind",None]=!="ClockCollapsedShape" || Length[grid]<2 || p<16 ||
      !AllTrue[Differences[grid],TrueQ[#>0]&] || !TrueQ[Equal@@Differences[grid]] ||
      !AllTrue[grid,TrueQ[#>=0] && FreeQ[#,_Real?MachineNumberQ] && Precision[#]>=p&],
    Return[Failure["ReducedClockShapeGrid",<|"MessageTemplate"->"Use an increasing, uniform exact/high-precision nonnegative log grid."|>]]];
  n=Length[grid]; sums=Subdivide[2 First[grid],2 Last[grid],2 n-2];
  logq=-N[Log[10],p] N[sums,p]-N[Log[4],p]; o=shape["Oscillators"];
  series=AssociationMap[Function[s,N[If[s=="t",N[10^-sums,p],1] Re[
    o[s]["ComplexAmplitude"] Exp[2 I N[o[s]["Nu"],p] logq]],p]],{"W","Z","t"}];
  total=Total[Values[series]];
  matrix=Table[Take[total,{j,j+n-1}],{j,n}];
  If[!MatrixQ[matrix,NumberQ] || !FreeQ[matrix,Indeterminate|DirectedInfinity[_]],
    Return[Failure["ReducedClockShapeEvaluation",<||>]]];
  Join[<|"SchemaVersion"->1,"CoefficientSourceHash"->coefficientSignature|>,shapeMetadata[shape,grid,p],
    <|"ComputedMeasures"->{"TOverE"},"GridDimensions"->{n,n},"Shape"->shape,
      "LogSumGrid"->sums,"SpeciesReducedSeries"->series,"TOverE"->matrix|>]
];
CachedReducedClockShapeData[path_String,shape_Association,grid_List,p_Integer,recompute_:False] :=
  cached[path,Join[shapeMetadata[shape,grid,p],<|"ComputedMeasures"->{"TOverE"}|>],
    Function[GenerateReducedClockShapeData[shape,grid,p]],recompute];
End[];
EndPackage[];
