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
Begin["`Private`"];
coefficientSignature=FileHash[FileNameJoin[{DirectoryName[$InputFileName],"coefficients.wl"}],"SHA256"];
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
End[];
EndPackage[];
