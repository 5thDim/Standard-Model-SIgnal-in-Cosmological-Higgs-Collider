BeginPackage["SignalPlots`AngularMaps`",{"SignalPlots`Data`","SignalPlots`Style`"}];
ClockAngularMapPlot::usage = "ClockAngularMapPlot[clockData,\"BandEdges\"->edges] plots actual tNL on theta1/theta3 axes. Supply common edges to compare phi3 slices.";
ClockShapeMapPlot::usage = "ClockShapeMapPlot[shapeData,\"T\"|\"TOverE\",options] plots a signed linear map on -log10(r1),-log10(r2), with a symmetric zero-centered color range.";
Options[ClockAngularMapPlot]={"BandEdges"->Automatic};
Options[ClockShapeMapPlot]={"BandCount"->20,"ColorRange"->Automatic,"ShowParameters"->True};
Begin["`Private`"];
ClockAngularMapPlot[data_Association,OptionsPattern[]] := Module[{rows,edges,range,grid,n,regular,input,extra},
  If[Lookup[data,"Kind",None]=!="ClockAngularGrid",
    Return[Failure["ClockAngularData",<||>]]];
  rows=AngularPlotData[data,"Magnitude"];
  range=MinMax[rows[[All,3]]];
  edges=Replace[OptionValue["BandEdges"],Automatic:>Subdivide[First[range],Last[range],16]];
  If[!TrueQ[First[edges]<=First[range] && Last[range]<=Last[edges]],
    Return[Failure["BandRange",<|"DataRange"->range,"BandEdges"->edges|>]]];
  grid=Lookup[data,"ThetaGrid",{}]; n=Length[grid];
  regular=n>=2 && Length[rows]==n^2 && TrueQ[Equal@@Differences[grid]];
  If[regular,
    input=Partition[rows[[All,3]],n];
    extra={DataRange->N[{{First[grid],Last[grid]},{First[grid],Last[grid]}},MachinePrecision]},
    input=rows; extra={}];
  (* Only the display copy uses machine numbers; the cached complex data keep
     their original precision. Matrix input preserves a regular theta grid. *)
  ListContourPlot[Developer`ToPackedArray[N[input,MachinePrecision]],
    Evaluate[Sequence@@extra],
    Evaluate[Sequence@@ClockMapPlotOptions[data["Species"],data["Nu"],data["Phi3"]]],
    Contours->N[Rest[Most[edges]],MachinePrecision],ContourShading->MapBandColors[edges],
    PlotLegends->ClockMapColorLegend[N[edges,MachinePrecision],data["Species"]]]
];
ClockShapeMapPlot[data_Association,measure_String,OptionsPattern[]] := Module[{values,grid,range,edges,n},
  If[Lookup[data,"Kind",None]=!="ClockShapeGrid" || !MemberQ[{"T","TOverE"},measure] || !KeyExistsQ[data,measure],
    Return[Failure["ClockShapeData",<||>]]];
  values=data[measure]; grid=data["LogGrid"]; n=OptionValue["BandCount"];
  If[!IntegerQ[n] || n<2 || OddQ[n],Return[Failure["ShapeBandCount",<||>]]];
  range=Replace[OptionValue["ColorRange"],Automatic:>With[{a=Max[Abs[Flatten[values]]]}, {-a,a}]];
  If[!MatchQ[range,{_?NumericQ,_?NumericQ}] ||
      !TrueQ[range[[1]]<range[[2]] && range[[1]]<=Min[values] && Max[values]<=range[[2]]],
    Return[Failure["ShapeColorRange",<||>]]];
  edges=N[Subdivide[First[range],Last[range],n],MachinePrecision];
  ListContourPlot[Developer`ToPackedArray[N[values,MachinePrecision]],
    DataRange->N[ConstantArray[{First[grid],Last[grid]},2],MachinePrecision],
    Evaluate[Sequence@@ClockShapePlotOptions[data["Angles"],data["Parameters"],OptionValue["ShowParameters"]]],
    Contours->Rest[Most[edges]],ContourShading->MapBandColors[edges],
    PlotLegends->ClockShapeColorLegend[edges,measure]]
];
End[];
EndPackage[];
