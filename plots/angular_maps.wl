BeginPackage["SignalPlots`AngularMaps`",{"SignalPlots`Data`","SignalPlots`Style`"}];
AngularMapPlot::usage = "AngularMapPlot[precomputedData,\"Magnitude\"|\"Phase\",options] plots stored complex angular rows.";
Options[AngularMapPlot]={"AmplitudeUnit"->1,"BandEdges"->Automatic};
Begin["`Private`"];
AngularMapPlot[data_Association,measure_String,OptionsPattern[]] := Module[{unit,rows,edges,range},
  unit=If[measure=="Magnitude",OptionValue["AmplitudeUnit"],1];
  If[!TrueQ[unit>0],Return[Failure["AmplitudeUnit",<|"Unit"->unit|>]]];
  rows=AngularPlotData[data,measure,unit]; If[FailureQ[rows],Return[rows]];
  range=MinMax[rows[[All,3]]];
  edges=Replace[OptionValue["BandEdges"],Automatic:>Subdivide[First[range],Last[range],16]];
  If[!TrueQ[First[edges]<=First[range] && Last[range]<=Last[edges]],
    Return[Failure["BandRange",<|"DataRange"->range,"BandEdges"->edges|>]]];
  ListContourPlot[Developer`ToPackedArray[N[rows,MachinePrecision]],
    Evaluate[Sequence@@MapPlotOptions[data["Nu"],data["Phi3"]]],
    Contours->Rest[Most[edges]],ContourShading->MapBandColors[edges],
    PlotLegends->MapColorLegend[edges,MapMeasureLabel[measure,data["Model"],unit]]]
];
End[];
EndPackage[];
