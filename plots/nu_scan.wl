BeginPackage["SignalPlots`NuScan`",{"SignalPlots`Data`","SignalPlots`Style`"}];
NuScanPlot::usage = "NuScanPlot[precomputedScan,\"LogMagnitude\"|\"Phase\"] plots stored coefficients without a physics evaluation.";
Begin["`Private`"];
NuScanPlot[scan_Association,measure_String] := Module[{d,styles,curves,curveStyles,phase,pos},
  d=NuPlotData[scan,measure]; If[FailureQ[d],Return[d]];
  styles=CurveStyles[Length[d["CoefficientNames"]],d["Model"]];
  curves=Join@@d["Segments"];
  curveStyles=Flatten[MapThread[ConstantArray[#1,Length[#2]]&,{styles,d["Segments"]}],1];
  phase=measure=="Phase";
  pos=$SignalPlotStyle[If[phase,"PhaseLegendPosition","MagnitudeLegendPosition"]];
  ListLinePlot[N[curves,MachinePrecision],Evaluate[Sequence@@ScanPlotOptions[]],
    FrameLabel->{Style["\[Nu]",Italic],If[phase,
      Row[{"arg ",Style["C",Italic]," (rad)"}],Row[{"ln |",Style["C",Italic],"|"}]]},
    PlotStyle->curveStyles,
    PlotLegends->CoefficientLegend[d["CoefficientNames"],styles,pos],
    PlotRange->If[phase,{All,{-Pi,Pi}},All]]
];
End[];
EndPackage[];
