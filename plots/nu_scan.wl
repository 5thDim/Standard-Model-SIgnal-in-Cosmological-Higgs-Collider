BeginPackage["SignalPlots`NuScan`",{"SignalPlots`Data`","SignalPlots`Style`"}];
NuScanPlot::usage = "NuScanPlot[precomputedScan,\"LogMagnitude\"|\"Phase\"] plots stored coefficients without a physics evaluation.";
ClockMassScanPlot::usage = "ClockMassScanPlot[speciesAssociation,\"Magnitude\"|\"Phase\"] plots curvature clock amplitudes on a logarithmic y axis or principal phases versus m/H, without evaluating physics.";
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
ClockMassScanPlot[scans_Association,measure_String] := Module[
    {species,data,styles,curves,curveStyles,phase,plot,range,angles,options},
  species=Keys[scans];
  If[species==={} || !AllTrue[Values[scans],AssociationQ],
    Return[Failure["ClockMassData",<||>]]];
  angles=First[Values[scans]]["Angles"];
  If[!AllTrue[Values[scans],#["Angles"]===angles&],
    Return[Failure["ClockMassAngles",<||>]]];
  data=ClockMassPlotData[scans[#],measure]& /@ species;
  If[!FreeQ[data,_Failure],Return[First[Cases[data,_Failure,Infinity]]]];
  styles=ClockMassCurveStyles[species];
  curves=Join@@Lookup[data,"Segments"];
  curveStyles=Flatten[MapThread[ConstantArray[#1,Length[#2["Segments"]]]&,{styles,data}],1];
  If[curves==={},Return[Failure["ClockMassEmpty",<||>]]];
  range=MinMax[Flatten[curves,1][[All,1]]];
  phase=measure=="Phase";
  plot=If[phase,ListLinePlot,ListLogPlot];
  options=If[phase,{FrameTicks->{{
    Table[{j Pi/2,TraditionalForm[j Pi/2]},{j,-2,2}],None},{Automatic,None}}},{}];
  plot[N[curves,MachinePrecision],Evaluate[Sequence@@ClockMassPlotOptions[angles]],
    Joined->True,PlotMarkers->None,
    FrameLabel->{Row[{Style["m",Italic],"/",Style["H",Italic]}],
      If[phase,Row[{Style["\[CurlyPhi]",Italic]," (rad)"}],
        Subscript[Style["t",Italic],Style["NL",Plain]]]},
    PlotStyle->curveStyles,PlotLegends->ClockMassLegend[species,styles,measure],
    PlotRange->If[phase,{range,{-Pi,Pi}},{range,All}],
    PlotRangePadding->{{0,0},{Scaled[0.02],Scaled[0.02]}},
    Evaluate[Sequence@@options]]
];
End[];
EndPackage[];
