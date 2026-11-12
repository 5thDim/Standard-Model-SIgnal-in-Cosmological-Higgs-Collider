(* All shared presentation choices. Editing this file needs only a replot. *)
BeginPackage["SignalPlots`Style`"];
$SignalPlotStyle::usage = "$SignalPlotStyle contains the shared font, sizes, line weights, padding, palettes and legend positions.";
ScanPlotOptions::usage = "ScanPlotOptions[] returns the original gauge summary plot geometry.";
CurveStyles::usage = "CurveStyles[n,model:Automatic] gives curve styles; Gauge uses seven spectrum colors and Top uses five.";
CoefficientLegend::usage = "CoefficientLegend[names,styles,position] is the framed two-column in-plot legend.";
MapPlotOptions::usage = "MapPlotOptions[nu,phi3] returns square angular plot options and parameter inset.";
MapBandColors::usage = "MapBandColors[edges] gives the discrete red-beige-blue fills.";
MapColorLegend::usage = "MapColorLegend[edges,label] gives a vertical discrete bar on the right.";
MapMeasureLabel::usage = "MapMeasureLabel[measure,model,unit] explicitly identifies the displayed quantity and units.";
Begin["`Private`"];
$SignalPlotStyle=<|
  "Font"->"Latin Modern Roman","TickSize"->10,"LabelSize"->12,"LegendSize"->10,
  "LineThickness"->1.4,"FrameThickness"->0.8,
  "ScanImageSize"->400,"ScanAspectRatio"->1/2,"ScanPadding"->{{55,12},{42,10}},
  "MapImageSize"->470,"MapPadding"->{{54,10},{44,12}},
  "ScanColors"->{RGBColor[123/255,44/255,191/255],RGBColor[36/255,85/255,195/255],
    RGBColor[35/255,139/255,69/255],RGBColor[168/255,139/255,0],
    RGBColor[201/255,42/255,42/255]},
  "GaugeScanColors"->{RGBColor[123/255,44/255,191/255],RGBColor[36/255,85/255,195/255],
    RGBColor[0,143/255,168/255],RGBColor[35/255,139/255,69/255],
    RGBColor[168/255,139/255,0],RGBColor[232/255,117/255,0],
    RGBColor[201/255,42/255,42/255]},
  "MapColors"->{RGBColor[0.46,0.13,0.20],RGBColor[0.66,0.29,0.31],
    RGBColor[0.82,0.56,0.46],RGBColor[0.90,0.85,0.67],RGBColor[0.68,0.81,0.80],
    RGBColor[0.31,0.59,0.70],RGBColor[0.15,0.31,0.54]},
  "MagnitudeLegendPosition"->{0.75,0.77},"PhaseLegendPosition"->{0.75,0.47},
  "ParameterPosition"->{0.72,0.93},"ColorBarSize"->{16,250},
  "ColorBarMargins"->{{8,40},{0,0}}|>;
fontStyle[size_] := {FontFamily->$SignalPlotStyle["Font"],FontSize->size};
CurveStyles[n_Integer,model_:Automatic] := (Directive[#,AbsoluteThickness[$SignalPlotStyle["LineThickness"]]]& /@
  (Blend[$SignalPlotStyle[If[model==="Gauge","GaugeScanColors","ScanColors"]],#]& /@ Subdivide[0,1,n-1]));
coefficientLabel["CBub"] := Superscript[Style["C",Italic],Style["Bub",Plain]];
coefficientLabel["CTri"] := Superscript[Style["C",Italic],Style["Tri",Plain]];
coefficientLabel[s_String] := If[StringStartsQ[s,"A"],
  Subscript[Style["A",Italic],StringDrop[s,1]],Style[s,Italic]];
legendLayout[pairs_] := Grid[Flatten[#,1]& /@
  Partition[If[OddQ[Length[pairs]],Append[pairs,{"",""}],pairs],2],
  Alignment->Left,Spacings->{0.6,0.3}];
CoefficientLegend[names_List,styles_List,position_] := Placed[
  LineLegend[styles,coefficientLabel /@ names,
    LabelStyle->Join[{Black},fontStyle[$SignalPlotStyle["LegendSize"]]],LegendLayout->legendLayout,
    LegendMarkerSize->{16,8},LegendMargins->0,
    LegendFunction->(Framed[#,Background->White,FrameStyle->Black,
      RoundingRadius->0,FrameMargins->5]&)],position];
ScanPlotOptions[] := {
  Frame->True,Axes->False,Background->White,
  FrameStyle->Directive[Black,AbsoluteThickness[$SignalPlotStyle["FrameThickness"]]],
  GridLines->Automatic,GridLinesStyle->Directive[GrayLevel[0.8],AbsoluteThickness[0.5]],
  BaseStyle->fontStyle[$SignalPlotStyle["LabelSize"]],
  FrameTicksStyle->Directive[Black,Sequence@@fontStyle[$SignalPlotStyle["TickSize"]]],
  LabelStyle->Directive[Black,Sequence@@fontStyle[$SignalPlotStyle["LabelSize"]]],
  ImageSize->$SignalPlotStyle["ScanImageSize"],AspectRatio->$SignalPlotStyle["ScanAspectRatio"],
  ImagePadding->$SignalPlotStyle["ScanPadding"]};
mapPalette[s_] := Blend[$SignalPlotStyle["MapColors"],s];
MapBandColors[edges_List] := Table[mapPalette[(j-1/2)/(Length[edges]-1)],{j,Length[edges]-1}];
MapPlotOptions[v_,phi_] := {
  Frame->True,Axes->False,AspectRatio->1,
  FrameStyle->Directive[Black,AbsoluteThickness[$SignalPlotStyle["FrameThickness"]]],
  FrameLabel->{Row[{"cos ",Subscript[Style["\[Theta]",Italic],1]}],
    Row[{"cos ",Subscript[Style["\[Theta]",Italic],3]}]},
  BaseStyle->fontStyle[$SignalPlotStyle["LabelSize"]],
  FrameTicksStyle->Directive[Black,Sequence@@fontStyle[$SignalPlotStyle["TickSize"]]],
  LabelStyle->Directive[Black,Sequence@@fontStyle[$SignalPlotStyle["LabelSize"]]],
  ContourStyle->Directive[GrayLevel[0.4],Opacity[0.55],AbsoluteThickness[0.4]],
  InterpolationOrder->1,PlotRange->All,PlotRangePadding->None,
  ImageSize->$SignalPlotStyle["MapImageSize"],ImagePadding->$SignalPlotStyle["MapPadding"],
  Epilog->Inset[Framed[Style[Row[{Style["\[Nu]",Italic]," = ",v,",  ",
    Subscript[Style["\[Phi]",Italic],3]," = ",phi}],
    Sequence@@fontStyle[$SignalPlotStyle["LabelSize"]]],
    Background->White,FrameStyle->Black,FrameMargins->5,RoundingRadius->0],
    Scaled[$SignalPlotStyle["ParameterPosition"]]]};
MapColorLegend[edges_List,label_] := Placed[
  BarLegend[{mapPalette,{First[edges],Last[edges]}},Rest[Most[edges]],
    LegendLabel->Placed[Style[label,Sequence@@fontStyle[$SignalPlotStyle["LabelSize"]]],Above],
    LabelStyle->fontStyle[$SignalPlotStyle["LegendSize"]],
    LegendMarkerSize->$SignalPlotStyle["ColorBarSize"],LegendMargins->$SignalPlotStyle["ColorBarMargins"]],Right];
MapMeasureLabel[measure_,model_,unit_] := Module[{c=If[model=="Top",
    Subscript[Style["C",Italic],Style["t",Italic]],
    Superscript[Style["C",Italic],Style["Box",Plain]]]},
  If[measure=="Phase",Row[{"arg ",c," (rad)"}],
    If[unit==1,Row[{"|",c,"|"}],Row[{"|",c,"| \[Times] ",
      If[unit==10^-6,Superscript[10,6],TraditionalForm[1/unit]]}]]]
];
End[];
EndPackage[];
