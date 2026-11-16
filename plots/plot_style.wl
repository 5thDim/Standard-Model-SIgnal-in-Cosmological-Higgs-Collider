(* All shared presentation choices. Editing this file needs only a replot. *)
BeginPackage["SignalPlots`Style`"];
$SignalPlotStyle::usage = "$SignalPlotStyle contains the shared font, sizes, line weights, padding, palettes and legend positions.";
ScanPlotOptions::usage = "ScanPlotOptions[] returns the original gauge summary plot geometry.";
CurveStyles::usage = "CurveStyles[n,model:Automatic] gives curve styles; Gauge uses seven spectrum colors and Top uses five.";
CoefficientLegend::usage = "CoefficientLegend[names,styles,position] is the framed two-column in-plot legend.";
MapPlotOptions::usage = "MapPlotOptions[nu,phi3] returns square angular plot options and parameter inset.";
MapBandColors::usage = "MapBandColors[edges] gives the discrete red-beige-blue fills.";
ClockMapPlotOptions::usage = "ClockMapPlotOptions[species,nu,phi3] gives theta1/theta3 map styling and a framed in-plot species/parameter label.";
ClockMapLegendLabel::usage = "ClockMapLegendLabel[species] labels the actual dimensionless tNL, without rescaling.";
ClockMapColorLegend::usage = "ClockMapColorLegend[edges,species] gives a white-background tNL legend with black labels.";
ClockMassPlotOptions::usage = "ClockMassPlotOptions[angles] gives the shared geometry for curvature mass scans without an angle title.";
ClockMassCurveStyles::usage = "ClockMassCurveStyles[species] gives consistent W, Z, and top colors.";
ClockMassLegend::usage = "ClockMassLegend[species,styles,measure] labels the amplitude or phase curves in a framed in-plot legend.";
ClockShapeLabel::usage = "ClockShapeLabel[\"T\"|\"TOverE\"] labels the signed nonlocal trispectrum or its explicitly de-enveloped value.";
ClockShapePlotOptions::usage = "ClockShapePlotOptions[angles,parameters,showParameters:True] gives log-ratio axes and optional fixed benchmark labels for the shape maps.";
ClockShapeColorLegend::usage = "ClockShapeColorLegend[edges,measure] gives a signed linear colorbar with scientific-notation ticks.";
Begin["`Private`"];
$SignalPlotStyle=<|
  "Font"->"Latin Modern Roman","TickSize"->10,"LabelSize"->12,"LegendSize"->10,
  "LineThickness"->1.4,"FrameThickness"->0.8,
  "ScanImageSize"->400,"ScanAspectRatio"->1/2,"ScanPadding"->{{55,12},{42,10}},
  "MapImageSize"->470,"MapPadding"->{{54,10},{44,12}},
  "ClockMassImageSize"->560,"ClockMassAspectRatio"->3/5,
  "ClockMassPadding"->{{68,16},{52,12}},"ClockMassLegendPosition"->{0.80,0.82},
  "ClockSpeciesColors"-><|"W"->RGBColor[201/255,42/255,42/255],
    "Z"->RGBColor[36/255,85/255,195/255],"t"->Black|>,
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
inPlotLegendBox[content_] := Framed[content,Background->White,FrameStyle->Black,
  RoundingRadius->0,FrameMargins->5];
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
ClockMapPlotOptions[species_,v_,phi_] := Module[{ticks,base},
  ticks=Table[{j Pi/4,TraditionalForm[j Pi/4]},{j,0,4}];
  base=DeleteCases[MapPlotOptions[v,phi],
    HoldPattern[(FrameLabel|FrameTicks|Epilog|ImagePadding|Background)->_]];
  Join[base,{
    Background->White,
    FrameLabel->{Subscript[Style["\[Theta]",Italic],1],Subscript[Style["\[Theta]",Italic],3]},
    FrameTicks->{{ticks,None},{ticks,None}},
    ImagePadding->{{50,14},{52,14}},PlotLabel->None,
    Epilog->Inset[inPlotLegendBox[Style[Row[{Style[If[species=="t","top",species],Italic],
      ":  ",Style["\[Nu]",Italic]," = ",v,",  ",
      Subscript[Style["\[Phi]",Italic],3]," = ",TraditionalForm[phi]}],
      Black,Sequence@@fontStyle[$SignalPlotStyle["LabelSize"]]]],
      Scaled[$SignalPlotStyle["ParameterPosition"]]]}]
];
ClockMapLegendLabel[species_] := Subsuperscript[
  Style["t",Italic],Style["NL",Plain],Style[species,Italic]];
ClockMapColorLegend[edges_List,species_] := Module[{ticks},
  ticks=Rest[Most[edges]][[1;;-1;;2]];
  Placed[BarLegend[{mapPalette,{First[edges],Last[edges]}},Rest[Most[edges]],
    Method->{"Ticks"->({#,ScientificForm[#,4]}& /@ ticks)},
    LegendLabel->Placed[Style[ClockMapLegendLabel[species],Black,
      Sequence@@fontStyle[$SignalPlotStyle["LabelSize"]]],Above],
    LabelStyle->Directive[Black,Sequence@@fontStyle[$SignalPlotStyle["LegendSize"]]],
    LegendMarkerSize->$SignalPlotStyle["ColorBarSize"],LegendMargins->$SignalPlotStyle["ColorBarMargins"],
    LegendFunction->(Framed[#,Background->White,FrameStyle->None,FrameMargins->0]&)],Right]
];
ClockMassPlotOptions[angles_List] := Join[
  DeleteCases[ScanPlotOptions[],HoldPattern[(ImageSize|AspectRatio|ImagePadding)->_]],{
    ImageSize->$SignalPlotStyle["ClockMassImageSize"],
    AspectRatio->$SignalPlotStyle["ClockMassAspectRatio"],
    ImagePadding->$SignalPlotStyle["ClockMassPadding"],PlotLabel->None}];
ClockMassCurveStyles[species_List] :=
  (Directive[$SignalPlotStyle["ClockSpeciesColors"][#],
    AbsoluteThickness[$SignalPlotStyle["LineThickness"]],
    Dashing[{}]]& /@ species);
ClockMassLegend[species_List,styles_List,measure_String] := Placed[
  LineLegend[styles,If[measure=="Phase",
    (Subscript[Style["\[CurlyPhi]",Italic],Style[#,Italic]]& /@ species),
    ClockMapLegendLabel /@ species],
    LabelStyle->Join[{Black},fontStyle[$SignalPlotStyle["LabelSize"]]],
    LegendLayout->"Column",LegendMarkerSize->{22,8},LegendMargins->0,
    LegendFunction->(inPlotLegendBox[#]&)],$SignalPlotStyle["ClockMassLegendPosition"]];
calligraphicT[] := Style[FromCharacterCode[16^^1D4AF],FontFamily->"Latin Modern Math",FontSlant->Plain];
ClockShapeLabel[measure_] := If[measure=="TOverE",
  Row[{Subscript[calligraphicT[],Style["NS",Plain]],"/",Style["E",Italic]}],
  Subscript[calligraphicT[],Style["NS",Plain]]];
ClockShapePlotOptions[angles_List,parameters_Association,showParameters_:True] := Join[
  DeleteCases[MapPlotOptions[0,0],HoldPattern[(FrameLabel|Epilog|ImagePadding)->_]],{
    Background->White,ImagePadding->{{56,14},{55,If[TrueQ[showParameters],50,12]}},
    FrameLabel->Table[Row[{"-",Subscript[Style["log",Plain],10]," ",Subscript[Style["r",Italic],j]}],{j,1,2}],
    PlotLabel->If[TrueQ[showParameters],Style[Column[{
      Row[{Style["\[Lambda]",Italic]," = ",N[parameters["lambda"]],",  ",
        Subscript[Style["h",Italic],0],"/",Style["H",Italic]," = ",parameters["h0OverH"]}],
      Row[{Subscript[Style["\[Theta]",Italic],1]," = ",TraditionalForm[angles[[1]]],",  ",
        Subscript[Style["\[Theta]",Italic],3]," = ",TraditionalForm[angles[[2]]],",  ",
        Subscript[Style["\[Phi]",Italic],3]," = ",TraditionalForm[angles[[3]]]}]},Alignment->Center],
      Black,Sequence@@fontStyle[$SignalPlotStyle["LabelSize"]]],None]}];
ClockShapeColorLegend[edges_List,measure_String] := Module[{ticks},
  ticks=Subdivide[First[edges],Last[edges],10][[2;;-2;;2]];
  Placed[BarLegend[{mapPalette,{First[edges],Last[edges]}},Rest[Most[edges]],
    Method->{"Ticks"->({#,If[Abs[#]<10^-14 Max[Abs[edges]],0,ScientificForm[#,3]]}& /@ ticks)},
    LegendLabel->Placed[Style[ClockShapeLabel[measure],Black,
      Sequence@@fontStyle[$SignalPlotStyle["LabelSize"]]],Above],
    LabelStyle->Directive[Black,Sequence@@fontStyle[$SignalPlotStyle["LegendSize"]]],
    LegendMarkerSize->$SignalPlotStyle["ColorBarSize"],LegendMargins->$SignalPlotStyle["ColorBarMargins"],
    LegendFunction->(Framed[#,Background->White,FrameStyle->None,FrameMargins->0]&)],Right]
];
End[];
EndPackage[];
