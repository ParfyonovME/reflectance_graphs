(* ::Package:: *)

ClearAll["Global`*"];

csvFile   = FileNameJoin[{NotebookDirectory[], "spline_grids_350-780_step0p1.csv"}];
wlName    = "Wavelength_nm";
qeName    = "QE";
lambdaMin = 362;
lambdaMax = 760;
stepNm    = 5;
dW        = 1.5;

baseW = <|
  "LED365"     -> 1,
  "LED385wide" -> 1,
  "LED430"     -> 1,
  "LED480"     -> 1,
  "LED535"     -> 1,
  "LED568wide" -> 1,
  "LED600"     -> 1,
  "LED635"     -> 1,
  "LED680"     -> 1,
  "LED700"     -> 1,
  "LED740"     -> 1,
  "LED765"     -> 1
|>;


raw = Import[csvFile, "CSV"];
If[raw === $Failed, Print["ERROR: cannot import CSV: ", csvFile]; Abort[]];
headerRaw = First[raw];
data = ToExpression /@ Rest[raw];
header = ToString /@ headerRaw;
colIndex = AssociationThread[header -> Range[Length[header]]];
If[!KeyExistsQ[colIndex, wlName],
  Print["ERROR: no column '", wlName, "'. Available: ", header]; Abort[];
];
If[!KeyExistsQ[colIndex, qeName],
  Print["ERROR: no column '", qeName, "'. Available: ", header]; Abort[];
];
wl = data[[All, colIndex[wlName]]];
ledNamesAll = DeleteCases[header, wlName | qeName];


ledCols = AssociationMap[data[[All, colIndex[#]]] &, ledNamesAll];
qeCol   = data[[All, colIndex[qeName]]];
fLED = AssociationMap[
  Interpolation[Transpose[{wl, ledCols[#]}], InterpolationOrder -> 3] &,
  ledNamesAll
];
qeInt = Interpolation[Transpose[{wl, qeCol}], InterpolationOrder -> 3];
qeFun[x_?NumericQ] := qeInt[x];


varLED = Intersection[ledNamesAll, Keys[baseW]];
If[Length[varLED] == 0,
  Print["ERROR: No overlap between LED columns and baseline weights keys."]; Abort[];
];
baseline = AssociationThread[varLED -> Lookup[baseW, varLED]];


fitGrid = Range[lambdaMin, lambdaMax, stepNm];
A = Transpose@Table[fLED[name] /@ fitGrid, {name, varLED}];
qeVec = qeInt /@ fitGrid;
nVar = Length[varLED];
vars = Array[x, nVar];
sigVec[wvec_?VectorQ] := (A . wvec) * qeVec;


epsFrac = 10^-6;  (* 1e-6 \:043e\:0442 \:0441\:0440\:0435\:0434\:043d\:0435\:0433\:043e (\:043f\:043e\:043b)*)
objVarVec[wvec_?VectorQ] := Module[{sig, mu},
  sig = sigVec[wvec];
  mu  = Mean[sig];
  If[mu <= 0, Return[10^9]];
  Total[(sig - mu)^2]/mu^2
];
obj := objVarVec[vars];


b0 = Lookup[baseline, varLED];
lb = Max[0., #] & /@ (b0 - dW);
ub = b0 + dW;
cons = Thread[lb <= vars <= ub];


sol = NMinimize[
  {obj, cons},
  vars,
  Method -> "DifferentialEvolution",
  MaxIterations -> 5000
];
bestVars = vars /. sol[[2]];
weightsOpt = AssociationThread[varLED -> bestVars];
Print["Emax/Emin (optimized) = ", sol[[1]]];
Print["weightsOpt = ", weightsOpt];


sigBase = sigVec[Lookup[baseline, varLED]];
sigOpt  = sigVec[Lookup[weightsOpt, varLED]];
ratioBase = Max[sigBase]/Max[Min[sigBase], epsFrac*Mean[sigBase]];
ratioOpt  = Max[sigOpt]/Max[Min[sigOpt],  epsFrac*Mean[sigOpt]];
Print["Emax/Emin baseline  = ", ratioBase];
Print["Emax/Emin optimized = ", ratioOpt];


Sdet[lam_?NumericQ, w_Association] :=
  (Total@Table[w[name]*fLED[name][lam], {name, Keys[w]}]) * qeFun[lam];
plotGrid = Range[lambdaMin, lambdaMax, 1];
ptsBase = Table[{lam, Sdet[lam, baseline]},   {lam, plotGrid}];
ptsOpt  = Table[{lam, Sdet[lam, weightsOpt]}, {lam, plotGrid}];

(*variance / CV / RMS deviation*)
muBase = Mean[sigBase]; 
muOpt  = Mean[sigOpt];
stdBase = StandardDeviation[sigBase];
stdOpt  = StandardDeviation[sigOpt];
cvBase = stdBase/muBase;   
cvOpt  = stdOpt/muOpt;
rmsRelBase = Sqrt[Mean[(sigBase/muBase - 1)^2]];
rmsRelOpt  = Sqrt[Mean[(sigOpt/muOpt - 1)^2]];

Print["Std baseline  = ", stdBase, " ; Std optimized  = ", stdOpt];
Print["CV baseline   = ", cvBase,  " (", 100*cvBase, "%)"];
Print["CV optimized  = ", cvOpt,   " (", 100*cvOpt, "%)"];
Print["RMS rel dev baseline  = ", rmsRelBase, " (", 100*rmsRelBase, "%)"];
Print["RMS rel dev optimized = ", rmsRelOpt,  " (", 100*rmsRelOpt,  "%)"];
peak2peakBase = (Max[sigBase]-Min[sigBase])/muBase;
peak2peakOpt  = (Max[sigOpt]-Min[sigOpt])/muOpt;

Print["(Emax-Emin)/mean baseline  = ", peak2peakBase, " (", 100*peak2peakBase, "%)"];
Print["(Emax-Emin)/mean optimized = ", peak2peakOpt,  " (", 100*peak2peakOpt,  "%)"];


baselineLabel  = Row[{Style["S", Italic], "\[Times]", Style["QE", Italic], "  (Baseline)"}];
optimizedLabel = Row[{Style["S", Italic], "\[Times]", Style["QE", Italic], "  (Optimized)"}];

comparePlotJournal =
  Legended[
    Show[
      ListLinePlot[ptsBase, PlotStyle -> Directive[Red,  AbsoluteThickness[4]]],  (* baseline = red *)
      ListLinePlot[ptsOpt,  PlotStyle -> Directive[Blue, AbsoluteThickness[4]]], (* optimized = blue *)
      Frame -> True,
      Background -> White,
      GridLines -> None,
      PlotRange -> {{lambdaMin, lambdaMax}, {0, All}},
      FrameLabel -> {"Wavelength (nm)", "S(\[Lambda]) \[Times] QE (a.u.)"},
      BaseStyle -> {FontFamily -> "Times", FontSize -> 22},
      ImageSize -> 1100,
      AspectRatio -> 1,
      FrameStyle -> Directive[Black, AbsoluteThickness[1.2]],
      TicksStyle -> Directive[Black, 18],
      PlotRangePadding -> Scaled[0.01]
    ],
    Placed[
      LineLegend[
        {Directive[Red, AbsoluteThickness[4]], Directive[Blue, AbsoluteThickness[4]]},
        {
          "Baseline (R=" <> ToString@NumberForm[ratioBase, {5, 3}] <> ")",
          "Optimized (R=" <> ToString@NumberForm[ratioOpt,  {5, 3}] <> ")"
        },
        LabelStyle -> {FontFamily -> "Times", FontSize -> 18},
        LegendLayout -> {"Column", 1},
        LegendMarkerSize -> {18, 12}
      ],
      {Right, Center}
    ]
  ];
