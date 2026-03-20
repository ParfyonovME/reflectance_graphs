(* ::Package:: *)

ClearAll["Global`*"];

csvFile   = FileNameJoin[{NotebookDirectory[], "spline_grids_362-760_step0p1.csv"}];
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

(*\:0412\:043e\:0441\:0441\:0442\:0430\:043d\:043e\:0432\:043b\:0435\:043d\:0438\:0435 \:0438\:043d\:0442\:0435\:0440\:043f\:043e\:043b\:044f\:0446\:0438\:0439 \:043f\:043e \:0434\:0430\:043d\:043d\:044b\:043c \:0438\:0437 .csv*)
ledCols = AssociationMap[data[[All, colIndex[#]]] &, ledNamesAll];
qeCol   = data[[All, colIndex[qeName]]];
fLED = AssociationMap[
  Interpolation[Transpose[{wl, ledCols[#]}], InterpolationOrder -> 3] &,
  ledNamesAll
];
qeInt = Interpolation[Transpose[{wl, qeCol}], InterpolationOrder -> 3];
qeFun[x_?NumericQ] := qeInt[x];

(*\:0417\:0430\:0434\:0430\:043d\:0438\:0435 \:043f\:0435\:0440\:0435\:043c\:0435\:043d\:043d\:044b\:0445 \:0444\:0443\:043d\:043a\:0446\:0438\:0438 \:043e\:043f\:0442\:0438\:043c\:0438\:0437*)
varLED = Intersection[ledNamesAll, Keys[baseW]];
If[Length[varLED] == 0,
  Print["ERROR: No overlap between LED columns and baseline weights keys."]; Abort[];
];
baseline = AssociationThread[varLED -> Lookup[baseW, varLED]];

(*\:0421\:0435\:0442\:043a\:0430 \:0438 \:043c\:0430\:0442\:0440\:0438\:0446\:0430*)
fitGrid = Range[lambdaMin, lambdaMax, stepNm];
A = Transpose@Table[fLED[name] /@ fitGrid, {name, varLED}];
qeVec = qeInt /@ fitGrid;
nVar = Length[varLED];
vars = Array[x, nVar];
sigVec[wvec_?VectorQ] := (A . wvec) * qeVec;

(*\:041e\:043f\:0442\:0438\:043c\:0438\:0437\:0430\:0446\:0438\:044f Emax/Emin*)
epsFrac = 10^-6;  (* 1e-6 \:043e\:0442 \:0441\:0440\:0435\:0434\:043d\:0435\:0433\:043e *)
objVarVec[wvec_?VectorQ] := Module[{sig, mu},
  sig = sigVec[wvec];
  mu  = Mean[sig];
  If[mu <= 0, Return[10^9]];
  Total[(sig - mu)^2]/mu^2
];
obj := objVarVec[vars];

(*\:043e\:0433\:0440\:0430\:043d\:0438\:0447\:0435\:043d\:0438\:044f baseline+-dw*)
b0 = Lookup[baseline, varLED];
lb = Max[0., #] & /@ (b0 - dW);
ub = b0 + dW;
cons = Thread[lb <= vars <= ub];

(*\:041e\:043f\:0442\:0438\:043c\:0438\:0437\:0430\:0446\:0438\:044f*)
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

(*Emax/Emin baseline*)
sigBase = sigVec[Lookup[baseline, varLED]];
sigOpt  = sigVec[Lookup[weightsOpt, varLED]];
ratioBase = Max[sigBase]/Max[Min[sigBase], epsFrac*Mean[sigBase]];
ratioOpt  = Max[sigOpt]/Max[Min[sigOpt],  epsFrac*Mean[sigOpt]];
Print["Emax/Emin baseline  = ", ratioBase];
Print["Emax/Emin optimized = ", ratioOpt];

(*\:0413\:0440\:0430\:0444\:0438\:043a \:0434\:043e/\:043f\:043e\:0441\:043b\:0435*)
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
cvBase = stdBase/muBase;   (*\:041a\:043e\:044d\:0444\:0444\:0438\:0446\:0438\:0435\:043d\:0442 \:0432\:0430\:0440\:0438\:0430\:0446\:0438\:0438*)
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
      ListLinePlot[ptsBase, PlotStyle -> Directive[Red,  AbsoluteThickness[4]]],
      ListLinePlot[ptsOpt,  PlotStyle -> Directive[Blue, AbsoluteThickness[4]]],
      Frame -> True,
      Background -> White,
      GridLines -> None,
      PlotRange -> {{lambdaMin, lambdaMax}, {0, All}},
      FrameLabel -> {"Wavelength (nm)", "S(\[Lambda]) \[Times] QE (a.u.)"},
      BaseStyle -> {FontFamily -> "Times", FontSize -> 32},
      ImageSize -> 1100,
      AspectRatio -> 1,
      FrameStyle -> Directive[Black, AbsoluteThickness[1.2]],
      TicksStyle -> Directive[Black, 18],
      PlotRangePadding -> Scaled[0.01]
    ],
    Placed[
      LineLegend[
        {Directive[Red, AbsoluteThickness[6]], 
         Directive[Blue, AbsoluteThickness[6]]},
        {
          "Baseline (R=" <> ToString@NumberForm[ratioBase, {5, 3}] <> ")",
          "Optimized (R=" <> ToString@NumberForm[ratioOpt,  {5, 3}] <> ")"
        },
        LabelStyle -> {FontFamily -> "Times", FontSize -> 26},   
        LegendLayout -> {"Column", 1},
        LegendMarkerSize -> {45, 30}                      
      ],
      {Right, Center}
    ]
  ];

comparePlotJournal

