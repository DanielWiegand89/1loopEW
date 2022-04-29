(* ::Package:: *)

(* -------------------------------------------------------------  *)
 (* Tools for 1-loop calculations                                 *)
 (* version 1.0 (April 28, 2022)                                  *)
 (* ------------------------------------------------------------- *)
 (* FAtoFCreplace.m                                               *)
 (* Daniel Wiegand (daniel.wiegand@northwestern.edu)              *)
 (* last revision 04/28/22                                        *)
 (* ------------------------------------------------------------- *)
 (* Mathematica file for cleaning up the mess after tensor        *)
 (* reduction                                                     *)
 (* ------------------------------------------------------------- *)


 FAtoFCtrans[exp_] :=
  Module[{FAchange, TadpoleReplaced, BubbleReplaced, TriangleReplaced,
    BoxReplaced},
   FAchange =
    ReplaceAll[
     FCPrepareFAAmp[exp,
      UndoChiralSplittings -> True], {DiracSpinor :> Spinor,
      SumOver[__] :> 1, GS :> gs, EL :> el, D :> $D,
      Pair[LorentzIndex[a_], Momentum[Polarization[b_, c_]]] :> 1,
      SUNT[a_, b_, c_] :> SUNTF[a, b, c],
      SUNF[a_, b_, c_, d_] :>
       SUNF[a, b, Index[Gluon, 111]] SUNF[c, d, Index[Gluon, 111]]}];
   TadpoleReplaced =
    ReplaceAll[
     FAchange, {FeynAmpDenominator[PropagatorDenominator[q1, m1_]] :>
       Pi^2*I*A0[m1^2]}];
   BubbleReplaced =
    ReplaceAll[
     TadpoleReplaced, {FeynAmpDenominator[
        PropagatorDenominator[q1, m1_],
        PropagatorDenominator[q1 + a_, m2_]] :>
       Pi^2*I*Bfct[a, m1^2, m2^2],
      FeynAmpDenominator[PropagatorDenominator[q1 + a_, m2_],
        PropagatorDenominator[q1, m1_]] :>
       Pi^2*I*Bfct[a, m1^2, m2^2]}];
   TriangleReplaced =
    ReplaceAll[
     BubbleReplaced, {FeynAmpDenominator[
        PropagatorDenominator[q1, m1_],
        PropagatorDenominator[q1 + a_, m2_],
        PropagatorDenominator[q1 + a_ + b_, m3_]] :>
       Pi^2*I*Cfct[a, b, m1^2, m2^2, m3^2],
      FeynAmpDenominator[PropagatorDenominator[q1, m1_],
        PropagatorDenominator[q1 + a_ + b_, m3_],
        PropagatorDenominator[q1 + a_, m2_]] :>
       Pi^2*I*Cfct[a, b, m1^2, m2^2, m3^2],
      FeynAmpDenominator[PropagatorDenominator[q1 + a_, m2_],
        PropagatorDenominator[q1, m1_],
        PropagatorDenominator[q1 + a_ + b_, m3_]] :>
       Pi^2*I*Cfct[a, b, m1^2, m2^2, m3^2],
      FeynAmpDenominator[PropagatorDenominator[q1 + a_ + b_, m3_],
        PropagatorDenominator[q1, m1_],
        PropagatorDenominator[q1 + a_, m2_]] :>
       Pi^2*I*Cfct[a, b, m1^2, m2^2, m3^2],
      FeynAmpDenominator[PropagatorDenominator[q1 + a_, m2_],
        PropagatorDenominator[q1 + a_ + b_, m3_],
        PropagatorDenominator[q1, m1_]] :>
       Pi^2*I*Cfct[a, b, m1^2, m2^2, m3^2],
      FeynAmpDenominator[PropagatorDenominator[q1 + a_ + b_, m3_],
        PropagatorDenominator[q1 + a_, m2_],
        PropagatorDenominator[q1, m1_]] :>
       Pi^2*I*Cfct[a, b, m1^2, m2^2, m3^2]}];
   BoxReplaced =
    OneLoopSum[
     PropagatorDenominatorExplicit[
      ReplaceAll[
       TriangleReplaced, {FeynAmpDenominator[
          PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1 + a_ + b_, m3_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1 + a_ + b_, m3_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1 + a_, m2_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1 + a_, m2_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1 + a_ + b_, m3_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1 + a_ + b_, m3_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1 + a_, m2_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1 + a_, m2_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_ + b_, m3_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_ + b_, m3_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_, m2_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1, m1_],
          PropagatorDenominator[q1 + a_, m2_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1, m1_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1, m1_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1, m1_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1, m1_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1, m1_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2],
        FeynAmpDenominator[
          PropagatorDenominator[q1 + a_ + b_ + c_, m4_],
          PropagatorDenominator[q1 + a_ + b_, m3_],
          PropagatorDenominator[q1 + a_, m2_],
          PropagatorDenominator[q1, m1_]] :>
         I*Pi^2*Dfct[a, b, c, m1^2, m2^2, m3^2, m4^2]}]]];
   Do[BoxReplaced = ReplaceAll[BoxReplaced, {HoldForm[a_] :> a}], {k,
     1, 10}]; Return[BoxReplaced]]
