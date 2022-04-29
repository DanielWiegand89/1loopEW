(* ::Package:: *)

(* -------------------------------------------------------------  *)
 (* Tools for 1-loop calculations                                 *)
 (* version 1.0 (April 28, 2022)                                  *)
 (* ------------------------------------------------------------- *)
 (* FinalTrans.m                                                  *)
 (* Daniel Wiegand (daniel.wiegand@northwestern.edu)              *)
 (* last revision 04/28/22                                        *)
 (* ------------------------------------------------------------- *)
 (* Mathematica file for cleaning up the mess after tensor        *)
 (* reduction                                                     *)
 (* ------------------------------------------------------------- *)

 CleanUp[x_] :=
  ReplaceAll[
   ScalarProductExpand[
    ReplaceAll[
     x, {Bfct[a_, b_, c_] :>
       B0PV[Pair[Momentum[a, $D], Momentum[a, $D]], b, c],
      Cfct[p1_, p2_, m12_, m22_, m32_] :>
       C0PV[Pair[Momentum[p1, $D], Momentum[p1, $D]],
        Pair[Momentum[p2, $D], Momentum[p2, $D]],
        Pair[Momentum[p1 + p2, $D], Momentum[p1 + p2, $D]], m12, m22,
        m32], Dfct[p1_, p2_, p3_, m12_, m22_, m32_, m42_] :>
       D0PV[Pair[Momentum[p1, $D], Momentum[p1, $D]],
        Pair[Momentum[p2, $D], Momentum[p2, $D]],
        Pair[Momentum[p3, $D], Momentum[p3, $D]],
        Pair[Momentum[p1 + p2 + p3, $D], Momentum[p1 + p2 + p3, $D]],
        Pair[Momentum[p1 + p2, $D], Momentum[p1 + p2, $D]],
        Pair[Momentum[p2 + p3, $D], Momentum[p2 + p3, $D]], m12, m22,
        m32, m42]}]], {B0PV[0, M_, M_] :> ($D - 2)*A0[M]/(2*M) /;
      M =!= 0,
    B0PV[M2_, M2_, 0] :> ($D - 2)/($D - 3)*A0[M2]/(2*M2) /; M2 =!= 0,
    B0PV[M2_, 0, M2_] :> ($D - 2)/($D - 3)*A0[M2]/(2*M2) /; M2 =!= 0,
    B0PV[0, M1_,
      M2_] :> (A0[M1] - A0[M2])/(M1 - M2) /; (M1 - M2) =!= 0}]
