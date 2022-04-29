(* ::Package:: *)

(* -------------------------------------------------------------   *)
 (* Tools for 1-loop calculations                             *)
 (* version 1.0 (April 28, 2022)                                  *)
 (* ------------------------------------------------------------- *)
 (* FctRename.m                                                   *)
 (* Daniel Wiegand (daniel.wiegand@northwestern.edu)              *)
 (* last revision 04/28/22                                        *)
 (* ------------------------------------------------------------- *)
 (* Mathematica file for translating the 1loop functions          *)
 (* into the appropriate format                                   *)
 (* ------------------------------------------------------------- *)

 MyRename[exp_, 1] :=
   Expand[ReplaceAll[
     ScalarProductExpand[
      exp], {Pair[Momentum[q1, $D], Momentum[q1, $D]]^2 :> q[2],
      Pair[Momentum[q1, $D], Momentum[q1, $D]] :> q[1],
      Pair[Momentum[q1, $D], Momentum[a_, $D]] :> Q[a],
      Pair[Momentum[a_, $D], Momentum[q1, $D]] :> Q[a]}]];

 MyRename[exp_, 2] :=
   Expand[ReplaceAll[
     exp, {Pair[Momentum[q2, $D], Momentum[q2, $D]]^2 :> Zq[2],
      Pair[Momentum[q2, $D], Momentum[q2, $D]] :> Zq[1],
      Pair[Momentum[q2, $D], Momentum[a_, $D]] :> ZQ[a],
      Pair[Momentum[a_, $D], Momentum[q2, $D]] :> ZQ[a]}]];

 MyRename[exp_, 3] :=
   Expand[ReplaceAll[
     exp, {Q[a_]^4 :> Q4[a], Q[a_]^3 :> Q3[a], Q[a_]^2 :> Q2[a],
      q[1]^2 :> q[2]}]];

 MyRename[exp_, 4] :=
   Expand[ReplaceAll[
     exp, {ZQ[a_]^4 :> ZQ4[a], ZQ[a_]^3 :> ZQ3[a], ZQ[a_]^2 :> ZQ2[a],
      Zq[1]^2 :> Zq[2]}]];

 MyRename[exp_, 5] :=
   ReplaceAll[
    exp, {Bfct[a_, m12_, m22_]*q[1] :> Bfctq[a, m12, m22],
     Bfct[a_, m12_, m22_]*Q[b_]*Q[c_] :> BfctB[b, c, a, m12, m22],
     Bfct[a_, m12_, m22_]*Q[b_] :> BfctA[b, a, m12, m22],
     Bfct[a_, m12_, m22_]*Q2[b_] :> BfctB[b, b, a, m12, m22]}];

 MyRename[exp_, 6] :=
   ReplaceAll[
    exp, {Cfct[a_, b_, m12_, m22_, m32_]*Q[p1_]*Q[p2_]*Q[p3_] :>
      CfctC[p1, p2, p3, a, b, m12, m22, m32],
     Cfct[a_, b_, m12_, m22_, m32_]*Q[p1_]*Q2[p2_] :>
      CfctC[p1, p2, p2, a, b, m12, m22, m32],
     Cfct[a_, b_, m12_, m22_, m32_]*Q3[p_] :>
      CfctC[p, p, p, a, b, m12, m22, m32],
     Cfct[a_, b_, m12_, m22_, m32_]*q[1]*Q[p_] :>
      CfctqA[p, a, b, m12, m22, m32],
     Cfct[a_, b_, m12_, m22_, m32_]*q[1] :> Cfctq[a, b, m12, m22, m32],
      Cfct[a_, b_, m12_, m22_, m32_]*Q[p1_]*Q[p2_] :>
      CfctB[p1, p2, a, b, m12, m22, m32],
     Cfct[a_, b_, m12_, m22_, m32_]*Q2[p_] :>
      CfctB[p, p, a, b, m12, m22, m32],
     Cfct[a_, b_, m12_, m22_, m32_]*Q[p_] :>
      CfctA[p, a, b, m12, m22, m32]}];

 MyRename[exp_, 7] :=
   ReplaceAll[
    exp, {Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*q[2] :>
      Dfctqq[a, b, c, m12, m22, m32, m42],
     Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*q[1]*Q[d_]*Q[e_] :>
      DfctqB[d, e, a, b, c, m12, m22, m32, m42],
     Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*q[1]*Q2[d_] :>
      DfctqB[d, d, a, b, c, m12, m22, m32, m42],
     Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*q[1]*Q[d_] :>
      DfctqA[d, a, b, c, m12, m22, m32, m42],
     Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*q[1] :>
      Dfctq[a, b, c, m12, m22, m32, m42],
     Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*Q[r_]*Q[p_]*Q[k_]*
       Q[l_] :> DfctD[r, p, k, l, a, b, c, m12, m22, m32, m42],
     Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*Q[k_]*Q[l_]*Q2[r_] :>
      DfctD[r, r, k, l, a, b, c, m12, m22, m32, m42],
     Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*Q2[p_]*Q2[k_] :>
      DfctD[p, p, k, k, a, b, c, m12, m22, m32, m42],
     Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*Q[p_]*Q3[k_] :>
      DfctD[p, k, k, k, a, b, c, m12, m22, m32, m42],
     Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*Q4[k_] :>
      DfctD[k, k, k, k, a, b, c, m12, m22, m32, m42],
     Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*Q[p_]*Q[k_]*Q[l_] :>
      DfctC[p, k, l, a, b, c, m12, m22, m32, m42],
     Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*Q[p_]*Q2[k_] :>
      DfctC[p, k, k, a, b, c, m12, m22, m32, m42],
     Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*Q3[p_] :>
      DfctC[p, p, p, a, b, c, m12, m22, m32, m42],
     Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*Q[p_]*Q[k_] :>
      DfctB[p, k, a, b, c, m12, m22, m32, m42],
     Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*Q2[p_] :>
      DfctB[p, p, a, b, c, m12, m22, m32, m42],
     Dfct[a_, b_, c_, m12_, m22_, m32_, m42_]*Q[p_] :>
      DfctA[p, a, b, c, m12, m22, m32, m42]}];

 RenameMomenta[x_] :=
  Module[{test1}, test1 = x; Do[test1 = MyRename[test1, k], {k, 1, 7}];
    Return[test1]]
