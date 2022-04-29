(* ::Package:: *)

(* -------------------------------------------------------------   *)
 (* Tools for 1-loop calculations                             *)
 (* version 1.0 (April 28, 2022)                                  *)
 (* ------------------------------------------------------------- *)
 (* DenCancel.m                                                   *)
 (* Daniel Wiegand (daniel.wiegand@northwestern.edu)              *)
 (* last revision 04/28/22                                        *)
 (* ------------------------------------------------------------- *)
 (* Mathematica file for canceling numerators and denonominators  *)
 (* of the 1loop functions                                        *)
 (* ------------------------------------------------------------- *)

 DList = {1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7,
    8, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 1, 2,
    3, 4, 5};

 CList = {1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7};

 BList = {1, 2, 3, 4};




 Dcancel[exp_, 1] :=
   ReplaceAll[
    exp, {Dfctqq[p1_, p2_, p3_, m12_, m22_, m32_, m42_] :>
      m12*Dfctq[p1, p2, p3, m12, m22, m32, m42] +
       Pair[Momentum[p1, $D], Momentum[p1, $D]]*
        Cfct[p2, p3, m22, m32,
         m42] + (m22*Cfct[p2, p3, m22, m32, m42] +
         Bfct[p3, m32, m42]) - 2*CfctA[p1, p2, p3, m22, m32, m42],
     DfctqB[a_, b_, p1_, p2_, p3_, m12_, m22_, m32_, m42_] :>
      m12*DfctB[a, b, p1, p2, p3, m12, m22, m32, m42] +
       CfctB[a, b, p2, p3, m22, m32, m42] -
       Pair[Momentum[p1, $D], Momentum[a, $D]]*
        CfctA[b, p2, p3, m22, m32, m42] -
       Pair[Momentum[p1, $D], Momentum[b, $D]]*
        CfctA[a, p2, p3, m22, m32, m42] +
       Pair[Momentum[p1, $D], Momentum[a, $D]]*
        Pair[Momentum[p1, $D], Momentum[b, $D]]*
        Cfct[p2, p3, m22, m32, m42],
     DfctqA[a_, p1_, p2_, p3_, m12_, m22_, m32_, m42_] :>
      m12*DfctA[a, p1, p2, p3, m12, m22, m32, m42] -
       Pair[Momentum[p1, $D], Momentum[a, $D]]*
        Cfct[p2, p3, m22, m32, m42] + CfctA[a, p2, p3, m22, m32, m42],
     Dfctq[p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> (Cfct[p2, p3, m22, m32, m42] +
        m12*Dfct[p1, p2, p3, m12, m22, m32, m42])}];

 Dcancel[exp_, 2] :=
   ReplaceAll[
    exp, {DfctD[p1_, a_, b_, c_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> ((m22 - m12 -
            Pair[Momentum[p1, $D], Momentum[p1, $D]])*
          DfctC[a, b, c, p1, p2, p3, m12, m22, m32, m42] +
         CfctC[a, b, c, p1 + p2, p3, m12, m32,
          m42] - (CfctC[a, b, c, p2, p3, m22, m32, m42] -
           Pair[Momentum[p1, $D], Momentum[a, $D]]*
            Pair[Momentum[p1, $D], Momentum[b, $D]]*
            Pair[Momentum[p1, $D], Momentum[c, $D]]*
            Cfct[p2, p3, m22, m32, m42] -
           Pair[Momentum[p1, $D], Momentum[c, $D]]*
            CfctB[a, b, p2, p3, m22, m32, m42] -
           Pair[Momentum[p1, $D], Momentum[a, $D]]*
            CfctB[b, c, p2, p3, m22, m32, m42] -
           
           Pair[Momentum[p1, $D], Momentum[b, $D]]*
            CfctB[a, c, p2, p3, m22, m32, m42] +
           Pair[Momentum[p1, $D], Momentum[a, $D]]*
            Pair[Momentum[p1, $D], Momentum[b, $D]]*
            CfctA[c, p2, p3, m22, m32, m42] +
           Pair[Momentum[p1, $D], Momentum[a, $D]]*
            Pair[Momentum[p1, $D], Momentum[c, $D]]*
            CfctA[b, p2, p3, m22, m32, m42] +
           Pair[Momentum[p1, $D], Momentum[b, $D]]*
            Pair[Momentum[p1, $D], Momentum[c, $D]]*
            CfctA[a, p2, p3, m22, m32, m42]))/2,
     DfctD[p2_, a_, b_, c_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> ((m32 - m22 -
            Pair[Momentum[p2, $D], Momentum[p2, $D]] -
            2*Pair[Momentum[p1, $D], Momentum[p2, $D]])*
          DfctC[a, b, c, p1, p2, p3, m12, m22, m32, m42] +
         CfctC[a, b, c, p1, p2 + p3, m12, m22, m42] -
         CfctC[a, b, c, p1 + p2, p3, m12, m32, m42])/2,
     DfctD[p3_, a_, b_, c_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> ((m42 - m32 -
            2*Pair[Momentum[p3, $D], Momentum[p1 + p2, $D]] -
            Pair[Momentum[p3, $D], Momentum[p3, $D]])*
          DfctC[a, b, c, p1, p2, p3, m12, m22, m32, m42] +
         CfctC[a, b, c, p1, p2, m12, m22, m32] -
         CfctC[a, b, c, p1, p2 + p3, m12, m22, m42])/2}];

 Dcancel[exp_, 3] :=
   ReplaceAll[
    exp, {DfctC[p1_, a_, b_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> ((m22 - m12 -
            Pair[Momentum[p1, $D], Momentum[p1, $D]])*
          DfctB[a, b, p1, p2, p3, m12, m22, m32, m42] +
         CfctB[a, b, p1 + p2, p3, m12, m32, m42] -
         Pair[Momentum[p1, $D], Momentum[a, $D]]*
          Pair[Momentum[p1, $D], Momentum[b, $D]]*
          Cfct[p2, p3, m22, m32, m42] -
         CfctB[a, b, p2, p3, m22, m32, m42] +
         Pair[Momentum[p1, $D], Momentum[a, $D]]*
          CfctA[b, p2, p3, m22, m32, m42] +
         Pair[Momentum[p1, $D], Momentum[b, $D]]*
          CfctA[a, p2, p3, m22, m32, m42])/2,
     DfctC[p2_, a_, b_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> ((m32 - m22 -
            Pair[Momentum[p2, $D], Momentum[p2, $D]] -
            2*Pair[Momentum[p1, $D], Momentum[p2, $D]])*
          DfctB[a, b, p1, p2, p3, m12, m22, m32, m42] -
         CfctB[a, b, p1 + p2, p3, m12, m32, m42] +
         CfctB[a, b, p1, p2 + p3, m12, m22, m42])/2,
     DfctC[p3_, a_, b_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> ((m42 - m32 -
            2*Pair[Momentum[p3, $D], Momentum[p1 + p2, $D]] -
            Pair[Momentum[p3, $D], Momentum[p3, $D]])*
          DfctB[a, b, p1, p2, p3, m12, m22, m32, m42] +
         CfctB[a, b, p1, p2, m12, m22, m32] -
         CfctB[a, b, p1, p2 + p3, m12, m22, m42])/2}];

 Dcancel[exp_, 4] :=
   ReplaceAll[
    exp, {DfctB[p1_, a_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> ((m22 - m12 -
            Pair[Momentum[p1, $D], Momentum[p1, $D]])*
          DfctA[a, p1, p2, p3, m12, m22, m32, m42] +
         CfctA[a, p1 + p2, p3, m12, m32, m42] +
         Pair[Momentum[p1, $D], Momentum[a, $D]]*
          Cfct[p2, p3, m22, m32, m42] -
         CfctA[a, p2, p3, m22, m32, m42])/2,
     DfctB[p2_, a_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> ((m32 - m22 -
            Pair[Momentum[p2, $D], Momentum[p2, $D]] -
            2*Pair[Momentum[p1, $D], Momentum[p2, $D]])*
          DfctA[a, p1, p2, p3, m12, m22, m32, m42] -
         CfctA[a, p1 + p2, p3, m12, m32, m42] +
         CfctA[a, p1, p2 + p3, m12, m22, m42])/2,
     DfctB[p3_, a_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> ((m42 - m32 -
            2*Pair[Momentum[p3, $D], Momentum[p1 + p2, $D]] -
            Pair[Momentum[p3, $D], Momentum[p3, $D]])*
          DfctA[a, p1, p2, p3, m12, m22, m32, m42] +
         CfctA[a, p1, p2, m12, m22, m32] -
         CfctA[a, p1, p2 + p3, m12, m22, m42])/2}];

 Dcancel[exp_, 5] :=
   ReplaceAll[
    exp, {DfctA[p1_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> (Cfct[p1 + p2, p3, m12, m32, m42] -
         Cfct[p2, p3, m22, m32,
          m42] + (m22 - m12 -
            Pair[Momentum[p1, $D], Momentum[p1, $D]])*
          Dfct[p1, p2, p3, m12, m22, m32, m42])/2,
     DfctA[p2_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> (Cfct[p1, p2 + p3, m12, m22, m42] -
         Cfct[p1 + p2, p3, m12, m32,
          m42] + (m32 - m22 -
            Pair[Momentum[p2, $D], Momentum[p2, $D]] -
            2*Pair[Momentum[p1, $D], Momentum[p2, $D]])*
          Dfct[p1, p2, p3, m12, m22, m32, m42])/2,
     DfctA[p3_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> (Cfct[p1, p2, m12, m22, m32] -
         Cfct[p1, p2 + p3, m12, m22,
          m42] + (m42 - m32 -
            Pair[Momentum[p3, $D], Momentum[p3, $D]] -
            2*Pair[Momentum[p3, $D], Momentum[p1 + p2, $D]])*
          Dfct[p1, p2, p3, m12, m22, m32, m42])/2}];

 Dcancel[exp_, 6] :=
   ReplaceAll[
    exp, {DfctD[a_, b_, c_, d_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> -DfctD[-a, b, c, d, p1, p2, p3, m12, m22, m32, m42],
     DfctC[a_, b_, c_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> -DfctC[-a, b, c, p1, p2, p3, m12, m22, m32, m42],
     DfctB[a_, b_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> -DfctB[-a, b, p1, p2, p3, m12, m22, m32, m42],
     DfctA[a_, p1_, p2_, p3_, m12_, m22_, m32_,
       m42_] :> -DfctA[-a, p1, p2, p3, m12, m22, m32, m42]}];

 Dcancel[exp_, 7] :=
   ReplaceAll[
    exp, {DfctD[k2, b_, c_, d_, p12_, p22_, p32_, m12_, m22_, m32_,
       m42_] :>
      DfctD[p1, b, c, d, p12, p22, p32, m12, m22, m32, m42] +
       DfctD[p2, b, c, d, p12, p22, p32, m12, m22, m32, m42] -
       DfctD[k1, b, c, d, p12, p22, p32, m12, m22, m32, m42],
     DfctC[k2, b_, c_, p12_, p22_, p32_, m12_, m22_, m32_, m42_] :>
      DfctC[p1, b, c, p12, p22, p32, m12, m22, m32, m42] +
       DfctC[p2, b, c, p12, p22, p32, m12, m22, m32, m42] -
       DfctC[k1, b, c, p12, p22, p32, m12, m22, m32, m42],
     DfctB[k2, b_, p12_, p22_, p32_, m12_, m22_, m32_, m42_] :>
      DfctB[p1, b, p12, p22, p32, m12, m22, m32, m42] +
       DfctB[p2, b, p12, p22, p32, m12, m22, m32, m42] -
       DfctB[k1, b, p12, p22, p32, m12, m22, m32, m42],
     DfctA[k2, p12_, p22_, p32_, m12_, m22_, m32_, m42_] :>
      DfctA[p1, p12, p22, p32, m12, m22, m32, m42] +
       DfctA[p2, p12, p22, p32, m12, m22, m32, m42] -
       DfctA[k1, p12, p22, p32, m12, m22, m32, m42],
     DfctA[k1, p12_, p22_, p32_, m12_, m22_, m32_, m42_] :>
      DfctA[p1, p12, p22, p32, m12, m22, m32, m42] +
       DfctA[p2, p12, p22, p32, m12, m22, m32, m42] -
       DfctA[k2, p12, p22, p32, m12, m22, m32, m42],
     DfctA[p1, p12_, p22_, p32_, m12_, m22_, m32_, m42_] :>
      DfctA[k1, p12, p22, p32, m12, m22, m32, m42] +
       DfctA[k2, p12, p22, p32, m12, m22, m32, m42] -
       DfctA[p2, p12, p22, p32, m12, m22, m32, m42]}];

 Dcancel[exp_, 8] :=
  ReplaceAll[
   exp, {DfctB[q1_, q2_, k1_, k2_, k3_, m12_, m22_, m32_, m42_] :>
     DfctB[q2, q1, k1, k2, k3, m12, m22, m32, m42]}]

 DRed[x_] :=
  Module[{test1}, test1 = x;
   Do[test1 = Dcancel[test1, DList[[k]]], {k, 1, Length[DList]}];
   Return[test1]]








 Ccancel[exp_, 1] :=
   ReplaceAll[
    exp, {CfctqA[a_, p1_, p2_, m12_, m22_,
       m32_] :> (m12*CfctA[a, p1, p2, m12, m22, m32] -
        Pair[Momentum[p1, $D], Momentum[a, $D]]*Bfct[p2, m22, m32] +
        BfctA[a, p2, m22, m32]),
     Cfctq[p1_, p2_, m12_, m22_,
       m32_] :> (Bfct[p2, m22, m32] +
        m12*Cfct[p1, p2, m12, m22, m32])}];

 Ccancel[exp_, 2] :=
   ReplaceAll[
    exp, {CfctC[p1_, a_, b_, p1_, p2_, m12_, m22_,
       m32_] :> ((m22 - m12 -
            Pair[Momentum[p1, $D], Momentum[p1, $D]])*
          CfctB[a, b, p1, p2, m12, m22, m32] +
         BfctB[a, b, p1 + p2, m12, m32] -
         Pair[Momentum[p1, $D], Momentum[a, $D]]*
          Pair[Momentum[p1, $D], Momentum[b, $D]]*Bfct[p2, m22, m32] -
         BfctB[a, b, p2, m22, m32] +
         Pair[Momentum[p1, $D], Momentum[a, $D]]*
          BfctA[b, p2, m22, m32] +
         Pair[Momentum[p1, $D], Momentum[b, $D]]*
          BfctA[a, p2, m22, m32])/2,
     CfctC[p2_, a_, b_, p1_, p2_, m12_, m22_,
       m32_] :> ((m32 - m22 -
            Pair[Momentum[p2, $D], Momentum[p2, $D]] -
            2*Pair[Momentum[p1, $D], Momentum[p2, $D]])*
          CfctB[a, b, p1, p2, m12, m22, m32] -
         BfctB[a, b, p1 + p2, m12, m32] + BfctB[a, b, p1, m12, m22])/
       2}];

 Ccancel[exp_, 3] :=
   ReplaceAll[
    exp, {CfctC[a_, p1_, b_, p1_, p2_, m12_, m22_,
       m32_] :> ((m22 - m12 -
            Pair[Momentum[p1, $D], Momentum[p1, $D]])*
          CfctB[a, b, p1, p2, m12, m22, m32] +
         BfctB[a, b, p1 + p2, m12, m32] -
         Pair[Momentum[p1, $D], Momentum[a, $D]]*
          Pair[Momentum[p1, $D], Momentum[b, $D]]*Bfct[p2, m22, m32] -
         BfctB[a, b, p2, m22, m32] +
         Pair[Momentum[p1, $D], Momentum[a, $D]]*
          BfctA[b, p2, m22, m32] +
         Pair[Momentum[p1, $D], Momentum[b, $D]]*
          BfctA[a, p2, m22, m32])/2,
     CfctC[a_, p2_, b_, p1_, p2_, m12_, m22_,
       m32_] :> ((m32 - m22 -
            Pair[Momentum[p2, $D], Momentum[p2, $D]] -
            2*Pair[Momentum[p1, $D], Momentum[p2, $D]])*
          CfctB[a, b, p1, p2, m12, m22, m32] -
         BfctB[a, b, p1 + p2, m12, m32] + BfctB[a, b, p1, m12, m22])/
       2}];

 Ccancel[exp_, 4] :=
   ReplaceAll[
    exp, {CfctC[a_, b_, p1_, p1_, p2_, m12_, m22_,
       m32_] :> ((m22 - m12 -
            Pair[Momentum[p1, $D], Momentum[p1, $D]])*
          CfctB[a, b, p1, p2, m12, m22, m32] +
         BfctB[a, b, p1 + p2, m12, m32] -
         Pair[Momentum[p1, $D], Momentum[a, $D]]*
          Pair[Momentum[p1, $D], Momentum[b, $D]]*Bfct[p2, m22, m32] -
         BfctB[a, b, p2, m22, m32] +
         Pair[Momentum[p1, $D], Momentum[a, $D]]*
          BfctA[b, p2, m22, m32] +
         Pair[Momentum[p1, $D], Momentum[b, $D]]*
          BfctA[a, p2, m22, m32])/2,
     CfctC[a_, b_, p2_, p1_, p2_, m12_, m22_,
       m32_] :> ((m32 - m22 -
            Pair[Momentum[p2, $D], Momentum[p2, $D]] -
            2*Pair[Momentum[p1, $D], Momentum[p2, $D]])*
          CfctB[a, b, p1, p2, m12, m22, m32] -
         BfctB[a, b, p1 + p2, m12, m32] + BfctB[a, b, p1, m12, m22])/
       2}];

 Ccancel[exp_, 5] :=
   ReplaceAll[
    exp, {CfctB[p1_, a_, p1_, p2_, m12_, m22_,
       m32_] :> ((m22 - m12 -
            Pair[Momentum[p1, $D], Momentum[p1, $D]])*
          CfctA[a, p1, p2, m12, m22, m32] +
         BfctA[a, p1 + p2, m12, m32] +
         Pair[Momentum[p1, $D], Momentum[a, $D]]*Bfct[p2, m22, m32] -
         BfctA[a, p2, m22, m32])/2,
     CfctB[p2_, a_, p1_, p2_, m12_, m22_,
       m32_] :> ((m32 - m22 -
            Pair[Momentum[p2, $D], Momentum[p2, $D]] -
            2*Pair[Momentum[p1, $D], Momentum[p2, $D]])*
          CfctA[a, p1, p2, m12, m22, m32] -
         BfctA[a, p1 + p2, m12, m32] + BfctA[a, p1, m12, m22])/2}];

 Ccancel[exp_, 6] :=
   ReplaceAll[
    exp, {CfctB[a_, p1_, p1_, p2_, m12_, m22_,
       m32_] :> ((m22 - m12 -
            Pair[Momentum[p1, $D], Momentum[p1, $D]])*
          CfctA[a, p1, p2, m12, m22, m32] +
         BfctA[a, p1 + p2, m12, m32] +
         Pair[Momentum[p1, $D], Momentum[a, $D]]*Bfct[p2, m22, m32] -
         BfctA[a, p2, m22, m32])/2,
     CfctB[a_, p2_, p1_, p2_, m12_, m22_,
       m32_] :> ((m32 - m22 -
            Pair[Momentum[p2, $D], Momentum[p2, $D]] -
            2*Pair[Momentum[p1, $D], Momentum[p2, $D]])*
          CfctA[a, p1, p2, m12, m22, m32] -
         BfctA[a, p1 + p2, m12, m32] + BfctA[a, p1, m12, m22])/2}];

 Ccancel[exp_, 7] :=
   ReplaceAll[
    exp, {CfctA[p1_, p1_, p2_, m12_, m22_,
       m32_] :> (Bfct[p1 + p2, m12, m32] -
         Bfct[p2, m22,
          m32] + (m22 - m12 -
            Pair[Momentum[p1, $D], Momentum[p1, $D]])*
          Cfct[p1, p2, m12, m22, m32])/2,
     CfctA[p2_, p1_, p2_, m12_, m22_,
       m32_] :> (Bfct[p1, m12, m22] -
         Bfct[p1 + p2, m12,
          m32] + (m32 - m22 -
            Pair[Momentum[p2, $D], Momentum[p2, $D]] -
            2*Pair[Momentum[p1, $D], Momentum[p2, $D]])*
          Cfct[p1, p2, m12, m22, m32])/2}];

 Ccancel[exp_, 8] :=
   ReplaceAll[
    exp, {CfctA[a_, p1_, p2_, m12_, m22_,
       m32_] :> -CfctA[-a, p1, p2, m12, m22, m32],
     CfctB[a_, b_, p1_, p2_, m12_, m22_, m32_] :>
      CfctB[-a, -b, p1, p2, m12, m22, m32],
     CfctC[a_, b_, c_, p1_, p2_, m12_, m22_,
       m32_] :> -CfctC[-a, -b, -c, p1, p2, m12, m22, m32]}];

 CRed[x_] :=
  Module[{test1}, test1 = x;
   Do[test1 = Ccancel[test1, CList[[k]]], {k, 1, Length[CList]}];
   Return[test1]]




 Bcancel[exp_, 1] :=
   ReplaceAll[exp,
    Bfctq[p1_, m12_, m22_] :> A0[m22] + m12*Bfct[p1, m12, m22]];

 Bcancel[exp_, 2] :=
   ReplaceAll[exp,
    BfctB[p1_, a_, p1_, m12_,
      m22_] :> ((m22 - m12 - Pair[Momentum[p1, $D], Momentum[p1, $D]])*
         BfctA[a, p1, m12, m22] +
        Pair[Momentum[p1, $D], Momentum[a, $D]]*A0[m22])/2];

 Bcancel[exp_, 3] :=
   ReplaceAll[exp,
    BfctB[a_, p1_, p1_, m12_,
      m22_] :> ((m22 - m12 - Pair[Momentum[p1, $D], Momentum[p1, $D]])*
         BfctA[a, p1, m12, m22] +
        Pair[Momentum[p1, $D], Momentum[a, $D]]*A0[m22])/2];

 Bcancel[exp_, 4] :=
   ReplaceAll[exp,
    BfctA[p1_, p1_, m12_,
      m22_] :> (A0[m12] -
        A0[m22] + (m22 - m12 -
           Pair[Momentum[p1, $D], Momentum[p1, $D]])*
         Bfct[p1, m12, m22])/2];

 BRed[x_] :=
  Module[{test1}, test1 = x;
   Do[test1 = Bcancel[test1, BList[[k]]], {k, 1, Length[BList]}];
   Return[test1]]




 AllCancel[exp_] := BRed[CRed[DRed[exp]]]
