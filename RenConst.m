(* ::Package:: *)

(* -------------------------------------------------------------  *)
 (* Tools for 1-loop calculations                                 *)
 (* version 1.0 (April 28, 2022)                                  *)
 (* ------------------------------------------------------------- *)
 (* RenConst.m                                                   *)
 (* Daniel Wiegand (daniel.wiegand@northwestern.edu)              *)
 (* last revision 04/28/22                                        *)
 (* ------------------------------------------------------------- *)
 (* Mathematica file for renormalization constants in the         *)
 (* conventions of Denner/FeynCalc   (on-Shell scheme)            *)
 (* ------------------------------------------------------------- *)


 
 RenConst[exp_] :=
   ReplaceAll[
    ReplaceAll[
     exp, {dSW1 :> CW^2/(2*SW) (dMZsq1/MZ^2 - dMWsq1/MW^2),
      dZe1 :> -1/2*dZAA1 - SW/(2*CW) dZZA1}], {dMZsq1 :>
      1/(48 CW^4 MZ^2 SW^2 (\[Pi] - \[Pi] $D))
        aem (3 CW^2 (MH^2 - MZ^2 (-1 + $D)) A0[MH^2] +
         6 CW^2 MZ^2 (2 CW^2 SW^2 (-2 + $D) - SW^4 (-2 + $D) +
            CW^4 (-14 + 15 $D - 4 $D^2)) A0[MW^2] -
         3 CW^2 (MH^2 + MZ^2 (-3 + $D)) A0[MZ^2] -
         3 (CW^2 (MH^4 - 4 MH^2 MZ^2) + 4 MW^2 MZ^2 (-1 + $D)) B0PV[
           MZ^2, MH^2, MZ^2] -
         CW^2 MZ^2 Nf (-2 (9 - 12 SW^2 + 8 SW^4) (-2 + $D) A0[MB^2] -
            2 (9 - 24 SW^2 + 32 SW^4) (-2 + $D) A0[MT^2] +
            2 MZ^2 (27 - 54 SW^2 + 76 SW^4) (-2 + $D) B0PV[MZ^2, 0,
              0] + (2 MB^2 (27 - 24 SW^2 + 16 SW^4 - 9 $D) +
               MZ^2 (9 - 12 SW^2 + 8 SW^4) (-2 + $D)) B0PV[MZ^2, MB^2,
              MB^2] + (2 MT^2 (27 - 48 SW^2 + 64 SW^4 - 9 $D) +
               MZ^2 (9 - 24 SW^2 + 32 SW^4) (-2 + $D)) B0PV[MZ^2, MT^2,
               MT^2]) -
         3 CW^2 MZ^2 (MZ^2 (-2 CW^2 SW^2 + SW^4 + CW^4 (9 - 12 $D)) +
            4 MW^2 (2 CW^2 SW^2 + CW^4 (9 - 6 $D) +
               SW^4 (-3 + 2 $D))) B0PV[MZ^2, MW^2, MW^2]),
     dZZZ1 :>
      1/(192 CW^4 \[Pi] SW^2 (-1 + $D))
        aem (1/(-4 MB^2 MZ^2 + MZ^4)
           4 CW^2 Nf (2 MB^2 (27 - 24 SW^2 + 16 SW^4 - 9 $D) +
             MZ^2 (9 - 12 SW^2 + 8 SW^4) (-2 + $D)) (-2 + $D) A0[
            MB^2] + (
         6 (CW^2 (-3 MH^6 + 12 MH^4 MZ^2 - 14 MH^2 MZ^4 + 2 MZ^6) -
            4 MW^2 MZ^2 (MH^2 - 2 MZ^2) (-1 + $D)) A0[MH^2])/(
         MZ^4 (MH^2 - MZ^2)^2) +
         1/(-4 MT^2 MZ^2 + MZ^4)
           4 CW^2 Nf (2 MT^2 (27 - 48 SW^2 + 64 SW^4 - 9 $D) +
             MZ^2 (9 - 24 SW^2 + 32 SW^4) (-2 + $D)) (-2 + $D) A0[
            MT^2] -
         1/(-4 MW^2 MZ^2 + MZ^4)
           12 CW^2 (-2 + $D) (2 CW^2 (-4 MW^2 + MZ^2) SW^2 -
            SW^4 (MZ^2 + 4 MW^2 (-3 + 2 $D)) +
            3 CW^4 (4 MW^2 (-3 + 2 $D) + MZ^2 (-3 + 4 $D))) A0[MW^2] -
         1/(MZ^4 (MH^2 - MZ^2)^2)
           3 (CW^2 (4 MZ^6 + MH^4 MZ^2 (34 - 5 $D) +
               4 MH^2 MZ^4 (-9 + $D) + MH^6 (-8 + $D)) +
            4 MW^2 MZ^2 (-MZ^2 (-6 + $D) +
               MH^2 (-4 + $D)) (-1 + $D)) A0[MZ^2] -
         4 CW^2 Nf (27 - 54 SW^2 + 76 SW^4) (-2 + $D)^2 B0PV[MZ^2, 0,
           0] -
         1/(-4 MB^2 MZ^2 + MZ^4)
           2 CW^2 Nf (8 MB^4 (27 - 24 SW^2 + 16 SW^4 - 9 $D) +
            MZ^4 (9 - 12 SW^2 + 8 SW^4) (-2 + $D)^2 -
            2 MB^2 MZ^2 (-48 SW^2 + 32 SW^4 +
               9 (8 - 5 $D + $D^2))) B0PV[MZ^2, MB^2, MB^2] + (
         6 (CW^2 (3 MH^4 - 8 MH^2 MZ^2) + 4 MW^2 MZ^2 (-1 + $D)) B0PV[
           MZ^2, MH^2, MZ^2])/MZ^4 -
         1/(-4 MT^2 MZ^2 + MZ^4)
           2 CW^2 Nf (8 MT^4 (27 - 48 SW^2 + 64 SW^4 - 9 $D) +
            MZ^4 (9 - 24 SW^2 + 32 SW^4) (-2 + $D)^2 -
            2 MT^2 MZ^2 (-96 SW^2 + 128 SW^4 +
               9 (8 - 5 $D + $D^2))) B0PV[MZ^2, MT^2, MT^2] +
         1/(-4 MW^2 MZ^2 + MZ^4)
           6 CW^2 (-2 CW^2 (4 MW^2 - MZ^2) SW^2 (4 MW^2 +
               MZ^2 (-2 + $D)) -
            SW^4 (MZ^4 (-2 + $D) + 16 MW^4 (-3 + 2 $D) +
               4 MW^2 MZ^2 (11 - 11 $D + 2 $D^2)) +
            3 CW^4 (16 MW^4 (-3 + 2 $D) +
               4 MW^2 MZ^2 (15 - 15 $D + 2 $D^2) +
               MZ^4 (6 - 11 $D + 4 $D^2))) B0PV[MZ^2, MW^2, MW^2]),
     dMWsq1 :>
      1/(16 CW^2 MW^2 SW^2 (\[Pi] - \[Pi] $D))
        aem (CW^2 (MH^2 - MW^2 (-1 + $D)) A0[MH^2] +
         CW^2 (-MH^2 + MZ^2 (-1 - 4 CW^2 (-2 + $D)) +
            
            2 MW^2 (1 + 3 $D - 2 $D^2 + 2 CW^2 (-5 + 3 $D) +
               2 SW^2 (-5 + 3 $D))) A0[MW^2] -
         CW^2 (1 + 4 CW^2 (-2 + $D)) (-MZ^2 + MW^2 (-1 + $D)) A0[
           MZ^2] + 16 CW^2 MW^4 SW^2 (-1 + $D) ((($D - 2)/($D - 3)*
             A0[MW^2]/(2 MW^2))) -
         6 CW^2 Nf ((MB^2 - MT^2 - MW^2 (-2 + $D)) A0[
              MB^2] - (MB^2 - MT^2 + MW^2 (-2 + $D)) A0[MT^2] +
            3 MW^4 (-2 + $D) B0PV[MW^2, 0,
              0] - ((MB^2 - MT^2)^2 + (MB^2 + MT^2) MW^2 (-3 + $D) -
               MW^4 (-2 + $D)) B0PV[MW^2, MB^2, MT^2]) -
         CW^2 (MH^4 - 4 MH^2 MW^2 + 4 MW^4 (-1 + $D)) B0PV[MW^2, MH^2,
           MW^2] - (CW^2 (-4 MW^2 MZ^2 + MZ^4) +
            4 MW^4 SW^4 (-1 + $D) -
            4 CW^4 (-MZ^4 (-2 + $D) + 5 MW^4 (-1 + $D) +
               MW^2 MZ^2 (-9 + 5 $D))) B0PV[MW^2, MW^2, MZ^2]),
     dZWW1 :>
      1/(64 \[Pi])
        aem ((12 Nf (3 MB^6 - 3 MT^6 + MT^2 MW^4 +
              MB^2 (9 MT^4 + 2 MT^2 MW^2 + MW^4 (5 - 2 $D)) +
              MB^4 (-9 MT^2 + MW^2 (-4 + $D)) - MT^4 MW^2 (-2 + $D) +
              MW^6 (-2 + $D)) A0[
             MB^2])/((MB^2 - MT^2)^2 MW^4 SW^2 (-1 + $D)) - (
         2 (3 MH^6 - 12 MH^4 MW^2 + 2 MW^6 (3 - 4 $D) +
            2 MH^2 MW^4 (5 + 2 $D)) A0[MH^2])/(
         MW^4 (MH^2 -
            MW^2)^2 SW^2 (-1 + $D)) + (6 Nf (6 MT^8 +
              MB^6 (-6 MT^2 + MW^2 (-2 + $D)) -
              MB^2 (18 MT^6 + MT^4 MW^2 (10 - 3 $D) +
                 2 MT^2 MW^4 (5 - 2 $D) + MW^6 (-2 + $D)^2) +
              MT^6 MW^2 (-2 + $D) + MT^2 MW^6 (8 - 6 $D + $D^2) -
              MT^4 MW^4 (8 - 5 $D + $D^2) +
              MB^4 (18 MT^4 + MT^2 MW^2 (14 - 5 $D) +
                 MW^4 (6 - 5 $D + $D^2))) A0[
             MT^2])/(MT^2 (MB^2 - MT^2)^2 MW^4 SW^2 (-1 + $D)) +
         1/MW^6 2 ((16 MW^4 (-2 + $D))/(-3 + $D) - (
            2 MW^2 (-MH^2 + MZ^2 (-1 - 4 CW^2 (-2 + $D)) +
               MW^2 (2 + 4 CW^2 (-2 + $D) + 4 SW^2 (-2 + $D))))/(
            SW^2 (-1 + $D)) + (
            MW^2 (MH^2 - 2 MW^2) (MH^4 - 4 MH^2 MW^2 +
               4 MW^4 (-1 + $D)))/((MH^2 -
               MW^2)^2 SW^2 (-1 + $D)) - ((MH^4 - 4 MH^2 MW^2 +
               4 MW^4 (-1 + $D)) (-2 + $D))/(
            2 SW^2 (-1 + $D)) + ((MH^2 - 2 MW^2) (MH^4 - 4 MH^2 MW^2 +
               4 MW^4 (-1 + $D)) (-2 + $D))/(
            2 (MH^2 -
               MW^2) SW^2 (-1 + $D)) - (MW^2 MZ^2 (CW^2 (4 MW^2 MZ^2 -
                    MZ^4) - 4 MW^4 SW^4 (-1 + $D) +
                 4 CW^4 (-MZ^4 (-2 + $D) + 5 MW^4 (-1 + $D) +
                    MW^2 MZ^2 (-9 + 5 $D))))/(CW^2 (MW^2 -
                 MZ^2)^2 SW^2 (-1 + $D))) A0[
           MW^2] + ((CW^2 (-6 MZ^8 + 4 MW^6 MZ^2 (-1 + $D) +
                 MW^2 MZ^6 (18 + $D) - MW^4 MZ^4 (2 + 5 $D)) +
              4 MW^4 SW^4 (-2 MZ^4 (-1 + $D) - MW^4 (2 - 3 $D + $D^2) +
                  MW^2 MZ^2 (2 - 3 $D + $D^2)) +
              4 CW^4 (-6 MZ^8 (-2 + $D) +
                 MW^4 MZ^4 (-8 + 21 $D - 6 $D^2) +
                 5 MW^8 (2 - 3 $D + $D^2) +
                 MW^2 MZ^6 (-38 + 18 $D + $D^2))) A0[
             MZ^2])/(CW^2 MW^4 MZ^2 (MW^2 - MZ^2)^2 SW^2 (-1 + $D)) - (
         16 (-2 + $D) A0IR[MW^2])/(MW^2 (-3 + $D)) - (
         36 Nf (-2 + $D)^2 B0PV[MW^2, 0, 0])/(SW^2 (-1 + $D)) -
         1/(MW^4 SW^2 (-1 + $D))
           12 Nf (3 MB^4 + 3 MT^4 + MB^2 (-6 MT^2 + MW^2 (-3 + $D)) +
            MT^2 MW^2 (-3 + $D) + MW^4 (-2 + $D)) B0PV[MW^2, MB^2,
           MT^2] + (
         2 (3 MH^4 - 8 MH^2 MW^2 + 4 MW^4 (-1 + $D)) B0PV[MW^2, MH^2,
           MW^2])/(MW^4 SW^2 (-1 + $D)) +
         1/(CW^2 MW^4 SW^2 (-1 + $D))
           2 (CW^2 (-8 MW^2 MZ^2 + 3 MZ^4) +
            
            4 CW^4 (MW^2 MZ^2 (17 - 9 $D) + 3 MZ^4 (-2 + $D) +
               3 MW^4 (-1 + $D)) + 4 MW^4 SW^4 (-1 + $D)) B0PV[MW^2,
           MW^2, MZ^2]),
     dZAZ1 :>
      1/(12 CW MZ^2 \[Pi] SW (-1 + $D))
        aem (6 (SW^2 + CW^2 (3 - 2 $D)) (2 - $D) A0[MW^2] +
         Nf (2 (-3 + 4 SW^2) (-2 + $D) A0[MB^2] +
            4 (-3 + 8 SW^2) (-2 + $D) A0[MT^2] -
            MZ^2 (-27 + 76 SW^2) (-2 + $D) B0PV[MZ^2, 0,
              0] + (3 - 4 SW^2) (4 MB^2 + MZ^2 (-2 + $D)) B0PV[MZ^2,
              MB^2, MB^2] +
            2 (3 - 8 SW^2) (4 MT^2 + MZ^2 (-2 + $D)) B0PV[MZ^2, MT^2,
              MT^2]) -
         3 (4 MW^2 (SW^2 (-2 + $D) + CW^2 (-4 + 3 $D)) +
            MZ^2 (SW^2 + CW^2 (-5 + 6 $D))) B0PV[MZ^2, MW^2, MW^2]),
     dZZA1 :> (aem*(CW^2 + SW^2)*(-2 + $D)*A0[MW^2])/(2*CW*MZ^2*Pi*
         SW),
     dZAA1 :> -((
        aem (-2 + $D) (16 MW^2 Nf A0[MT^2] +
           3 MT^2 (-13 + $D) A0[MW^2]))/(72 MT^2 MW^2 \[Pi])) -
       Nf*(DeltaAlpha + (1/(6 MZ^2 \[Pi] (-1 + $D))
             aem (-2 (-2 + $D) A0[MB^2] +
              19 MZ^2 (-2 + $D) B0PV[MZ^2, 0,
                0] + (4 MB^2 + MZ^2 (-2 + $D)) B0PV[MZ^2, MB^2,
                MB^2]))),
     dZfL1[1, 1, 1] :> -((
       aem (-2 + $D)^2 (2 CW^2 MZ^2 A0[MW^2] + MW^2 A0[MZ^2]))/(
       16 CW^2 MW^2 MZ^2 \[Pi] SW^2 $D)), dZfR1[1, 1, 1] :> 0,
     dZfL1[2, 1,
       1] :> -((
        aem (-2 + $D)^2 (2 CW^2 MZ^2 A0[MW^2] +
           MW^2 (1 - 2 SW^2)^2 A0[MZ^2]))/(
        16 CW^2 MW^2 MZ^2 \[Pi] SW^2 $D)) + (-((
          aem (-2 + $D) B0PV[0, 0, 0])/(8 \[Pi])))*aqed,
     dZfR1[2, 1,
       1] :> -((aem SW^2 (-2 + $D)^2 A0[MZ^2])/(
        4 CW^2 MZ^2 \[Pi] $D)) + (-((aem (-2 + $D) B0PV[0, 0, 0])/(
          8 \[Pi])))*aqed,
     dZfL1[3, 1,
       1] :> -((
        aem (-2 + $D)^2 (18 CW^2 MZ^2 A0[MW^2] +
           MW^2 (3 - 4 SW^2)^2 A0[MZ^2]))/(
        144 CW^2 MW^2 MZ^2 \[Pi] SW^2 $D)) + (-((
          aem (-2 + $D) B0PV[0, 0, 0])/(18 \[Pi])))*aqed,
     dZfR1[3, 1,
       1] :> -((aem SW^2 (-2 + $D)^2 A0[MZ^2])/(
        9 CW^2 MZ^2 \[Pi] $D)) + (-((aem (-2 + $D) B0PV[0, 0, 0])/(
          18 \[Pi])))*aqed,
     dZfL1[4, 1,
       1] :> -((
        aem (-2 + $D)^2 (18 CW^2 MZ^2 A0[MW^2] +
           MW^2 (3 - 2 SW^2)^2 A0[MZ^2]))/(
        144 CW^2 MW^2 MZ^2 \[Pi] SW^2 $D)) + (-((
          aem (-2 + $D) B0PV[0, 0, 0])/(72 \[Pi])))*aqed,
     dZfR1[4, 1,
       1] :> -((aem SW^2 (-2 + $D)^2 A0[MZ^2])/(
        36 CW^2 MZ^2 \[Pi] $D)) + (-((aem (-2 + $D) B0PV[0, 0, 0])/(
          72 \[Pi])))*aqed,
     dMHsq1 :>
      1/(32 CW^4 MW^2 \[Pi] SW^2)
        aem (3 CW^4 MH^2 A0[MH^2] +
         2 CW^4 (MH^2 + 2 MW^2 (-1 + $D)) A0[MW^2] +
         CW^2 (CW^2 MH^2 + 2 MW^2 (-1 + $D)) A0[MZ^2] +
         9 CW^4 MH^4 B0PV[MH^2, MH^2, MH^2] -
         12 CW^4 MT^2 Nf (2 A0[MT^2] - (MH^2 - 4 MT^2) B0PV[MH^2, MT^2,
               MT^2]) +
         2 CW^4 (MH^4 - 4 MH^2 MW^2 + 4 MW^4 (-1 + $D)) B0PV[MH^2,
           MW^2, MW^2] + (CW^4 MH^4 - 4 CW^2 MW^2 (MH^2 + MZ^2) +
            4 MW^4 $D) B0PV[MH^2, MZ^2, MZ^2]),
     dZH1 :> 1/(64 MW^2 \[Pi] SW^2)
        aem (-6 (-2 + $D) A0[MH^2] + (24 MT^2 Nf (-2 + $D) A0[MT^2])/
         MH^2 + (4 (MH^4 - 4 MH^2 MW^2 +
            4 MW^4 (-1 + $D)) (-2 + $D) A0[MW^2])/(
         MH^4 - 4 MH^2 MW^2) + (
         2 (-2 + $D) (CW^4 MH^4 - 4 CW^2 MW^2 (MH^2 + MZ^2) +
            4 MW^4 $D) A0[MZ^2])/(CW^4 MH^2 (MH^2 - 4 MZ^2)) +
         3 MH^2 $D B0PV[MH^2, MH^2, MH^2] - (
         12 MT^2 Nf (4 MT^2 + MH^2 (-2 + $D)) B0PV[MH^2, MT^2, MT^2])/
         MH^2 - 1/(MH^4 - 4 MH^2 MW^2)
           2 (MH^6 (-4 + $D) - 4 MH^4 MW^2 (-3 + $D) +
            16 MW^6 (-1 + $D) + 4 MH^2 MW^4 (8 - 5 $D + $D^2)) B0PV[
           MH^2, MW^2, MW^2] -
         1/(CW^4 MH^2 (MH^2 -
             4 MZ^2)) (CW^4 (4 MH^4 MZ^2 + MH^6 (-4 + $D)) -
            4 CW^2 MW^2 (4 MZ^4 + MH^2 MZ^2 (-8 + $D) +
               MH^4 (-2 + $D)) +
            4 MW^4 (4 MZ^2 + MH^2 (-4 + $D)) $D) B0PV[MH^2, MZ^2,
           MZ^2]), dMf1[3, 3, 3] :>
      1/(576 CW^2 MT MW^2 \[Pi] SW^2)
        aem (-18 CW^2 (MB^2 + MT^2 + MW^2 (-2 + $D)) A0[MB^2] +
         18 CW^2 MT^2 A0[MH^2] - 36 CW^2 MT^2 A0[MT^2] +
         18 MW^2 A0[MT^2] - 48 MW^2 SW^2 A0[MT^2] +
         64 CW^2 MW^2 SW^2 A0[MT^2] + 64 MW^2 SW^4 A0[MT^2] -
         9 MW^2 $D A0[MT^2] + 24 MW^2 SW^2 $D A0[MT^2] -
         32 CW^2 MW^2 SW^2 $D A0[MT^2] - 32 MW^2 SW^4 $D A0[MT^2] +
         18 CW^2 MB^2 A0[MW^2] + 18 CW^2 MT^2 A0[MW^2] -
         36 CW^2 MW^2 A0[MW^2] + 18 CW^2 MW^2 $D A0[MW^2] +
         18 CW^2 MT^2 A0[MZ^2] - 18 MW^2 A0[MZ^2] +
         48 MW^2 SW^2 A0[MZ^2] - 64 MW^2 SW^4 A0[MZ^2] +
         9 MW^2 $D A0[MZ^2] - 24 MW^2 SW^2 $D A0[MZ^2] +
         32 MW^2 SW^4 $D A0[MZ^2] -
         128 CW^2 MT^2 MW^2 SW^2 B0PV[MT^2, 0, MT^2] +
         18 CW^2 MB^4 B0PV[MT^2, MB^2, MW^2] -
         36 CW^2 MB^2 MT^2 B0PV[MT^2, MB^2, MW^2] +
         18 CW^2 MT^4 B0PV[MT^2, MB^2, MW^2] -
         54 CW^2 MB^2 MW^2 B0PV[MT^2, MB^2, MW^2] -
         54 CW^2 MT^2 MW^2 B0PV[MT^2, MB^2, MW^2] +
         36 CW^2 MW^4 B0PV[MT^2, MB^2, MW^2] +
         18 CW^2 MB^2 MW^2 $D B0PV[MT^2, MB^2, MW^2] +
         18 CW^2 MT^2 MW^2 $D B0PV[MT^2, MB^2, MW^2] -
         18 CW^2 MW^4 $D B0PV[MT^2, MB^2, MW^2] -
         18 CW^2 MH^2 MT^2 B0PV[MT^2, MH^2, MT^2] +
         72 CW^2 MT^4 B0PV[MT^2, MH^2, MT^2] -
         36 MT^2 MW^2 B0PV[MT^2, MT^2, MZ^2] -
         18 CW^2 MT^2 MZ^2 B0PV[MT^2, MT^2, MZ^2] +
         18 MW^2 MZ^2 B0PV[MT^2, MT^2, MZ^2] +
         96 MT^2 MW^2 SW^2 B0PV[MT^2, MT^2, MZ^2] -
         48 MW^2 MZ^2 SW^2 B0PV[MT^2, MT^2, MZ^2] -
         128 MT^2 MW^2 SW^4 B0PV[MT^2, MT^2, MZ^2] +
         64 MW^2 MZ^2 SW^4 B0PV[MT^2, MT^2, MZ^2] +
         18 MT^2 MW^2 $D B0PV[MT^2, MT^2, MZ^2] -
         9 MW^2 MZ^2 $D B0PV[MT^2, MT^2, MZ^2] +
         24 MW^2 MZ^2 SW^2 $D B0PV[MT^2, MT^2, MZ^2] -
         32 MW^2 MZ^2 SW^4 $D B0PV[MT^2, MT^2, MZ^2]),
     dMf1[3, 3, 3] :>
      1/(576 CW^2 MB MW^2 \[Pi] SW^2)
        aem ((-36 CW^2 MB^2 -
            MW^2 (9 - 12 SW^2 + 8 SW^4) (-2 + $D)) A0[MB^2] +
         18 CW^2 MB^2 A0[MH^2] - 18 CW^2 MB^2 A0[MT^2] -
         18 CW^2 MT^2 A0[MT^2] + 36 CW^2 MW^2 A0[MT^2] -
         18 CW^2 MW^2 $D A0[MT^2] + 18 CW^2 MB^2 A0[MW^2] +
         18 CW^2 MT^2 A0[MW^2] - 36 CW^2 MW^2 A0[MW^2] +
         18 CW^2 MW^2 $D A0[MW^2] + 18 CW^2 MB^2 A0[MZ^2] -
         18 MW^2 A0[MZ^2] + 24 MW^2 SW^2 A0[MZ^2] -
         16 MW^2 SW^4 A0[MZ^2] + 9 MW^2 $D A0[MZ^2] -
         12 MW^2 SW^2 $D A0[MZ^2] + 8 MW^2 SW^4 $D A0[MZ^2] +
         72 CW^2 MB^4 B0PV[MB^2, MB^2, MH^2] -
         18 CW^2 MB^2 MH^2 B0PV[MB^2, MB^2, MH^2] -
         36 MB^2 MW^2 B0PV[MB^2, MB^2, MZ^2] -
         18 CW^2 MB^2 MZ^2 B0PV[MB^2, MB^2, MZ^2] +
         18 MW^2 MZ^2 B0PV[MB^2, MB^2, MZ^2] +
         48 MB^2 MW^2 SW^2 B0PV[MB^2, MB^2, MZ^2] -
         24 MW^2 MZ^2 SW^2 B0PV[MB^2, MB^2, MZ^2] -
         32 MB^2 MW^2 SW^4 B0PV[MB^2, MB^2, MZ^2] +
         16 MW^2 MZ^2 SW^4 B0PV[MB^2, MB^2, MZ^2] +
         18 MB^2 MW^2 $D B0PV[MB^2, MB^2, MZ^2] -
         9 MW^2 MZ^2 $D B0PV[MB^2, MB^2, MZ^2] +
         12 MW^2 MZ^2 SW^2 $D B0PV[MB^2, MB^2, MZ^2] -
         8 MW^2 MZ^2 SW^4 $D B0PV[MB^2, MB^2, MZ^2] +
         18 CW^2 MB^4 B0PV[MB^2, MT^2, MW^2] -
         36 CW^2 MB^2 MT^2 B0PV[MB^2, MT^2, MW^2] +
         18 CW^2 MT^4 B0PV[MB^2, MT^2, MW^2] -
         54 CW^2 MB^2 MW^2 B0PV[MB^2, MT^2, MW^2] -
         54 CW^2 MT^2 MW^2 B0PV[MB^2, MT^2, MW^2] +
         36 CW^2 MW^4 B0PV[MB^2, MT^2, MW^2] +
         18 CW^2 MB^2 MW^2 $D B0PV[MB^2, MT^2, MW^2] +
         18 CW^2 MT^2 MW^2 $D B0PV[MB^2, MT^2, MW^2] -
         18 CW^2 MW^4 $D B0PV[MB^2, MT^2, MW^2]),
     dZW1 :> ((3*
           aem*(3*MB^6 - 3*MT^6 + MT^2*MW^4 +
             MB^2*(9*MT^4 + 2*MT^2*MW^2 + MW^4*(5 - 2*$D)) +
             MB^4*(-9*MT^2 + MW^2*(-4 + $D)) - MT^4*MW^2*(-2 + $D) +
             MW^6*(-2 + $D))*A0[MB^2]*Nf)/(16*(MB^2 - MT^2)^2*MW^4*Pi*
           SW^2*(-1 + $D)) - (aem*(3*MH^6 - 12*MH^4*MW^2 +
             2*MW^6*(3 - 4*$D) + 2*MH^2*MW^4*(5 + 2*$D))*A0[MH^2])/(32*
           MW^4*(MH^2 - MW^2)^2*Pi*SW^2*(-1 + $D)) + (3*
           aem*(6*MT^8 + MB^6*(-6*MT^2 + MW^2*(-2 + $D)) -
             MB^2*(18*MT^6 + MT^4*MW^2*(10 - 3*$D) +
                2*MT^2*MW^4*(5 - 2*$D) + MW^6*(-2 + $D)^2) +
             MT^6*MW^2*(-2 + $D) + MT^2*MW^6*(8 - 6*$D + $D^2) -
             MT^4*MW^4*(8 - 5*$D + $D^2) +
             MB^4*(18*MT^4 + MT^2*MW^2*(14 - 5*$D) +
                MW^4*(6 - 5*$D + $D^2)))*A0[MT^2]*Nf)/(32*
           MT^2*(MB^2 - MT^2)^2*MW^4*Pi*
           SW^2*(-1 + $D)) + (aem*((MW^2*(MH^2 - 2*MW^2)*(MH^4 -
                  4*MH^2*MW^2 + 4*MW^4*(-1 + $D)))/((MH^2 - MW^2)^2*
                SW^2) - ((MH^4 - 4*MH^2*MW^2 +
                  4*MW^4*(-1 + $D))*(-2 + $D))/(2*
                SW^2) + ((MH^2 - 2*MW^2)*(MH^4 - 4*MH^2*MW^2 +
                  4*MW^4*(-1 + $D))*(-2 + $D))/(2*(MH^2 - MW^2)*
                SW^2) + (16*MW^4*(-2 + $D)*(-1 + $D))/(-3 + $D) + (4*
                MW^4*(2*CW^2 + 2*SW^2 - $D)*(-3 + 2*$D))/
              SW^2 - (2*
                MW^2*(-MH^2 + MZ^2*(-1 - 4*CW^2*(-2 + $D)) +
                  2*MW^2*(1 + 3*$D - 2*$D^2 + 2*CW^2*(-5 + 3*$D) +
                     2*SW^2*(-5 + 3*$D))))/
              
              SW^2 - (MW^2*
                MZ^2*(CW^2*(4*MW^2*MZ^2 - MZ^4) -
                  4*MW^4*SW^4*(-1 + $D) +
                  4*CW^4*(-(MZ^4*(-2 + $D)) + 5*MW^4*(-1 + $D) +
                     MW^2*MZ^2*(-9 + 5*$D))))/(CW^2*(MW^2 - MZ^2)^2*
                SW^2))*A0[MW^2])/(32*MW^6*
           Pi*(-1 + $D)) + (aem*(CW^2*(-6*MZ^8 +
                4*MW^6*MZ^2*(-1 + $D) + MW^2*MZ^6*(18 + $D) -
                MW^4*MZ^4*(2 + 5*$D)) +
             4*MW^4*SW^4*(-2*MZ^4*(-1 + $D) - MW^4*(2 - 3*$D + $D^2) +
                MW^2*MZ^2*(2 - 3*$D + $D^2)) +
             4*CW^4*(-6*MZ^8*(-2 + $D) +
                MW^4*MZ^4*(-8 + 21*$D - 6*$D^2) +
                5*MW^8*(2 - 3*$D + $D^2) +
                MW^2*MZ^6*(-38 + 18*$D + $D^2)))*A0[MZ^2])/(64*CW^2*
           MW^4*MZ^2*(MW^2 - MZ^2)^2*Pi*
           SW^2*(-1 + $D)) - (aem*(-2 + $D)*A0IR[MW^2])/(4*MW^2*
           Pi*(-3 + $D)) - (9*aem*(-2 + $D)^2*B0PV[MW^2, 0, 0]*Nf)/(16*
           Pi*SW^2*(-1 + $D)) - (3*
           aem*(3*MB^4 + 3*MT^4 + MB^2*(-6*MT^2 + MW^2*(-3 + $D)) +
             MT^2*MW^2*(-3 + $D) + MW^4*(-2 + $D))*
           B0PV[MW^2, MB^2, MT^2]*Nf)/(16*MW^4*Pi*
           SW^2*(-1 + $D)) + (aem*(3*MH^4 - 8*MH^2*MW^2 +
             4*MW^4*(-1 + $D))*B0PV[MW^2, MH^2, MW^2])/(32*MW^4*Pi*
           
           SW^2*(-1 + $D)) + (aem*(CW^2*(-8*MW^2*MZ^2 + 3*MZ^4) +
             4*CW^4*(MW^2*MZ^2*(17 - 9*$D) + 3*MZ^4*(-2 + $D) +
                3*MW^4*(-1 + $D)) + 4*MW^4*SW^4*(-1 + $D))*
           B0PV[MW^2, MW^2, MZ^2])/(32*CW^2*MW^4*Pi*SW^2*(-1 + $D)))}];
