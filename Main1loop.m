(* ::Package:: *)

 (* ------------------------------------------------------------- *)
 (* Tools for 1-loop calculations                                 *)
 (* version 1.0 (April 28, 2022)                                  *)
 (* ------------------------------------------------------------- *)
 (* Main1loop.m                                                   *)
 (* Daniel Wiegand (daniel.wiegand@northwestern.edu)              *)
 (* last revision 04/28/22                                        *)
 (* ------------------------------------------------------------- *)
 (* Mathematica file for loading all files necessary for the      *)
 (* 1loop calculation                                             *)
 (* ------------------------------------------------------------- *)

CurrentPath = NotebookDirectory[];

 <<(CurrentPath<>"FAtoFCreplace.m");
 <<(CurrentPath<>"FctRename.m");
 <<(CurrentPath<>"DenCancel.m");
 <<(CurrentPath<>"PVred.m");
 <<(CurrentPath<>"FinalTrans.m");
 
 EndPackage[]
