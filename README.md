# 1loopEW
Mathematica package for automating 1-loop calculations based on FeynArts/FeynCalc

This was part of my PhD thesis - it's a little clunky - proper manual is coming. Only works with FeynCalc already running.

commands/functions:

FAtoFCtrans[--] - translates FeynArts Born/1-loop output into useable FeynCalc input

RenameMomenta[--] - after squaring the amplitude, translates and collects all products of loop momenta

AllCancel[--] - Cancels all numerators vs denominators as much as possible

AllPV[--] - performs a Passarino-Veltman reduction to get rid of the remaining loop-momenta in numerator

CleanUp[--] - translates and cleans up the final reduced expressions in terms of basis integrals

RenConst[--] - 1-loop SM counterterms, "aqed" multiplies photon terms of external fermion wavefunctions (proportional to B0[0,0,0])

The packages can be loaded into a mathematica session separately or all together using Main1loop.m

If using/modifying make sure to cite the proper sources.

Daniel Wiegand (04/28/2022)

