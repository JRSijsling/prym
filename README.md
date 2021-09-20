Description
--

This repository contains Magma code to verify the statements in the article

Nicolas Mascot, Jeroen Sijsling, and John Voight  
*A Prym variety with everywhere good reduction over QQ (sqrt (61))*
Preprint available at [arXiv:1908.00421](https://arxiv.org/abs/1908.00421)

The file `PariPeriodsOutput.txt` and `Matrices.m` contain the output of the period matrix calculations in Section 2. The files whose names start with `FindEquation` perform the isolation of the building block in Section 3 and the recovery of a first equation for an algebraic curve in Section 4. The simplification and verification in Section 5 takes place in the places starting with `Simplify`, `Certify` (verification of endomorphisms), and `Verify` (verification of the curve and its modularity). These algorithms also require some very summary routines for invariants of binary quartics `binary_quartics.m`.

Prerequisites
--

An installation of Magma and the dependency [`JRSijsling/curve_reconstruction`](https://github.com/JRSijsling/curve_reconstruction).

