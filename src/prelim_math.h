/*****************************************************************

MYDFT: Basic DFT-solver for simple molecules using STO-NG

Authors: B. Klumpers
         I.A.W. Filot

Published under GNU General Public License 3.0
Source code available at: https://github.com/BKlumpers/dft

Allows for SCF-computation of molecular energies for simple molecules
using both HF and DFT. Supports both closed and open shell systems.
Includes testcases for: H, He, H2, HeH+, HeH, He2, CO, and H2O.

*****************************************************************/

#ifndef _PRELIM_MATH_H
#define _PRELIM_MATH_H

#include <iostream>
#include <cmath>

using namespace std;

const bool EFLAG = false; //turn mathematical errors off

//define function calls
int factorial(int n);
int doublefactorial_odd(int k);
int binomialCoeff(int n, int k);

#endif //_PRELIM_MATH_H