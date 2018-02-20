/*****************************************************************

Basic closed-shell spin-restricted HF-solver for simple molecules using STO-NG

Author: B. Klumpers (bartkl@live.nl)

Published under GNU General Public License 3.0

Allows for SCF-computation of molecular energies for simple molecules.
Testcases for H, He, H2, HeH+ and He2 are included.

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