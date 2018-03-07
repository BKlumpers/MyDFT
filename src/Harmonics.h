/*****************************************************************

Basic closed-shell spin-restricted DFT-solver for simple molecules using STO-NG

Authors: B. Klumpers
		 I.A.W. Filot

Published under GNU General Public License 3.0
Source code available at: https://github.com/BKlumpers/dft

Allows for SCF-computation of molecular energies for simple molecules.
Includes testcases for: H, He, H2, HeH+, He2, CO, and H2O.

*****************************************************************/

#ifndef _HARMONICS_H
#define _HARMONICS_H

#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include <boost/math/special_functions/spherical_harmonic.hpp>

typedef Eigen::Vector3d vec3; //define vec3 as the Eigen::Vector3d-object

using namespace std;

double HarmonicReal(int l, int m, vec3 gridpoint);

#endif //_HARMONICS_H