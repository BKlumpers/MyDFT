/*****************************************************************

Basic closed-shell spin-restricted DFT-solver for simple molecules using STO-NG

Authors: B. Klumpers
		 I.A.W. Filot

Published under GNU General Public License 3.0
Source code available at: https://github.com/BKlumpers/dft

Allows for SCF-computation of molecular energies for simple molecules.
Includes testcases for: H, He, H2, HeH+, He2, CO, and H2O.

*****************************************************************/

#ifndef _DFT_H
#define _DFT_H

#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "integrals.h"
#include "orbitals.h"
#include "solvers.h"
#include "quadrature.h"

typedef Eigen::Vector3d vec3; //define vec3 as the Eigen::Vector3d-object

using namespace std;

const bool Quadrature = true; //switch between Quadrature and Cartesian integration

//define function calls
SCF_E SCF_DFT_energy(vector<CGF> AO_list, const vector<vec3>& pos_list, const vector<double>& charge_list, const vector<int>& nelec_list, const vector<int>& atnum_list);
double density(const vec3& pos, const Eigen::MatrixXd& Pmatrix, vector<CGF> AO_list);
double density(int atom, int rpnt, int apnt, const Eigen::MatrixXd& Pmatrix, grid Grid);
double exchange_Dirac(double Density);
double correlation_VWN(double Density);
double VWN_potential(double Density);

#endif //_DFT_H