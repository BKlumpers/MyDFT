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

#ifndef _GEO_H
#define _GEO_H

#include <iostream>
//#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "orbitals.h"
#include "hf.h"
#include "dft.h"
#include "quadrature.h"

typedef Eigen::Vector3d vec3; //define vec3 as the Eigen::Vector3d-object

using namespace std;

//perform geometry optimisation through numerical computation of 1st derivatives

class GEO{
public:
	//contain output info -> optim geometry parameters, optim_energy
	double energy;
	vector<vec3> pos;
};

GEO GEO_HF_nograd(vector<CGF> AO_list, const vector<vec3>& pos_init, const vector<double>& charge_list, const vector<int>& nelec_list);
GEO GEO_HF_numerical(vector<CGF> AO_list, const vector<vec3>& pos_init, const vector<double>& charge_list, const vector<int>& nelec_list);

#endif //_GEO_H