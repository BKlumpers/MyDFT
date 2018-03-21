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

#ifndef _SOLVERS_H
#define _SOLVERS_H

#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "integrals.h"
#include "orbitals.h"

typedef Eigen::Vector3d vec3; //define vec3 as the Eigen::Vector3d-object

using namespace std;

//parse output of SCF energy minimisation
class SCF_E {
public:
    double energy; //total energy of the system
    Eigen::VectorXd orbital_energies; //energies of the molecular orbitals
    Eigen::MatrixXd C_vector; //coefficient matrix
    void SCF_result(double Result_energy, const Eigen::VectorXd& Result_orbital_energies, const Eigen::MatrixXd& Result_C_vector)
    {
        energy = Result_energy;
        orbital_energies = Result_orbital_energies;
        C_vector = Result_C_vector;
    }
};

//define function calls
SCF_E SCF_HF_energy(vector<CGF> AO_list, const vector<vec3>& pos_list, const vector<double>& charge_list, const vector<int>& nelec_list);

//function call to spin unrestricted: (part of uhf.cpp)
SCF_E SCF_UHF_energy(vector<CGF> AO_list, const vector<vec3>& pos_list, const vector<double>& charge_list, const vector<int>& nelec_list);

#endif //_SOLVERS_H