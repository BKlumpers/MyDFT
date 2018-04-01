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

#ifndef _DFT_H
#define _DFT_H

#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "integrals.h"
#include "orbitals.h"
#include "hf.h"
#include "quadrature.h"

typedef Eigen::Vector3d vec3; //define vec3 as the Eigen::Vector3d-object

using namespace std;

const bool Quadrature = true; //switch between Quadrature and Cartesian integration
const bool Poisson = false; //switch between Poisson potential and repulsion matrix
const bool GGA = true;

//define function calls
SCF_E SCF_DFT_energy(vector<CGF> AO_list, const vector<vec3>& pos_list, const vector<double>& charge_list, const vector<int>& nelec_list, const vector<int>& atnum_list);
double density(const vec3& pos, const Eigen::MatrixXd& Pmatrix, vector<CGF> AO_list);
double density(int atom, int rpnt, int apnt, const Eigen::MatrixXd& Pmatrix, grid Grid);
double exchange_Dirac(double Density);
double correlation_VWN(double Density);
double VWN_potential(double Density);

//additional function calls for udft:
SCF_E SCF_UDFT_energy(vector<CGF> AO_list, const vector<vec3>& pos_list, const vector<double>& charge_list, const vector<int>& nelec_list, const vector<int>& atnum_list);
double exchange_Dirac(double Density_alpha, double Density_beta);
double exchange_potential(double Density_alpha, double Density_beta, bool ab);
double correlation_VWN(double Density_alpha, double Density_beta);
double VWN_potential(double Density_alpha, double Density_beta, bool ab);
Eigen::VectorXd exchange_correlation(double Density_alpha, double Density_beta);
Eigen::VectorXd PBE(double Density, double gradient);
Eigen::VectorXd PBE(double Density_alpha, double Density_beta, Eigen::VectorXd gradienta, Eigen::VectorXd gradientb, Eigen::MatrixXd Hessiana, Eigen::MatrixXd Hessianb);

#endif //_DFT_H