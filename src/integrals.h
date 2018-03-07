/*****************************************************************

Basic closed-shell spin-restricted DFT-solver for simple molecules using STO-NG

Authors: B. Klumpers
		 I.A.W. Filot

Published under GNU General Public License 3.0
Source code available at: https://github.com/BKlumpers/dft

Allows for SCF-computation of molecular energies for simple molecules.
Includes testcases for: H, He, H2, HeH+, He2, CO, and H2O.

*****************************************************************/

#ifndef _INTEGRALS_H
#define _INTEGRALS_H

#include <iostream>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <Eigen/Dense>
#include "prelim_math.h"
#include "orbitals.h"

typedef Eigen::Vector3d vec3; //define vec3 as the Eigen::Vector3d-object

using namespace std;

double overlapGTO(GTO GTO1, GTO GTO2);
double binomialComposite(int k, int a1, int a2, double r1, double r2);
double binomAlt(int j, int a1, int a2, double r1, double r2);
double kineticGTO(GTO GTO1, GTO GTO2);
double nuclearGTO(GTO GTO1, GTO GTO2, double charge, const vec3& nucleus);
double nuclearComponent(int l, int r, int i, int a1, int a2, double PA, double PB, double PC, double alpha_c);
double fmch(int nu, double x);
double two_electronGTO(GTO GTO1, GTO GTO2, GTO GTO3, GTO GTO4);
double electricComponent(int l1, int l2, int r1, int r2, int i, int a1, int a2, int a3, int a4, double AX, double BX, double CX, double DX, double PX, double QX, double g12, double g34, double delta);
double electricFourier(int l, int a1, int a2, double PA, double PB, double r1, double gamma);
double overlapCGF(CGF CGF1, CGF CGF2);
double kineticCGF(CGF CGF1, CGF CGF2);
double nuclearCGF(CGF CGF1, CGF CGF2, double charge, const vec3& pos);
double two_electronCGF(CGF CGF1, CGF CGF2, CGF CGF3, CGF CGF4);
int two_electronSort(int i, int j, int k, int l);

#endif //_INTEGRALS_H