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

#ifndef _QUADRATURE_H
#define _QUADRATURE_H

#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "lebedev.h"
#include "orbitals.h"
#include "Harmonics.h"
#include <boost/math/interpolators/cubic_b_spline.hpp>

typedef Eigen::Vector3d vec3; //define vec3 as the Eigen::Vector3d-object

using namespace std;

const bool Bragg = false; //set whether to include scaling of the radial grid according to modified Slater-Bragg radii

class Hessian{
private:
	Eigen::MatrixXd Hess = Eigen::MatrixXd::Zero(3,3);
public:
	void set_Hess(const Eigen::MatrixXd& H){
		Hess = H;
	}
	Eigen::MatrixXd get_Hess(){ //function overload
		return Hess;
	}
	double get_Hess(int i, int j){
		return Hess(i,j);
	}
};

//Hessian of wavefunction element ij, for each of the gridpoints
class waveHess{
private:
	vector<Hessian> Hessians;
public:
	waveHess(int pts){
		Hessian Hess;
		for(int i=0; i<pts; i++){
			Hessians.push_back(Hess);
		}
	}
	void set_Hess(int idx, const Eigen::MatrixXd& H){
		Hessians[idx].set_Hess(H);
	}
	double get_Hess(int idx, int i, int j){ //function overload
		return Hessians[idx].get_Hess(i,j);
	}
	Eigen::MatrixXd get_Hess(int idx){
		return Hessians[idx].get_Hess();
	}
};

//object for atomic integration grid
//{gridpoints} matrix[3,m] with each column m representing a gridpoint with its {x,y,z}-coordinates stored in the 3 rows
//{amps} matrix[N+1,m] with each row containing the amplitudes of 1 basis function at each of the aforementioned m gridpoints
//-------------------->the N+1th row is used for storing the total local electron density as computed by grid::set_density
//{Laplace} contains the Laplace coefficients for the local Poisson potential; matrix[rpnts, lmttl], each value of {l,m} has values at each radial point
//{Ulm} contains the radial potential function, hence dimensions are equal to Laplace
//{potential} contains the local and external Poisson potential at each radial point; matrix[2, radial_points], 1st row local, 2nd row external
class atomicgrid {
public:
	Eigen::MatrixXd gridpoints; //3xPts matrix contain {x,y,z} of each point on the grid
	Eigen::MatrixXd amps; //N+1 x Pts matrix containing tha amplitudes of each of the N wavefunctions at each point, final row stores densities
	Eigen::MatrixXd spin; //2xPts matrix containing alpha and beta densities at each point
	Eigen::MatrixXd Laplace; //rpnt x ns contains the values for the spherical decompositions of the radial density for each radial shell
	Eigen::MatrixXd Ulm; //contains the radial part of the internal Poisson potential
	Eigen::MatrixXd potential; //contains the internal and external Poisson potentials
	Eigen::MatrixXd gradients_x; //contains the value of the wavefunction derivatives at each point, for each basis function couple ij: N*N x pts
	Eigen::MatrixXd gradients_y;
	Eigen::MatrixXd gradients_z;
	Eigen::MatrixXd gradients_norm;
	Eigen::MatrixXd PBE; //contains the value of the density gradient at each point
	Eigen::MatrixXd PBEC; //contains the value of the alpha and beta gradients at each point
	vector<Hessian> rho_Hess; //contain the Hessian of the density at each point
	vector<Hessian> alpha_Hess;
	vector<Hessian> beta_Hess;
	vector<waveHess> wave_Hess; //contain a vector of ij Hessians at each point
	double point_charge;
	double Slater;
	//class constructor
	atomicgrid(int radial_points, int angular_points, int Basis, int lmttl){
		int pts = radial_points*angular_points;
		gridpoints = Eigen::MatrixXd::Zero(3, pts); //create empty atomic grid
		amps = Eigen::MatrixXd::Zero(Basis+1, pts); //create empty amplitude and density maps
		spin = Eigen::MatrixXd::Zero(2, pts); //store alpha and beta densities
		//make Laplace contain 2nd index for each l,m couple; first index along radii
		Laplace = Eigen::MatrixXd::Zero(radial_points, lmttl);
		//make potential contain internal + external potentials, for each gridpoint
		potential = Eigen::MatrixXd::Zero(2, pts);
		//Ulm has the same dimensions as Laplace
		Ulm = Laplace; //same initialisation
		gradients_x = Eigen::MatrixXd::Zero(Basis*Basis,pts);
		gradients_y = Eigen::MatrixXd::Zero(Basis*Basis,pts);
		gradients_z = Eigen::MatrixXd::Zero(Basis*Basis,pts);
		gradients_norm = Eigen::MatrixXd::Zero(Basis*Basis,pts);
		PBE = Eigen::MatrixXd::Zero(4,pts); //first 3 rows contain x,y,z; final row contains norm
		PBEC = Eigen::MatrixXd::Zero(8,pts); //first 4 rows contain alpha, second block contains beta
		Hessian Hess;
		for(int i=0; i<pts; i++){
			rho_Hess.push_back(Hess); //store 1 Hessian per point
			alpha_Hess.push_back(Hess);
			beta_Hess.push_back(Hess);
		}
		waveHess wave(pts);
		for(int i=0; i<Basis*Basis; i++){
			wave_Hess.push_back(wave); //store i*j Hessians per point
		}
	}
};

//define output structure containing the grid
//
//comprises:
//
//matrix[n,m] of vec3 elements containing the mesh-centra
//matrix[n,m] containing the weights of the mesh-units
//where n is the index of the nuclei and m are the mesh-centra and gridpoints respectively
//
//{lmttl} is the total number of real spherical harmonics in the basis
//{lebedev_order} is the order of the Lebedev quadrature, which relates to the highest angular momentum in the harmonics expansion as: lmax = lebedev_order/2
class grid {
private:
	static constexpr double pi = 3.14159265359;
	int lmttl, lmax, angular_index;
	int radial_points, angular_points; //store number of radial and angular gridpoints respectively
	Eigen::MatrixXd weights; //Matrix containing the weights, first index is the nucleic label, second is the gridpoint-label
	Eigen::MatrixXd Gauss;
	vector<atomicgrid> gridpoints; //vector containing the atomic grid, index runs across the nucleic labels
	Eigen::MatrixXd U7; //heptadiagonal matrix storage for computation of the radial Poisson potential
	vector<vec3> pos_list; //need the nucleic coordinates for spline interpolation of the external Poisson potentials
public:
	//initialise object with predefined gridsizes and known number of atoms; //----------> pos_list to replace int atoms
	grid(const vector<vec3>& Pos_list, int radial, int angular, int Basis, int Lebedev_order, int Index) {
		radial_points = radial; //store number of radial gridpoints inside the class
		angular_points = angular; //store number of angular gridpoints inside the class
		angular_index = Index;
		pos_list = Pos_list; //store nucleic positions
		//Poisson:
		lmax = Lebedev_order + 1; //store this parameter, since Lebedev::lebedev_num = n-1 gives Lebedev_order 2n+1, +1 is to correct for n starting at 0
		int lmttl = 0;
		for(int i=0; i<(1+lmax); i++){ //lmax = lebedev_order/2;
			lmttl += 2*i+1; //add the number of m-values for degree l, sum over all degrees to get the total number of harmonics
		}
		atomicgrid emptygrid(radial_points, angular_points, Basis, lmttl); //create empty atomic grid
		weights = Eigen::MatrixXd::Zero(pos_list.size(), radial_points*angular_points); //create empty matrix for the weights
		Gauss = weights;
		for(int i=0; i<pos_list.size(); i++){
			gridpoints.push_back(emptygrid); //initialise an empty atomic grid for each atom
		}
		U7 = Eigen::MatrixXd::Zero(radial_points, radial_points); //heptadiagonal matrix for finite difference solution of radial Poisson potential
	}
	//obtain coordinates of a particular gridpoint
	vec3 get_gridpoint(int atomnr, int radial, int angular){
		static vec3 point;
		point << gridpoints[atomnr].gridpoints.col(radial*angular_points + angular);
		return point;
	}
	//assign value to gridpoint
	void set_gridpoint(int atomnr, int radial, int angular, const vec3& point){
		static int index;
		index = radial*angular_points + angular;
		gridpoints[atomnr].gridpoints(0,index) = point[0];
		gridpoints[atomnr].gridpoints(1,index) = point[1];
		gridpoints[atomnr].gridpoints(2,index) = point[2];
	}
	//obtain value of specific weight
	double get_weight(int atomnr, int radial, int angular){
		static int index;
		index = radial*angular_points + angular;
		return weights(atomnr, index);
	}
	//assign value to specific weight
	void set_weight(int atomnr, int radial, int angular, double weight){
		static int index;
		index = radial*angular_points + angular;
		weights(atomnr, index) = weight;
	}
	//apply Gaussian weight
	void Gauss_weight(int atomnr, int radial, int angular, double weight){
		static int index;
		index = radial*angular_points + angular;
		weights(atomnr, index) *= weight; //apply Gauss weight
		Gauss(atomnr, index) = weight; //store Gauss weight
	}
	//write all the wave-function amplitudes to atomicgrid
	void write_amps(vector<CGF> AO_list);

	//get specific wavefunction amplitude
	double get_amps(int atomnr, int radial, int angular, int Basis){
		return gridpoints[atomnr].amps(Basis,radial*angular_points+angular);
	}
	//write values for the electron density
	void set_density(const Eigen::MatrixXd& Pmatrix);

	//write values in spin-polarised systems
	void set_density(const Eigen::MatrixXd& P_alpha, const Eigen::MatrixXd& P_beta); //function overload

	//retrieve value of the density at specified gridpoint
	double get_density(int atomnr, int radial, int angular){
		return gridpoints[atomnr].amps(gridpoints[atomnr].amps.rows()-1, radial*angular_points + angular);
	}
	//retrieve alpha/beta density in case of spin polarisation
	double get_density(int atomnr, int radial, int angular, bool ab){ //function overload
		static int id = 0;
		if(ab){
			id = 0; //get alpha density
		}
		if(!ab){
			id = 1; //get beta density
		}
		return gridpoints[atomnr].spin(id, radial*angular_points + angular);
	}
	//retrieve point-charge potential of specified atomic grid
	double get_point_charge(int atomnr){
		return gridpoints[atomnr].point_charge;
	}
	//retrieve number of radial points
	int get_rpnt(){
		return radial_points;
	}
	//retrieve number of angular points
	int get_apnt(){
		return angular_points;
	}
	//retrieve number of atomic grids
	int get_atomnr(){
		return gridpoints.size();
	}
	//retrieve Gauss weight of the gridpoint
	double get_Gauss(int atnum, int rpnt, int apnt){
		static int index;
		index = rpnt*angular_points + apnt;
		return Gauss(atnum, index);
	}
	//store value of the SlaterBragg radius
	void set_Slater(int atomnr, double SlaterBragg){
		gridpoints[atomnr].Slater = SlaterBragg;
	}
	//get Slater-Bragg radius belonging to the nuclei {atomnr}
	double get_Slater(int atomnr){
		return gridpoints[atomnr].Slater;
	}
	//compute the wavefunction gradients and Hessians at each point
	void write_gradient(vector<CGF> AO_list);

	//compute the density gradient and Hessian at each point
	void set_gradient(const Eigen::MatrixXd& Pmatrix);

	//comput the alpha and beta gradients and Hessians at each point
	void set_gradient(const Eigen::MatrixXd& P_alpha, const Eigen::MatrixXd& P_beta); //function overload

	//get density gradient
	//{element} = {0,1,2,3} <> {x,y,z,norm}
	double get_gradient(int atom, int rpnt, int apnt, int element){
		static double index;
		index = rpnt*angular_points + apnt;
		return gridpoints[atom].PBE(element, index);
	}
	//get alpha or beta density gradient
	//{element} = {0,1,2,3} <> {x,y,z,norm}
	double get_gradient(int atom, int rpnt, int apnt, int element, bool ab){ //function overload
		static int id = 0;
		if(ab){
			id = 0; //get alpha density
		}
		if(!ab){
			id = 4; //get beta density
		}
		return gridpoints[atom].PBEC(id+element, rpnt*angular_points + apnt);
	}
	Eigen::VectorXd get_gradient(int atom, int rpnt, int apnt){
		static double index;
		index = rpnt*angular_points + apnt;
		return gridpoints[atom].PBE.col(index);
	}
	Eigen::VectorXd get_gradient(int atom, int rpnt, int apnt, bool ab){
		Eigen::VectorXd id = Eigen::VectorXd::Zero(4);
		static int index;
		index = rpnt*angular_points + apnt;
		if(ab){
			id << gridpoints[atom].PBEC(0,index),gridpoints[atom].PBEC(1,index),gridpoints[atom].PBEC(2,index),gridpoints[atom].PBEC(3,index); //get alpha density
		}
		if(!ab){
			id << gridpoints[atom].PBEC(4,index),gridpoints[atom].PBEC(5,index),gridpoints[atom].PBEC(6,index),gridpoints[atom].PBEC(7,index); //get beta density
		}
		return id;
	}
	//get density Hessian
	//{element} = {0,1,2} <> {x,y,z}
	//e.g. {element1 = 0, element2 = 1} <> {xy}
	double get_hessian(int atom, int rpnt, int apnt, int element1, int element2){
		static double index;
		index = rpnt*angular_points + apnt;
		return gridpoints[atom].rho_Hess[index].get_Hess(element1, element2);
	}
	//get alpha or beta density Hessian
	//{element} = {0,1,2,3} <> {x,y,z,norm}
	double get_hessian(int atom, int rpnt, int apnt, int element1, int element2, bool ab){ //function overload
		static double index;
		index = rpnt*angular_points + apnt;
		if(ab){
			return gridpoints[atom].alpha_Hess[index].get_Hess(element1, element2); //get alpha density
		}
		if(!ab){
			return gridpoints[atom].beta_Hess[index].get_Hess(element1, element2); //get beta density
		}
	}
	Eigen::MatrixXd get_hessian(int atom, int rpnt, int apnt){ //function overload
		static double index;
		index = rpnt*angular_points + apnt;
		return gridpoints[atom].rho_Hess[index].get_Hess();
	}
	Eigen::MatrixXd get_hessian(int atom, int rpnt, int apnt, bool ab){ //function overload
		static double index;
		index = rpnt*angular_points + apnt;
		if(ab){
			return gridpoints[atom].alpha_Hess[index].get_Hess(); //get alpha density
		}
		if(!ab){
			return gridpoints[atom].beta_Hess[index].get_Hess(); //get beta density
		}
	}

	//*******************************************
	//Poisson potential:

	//compute the Laplace coefficients
	void set_Laplace();

	//confirm proper spherical harmonic decomposition according to identity
	void PC_check(int atomnr);

	//get index for mapping all {l,m} values to 1D
	//example: [(0,0), (1,-1), (1,0), (1,1), (2,-2), (2,-1), (2,0), (2,1), (2,2), (3,-3), ...]
	int lm_index(int l, int m){
		static int sum;
		sum = 0; //reset sum
		for(int i=0; i<l+1; i++){
			sum += 2*i+1; //shift the index to skip all sets {l,m} for l < current l
		}
		return sum + m - 1 - l; //shift the index to match current value of m
	}
	//compute both local and external Poisson potentials
	void set_Poisson(){
		set_Laplace(); //compute the Laplace coefficients
		set_Poisson_L(); //set local Poisson potential
		set_Poisson_E(); //set external Poisson potential
	}
	//retrieve value of the Poisson potential at the point specified by atom, rpnt
	double get_Poisson(int atomnr, int rpnt, int apnt){
		return gridpoints[atomnr].potential(0,rpnt*angular_points + apnt) + gridpoints[atomnr].potential(1,rpnt*angular_points + apnt); //return sum of local and external potentials
	}
	//compute the local Poisson potential at all points using finite differences
	void set_Poisson_L();

	//compute the external Poisson potential at all points using spline interpolation
	void set_Poisson_E();

	//generate the heptadiagonal matrix
	//{l} is the degree of the associated real spherical harmonic
	void set_hepta(int l, int atomnr);

	//calculate the non-constant ODE coefficients
	double hepta_coeff(int tag, double z, double Slater);

	//update the harmonic dependency of the heptadiagonal matrix only
	//when changing l, only a small part of U7 changes
	//update only these terms -> diagonals
	void update_hepta(int l, int atomnr);

	//solve the matrix equation for the transformed Poisson potential
	//{U7} is the heptadiagonal matrix from the modified Bickley finite difference scheme with dimensions rpnt x rpnt (15x15)
	//{rho} is the vector of Laplace coefficients for a specific {l,m} at each radial point, so a vector[rpnt] | {l,m}
	//
	void Ulm_solve(int atom, int index);

	//compute the contribution to the external Poisson potential through cubic B-spline interpolation
	double Poisson_spline(int atomnr, double distance, int lm);
	//alternative: linear interpolation
	double lazy_spline(int atomnr, double distance, int lm);

};

//define function calls
grid make_grid(const vector<vec3>& pos_list, const vector<int>& atnum_list, int Basis);
double Gauss_weight(int atomnr, const vector<vec3>& pos_list, const vector<int>& atnum_list, const vec3& point);
double fuzzy_cutoff(double mu);
double Voronoi_Hetero(int atnum1, int atnum2, double mu);
double SlaterBragg(int atnum);

#endif //_QUADRATURE_H