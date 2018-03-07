/*****************************************************************

Basic closed-shell spin-restricted DFT-solver for simple molecules using STO-NG

Authors: B. Klumpers
		 I.A.W. Filot

Published under GNU General Public License 3.0
Source code available at: https://github.com/BKlumpers/dft

Allows for SCF-computation of molecular energies for simple molecules.
Includes testcases for: H, He, H2, HeH+, He2, CO, and H2O.

*****************************************************************/

#ifndef _QUADRATURE_H
#define _QUADRATURE_H

//Poisson potential should be per unique 2e-integral ---> how to do the storage?   ----> NOT NEEDED! alleviate {k,l} summation through commutivity of densities

#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "lebedev.h"
#include "orbitals.h"
#include "Harmonics.h"

typedef Eigen::Vector3d vec3; //define vec3 as the Eigen::Vector3d-object

using namespace std;

//object for atomic integration grid
//{gridpoints} matrix[3,m] with each column m representing a gridpoint with its {x,y,z}-coordinates stored in the 3 rows
//{amps} matrix[N+1,m] with each row containing the amplitudes of 1 basis function at each of the aforementioned m gridpoints
//-------------------->the N+1th row is used for storing the total local electron density as computed by grid::set_density
//{Laplace} contains the Laplace coefficients for the local Poisson potential; matrix[rpnts, lmttl], each value of {l,m} has values at each radial point
//{Ulm} contains the radial potential function, hence dimensions are equal to Laplace
//{potential} contains the local and external Poisson potential at each radial point; matrix[2, radial_points], 1st row local, 2nd row external
//{lmttl} is the total number of real spherical harmonics in the basis
//{lebedev_order} is the order of the Lebedev quadrature, which relates to the highest angular momentum in the harmonics expansion as: lmax = lebedev_order/2
class atomicgrid {
public:
	//int lmttl, lmax;
	Eigen::MatrixXd gridpoints;
	Eigen::MatrixXd amps;
	Eigen::MatrixXd Laplace;
	Eigen::MatrixXd Ulm;
	Eigen::MatrixXd potential;
	//class constructor
	atomicgrid(int radial_points, int angular_points, int Basis, int lmttl){
		gridpoints = Eigen::MatrixXd::Zero(3, radial_points*angular_points); //create empty atomic grid
		amps = Eigen::MatrixXd::Zero(Basis+1, radial_points*angular_points); //create empty amplitude and density maps
		//Poisson:
		// lmax = Lebedev_order + 1; //store this parameter, since Lebedev::lebedev_num = n-1 gives Lebedev_order 2n+1, +1 is to correct for n starting at 0
		// int lmttl = 0;
		// for(int i=0; i<(1+lmax); i++){ //lmax = lebedev_order/2;
		// 	lmttl += 2*i+1; //add the number of m-values for degree l, sum over all degrees to get the total number of harmonics
		// }
		//make Laplace contain first index for each l,m couple; second index along radii: atomnr*radial + rpnt
		Laplace = Eigen::MatrixXd::Zero(radial_points, lmttl);
		//make potential contain internal + external potentials, for each gridpoint
		potential = Eigen::MatrixXd::Zero(2, radial_points*angular_points);
		//Ulm has the same dimensions as Laplace
		Ulm = Laplace; //same initialisation
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
class grid {
private:
	static constexpr double pi = 3.14159265359;
	int lmttl, lmax;
	int radial_points, angular_points; //store number of radial and angular gridpoints respectively
	Eigen::MatrixXd weights; //Matrix containing the weights, first index is the nucleic label, second is the gridpoint-label
	vector<atomicgrid> gridpoints; //vector containing the atomic grid, index runs across the nucleic labels
	Eigen::MatrixXd U7; //heptadiagonal matrix storage for computation of the radial Poisson potential
	vector<vec3> pos_list; //need the nucleic coordinates for spline interpolation of the external Poisson potentials
public:
	//initialise object with predefined gridsizes and known number of atoms; //----------> pos_list to replace int atoms
	grid(vector<vec3> Pos_list, int radial, int angular, int Basis, int Lebedev_order) {
		radial_points = radial; //store number of radial gridpoints inside the class
		angular_points = angular; //store number of angular gridpoints inside the class
		pos_list = Pos_list;
		//Poisson:
		lmax = Lebedev_order + 1; //store this parameter, since Lebedev::lebedev_num = n-1 gives Lebedev_order 2n+1, +1 is to correct for n starting at 0
		int lmttl = 0;
		for(int i=0; i<(1+lmax); i++){ //lmax = lebedev_order/2;
			lmttl += 2*i+1; //add the number of m-values for degree l, sum over all degrees to get the total number of harmonics
		}
		atomicgrid emptygrid(radial_points, angular_points, Basis, lmttl); //create empty atomic grid
		weights = Eigen::MatrixXd::Zero(pos_list.size(), radial_points*angular_points); //create empty matrix for the weights
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
		weights(atomnr, index) *= weight;
	}
	//write all the wave-function amplitudes to atomicgrid
	void write_amps(vector<CGF> AO_list){
		for(int i=0; i<gridpoints.size(); i++){ //loop across all atomic grids
			for(int rpnt=0; rpnt<radial_points; rpnt++){ //loop across all radial points
				for(int apnt=0; apnt<angular_points; apnt++){ //loop across all angular points
					static int index;
					index = rpnt*angular_points + apnt;
					for(int j=0; j<AO_list.size(); j++){ //loop across all basis functions
						gridpoints[i].amps(j,index) = AO_list[j].getvalue(get_gridpoint(i,rpnt,apnt));
					}
				}
			}
		}
	}
	//get specific wavefunction amplitude
	double get_amps(int atomnr, int radial, int angular, int Basis){
		return gridpoints[atomnr].amps(Basis,radial*angular_points+angular);
	}
	//write values for the electron density
	void set_density(const Eigen::MatrixXd& Pmatrix){
		//write all the densities by calling get_amps
		static double density_value;
		//loop across all gridpoints and all basis functions
		for(int atomnr=0; atomnr<gridpoints.size(); atomnr++){ //loop across all atomic grids
			for(int rpnt=0; rpnt<radial_points; rpnt++){ //loop across all radial points
				for(int apnt=0; apnt<angular_points; apnt++){ //loop across all angular points
					density_value = 0.0; //reset density storage
					//loop across all basis functions
					for(int i=0; i<Pmatrix.rows(); i++){
						for(int j=0; j<Pmatrix.rows(); j++){
							density_value += Pmatrix(i,j)*get_amps(atomnr,rpnt,apnt,i)*get_amps(atomnr,rpnt,apnt,j);
						}
					}
					gridpoints[atomnr].amps(gridpoints[atomnr].amps.rows()-1, rpnt*angular_points+apnt) = density_value; //write total local electron density to atomic grid
				}
			}
		}
	}
	//retrieve value of the density at specified gridpoint
	double get_density(int atomnr, int radial, int angular){
		return gridpoints[atomnr].amps(gridpoints[atomnr].amps.rows()-1, radial*angular_points + angular);
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

	//*******************************************
	//Poisson potential:

	//compute the Laplace coefficients
	void set_Laplace(){
		static const double omega = 4.0*pi/(angular_points);

		cout << "confirm static omega: " << omega << endl;

		//static vec3 point;
		static double delta;
		static int index;
		//for each l,m:
		for(int atom=0; atom<gridpoints.size(); atom++){ //compute them for each atomic grid
			for(int rpnt=0; rpnt<radial_points; rpnt++){ //compute each radial point
				delta = omega * (get_gridpoint(atom, rpnt, 0) - pos_list[atom]).squaredNorm(); //r^2 * omega; invariant for apnt since all apnt are positioned on a sphere with set radius
				//compute all {l,m} components
				for(int l=0; l<lmax+1; l++){
					for(int m=-l; m<l+1; m++){
						index = lm_index(l,m);
						gridpoints[atom].Laplace(rpnt, index) = 0.0; //reset coefficient
						//sum across all angular points:
						for(int apnt=0; apnt<angular_points; apnt++){
							//point << (get_gridpoint(atom, rpnt, 0) - pos_list[atom]); //shift coordinate system to position centre of the atomic grid at the origin
							gridpoints[atom].Laplace(rpnt, index) += get_weight(atom, rpnt, apnt)*get_density(atom, rpnt, apnt) * HarmonicReal(l, m, get_gridpoint(atom, rpnt, apnt) - pos_list[atom]) * delta;
						}
					}
				}
			}
		}
	}
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
	//compute the local Poisson potential at all points using finite differences
	void set_Poisson_L(){

		//reset potentials

		static int index;
		static double radius;
		static vec3 point;
		set_hepta(0); //move inside the atomic loop if scaling of r through rm is implemented since rm = f(atom)
		for(int atom=0; atom<gridpoints.size(); atom++){ //evaluate the potential on each atomic grid
			//for each rpnt, sum across all {l,m}
			for(int l=0; l<lmax+1; l++){
				update_hepta(l);
				for(int m=-l; m<l+1; m++){
					index = lm_index(l,m);
					//get a vector of rho_lm(r), create a matrix of these vectors for each m value given the current l
					//call solver for all of these?
					//->bypass U00 by adjusting for BC in this function already

					//Ulm_solve(atom, rho_lm); //fill gridpoints[atom].Ulm given Laplace + Boundary conditions in rho_lm -> or simply call gridpoint[atom].Laplace and only parse BC
					for(int rpnt=0; rpnt<radial_points; rpnt++){
						point << (get_gridpoint(atom, rpnt, 0) - pos_list[atom]);
						radius = point.norm(); //radius is invariant under apnt
						for(int apnt=0; apnt<angular_points; apnt++){
							gridpoints[atom].potential(1,rpnt*angular_points+apnt) += gridpoints[atom].Ulm(rpnt, index) * HarmonicReal(l, m, get_gridpoint(atom, rpnt, apnt) - pos_list[atom]) / radius;
						}
					}
				}
			}
			//
			//for given l, construct U7
			//simultaneously solve for all values of m, since U7 is the same for each (only depends on l)
			//do recalculate U7 for each atom since formally the grids have different scaling for different atoms, hence values of r in U7 are different as well
		}
	}
	//compute the external Poisson potential at all points using spline interpolation
	void set_Poisson_Eeee(){
		static double distance;

		for(int atom1=0; atom1<gridpoints.size(); atom1++){
			for(int rpnt=0; rpnt<radial_points; rpnt++){ //loop across all the radial points of atom1 to get external potential at each of these points
				gridpoints[atom1].potential(2,rpnt) = 0.0; //reset external potential
				for(int atom2=0; atom2<gridpoints.size(); atom2++){
					//compute the external potential for all other atoms
					if(atom1 != atom2){
						//total distance between nuclei - distance between gridpoint and atom1 = distance between gridpoint and atom2
						distance = (pos_list[atom1] - pos_list[atom2]).norm() - (get_gridpoint(atom1,rpnt,0) - pos_list[atom1]).norm();
						//call cubic (2nd order polynomial) spline interpolation function to get the external potential due to atom2 at distance from its origin
						gridpoints[atom1].potential(2,rpnt) += spline(atom2, distance);
					}
				}
			}
		}
	}
	//generate the heptadiagonal matrix
	//{l} is the degree of the associated real spherical harmonic
	void set_hepta(int l){
		static const double Np1 = 1.0 / double(radial_points + 1);
		static const double jB = pow(Np1, 2); //Bickley coefficient for the 2nd derivative
		static const double i1 = 1.0 / (60.0 * Np1);
		static const double i2 = 1.0 / (180.0 * jB);
		static double z;
		static double r;
		static double L;
		L = double(l*(l+1)); 

		cout << "confirm static Np1: " << Np1 << " for " << (radial_points+1) << endl;

		for(int row=3; row<radial_points-3; row++){ //run from 4th row till N-4th row
			//write coefficients of 7 elements
			z = Np1 * (row-2);
			U7(row,row-3) = 2.0 * hepta_coeff(1, z) * i2 - hepta_coeff(2, z) * i1;
			z = Np1 * (row-1);
			U7(row,row-2) = -27.0 * hepta_coeff(1, z) * i2 + 9.0 * hepta_coeff(2, z) * i1;
			z = Np1 * row;
			U7(row,row-1) = 270.0 * hepta_coeff(1, z) * i2 - 45.0 * hepta_coeff(2, z) * i1;
			z = Np1 * (row+1);
			r = (1.0 + cos(pi*z)) / (1 - cos(pi*z));
			U7(row,row) = -490.0 * hepta_coeff(1, z) * i2 - L / pow(r,3);
			z = Np1 * (row+2);
			U7(row,row+1) = 270.0 * hepta_coeff(1, z) * i2 + 45.0 * hepta_coeff(2, z) * i1;
			z = Np1 * (row+3);
			U7(row,row+2) = -27.0 * hepta_coeff(1, z) * i2 - 9.0 * hepta_coeff(2, z) * i1;
			z = Np1 * (row+4);
			U7(row,row+3) = 2.0 * hepta_coeff(1, z) * i2 + hepta_coeff(2, z) * i1;
		}

		//do rows 1,2 and N,N-1 separately as they use different finite difference expansions, 3 and N-2 to include the boundary conditions
		
		//row 1
		z = Np1;
		r = (1.0 + cos(pi*z)) / (1 - cos(pi*z));
		U7(0,0) = (-26.0/3.0) * jB * hepta_coeff(1, z) + 4.0 * Np1 * hepta_coeff(2, z) - L / pow(r,3);
		z = 2.0 * Np1;
		U7(0,1) = 9.5 * jB * hepta_coeff(1, z) - 3.0 * hepta_coeff(2, z) * Np1;
		z = 3.0 * Np1;
		U7(0,2) = (-14.0/3.0) * jB * hepta_coeff(1, z) + (4.0/3.0) * hepta_coeff(2, z) * Np1;
		z = 4.0 * Np1;
		U7(0,3) = (11.0/12.0) * jB * hepta_coeff(1, z) - 0.25 * hepta_coeff(2, z) * Np1;

		//row 2
		z = 1.0 * Np1;
		U7(1,0) = (-77.0/6.0) * jB * hepta_coeff(1, z) + 5.0 * Np1 * hepta_coeff(2, z);
		z = 2.0 * Np1;
		r = (1.0 + cos(pi*z)) / (1 - cos(pi*z));
		U7(1,1) = (107.0/6.0) * jB * hepta_coeff(1, z) - 5.0 * hepta_coeff(2, z) * Np1 - L / pow(r,3);
		z = 3.0 * Np1;
		U7(1,2) = -13.0 * jB * hepta_coeff(1, z) + (10.0/3.0) * Np1 * hepta_coeff(2, z);
		z = 4.0 * Np1;
		U7(1,3) = (61.0/12.0) * jB * hepta_coeff(1, z) - 1.25 * Np1 * hepta_coeff(2, z);
		z = 5.0 * Np1;
		U7(1,4) = (-5.0/6.0) * jB * hepta_coeff(1, z) + 0.2 * Np1 * hepta_coeff(2, z);

		//row 3
		z = 1.0 * Np1;
		U7(2,0) = -27.0 * hepta_coeff(1, z) * i2 + 9.0 * hepta_coeff(2, z) * i1;
		z = 2.0 * Np1;
		U7(2,1) = 270.0 * hepta_coeff(1, z) * i2 - 45.0 * hepta_coeff(2, z) * i1;
		z = 3.0 * Np1;
		r = (1.0 + cos(pi*z)) / (1 - cos(pi*z));
		U7(2,2) = -490.0 * hepta_coeff(1, z) * i2 - L / pow(r,3);
		z = 4.0 * Np1;
		U7(2,3) = 270.0 * hepta_coeff(1, z) * i2 + 45.0 * hepta_coeff(2, z) * i1;
		z = 5.0 * Np1;
		U7(2,4) = -27.0 * hepta_coeff(1, z) * i2 - 9.0 * hepta_coeff(2, z) * i1;
		z = 6.0 * Np1;
		U7(2,5) = 2.0 * hepta_coeff(1, z) * i2 + hepta_coeff(2, z) * i1;

		//row N-2
		z = double(radial_points - 5) / Np1;
		U7(radial_points - 3,radial_points-6) = 2.0 * hepta_coeff(1, z) * i2 - hepta_coeff(2, z) * i1;
		z = double(radial_points - 4) / Np1;
		U7(radial_points - 3,radial_points-5) = -27.0 * hepta_coeff(1, z) * i2 + 9.0 * hepta_coeff(2, z) * i1;
		z = double(radial_points - 3) / Np1;
		U7(radial_points - 3,radial_points-4) = 270.0 * hepta_coeff(1, z) * i2 - 45.0 * hepta_coeff(2, z) * i1;
		z = double(radial_points - 2) / Np1;
		r = (1.0 + cos(pi*z)) / (1 - cos(pi*z));
		U7(radial_points - 3,radial_points-3) = -490.0 * hepta_coeff(1, z) * i2 - L / pow(r,3);
		z = double(radial_points - 1) / Np1;
		U7(radial_points - 3,radial_points-2) = 270.0 * hepta_coeff(1, z) * i2 + 45.0 * hepta_coeff(2, z) * i1;
		z = double(radial_points) / Np1;
		U7(radial_points - 3,radial_points-1) = -27.0 * hepta_coeff(1, z) * i2 - 9.0 * hepta_coeff(2, z) * i1;

		//row N-1
		z = double(radial_points - 4) / Np1;
		U7(radial_points - 2,radial_points-5) = (-5.0/6.0) * jB * hepta_coeff(1, z) - 0.2 * Np1 * hepta_coeff(2, z);
		z = double(radial_points - 3) / Np1;
		U7(radial_points - 2,radial_points-4) = (61.0/12.0) * jB * hepta_coeff(1, z) + 1.25 * Np1 * hepta_coeff(2, z);
		z = double(radial_points - 2) / Np1;
		U7(radial_points - 2,radial_points-3) = -13.0 * jB * hepta_coeff(1, z) - (10.0/3.0) * Np1 * hepta_coeff(2, z);
		z = double(radial_points - 1) / Np1;
		r = (1.0 + cos(pi*z)) / (1 - cos(pi*z));
		U7(radial_points - 2,radial_points-2) = (107.0/6.0) * jB * hepta_coeff(1, z) + 5.0 * hepta_coeff(2, z) * Np1 - L / pow(r,3);
		z = double(radial_points) / Np1;
		U7(radial_points - 2,radial_points-1) = (-77.0/6.0) * jB * hepta_coeff(1, z) - 5.0 * Np1 * hepta_coeff(2, z);

		//row N
		z = double(radial_points - 3) / Np1;
		U7(radial_points - 1,radial_points - 4) = (11.0/12.0) * jB * hepta_coeff(1, z) + 0.25 * hepta_coeff(2, z) * Np1;
		z = double(radial_points - 2) / Np1;
		U7(radial_points - 1,radial_points - 3) = (-14.0/3.0) * jB * hepta_coeff(1, z) - (4.0/3.0) * hepta_coeff(2, z) * Np1;
		z = double(radial_points - 1) / Np1;
		U7(radial_points - 1,radial_points - 2) = 9.5 * jB * hepta_coeff(1, z) + 3.0 * hepta_coeff(2, z) * Np1;
		z = double(radial_points) / Np1;
		r = (1.0 + cos(pi*z)) / (1 - cos(pi*z));
		U7(radial_points - 1,radial_points - 1) = (-26.0/3.0) * jB * hepta_coeff(1, z) - 4.0 * Np1 * hepta_coeff(2, z) - L / pow(r,3);
	}
	//calculate the non-constant ODE coefficients
	double hepta_coeff(int tag, double z){
		static double r;
		r = (1.0 + cos(pi*z)) / (1 - cos(pi*z)); //transform z back to r
		switch(tag){
			case 1:
			return pow(hepta_aux(1,z), 2) / r;
			break;
			case 2:
			return hepta_aux(1,z) * hepta_aux(2,z) / r;
			break;
			// case 3:
			// return double(l*(l+1)) / pow(r,3);
			// break;
			default:
			cout << "Unknown tag for ODE-coefficient function: tag = " << tag << endl;
			cout << "Defaulting to: tag = 1 --> [ dz/dr ]" << endl;
			return hepta_coeff(1, z);
			break;
		}
	}
	//auxiliary functions for the ODE coefficients
	double hepta_aux(int tag, double z){
		switch(tag){
			case 1: //function a(z) = dz/dr
			return pow(1.0 - cos(pi*z), 2) / (-2.0*pi*sin(pi*z));
			break;
			case 2: //function b(z) = d2z/dr2 / dz/dr
			return (1.0 - cos(pi*z)) * (1.0 + 0.5*cos(pi*z)/pow(sin(pi*z),2));
			break;
			default:
			cout << "Unknown tag for ODE-coefficient auxiliary function: tag = " << tag << endl;
			cout << "Defaulting to: tag = 1 --> [ a(z) ]" << endl;
			return hepta_aux(2, z);
			break;
		}
	}
	//update the harmonic dependency of the heptadiagonal matrix only
	//when changing l, only a small part of U7 changes
	//updat only these terms -> diagonals
	void update_hepta(int l){
		static const double Np1 = 1.0 / double(radial_points + 1);
		static const double jB = pow(Np1, 2); //Bickley coefficient for the 2nd derivative
		static const double i2 = 1.0 / (180.0 * jB);
		static double z;
		static double r;
		static double L;
		L = double(l*(l+1));

		for(int row=2; row<radial_points-2; row++){
			z = Np1 * (row+1);
			r = (1.0 + cos(pi*z)) / (1 - cos(pi*z));
			U7(row,row) = -490.0 * hepta_coeff(1, z) * i2 - L / pow(r,3);
		}
		//handle different finite difference schemes for first and last 2 rows:
		//row 1
		z = Np1;
		r = (1.0 + cos(pi*z)) / (1 - cos(pi*z));
		U7(0,0) = (-26.0/3.0) * jB * hepta_coeff(1, z) + 4.0 * Np1 * hepta_coeff(2, z) - L / pow(r,3);
		//row 2
		z = 2.0 * Np1;
		r = (1.0 + cos(pi*z)) / (1 - cos(pi*z));
		U7(1,1) = (107.0/6.0) * jB * hepta_coeff(1, z) - 5.0 * hepta_coeff(2, z) * Np1 - L / pow(r,3);
		//row N-1
		z = double(radial_points - 1) / Np1;
		r = (1.0 + cos(pi*z)) / (1 - cos(pi*z));
		U7(radial_points - 2,radial_points-2) = (107.0/6.0) * jB * hepta_coeff(1, z) + 5.0 * hepta_coeff(2, z) * Np1 - L / pow(r,3);
		//row N
		z = double(radial_points) / Np1;
		r = (1.0 + cos(pi*z)) / (1 - cos(pi*z));
		U7(radial_points - 1,radial_points - 1) = (-26.0/3.0) * jB * hepta_coeff(1, z) - 4.0 * Np1 * hepta_coeff(2, z) - L / pow(r,3);
	}
	//solve the matrix equation for the transformed Poisson potential
	//{U7} is the heptadiagonal matrix from the modified Bickley finite difference scheme with dimensions rpnt x rpnt (15x15)
	//{rho} is the vector of Laplace coefficients for a specific {l,m} at each radial point, so a vector[rpnt] | {l,m}
	//
	//scheme will be made more efficient through simultaneously solving the equation for fixed l, for each m, since U7 != f(m)
	//
	void Ulm_solve(int atom, Eigen::VectorXd rho){
		//solve matrix equation: U7 * Ui = -4pi * rho
		static const double pre = -4.0 * pi;
		static Eigen::VectorXd Ui;
		Ui = U7.colPivHouseholderQr().solve(pre*rho);
		//return Ui; ---> write to atomicgrid.Ulm
	}
	//seperate solver for Ulm | {l,m} = {0,0} --------> merge with Ulm_solve by incorporating BC in set_Poisson_L
	void U00_solve(int atom, Eigen::VectorXd rho){
		//solve matrix equation: U7 * Ui = -4pi * rho
		static const double pre = -4.0 * pi;
		static Eigen::VectorXd Ui;
		static Eigen::VectorXd boundary; //implement correction for boundary conditions by modifying rho

		Ui = U7.colPivHouseholderQr().solve(pre*rho);
		//return Ui; ---> write to atomicgrid.Ulm
	}
	//compute the contribution to the external Poisson potential through cubic spline interpolation
	double spline(int atomnr, double distance){
		static double Poisson_Ext;
		Poisson_Ext = 0.0;
		//interpolate the Poisson potential to a point at {distance} away from the nucleus indexed by {atomnr}
		// #include <boost/math/interpolators/cubic_b_spline.hpp>
		return Poisson_Ext;
	}
};

//define function calls
grid make_grid(const vector<vec3>& pos_list, const vector<int>& atnum_list, int Basis);
double Gauss_weight(int atomnr, const vector<vec3>& pos_list, const vector<int>& atnum_list, const vec3& point);
double fuzzy_cutoff(double mu);
double Voronoi_Hetero(int atnum1, int atnum2, double mu);
double SlaterBragg(int atnum);

#endif //_QUADRATURE_H