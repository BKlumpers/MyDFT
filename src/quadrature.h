/*****************************************************************

Basic closed-shell spin-restricted DFT-solver for simple molecules using STO-NG

Authors: B. Klumpers (bartkl@live.nl)
		 I.A.W. Filot

Published under GNU General Public License 3.0

Allows for SCF-computation of molecular energies for simple molecules.
Testcases for H, He, H2, HeH+, He2 and CO are included.

*****************************************************************/

#ifndef _QUADRATURE_H
#define _QUADRATURE_H

#include <iostream>
#include <cmath>
#include <vector>
#include <Eigen/Core>
#include "lebedev.h"
#include "orbitals.h"

typedef Eigen::Vector3d vec3; //define vec3 as the Eigen::Vector3d-object

using namespace std;

//object for atomic integration grid
//{gridpoints} matrix[3,m] with each column m representing a gridpoint with its {x,y,z}-coordinates stored in the 3 rows
//{amps} matrix[N,m] with each row containing the amplitudes of 1 basis function at each of the aforementioned m gridpoints
class atomicgrid {
public:
	Eigen::MatrixXd gridpoints;
	Eigen::MatrixXd amps;
	//class constructor
	atomicgrid(int radial_points, int angular_points, int Basis){
		gridpoints = Eigen::MatrixXd::Zero(3, radial_points*angular_points); //create empty atomic grid
		amps = Eigen::MatrixXd::Zero(Basis, radial_points*angular_points); //create empty amplitude map
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
	int radial_points, angular_points; //store number of radial and angular gridpoints respectively
	Eigen::MatrixXd weights; //Matrix containing the weights, first index is the nucleic label, second is the gridpoint-label
	vector<atomicgrid> gridpoints; //vector containing the atomic grid, index runs across the nucleic labels
public:
	//initialise object with predefined gridsizes and known number of atoms
	grid(int atoms, int radial, int angular, int Basis) {
		radial_points = radial; //store number of radial gridpoints inside the class
		angular_points = angular; //store number of angular gridpoints inside the class
		static const atomicgrid emptygrid(radial_points, angular_points, Basis); //create empty atomic grid
		weights = Eigen::MatrixXd::Zero(atoms, radial_points*angular_points); //create empty matrix for the weights
		for(int i=0; i<atoms; i++){
			gridpoints.push_back(emptygrid); //initialise an empty atomic grid for each atom
		}
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
};

//define function calls
grid make_grid(const vector<vec3>& pos_list, const vector<int>& atnum_list, int Basis);
double Gauss_weight(int atomnr, const vector<vec3>& pos_list, const vector<int>& atnum_list, const vec3& point);
double fuzzy_cutoff(double mu);
double Voronoi_Hetero(int atnum1, int atnum2, double mu);
double SlaterBragg(int atnum);

#endif //_QUADRATURE_H