/*****************************************************************

Basic closed-shell spin-restricted DFT-solver for simple molecules using STO-NG

Authors: B. Klumpers
		 I.A.W. Filot

Published under GNU General Public License 3.0
Source code available at: https://github.com/BKlumpers/dft

Allows for SCF-computation of molecular energies for simple molecules.
Includes testcases for: H, He, H2, HeH+, He2, CO, and H2O.

*****************************************************************/

#include "quadrature.h"


//****************************************************************
//begin main codeblock for grid::member functions


//write all the wave-function amplitudes to atomicgrid
void grid::write_amps(vector<CGF> AO_list){
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

//write values for the electron density
void grid::set_density(const Eigen::MatrixXd& Pmatrix){
	//write all the densities by calling get_amps
	static double density_value;
	//loop across all gridpoints and all basis functions
	for(int atomnr=0; atomnr<gridpoints.size(); atomnr++){ //loop across all atomic grids
		gridpoints[atomnr].point_charge = 0.0; //reset magnitude of point-charge potential
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
				gridpoints[atomnr].point_charge += get_weight(atomnr,rpnt,apnt) * density_value;
			}
		}
	}
	//set_Poisson(); //update the Poisson potentials as well
}

//*******************************************
//Poisson potential:

//compute the Laplace coefficients
void grid::set_Laplace(){
	static const double omega = 4.0*pi;//  /(angular_points);

	//replace 1/apnt in omega with appropriate lebedev coefficients --> store Lebedev offset from make_grid

	//static vec3 point;
	static double delta;
	static int index;
	//for each l,m:
	for(int atom=0; atom<gridpoints.size(); atom++){ //compute them for each atomic grid
		for(int rpnt=0; rpnt<radial_points; rpnt++){ //compute each radial point
			//delta = omega * (get_gridpoint(atom, rpnt, 0) - pos_list[atom]).squaredNorm(); //r^2 * omega; invariant for apnt since all apnt are positioned on a sphere with set radius
			//compute all {l,m} components
			for(int l=0; l<lmax+1; l++){
				for(int m=-l; m<l+1; m++){
					index = lm_index(l,m);
					gridpoints[atom].Laplace(rpnt, index) = 0.0; //reset coefficient
					//sum across all angular points:
					for(int apnt=0; apnt<angular_points; apnt++){
						//point << (get_gridpoint(atom, rpnt, 0) - pos_list[atom]); //shift coordinate system to position centre of the atomic grid at the origin
						gridpoints[atom].Laplace(rpnt, index) += get_weight(atom, rpnt, apnt)*get_density(atom, rpnt, apnt) * HarmonicReal(l, m, get_gridpoint(atom, rpnt, apnt) - pos_list[atom]) * omega * Lebedev::lebedev_coefficients[apnt+angular_index][3];//delta;
						//cout << "harmonictest: " << HarmonicReal(l, m, get_gridpoint(atom, rpnt, apnt) - pos_list[atom]) << " at " << gridpoints[atom].Laplace(rpnt, index) << endl;
					}
				}
			}					//	get_weight(atom, rpnt, apnt)* ->>should be part of above formula
		}
	}
	//std::ofstream file("Laplace.txt"); //set file as output stream
	//file << gridpoints[0].Laplace << '\n'; //write image to file
	//cout << "rho_lm[0]: " << endl << gridpoints[0].Laplace << endl;
	//cout << "rho_lm[1]: " << endl << gridpoints[1].Laplace << endl;
}
//confirm proper spherical harmonic decomposition according to identity
void grid::PC_check(int atomnr){
	static int index;
	double PC = 0.0;
	for(int rpnt=0; rpnt<radial_points; rpnt++){
		for(int apnt=0; apnt<angular_points; apnt++){
			for(int l=0; l<lmax+1; l++){
				for(int m=-l; m<l+1; m++){
					index = lm_index(l,m);
					PC += gridpoints[atomnr].Laplace(rpnt, index) * HarmonicReal(l, m, get_gridpoint(atomnr, rpnt, apnt) - pos_list[atomnr]);
				}
			}
		}
	}
	cout << "point charge for atom " << atomnr+1 << " = " << PC << endl;
}

//compute the local Poisson potential at all points using finite differences
void grid::set_Poisson_L(){

	static int index;
	static double radius;

	//do recalculate U7 for each atom since formally the grids have different scaling for different atoms, hence values of r in U7 are different as well
	set_hepta(0); //move inside the atomic loop if scaling of r through rm is implemented since rm = f(atom)

	for(int atom=0; atom<gridpoints.size(); atom++){ //evaluate the potential on each atomic grid
		gridpoints[atom].potential.row(0) = Eigen::VectorXd::Zero(gridpoints[atom].potential.cols()); //reset potential
		//for each rpnt, sum across all {l,m}
		for(int l=0; l<lmax+1; l++){
			update_hepta(l);			//for given l, construct U7
			for(int m=-l; m<l+1; m++){
				index = lm_index(l,m);
				Ulm_solve(atom, index); //fill gridpoints[atom].Ulm given Laplace + Boundary conditions


				//add loop for external here, no need to parse Ulm since it is private
				//loop over all other atoms and interpolate over their atomic grids the radial potential
				//then multiply this with the HarmonicReal computed for their Cartesian coord. - atompos(atom2) to get the Ylm term -> finally divide by distance
				//---> do all of this in Poisson_E -> call this function here
				//---> is this necessary since Ulm are all stored for each grid anyway...


				for(int rpnt=0; rpnt<radial_points; rpnt++){
					radius = gridpoints[atom].Ulm(rpnt, index) / ( (get_gridpoint(atom, rpnt, 0) - pos_list[atom]).norm() ); //radius is invariant under apnt
					for(int apnt=0; apnt<angular_points; apnt++){
						gridpoints[atom].potential(0,rpnt*angular_points+apnt) +=  radius * HarmonicReal(l, m, get_gridpoint(atom, rpnt, apnt) - pos_list[atom]);
					}
				}
			}
		}
	}
	//cout << "Ulm[0]: " << endl << gridpoints[0].Ulm << endl;
	//cout << "Ulm[1]: " << endl << gridpoints[1].Ulm << endl;
}
//compute the external Poisson potential at all points using spline interpolation
void grid::set_Poisson_E(){
	static double distance;
	static int lm;

	for(int atom1=0; atom1<gridpoints.size(); atom1++){
		for(int rpnt=0; rpnt<radial_points; rpnt++){ //loop across all the radial points of atom1 to get external potential at each of these points
			for(int apnt=0; apnt<angular_points; apnt++){
				gridpoints[atom1].potential(1,rpnt*angular_points+apnt) = 0.0; //reset external potential
				for(int atom2=0; atom2<gridpoints.size(); atom2++){
					//compute the external potential for all other atoms
					if(atom1 != atom2){
						//total distance between nuclei - distance between gridpoint and atom1 = distance between gridpoint and atom2
						//distance = abs((pos_list[atom1] - pos_list[atom2]).norm() - (get_gridpoint(atom1,rpnt,0) - pos_list[atom1]).norm()); //take abs as gridpoint may exceed internuclear distance
						distance = (get_gridpoint(atom1,rpnt,apnt) - pos_list[atom2]).norm(); //norm between gridpoint and origin of atom2
						//cout << "dist: " << (pos_list[atom1] - pos_list[atom2]).norm() << "  rad:  " << (get_gridpoint(atom1,rpnt,0) - pos_list[atom1]).norm() << endl;
						//call cubic (2nd order polynomial) spline interpolation function to get the external potential due to atom2 at distance from its origin
						for(int l=0; l<lmax+1; l++){
							for(int m=-l; m<l+1; m++){
								lm = lm_index(l,m);
								gridpoints[atom1].potential(1,rpnt*angular_points+apnt) += lazy_spline(atom2, distance, lm) * HarmonicReal(l,m,get_gridpoint(atom1,rpnt,apnt) - pos_list[atom2]) / distance;
							}
						}
						//gridpoints[atom1].potential(1,rpnt*angular_points+apnt) += Poisson_spline(atom2, distance);
						//gridpoints[atom1].potential(1,rpnt*angular_points+apnt) += lazy_spline(atom2, distance);
					}
				}
			}
		}
	}
}
void grid::set_hepta(int l){
	static const double Np1 = 1.0 / double(radial_points + 1); // = h
	static const double jB = Np1 * Np1; //Bickley coefficient for the 2nd derivative
	static const double i1 = 1.0 / (60.0 * Np1);
	static const double i2 = 1.0 / (180.0 * jB);
	static double z;
	static double r;
	static double L;
	L = double(l*(l+1));
	static double c1;
	static double c2;

	for(int row=3; row<radial_points-3; row++){ //run from 4th row till N-4th row
		//write coefficients of 7 elements
		z = Np1 * (row+1);
		r = (1.0 + cos(pi*z)) / (1 - cos(pi*z));
		c1 = hepta_coeff(1, z) * i2;
		c2 = hepta_coeff(2, z) * i1;
		U7(row,row-3) = 2.0 * c1 - c2;
		U7(row,row-2) = -27.0 * c1 + 9.0 * c2;
		U7(row,row-1) = 270.0 * c1 - 45.0 * c2;
		U7(row,row) = -490.0 * c1 - (L / (r*r*r));
		U7(row,row+1) = 270.0 * c1 + 45.0 * c2;
		U7(row,row+2) = -27.0 * c1 - 9.0 * c2;
		U7(row,row+3) = 2.0 * c1 + c2;
	}
	//do rows 1,2 and N,N-1 separately as they use different finite difference expansions, 3 and N-2 to include the boundary conditions

	//row 1
	z = Np1;
	r = (1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z) / (12.0 * jB);
	c2 = hepta_coeff(2, z) / (12.0 * Np1);
	U7(0,0) = -20.0 * c1 - 10.0 * c2 - (L / (r*r*r));
	U7(0,1) = 6.0 * c1 + 18.0 * c2;
	U7(0,2) = 4.0 * c1 - 6.0 * c2;
	U7(0,3) = -c1 + c2;

	//row 2
	z = 2.0 * Np1;
	r = (1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z) / (12.0 * jB);
	c2 = hepta_coeff(2, z) / (60.0 * Np1);
	U7(1,0) = 16.0* c1 - 30.0 * c2;
	U7(1,1) = -30.0 * c1 - 20.0 * c2 - (L / (r*r*r));
	U7(1,2) = 16.0 * c1 + 60 * c2;
	U7(1,3) = -1.0 * c1 - 15.0 * c2;
	U7(1,4) = 2.0 * c2;

	//row 3
	z = 3.0 * Np1;
	r = (1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z) * i2;
	c2 = hepta_coeff(2, z) * i1;
	U7(2,0) = -27.0 * c1 + 9.0 * c2;
	U7(2,1) = 270.0 * c1 - 45.0 * c2;
	U7(2,2) = -490.0 * c1 - (L / (r*r*r));
	U7(2,3) = 270.0 * c1 + 45.0 * c2;
	U7(2,4) = -27.0 * c1 - 9.0 * c2;
	U7(2,5) = 2.0 * c1 + c2;

	//row N-2
	z = double(radial_points - 2) * Np1;
	r = (1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z) * i2;
	c2 = hepta_coeff(2, z) * i1;
	U7(radial_points - 3,radial_points-6) = 2.0 * c1 - c2;
	U7(radial_points - 3,radial_points-5) = -27.0 * c1 + 9.0 * c2;
	U7(radial_points - 3,radial_points-4) = 270.0 * c1 - 45.0 * c2;
	U7(radial_points - 3,radial_points-3) = -490.0 * c1 - (L / (r*r*r));
	U7(radial_points - 3,radial_points-2) = 270.0 * c1 + 45.0 * c2;
	U7(radial_points - 3,radial_points-1) = -27.0 * c1 - 9.0 * c2;

	//row N-1
	z = double(radial_points - 1) * Np1;
	r = (1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z) / (12.0 * jB);
	c2 = hepta_coeff(2, z) / (60.0 * Np1);
	U7(radial_points - 2,radial_points-5) =  -2.0 * c2;
	U7(radial_points - 2,radial_points-4) = -c1 + 15.0* c2;
	U7(radial_points - 2,radial_points-3) = 16.0 * c1 - 60.0 * c2;
	U7(radial_points - 2,radial_points-2) = -30.0 * c1 + 20.0 * c2 - (L / (r*r*r));
	U7(radial_points - 2,radial_points-1) = 16.0 * c1 + 30.0 * c2;

	//row N
	z = double(radial_points) * Np1;
	r = (1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z) / (12.0 * jB);
	c2 = hepta_coeff(2, z) / (12.0 * Np1);
	U7(radial_points - 1,radial_points - 4) = -c1 - c2;
	U7(radial_points - 1,radial_points - 3) = 4.0* c1 + 6.0 * c2;
	U7(radial_points - 1,radial_points - 2) = 6.0 * c1 - 18.0 * c2;
	U7(radial_points - 1,radial_points - 1) = -20.0 * c1 + 10.0 * c2 - (L / (r*r*r));

	//cout << "U7: " << endl << U7 << endl;

}

//calculate the non-constant ODE coefficients
double grid::hepta_coeff(int tag, double z){
	static double r;
	static const double H_pre = 1.0/(2.0*pi);
	r = (1.0 + cos(pi*z)) / (1.0 - cos(pi*z)); //transform z back to r
	switch(tag){
		case 1:
		return pow(hepta_aux(1,z), 2) / r;
		break;
		case 2:
		return H_pre*(3.0*r+1.0)/(r*r*sqrt(r)*(r+1.0)*(r+1.0));
		break;
		default:
		cout << "Unknown tag for ODE-coefficient function: tag = " << tag << endl;
		cout << "Defaulting to: tag = 1 --> [ dz/dr ]" << endl;
		return hepta_coeff(1, z);
		break;
	}
}
//auxiliary functions for the ODE coefficients 						--------------> can probably be merged into hepta_coeff tbh
double grid::hepta_aux(int tag, double z){
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
		return hepta_aux(1, z);
		break;
	}
}
//update the harmonic dependency of the heptadiagonal matrix only
//when changing l, only a small part of U7 changes
//updat only these terms -> diagonals
void grid::update_hepta(int l){
	static const double Np1 = 1.0 / double(radial_points + 1); //static const 
	static const double jB = Np1 * Np1; //Bickley coefficient for the 2nd derivative static const 
	static const double i2 = 1.0 / (180.0 * jB);
	static double z;
	static double r;
	static double L;
	L = double(l*(l+1));
	static double c1;
	static double c2;

	for(int row=2; row<radial_points-2; row++){
		z = Np1 * (row+1);
		r = (1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
		U7(row,row) = -490.0 * hepta_coeff(1, z) * i2 - L / pow(r,3);
	}
	//handle different finite difference schemes for first and last 2 rows:

	//row 1
	z = Np1;
	r = (1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z) / (12.0 * jB);
	c2 = hepta_coeff(2, z) / (12.0 * Np1);
	U7(0,0) = -20.0 * c1 - 10.0 * c2 - (L / (r*r*r));
	//row 2
	z = 2.0 * Np1;
	r = (1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z) / (12.0 * jB);
	c2 = hepta_coeff(2, z) / (60.0 * Np1);
	U7(1,1) = -30.0 * c1 - 20.0 * c2 - (L / (r*r*r));
	//row N-1
	z = double(radial_points - 1) * Np1;
	r = (1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z) / (12.0 * jB);
	c2 = hepta_coeff(2, z) / (60.0 * Np1);
	U7(radial_points - 2,radial_points-2) = -30.0 * c1 + 20.0 * c2 - (L / (r*r*r));
	//row N
	z = double(radial_points) * Np1;
	r = (1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z) / (12.0 * jB);
	c2 = hepta_coeff(2, z) / (12.0 * Np1);
	U7(radial_points - 1,radial_points - 1) = -20.0 * c1 + 10.0 * c2 - (L / (r*r*r));

}
//solve the matrix equation for the transformed Poisson potential
//{U7} is the heptadiagonal matrix from the modified Bickley finite difference scheme with dimensions rpnt x rpnt (15x15)
//{rho} is the vector of Laplace coefficients for a specific {l,m} at each radial point, so a vector[rpnt] | {l,m}
//
//scheme will be made more efficient through simultaneously solving the equation for fixed l, for each m, since U7 != f(m)
//
void grid::Ulm_solve(int atom, int index){
	//solve matrix equation: U7 * Ui = -4pi * rho
	static const double pre = -4.0 * pi;
	static const double root = sqrt(-pre);
	static const double h = 1.0 / double(radial_points + 1);
	static const double z = radial_points * h;
	// static const double BC1 = (-11.0 * hepta_coeff(1, z) / (12.0 * h * h)) - 3.0 * (hepta_coeff(2, z) / (12.0 * h)); //U_N+1 term of U_N
	// static const double BC2 = (hepta_coeff(1, z - h) / (12.0 * h * h)) + 3.0 * (hepta_coeff(2, z - h) / (60.0 * h)); //U_N+1 term of U_N-1
	// static const double BC3 = -2.0 * (hepta_coeff(1, z - 2.0*h) / (180.0 * h * h)) - (hepta_coeff(2, z - 2.0*h) / (60 * h)); //U_N+1 term of U_N-2
	static const double BC1 = (-11.0 * hepta_coeff(1, h) / (12.0 * h * h)) + 3.0 * (hepta_coeff(2, h) / (12.0 * h)); //U_N+1 term of U_N
	static const double BC2 = (hepta_coeff(1, 2.0*h) / (12.0 * h * h)) - 3.0 * (hepta_coeff(2, 2.0*h) / (60.0 * h)); //U_N+1 term of U_N-1
	static const double BC3 = -2.0 * (hepta_coeff(1, 3.0*h) / (180.0 * h * h)) + (hepta_coeff(2, 3.0*h) / (60 * h)); //U_N+1 term of U_N-2
	//Boundary condition:
	static Eigen::VectorXd BC = Eigen::VectorXd::Zero(radial_points); //VectorXd initialises as a column-vector

	if(index == 0){
		//BC(radial_points-1) = BC1;
		//BC(radial_points-2) = BC2;
		//BC(radial_points-3) = BC3;
		BC(0) = BC1; //apply BCs in reverse since r->inf :: z->0 /\ r->0 :: z -> 1
		BC(1) = BC2;
		BC(2) = BC3;
		gridpoints[atom].Ulm.col(index) = U7.colPivHouseholderQr().solve(pre * gridpoints[atom].Laplace.col(index) + BC * root * get_point_charge(atom));
		//cout << "rhovec " << endl << pre * gridpoints[atom].Laplace.col(index) << endl;
		//cout << endl << "b: " << endl << gridpoints[atom].Ulm.col(index) << endl;
		//cout << endl << "rev " << endl << (BC * root * get_point_charge(atom)).colwise().reverse() << endl;
	}
	else{
		gridpoints[atom].Ulm.col(index) = U7.colPivHouseholderQr().solve(pre * gridpoints[atom].Laplace.col(index));
	}
	// cout << "rhovec " << endl << pre * gridpoints[atom].Laplace.col(index) << endl;
	// cout << endl << "b: " << endl << gridpoints[atom].Ulm.col(index) << endl;
	// cout << endl << "BC " << endl << BC * root * get_point_charge(atom) << endl;
	// cout << endl << "corr " << endl << pre * gridpoints[atom].Laplace.col(index) + BC * root * get_point_charge(atom) << endl;
	// cout << endl << "check " << endl << U7 * gridpoints[atom].Ulm.col(index) << endl << endl;
}

//compute the contribution to the external Poisson potential through cubic B-spline interpolation
double grid::Poisson_spline(int atomnr, double distance){
	//interpolate the Poisson potential to a point at {distance} away from the nucleus indexed by {atomnr}
	//
	//std::vector<double> f{0.01, -0.02, 0.3, 0.8, 1.9, -8.78, -22.6}; 
	//double t0 = 1; //variable denoting the starting point, so setting z = 0
	//static Eigen::VectorXd localdata; //local storage to assure proper communication of row pointer and iterator for boost::spline


	static double z_transform;
	Eigen::VectorXd localdata = Eigen::VectorXd::Zero(gridpoints[atomnr].potential.cols());
	localdata = gridpoints[atomnr].potential.row(0); //vector of data
	static const double h = 1.0/double(radial_points+1); //stepsize -------> fit the B-spline in z-domain since there the points are equispaced
	//cout << localdata << endl;
	boost::math::cubic_b_spline<double> spline(localdata.data(), localdata.size(), h, h); //fit the b-spline
	z_transform = acos((distance - 1.0)/(distance + 1.0)) / pi; //transform distance to the appropriate fractional z-value since the spline = f(z)


	//cout << z_transform << "   " << distance << "   " << (distance - 1.0)/(distance + 1.0) << "   " << acos((distance - 1.0)/(distance + 1.0)) << endl;
	//cout << "fit spline" << endl;
	//cout << spline(z_transform) << "  for  " << z_transform << " at " << distance << endl;

	return spline(z_transform*angular_points);
}

// double grid::lazy_spline(int atomnr, double distance){
// 	static const double h = 1.0 / double(radial_points + 1);
// 	static double z_transform;
// 	static double z_index;
// 	static int index1;
// 	static int index2;
// 	static double avg;

// 	z_transform = acos((distance - 1.0)/(distance + 1.0)) / pi; //transform r to z-coordinate
// 	z_index = z_transform / h; //get the number of integer steps that need to be taken to get to this point -> extract j
// 	index1 = int(std::floor(z_index)); //get the point just below the radial point
// 	index2 = int(std::ceil(z_index));  //get the point just above the radial point
// 	if(index1 < 0){
// 		index1 == 0; 
// 	}
// 	// if(index2 < index1){
// 	// 	index2 == index1 + 1;
// 	// }

// 	avg =  0.5 * (gridpoints[atomnr].potential(0,index1*angular_points) + gridpoints[atomnr].potential(0, index2*angular_points)); //just take the average...
// 	cout << "spline: " << index1 << "   " << index2 << "    " << z_index << "    " << avg << endl;
// 	return avg;
// }
double grid::lazy_spline(int atomnr, double distance, int lm){
	static const double h = 1.0 / double(radial_points + 1);
	static double z_transform;
	static double z_index;
	static int index1;
	static int index2;
	static double avg;

	z_transform = acos((distance - 1.0)/(distance + 1.0)) / pi; //transform r to z-coordinate
	z_index = z_transform / h; //get the number of integer steps that need to be taken to get to this point -> extract j
	index1 = int(std::floor(z_index)); //get the point just below the radial point
	index2 = int(std::ceil(z_index));  //get the point just above the radial point
	if(index1 < 0){
		index1 == 0;
	}
	// if(index2 < index1){
	// 	index2 == index1 + 1;
	// }

	avg =  0.5 * (gridpoints[atomnr].Ulm(index1,lm) + gridpoints[atomnr].Ulm(index2, lm)); //just take the average...
	//cout << "spline: " << index1 << "   " << index2 << "    " << z_index << "    " << avg << endl;
	return avg;
}



//****************************************************************
//begin main codeblock


//function generating the grid-object
//
//{pos_list} is a vector containing the vec3 coordinates of the nucleic centra
//{atnum_list} is a vector containing the atomic numbers of the nuclei
//{Basis} is the number of basis functions
//
grid make_grid(const vector<vec3>& pos_list, const vector<int>& atnum_list, int Basis)
{

	static const double pi = 3.14159265359;
	static const int radial_points = 15;
	static const int lebedev_order = Lebedev::LEBEDEV_110;
	static const int angular_points = Lebedev::num_lebedev_points[lebedev_order];

	// The Lebedev coefficients and radial points are stored in a large matrix.
    // Given a particular order for the integration, we have to start reading
    // this matrix from a specific position. Here, we calculate that position.
    int angular_index = 0;
    for(int i=0; i<lebedev_order; i++) {
        angular_index += Lebedev::num_lebedev_points[i];
    }

	//initialise grid-object
	grid Grid(pos_list, radial_points, angular_points, Basis, lebedev_order, angular_index);

	if(pos_list.size() != atnum_list.size()){
		cout << "List of atomic numbers does not match number of nuclei provided." << endl;
		return Grid;
	}

	int lidx;
	double x_std, r_point;
	double Bragg_factor, Chebyshev;
	vec3 point;

	//use Gauss-Chebyshev quadrature to generate the gridpoints
	for(int atom=0; atom<pos_list.size(); atom++){ //loop across all atoms
		
		//****************************************************************
		//generate atomic grid:


		Bragg_factor = 0.5*SlaterBragg(atnum_list[atom]); //reset Slater-Bragg radius
		if(atnum_list[atom]==1){
			Bragg_factor = 2.0*Bragg_factor; //do not use the 0.5 prefactor for hydrogen
		}

		//loop across all radial points
		for(int j=1; j<radial_points+1; j++){ //start index at 1 to avoid singularities in the radial coordinate
			
			//z_std = j/(radial_points+1); //get z-coordinate of the gridpoint in the transformed Gauss-Chebyshev domain
			x_std = cos(pi*j/(radial_points+1)); //get x-coordinate of the gridpoint in the standard Gauss-Chebyshev domain
			r_point = (1.0 + x_std)/(1.0 - x_std); //transform x-coordinate back to radial coordinate ------------- Bragg_factor* ?
			//compute Chebyshev-weight:
			Chebyshev = pow(sin(pi*j/(radial_points+1)),2.0) * (pi/(radial_points+1)); //get weight associated with the Gauss-Chebyshev discretisation
			//convert functional form of the Gauss-Chebyshev integral of the second kind:
			Chebyshev *= 1.0 / sqrt( (1-pow(x_std, 2.0)) );
			//correct for transforming the integration variable due to mapping x -> r
			Chebyshev *= 2.0 / pow( (1-x_std), 2.0);

			//get all angular points
			for(int a=0; a<angular_points; a++){
				lidx = a + angular_index; //transform index of the angular points to index for storage of the Lebedev coefficients
				//get coordinate of this angular point:
				//Lebedev coefficients are stored as an Nx4 matrix with x_angle, y_angle, z_angle, coeff
				//the angular point is thus pos_nucleus + coeff * [radius * x_angle; radius * y_angle; radius * z_angle]
				point << pos_list[atom] + vec3(Lebedev::lebedev_coefficients[lidx][0], Lebedev::lebedev_coefficients[lidx][1], Lebedev::lebedev_coefficients[lidx][2]) * r_point;
				Grid.set_gridpoint(atom, j-1, a, point);

				//write total weight as product of Chebyshev, transforms, Lebedev = 4*pi*coeff,
				//and add r^2 from the Jacobian for using solid angle spherical coordinate system
				//note that the Gaussian weights are still missing
				Grid.set_weight(atom, j-1, a, Chebyshev * 4 * pi * Lebedev::lebedev_coefficients[lidx][3] * ((point - pos_list[atom]).squaredNorm()) );
			}

		}
		double Gauss = 1.0; //dummy variable for holding the value of each Gaussian weight
		double nucleic_weight = 1.0; //weight of the nucleus under consideration;
		double norm = 0.0; //normalisation factor for the Gaussian weights
		cout << "commencing Gauss: " << endl;
		//loop across all gridpoints, for every nuclei, to get the Gaussian weights:
	    for(int rpnt=0; rpnt<radial_points; rpnt++){ //loop across radial points
	    	for(int apnt=0; apnt<angular_points; apnt++){ //loop across angular points
	    		norm = 0.0; //reset norm
	    		for(int i=0; i<atnum_list.size(); i++){ //loop across all atoms
					Gauss = Gauss_weight(i, pos_list, atnum_list, Grid.get_gridpoint(atom, rpnt, apnt)); //compute term of the Gaussian weight -> P_A
					norm += Gauss; //sum over all weights to get normalisation factor
					if(i == atom){
						nucleic_weight = Gauss; //retrieve Gaussian weight associated with the current atom -> P_n
					}
	    		}
	    		
	    		//map this weight to the grid
	    		if(abs(norm) > 1e-5){ //avoid singularities
	    			Grid.Gauss_weight(atom, rpnt, apnt, nucleic_weight/norm);
	    		}
	    	}
	    }

	    //Atomic grid end;
		//****************************************************************
	}

	return Grid;
    //Total grid end;
	//****************************************************************

}

//end main codeblock
//****************************************************************


//compute the Gaussian weights associated with the fuzzy Voronoi-cell
//
//{atomnr} is the index of the atom corresponding to the current atomic grid
//{pos_list} is a vector with the positions of the atoms in the system
//{atnum_list} is a vector with the atomic numbers of the atoms in the system
//{point} is the position vector of the mesh-centre of the current gridpoint
//
double Gauss_weight(int atomnr, const vector<vec3>& pos_list, const vector<int>& atnum_list, const vec3& point)
{
	double weight = 1.0; //the current weight: omega_n = P_n|fuzzy
	double mu = 0.0; //initialise diatomic confocal elliptical coordinate
	for(int i=0; i<pos_list.size(); i++){ //loop across all *other* atoms
		if(i!=atomnr){
			//get diatomic confocal elliptical coordinate
			mu = (   (pos_list[atomnr] - point).norm() - (pos_list[i] - point).norm()   )/(   (pos_list[atomnr] - pos_list[i]).norm()   );
			mu = Voronoi_Hetero(atnum_list[atomnr], atnum_list[i], mu); //correct for hetero-atomics
			weight *= fuzzy_cutoff(mu); //add to Voronoi-cell product-function
		}
	}
	//do not write weights to atomic grid here, do this in the main function
	return weight;
}

//smoothened cut-off function for the fuzzy Voronoi-cell boundary
//using recursion order 3 as recommended by Becke (A.D. Becke, J.Chem.Phys., Vol.88, p.2547, 1988)
//
//{mu} is the diatomic confocal elliptical coordinate (assumed to have been corrected in case of hetero-atomics)
//
double fuzzy_cutoff(double mu)
{
	double p_mu = mu; //initialisation
	for(int k = 0; k < 3; k++){
		p_mu = 1.5*p_mu - 0.5*p_mu*p_mu*p_mu; //apply recursion relation for the transformed, smoothened step-function
	}
	return 0.5-0.5*p_mu; //return cut-off function s(mu)
}

//compute the correction to the cut-off radius for hetero-atomics
//based on the procedure from A.D. Becke, J.Chem.Phys., Vol.88, p.2547, 1988
//
//{atnum1} is the atomic number of the first atom
//{atnum2} is the atomic number of the second atom
//{mu} is the non-corrected diatomic confocal elliptical coordinate
//
double Voronoi_Hetero(int atnum1, int atnum2, double mu)
{
    double chi = SlaterBragg(atnum1)/SlaterBragg(atnum2); //ratio of Slater-Bragg radii
    double u = (chi - 1.0)/(chi + 1.0); //reduced ratio
    double a = u/(u*u - 1.0); //correction factor for heteronuclearity
    //impose limiting values of the cut-off radii
    if(a > 0.5){
    	a = 0.5;
    }
    else if(a < -0.5){
    	a = -0.5;
    }

    return mu + a * (1.0 - mu*mu); //apply correction factor to the cut-off boundary
}

//get Slater-Bragg radius for atom with atomic number {atnum}
//Noble gas radii assumed equal to the radius of {atnum-1}
//Slater-Bragg radii are adapted from J.C. Slater, J. Chem. Phys., Vol.41, p.3199, 1964;
double SlaterBragg(int atnum)
{
	double radius; //initialise output: atomic radius
	//switch for assigning atomic radius to corresponding atomic number
	switch(atnum){
		case 1: //H
		radius = 0.35; //Becke recommends 0.35 ipv 0.25, evaluate performance before updating case 1 and default
		break;
		case 2: //He
		radius = 0.28;
		break;
		case 3: //Li
		radius = 1.45;
		break;
		case 4: //Be
		radius = 1.05;
		break;
		case 5: //B
		radius = 0.85;
		break;
		case 6: //C
		radius = 0.70;
		break;
		case 7: //N
		radius = 0.65;
		break;
		case 8: //O
		radius = 0.60;
		break;
		case 9: //F
		radius = 0.50;
		break;
		case 10: //Ne
		radius = 0.50;
		break;
		default:
		cout << "Atomic Number: " << atnum << " , currently not implemented or unknown." << endl;
		cout << "Defaulting to Hydrogen." << endl;
		radius = 0.35;
		break;
	}
	return radius;
}

//End of file
//****************************************************************

