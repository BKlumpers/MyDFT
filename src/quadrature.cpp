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

//write values in spin-polarised systems
void grid::set_density(const Eigen::MatrixXd& P_alpha, const Eigen::MatrixXd& P_beta){
	//write total density

	//write alpha and beta densities by calling get_amps
	static double density_value_beta;
	static double density_value_alpha;
	//loop across all gridpoints and all basis functions
	for(int atomnr=0; atomnr<gridpoints.size(); atomnr++){ //loop across all atomic grids
		gridpoints[atomnr].point_charge = 0.0; //reset magnitude of point-charge potential
		for(int rpnt=0; rpnt<radial_points; rpnt++){ //loop across all radial points
			for(int apnt=0; apnt<angular_points; apnt++){ //loop across all angular points
				density_value_alpha = 0.0; //reset density storage
				density_value_beta = 0.0;
				//loop across all basis functions
				for(int i=0; i<P_alpha.rows(); i++){
					for(int j=0; j<P_alpha.rows(); j++){
						density_value_alpha += P_alpha(i,j)*get_amps(atomnr,rpnt,apnt,i)*get_amps(atomnr,rpnt,apnt,j);
						density_value_beta += P_beta(i,j)*get_amps(atomnr,rpnt,apnt,i)*get_amps(atomnr,rpnt,apnt,j);
					}
				}
				gridpoints[atomnr].spin(0, rpnt*angular_points+apnt) = density_value_alpha; //write total local electron density to atomic grid
				gridpoints[atomnr].spin(1, rpnt*angular_points+apnt) = density_value_beta;
				gridpoints[atomnr].amps(gridpoints[atomnr].amps.rows()-1, rpnt*angular_points+apnt) = density_value_alpha + density_value_beta;
				gridpoints[atomnr].point_charge += get_weight(atomnr,rpnt,apnt) * (density_value_alpha + density_value_beta);
			}
		}
	}
}

//compute the wavefunction gradients at each point
void grid::write_gradient(vector<CGF> AO_list){
	//get element ij for each point
	//store in matrix with 1 column per point, one row per element of ij
	vec3 point;
	vec3 gradi;
	vec3 gradj;
	vec3 wave;
	Eigen::MatrixXd Hessi;
	Eigen::MatrixXd Hessj;
	Eigen::MatrixXd H;
	//assume object gradients, same size as amps
	for(int atom=0; atom<gridpoints.size(); atom++){ //loop across all atomic grids
		for(int rpnt=0; rpnt<radial_points; rpnt++){ //loop across all radial points
			for(int apnt=0; apnt<angular_points; apnt++){ //loop across all angular points
				static int index;
				index = rpnt*angular_points + apnt;
				int rowidx = 0; //allows indexing of each element ij
				point = get_gridpoint(atom,rpnt,apnt);
				for(int i=0; i<AO_list.size(); i++){
					for(int j=0; j<AO_list.size(); j++){ //loop across all basis functions
						//get grad of i
						gradi = AO_list[i].getderiv(point);
						//get grad of j
						gradj = AO_list[j].getderiv(point);
						//get these grads at each point -> store above values
						//get total grad of element ij
						//store each of these ij:
						//gridpoints[atom].gradients[rowidx,index] = ( gridpoints[atom].amps(j,index) * gradi + gridpoints[atom].amps(i,index) * gradj ).norm;
						wave = gridpoints[atom].amps(j,index) * gradi + gridpoints[atom].amps(i,index) * gradj;
						gridpoints[atom].gradients_x(rowidx,index) = wave[0];
						gridpoints[atom].gradients_y(rowidx,index) = wave[1];
						gridpoints[atom].gradients_z(rowidx,index) = wave[2];
						gridpoints[atom].gradients_norm(rowidx,index) = wave.norm();
						//Hessian:
						Hessi = AO_list[i].getHessian(point);
						Hessj = AO_list[j].getHessian(point);

						H = Eigen::MatrixXd::Zero(3,3);
						H(0,0) = Hessi(0,0)*gridpoints[atom].amps(j,index) + 2.0*gradi[0]*gradj[0] + Hessj(0,0)*gridpoints[atom].amps(i,index);
						H(1,0) = H(0,1) = Hessi(1,0)*gridpoints[atom].amps(j,index) + gradi[1]*gradj[0] + gradi[0]*gradj[1] + Hessj(1,0)*gridpoints[atom].amps(i,index);
						H(2,0) = H(0,2) = Hessi(2,0)*gridpoints[atom].amps(j,index) + gradi[2]*gradj[0] + gradi[0]*gradj[2] + Hessj(2,0)*gridpoints[atom].amps(i,index);
						H(1,1) = Hessi(1,1)*gridpoints[atom].amps(j,index) + 2.0*gradi[1]*gradj[1] + Hessj(1,1)*gridpoints[atom].amps(i,index);
						H(1,2) = H(2,1) = Hessi(1,2)*gridpoints[atom].amps(j,index) + gradi[1]*gradj[2] + gradi[2]*gradj[1] + Hessj(1,2)*gridpoints[atom].amps(i,index);
						H(2,2) = Hessi(2,2)*gridpoints[atom].amps(j,index) + 2.0*gradi[2]*gradj[2] + Hessj(2,2)*gridpoints[atom].amps(i,index);

						gridpoints[atom].wave_Hess[rowidx].set_Hess(index,H);
						rowidx++;
					}
				}
			}
		}
	}
}

//compute the density gradient at each point
void grid::set_gradient(const Eigen::MatrixXd& Pmatrix){
	static double gradient;
	static double x;
	static double y;
	static double z;
	static Eigen::MatrixXd Hessian;
	for(int atom=0; atom<gridpoints.size(); atom++){
		for(int rpnt=0; rpnt<radial_points; rpnt++){
			for(int apnt=0; apnt<angular_points; apnt++){
				gradient = 0.0;
				x = 0.0;
				y = 0.0;
				z = 0.0;
				Hessian = Eigen::MatrixXd::Zero(3,3);
				static int index;
				index = rpnt*angular_points + apnt;
				int rowidx = 0;
				for(int i=0; i<Pmatrix.rows(); i++){
					for(int j=0; j<Pmatrix.rows(); j++){
						gradient += Pmatrix(i,j)*gridpoints[atom].gradients_norm(rowidx,index);
						x += Pmatrix(i,j)*gridpoints[atom].gradients_x(rowidx,index);
						y += Pmatrix(i,j)*gridpoints[atom].gradients_y(rowidx,index);
						z += Pmatrix(i,j)*gridpoints[atom].gradients_z(rowidx,index);
						//Hessian:
						Hessian += Pmatrix(i,j)*gridpoints[atom].wave_Hess[rowidx].get_Hess(index);

						rowidx++;
					}
				}
				gridpoints[atom].PBE(3,index) = gradient; //write density gradient to grid
				gridpoints[atom].PBE(0,index) = x; //write x,y,z elements:
				gridpoints[atom].PBE(1,index) = y;
				gridpoints[atom].PBE(2,index) = z;
				gridpoints[atom].rho_Hess[index].set_Hess(Hessian);
			}
		}
	}
}

//comput the alpha and beta gradients at each point
void grid::set_gradient(const Eigen::MatrixXd& P_alpha, const Eigen::MatrixXd& P_beta){ //function overload
	static double gradienta;
	static double gradientb;
	static double xa;
	static double ya;
	static double za;
	static double xb;
	static double yb;
	static double zb;
	static Eigen::MatrixXd Hessiana;
	static Eigen::MatrixXd Hessianb;
	for(int atom=0; atom<gridpoints.size(); atom++){
		for(int rpnt=0; rpnt<radial_points; rpnt++){
			for(int apnt=0; apnt<angular_points; apnt++){
				gradienta = 0.0;
				gradientb = 0.0;
				xa = 0.0;
				ya = 0.0;
				za = 0.0;
				xb = 0.0;
				yb = 0.0;
				zb = 0.0;
				Hessiana = Eigen::MatrixXd::Zero(3,3);
				Hessianb = Eigen::MatrixXd::Zero(3,3);
				static int index;
				index = rpnt*angular_points + apnt;
				int rowidx = 0;
				for(int i=0; i<P_alpha.rows(); i++){
					for(int j=0; j<P_alpha.rows(); j++){
						gradienta += P_alpha(i,j)*gridpoints[atom].gradients_norm(rowidx,index);
						gradientb += P_beta(i,j)*gridpoints[atom].gradients_norm(rowidx,index);
						
						xa += P_alpha(i,j)*gridpoints[atom].gradients_x(rowidx,index);
						ya += P_alpha(i,j)*gridpoints[atom].gradients_y(rowidx,index);
						za += P_alpha(i,j)*gridpoints[atom].gradients_z(rowidx,index);
						xb += P_beta(i,j)*gridpoints[atom].gradients_x(rowidx,index);
						yb += P_beta(i,j)*gridpoints[atom].gradients_y(rowidx,index);
						zb += P_beta(i,j)*gridpoints[atom].gradients_z(rowidx,index);
						
						//Hessian:
						Hessiana += P_alpha(i,j)*gridpoints[atom].wave_Hess[rowidx].get_Hess(index);
						Hessianb += P_beta(i,j)*gridpoints[atom].wave_Hess[rowidx].get_Hess(index);
						
						rowidx++;
					}
				}
				vec3 temp;
				temp << xa+xb,ya+yb,za+zb;
				
				gridpoints[atom].PBEC(3,index) = gradienta; //write density gradient to grid
				gridpoints[atom].PBEC(7,index) = gradientb;
				gridpoints[atom].PBEC(0,index) = xa; //write x,y,z elements:
				gridpoints[atom].PBEC(1,index) = ya;
				gridpoints[atom].PBEC(2,index) = za;
				gridpoints[atom].PBEC(4,index) = xb; //write x,y,z elements:
				gridpoints[atom].PBEC(5,index) = yb;
				gridpoints[atom].PBEC(6,index) = zb;
				
				gridpoints[atom].alpha_Hess[index].set_Hess(Hessiana);
				gridpoints[atom].beta_Hess[index].set_Hess(Hessianb);
				gridpoints[atom].rho_Hess[index].set_Hess(Hessiana+Hessianb);
				
				gridpoints[atom].PBE(3,index) = temp.norm(); //write density gradient to grid
				gridpoints[atom].PBE(0,index) = xa+xb; //write x,y,z elements:
				gridpoints[atom].PBE(1,index) = ya+yb;
				gridpoints[atom].PBE(2,index) = za+zb;
			}
		}
	}
}


//*******************************************
//Poisson potential:

//compute the Laplace coefficients
void grid::set_Laplace(){
	static const double omega = 4.0*pi;

	static int index;
	//for each l,m:
	for(int atom=0; atom<gridpoints.size(); atom++){ //compute them for each atomic grid
		for(int rpnt=0; rpnt<radial_points; rpnt++){ //compute each radial point
			//compute all {l,m} components
			for(int l=0; l<lmax+1; l++){
				for(int m=-l; m<l+1; m++){
					index = lm_index(l,m);
					gridpoints[atom].Laplace(rpnt, index) = 0.0; //reset coefficient
					//sum across all angular points:
					for(int apnt=0; apnt<angular_points; apnt++){
						//point << (get_gridpoint(atom, rpnt, 0) - pos_list[atom]); //shift coordinate system to position centre of the atomic grid at the origin
						gridpoints[atom].Laplace(rpnt, index) += get_Gauss(atom, rpnt, apnt)*get_density(atom, rpnt, apnt) * HarmonicReal(l, m, get_gridpoint(atom, rpnt, apnt) - pos_list[atom]) * omega * Lebedev::lebedev_coefficients[apnt+angular_index][3];//delta;
					}
				}
			}
		}
	}
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
					if(isnan(get_weight(atomnr,rpnt,apnt)/get_Gauss(atomnr, rpnt, apnt))){
						continue; //skip entries where get_Gauss ~= zero
					}
					PC += (get_weight(atomnr,rpnt,apnt)/get_Gauss(atomnr, rpnt, apnt))*gridpoints[atomnr].Laplace(rpnt, index) * HarmonicReal(l, m, get_gridpoint(atomnr, rpnt, apnt) - pos_list[atomnr]);
				}	//divide by Gauss coefficient to avoid double counting
			}		//Gauss coeff. in Laplace to account for atomic density, but this is not the density within the grid-volume, rather it is the "amplitude" of the density here!
		}			//PC sums over the atomic densities already, so the Gauss coeff. which is included in get_weight, should not be counted twice
	}
	cout << "point charge for atom " << atomnr+1 << " = " << PC << endl;
}

//compute the local Poisson potential at all points using finite differences
void grid::set_Poisson_L(){

	static int index;
	static double radius;

	//do recalculate U7 for each atom since formally the grids have different scaling for different atoms, hence values of r in U7 are different as well
	if(Bragg == false){
		set_hepta(0,0); //move inside the atomic loop if scaling of r through rm is implemented since rm = f(atom)
	}

	for(int atom=0; atom<gridpoints.size(); atom++){ //evaluate the potential on each atomic grid
		gridpoints[atom].potential.row(0) = Eigen::VectorXd::Zero(gridpoints[atom].potential.cols()); //reset potential
		if(Bragg == true){
			set_hepta(0,atom);
		}
		//for each rpnt, sum across all {l,m}
		for(int l=0; l<lmax+1; l++){
			update_hepta(l, atom);		//for given l, construct U7
			for(int m=-l; m<l+1; m++){
				index = lm_index(l,m);
				Ulm_solve(atom, index); //fill gridpoints[atom].Ulm given Laplace + Boundary conditions
				for(int rpnt=0; rpnt<radial_points; rpnt++){
					radius = gridpoints[atom].Ulm(rpnt, index) / ( (get_gridpoint(atom, rpnt, 0) - pos_list[atom]).norm() ); //radius is invariant under apnt
					for(int apnt=0; apnt<angular_points; apnt++){
						gridpoints[atom].potential(0,rpnt*angular_points+apnt) +=  radius * HarmonicReal(l, m, get_gridpoint(atom, rpnt, apnt) - pos_list[atom]);
					}
				}
			}
		}
	}
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
						distance = (get_gridpoint(atom1,rpnt,apnt) - pos_list[atom2]).norm(); //norm between gridpoint and origin of atom2
						//call cubic (2nd order polynomial) spline interpolation function to get the external potential due to atom2 at distance from its origin
						for(int l=0; l<lmax+1; l++){
							for(int m=-l; m<l+1; m++){
								lm = lm_index(l,m);
								//gridpoints[atom1].potential(1,rpnt*angular_points+apnt) += lazy_spline(atom2, distance, lm) * HarmonicReal(l,m,get_gridpoint(atom1,rpnt,apnt) - pos_list[atom2]) / distance;
								gridpoints[atom1].potential(1,rpnt*angular_points+apnt) += Poisson_spline(atom2, distance, lm) * HarmonicReal(l,m,get_gridpoint(atom1,rpnt,apnt) - pos_list[atom2]) / distance;
							}
						}
					}
				}
			}
		}
	}
}
void grid::set_hepta(int l, int atomnr){
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
	static double Slater;
	Slater = gridpoints[atomnr].Slater;

	for(int row=3; row<radial_points-3; row++){ //run from 4th row till N-4th row
		//write coefficients of 7 elements
		z = Np1 * (row+1);
		r = Slater*(1.0 + cos(pi*z)) / (1 - cos(pi*z));
		c1 = hepta_coeff(1, z, Slater) * i2;
		c2 = hepta_coeff(2, z, Slater) * i1;
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
	r = Slater*(1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z, Slater) / (12.0 * jB);
	c2 = hepta_coeff(2, z, Slater) / (12.0 * Np1);
	U7(0,0) = -20.0 * c1 - 10.0 * c2 - (L / (r*r*r));
	U7(0,1) = 6.0 * c1 + 18.0 * c2;
	U7(0,2) = 4.0 * c1 - 6.0 * c2;
	U7(0,3) = -c1 + c2;

	//row 2
	z = 2.0 * Np1;
	r = Slater*(1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z, Slater) / (12.0 * jB);
	c2 = hepta_coeff(2, z, Slater) / (60.0 * Np1);
	U7(1,0) = 16.0* c1 - 30.0 * c2;
	U7(1,1) = -30.0 * c1 - 20.0 * c2 - (L / (r*r*r));
	U7(1,2) = 16.0 * c1 + 60 * c2;
	U7(1,3) = -1.0 * c1 - 15.0 * c2;
	U7(1,4) = 2.0 * c2;

	//row 3
	z = 3.0 * Np1;
	r = Slater*(1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z, Slater) * i2;
	c2 = hepta_coeff(2, z, Slater) * i1;
	U7(2,0) = -27.0 * c1 + 9.0 * c2;
	U7(2,1) = 270.0 * c1 - 45.0 * c2;
	U7(2,2) = -490.0 * c1 - (L / (r*r*r));
	U7(2,3) = 270.0 * c1 + 45.0 * c2;
	U7(2,4) = -27.0 * c1 - 9.0 * c2;
	U7(2,5) = 2.0 * c1 + c2;

	//row N-2
	z = double(radial_points - 2) * Np1;
	r = Slater*(1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z, Slater) * i2;
	c2 = hepta_coeff(2, z, Slater) * i1;
	U7(radial_points - 3,radial_points-6) = 2.0 * c1 - c2;
	U7(radial_points - 3,radial_points-5) = -27.0 * c1 + 9.0 * c2;
	U7(radial_points - 3,radial_points-4) = 270.0 * c1 - 45.0 * c2;
	U7(radial_points - 3,radial_points-3) = -490.0 * c1 - (L / (r*r*r));
	U7(radial_points - 3,radial_points-2) = 270.0 * c1 + 45.0 * c2;
	U7(radial_points - 3,radial_points-1) = -27.0 * c1 - 9.0 * c2;

	//row N-1
	z = double(radial_points - 1) * Np1;
	r = Slater*(1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z, Slater) / (12.0 * jB);
	c2 = hepta_coeff(2, z, Slater) / (60.0 * Np1);
	U7(radial_points - 2,radial_points-5) =  -2.0 * c2;
	U7(radial_points - 2,radial_points-4) = -c1 + 15.0* c2;
	U7(radial_points - 2,radial_points-3) = 16.0 * c1 - 60.0 * c2;
	U7(radial_points - 2,radial_points-2) = -30.0 * c1 + 20.0 * c2 - (L / (r*r*r));
	U7(radial_points - 2,radial_points-1) = 16.0 * c1 + 30.0 * c2;

	//row N
	z = double(radial_points) * Np1;
	r = Slater*(1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z, Slater) / (12.0 * jB);
	c2 = hepta_coeff(2, z, Slater) / (12.0 * Np1);
	U7(radial_points - 1,radial_points - 4) = -c1 - c2;
	U7(radial_points - 1,radial_points - 3) = 4.0* c1 + 6.0 * c2;
	U7(radial_points - 1,radial_points - 2) = 6.0 * c1 - 18.0 * c2;
	U7(radial_points - 1,radial_points - 1) = -20.0 * c1 + 10.0 * c2 - (L / (r*r*r));

}

//calculate the non-constant ODE coefficients
double grid::hepta_coeff(int tag, double z, double Slater){
	static double r;
	static const double H_pre = 1.0/(2.0*pi);
	static const double pre = 1.0/(pi*pi);
	r = Slater * (1.0 + cos(pi*z)) / (1.0 - cos(pi*z)); //transform z back to r
	switch(tag){
		case 1:
		return Slater*pre/(r*r*(r+Slater)*(r+Slater));
		break;
		case 2:
		return Slater*Slater*(3.0*r+Slater)*H_pre/(r*r*Slater*(Slater+r)*(Slater+r)*sqrt(r*Slater));
		break;
		default:
		cout << "Unknown tag for ODE-coefficient function: tag = " << tag << endl;
		cout << "Defaulting to: tag = 1 --> [ dz/dr ]" << endl;
		return hepta_coeff(1, z, Slater);
		break;
	}
}

//update the harmonic dependency of the heptadiagonal matrix only
//when changing l, only a small part of U7 changes
//updat only these terms -> diagonals
void grid::update_hepta(int l, int atomnr){
	static const double Np1 = 1.0 / double(radial_points + 1); //static const 
	static const double jB = Np1 * Np1; //Bickley coefficient for the 2nd derivative static const 
	static const double i2 = 1.0 / (180.0 * jB);
	static double z;
	static double r;
	static double L;
	L = double(l*(l+1));
	static double c1;
	static double c2;
	static double Slater;
	Slater = gridpoints[atomnr].Slater;

	for(int row=2; row<radial_points-2; row++){
		z = Np1 * (row+1);
		r = Slater*(1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
		U7(row,row) = -490.0 * hepta_coeff(1, z, Slater) * i2 - L / pow(r,3);
	}
	//handle different finite difference schemes for first and last 2 rows:

	//row 1
	z = Np1;
	r = Slater*(1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z, Slater) / (12.0 * jB);
	c2 = hepta_coeff(2, z, Slater) / (12.0 * Np1);
	U7(0,0) = -20.0 * c1 - 10.0 * c2 - (L / (r*r*r));
	//row 2
	z = 2.0 * Np1;
	r = Slater*(1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z, Slater) / (12.0 * jB);
	c2 = hepta_coeff(2, z, Slater) / (60.0 * Np1);
	U7(1,1) = -30.0 * c1 - 20.0 * c2 - (L / (r*r*r));
	//row N-1
	z = double(radial_points - 1) * Np1;
	r = Slater*(1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z, Slater) / (12.0 * jB);
	c2 = hepta_coeff(2, z, Slater) / (60.0 * Np1);
	U7(radial_points - 2,radial_points-2) = -30.0 * c1 + 20.0 * c2 - (L / (r*r*r));
	//row N
	z = double(radial_points) * Np1;
	r = Slater*(1.0 + cos(pi*z)) / (1.0 - cos(pi*z));
	c1 = hepta_coeff(1, z, Slater) / (12.0 * jB);
	c2 = hepta_coeff(2, z, Slater) / (12.0 * Np1);
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
	static double Slater;
	Slater = gridpoints[atom].Slater;
	// static const double BC1 = (-11.0 * hepta_coeff(1, z) / (12.0 * h * h)) - 3.0 * (hepta_coeff(2, z) / (12.0 * h)); //U_N+1 term of U_N
	// static const double BC2 = (hepta_coeff(1, z - h) / (12.0 * h * h)) + 3.0 * (hepta_coeff(2, z - h) / (60.0 * h)); //U_N+1 term of U_N-1
	// static const double BC3 = -2.0 * (hepta_coeff(1, z - 2.0*h) / (180.0 * h * h)) - (hepta_coeff(2, z - 2.0*h) / (60 * h)); //U_N+1 term of U_N-2
	static const double BC1 = (-11.0 * hepta_coeff(1, h, Slater) / (12.0 * h * h)) + 3.0 * (hepta_coeff(2, h, Slater) / (12.0 * h)); //U_N+1 term of U_N
	static const double BC2 = (hepta_coeff(1, 2.0*h, Slater) / (12.0 * h * h)) - 3.0 * (hepta_coeff(2, 2.0*h, Slater) / (60.0 * h)); //U_N+1 term of U_N-1
	static const double BC3 = -2.0 * (hepta_coeff(1, 3.0*h, Slater) / (180.0 * h * h)) + (hepta_coeff(2, 3.0*h, Slater) / (60 * h)); //U_N+1 term of U_N-2
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
	}
	else{
		gridpoints[atom].Ulm.col(index) = U7.colPivHouseholderQr().solve(pre * gridpoints[atom].Laplace.col(index));
	}
}

//compute the contribution to the external Poisson potential through cubic B-spline interpolation
double grid::Poisson_spline(int atomnr, double distance, int lm){
	//interpolate the Poisson potential to a point at {distance} away from the nucleus indexed by {atomnr}
	//
	//std::vector<double> f{0.01, -0.02, 0.3, 0.8, 1.9, -8.78, -22.6}; 
	//double t0 = 1; //variable denoting the starting point, so setting z = 0
	//static Eigen::VectorXd localdata; //local storage to assure proper communication of row pointer and iterator for boost::spline

	static double z_transform;
	Eigen::VectorXd localdata = Eigen::VectorXd::Zero(gridpoints[atomnr].Ulm.rows());
	localdata = gridpoints[atomnr].Ulm.col(lm); //vector of data
	static const double h = 1.0/double(radial_points+1); //stepsize -------> fit the B-spline in z-domain since there the points are equispaced
	boost::math::cubic_b_spline<double> spline(localdata.data(), localdata.size(), h, h); //fit the b-spline
	z_transform = acos(((distance/gridpoints[atomnr].Slater) - 1.0)/((distance/gridpoints[atomnr].Slater) + 1.0)) / pi; //transform distance to the appropriate fractional z-value since the spline = f(z)

	return spline(z_transform);
}

double grid::lazy_spline(int atomnr, double distance, int lm){
	static const double h = 1.0 / double(radial_points + 1);
	static double z_transform;
	static double z_index;
	static int index1;
	static int index2;
	static double avg;

	z_transform = acos(((distance/gridpoints[atomnr].Slater) - 1.0)/((distance/gridpoints[atomnr].Slater) + 1.0)) / pi; //transform r to z-coordinate
	z_index = z_transform / h; //get the number of integer steps that need to be taken to get to this point -> extract j
	index1 = int(std::floor(z_index)); //get the point just below the radial point
	index2 = int(std::ceil(z_index));  //get the point just above the radial point
	if(index1 < 0){
		index1 == 0;
	}

	avg =  0.5 * (gridpoints[atomnr].Ulm(index1,lm) + gridpoints[atomnr].Ulm(index2, lm)); //just take the average...
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
	static const int setting = 0;
	int radial_points, lebedev_order;
	//retrieve preset grids
	switch(setting){
		case 0: //medium grid
		radial_points = 15;
		lebedev_order = Lebedev::LEBEDEV_110;
		break;
		case 1: //ultrafine grid
		radial_points = 30;
		lebedev_order = Lebedev::LEBEDEV_194;
		break;
		case 2: //fine grid
		radial_points = 20;
		lebedev_order = Lebedev::LEBEDEV_146;
		break;
		case 3: //coarse grid
		radial_points = 10;
		lebedev_order = Lebedev::LEBEDEV_50;
		break;
		default:
		radial_points = 15;
		lebedev_order = Lebedev::LEBEDEV_110;
		break;
	}
	static const int angular_points = Lebedev::num_lebedev_points[lebedev_order];

	//get the starting position from which to start reading the Lebedev coefficients
	//corresponding to the defined order
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

		if(Bragg == true){
			Bragg_factor = 0.5*SlaterBragg(atnum_list[atom]); //reset Slater-Bragg radius
			if(atnum_list[atom]==1 || atnum_list[atom]==2){
				Bragg_factor = 1.0; //do not use the 0.5 prefactor for hydrogen, Bragg_factors of 1 are found to work better for H/He
			}
		}
		else{
			Bragg_factor = 1.0;
		}
		Grid.set_Slater(atom, Bragg_factor); //store Slater-Bragg radius in the atomic grid for computation of the Poisson potential

		//loop across all radial points
		for(int j=1; j<radial_points+1; j++){ //start index at 1 to avoid singularities in the radial coordinate
			
			//z_std = j/(radial_points+1); //get z-coordinate of the gridpoint in the transformed Gauss-Chebyshev domain
			x_std = cos(pi*j/(radial_points+1)); //get x-coordinate of the gridpoint in the standard Gauss-Chebyshev domain
			r_point = Bragg_factor*(1.0 + x_std)/(1.0 - x_std); //transform x-coordinate back to radial coordinate ------------- Bragg_factor* ?
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
//Noble gas radii assumed equal to the radius of {atnum-1} --> (reasonable assumption)
//Slater-Bragg radii are adapted from J.C. Slater, J. Chem. Phys., Vol.41, p.3199, 1964;
//values for H and He have been modified to assure reasonable calculation results
double SlaterBragg(int atnum)
{
	double radius = 1.88973; //initialise output: atomic radius, factor 1.88973 is to convert radii from Angstrom to Bohr
	//switch for assigning atomic radius to corresponding atomic number
	switch(atnum){
		case 1: //H
		radius *= 0.35; //Becke recommends 0.35 ipv 0.25, evaluate performance before updating case 1 and default, additional factor 2 since the 0.5 prefactor should not be used for hydrogen
		break;
		case 2: //He
		radius *= 0.25;
		break;
		case 3: //Li
		radius *= 1.45;
		break;
		case 4: //Be
		radius *= 1.05;
		break;
		case 5: //B
		radius *= 0.85;
		break;
		case 6: //C
		radius *= 0.70;
		break;
		case 7: //N
		radius *= 0.65;
		break;
		case 8: //O
		radius *= 0.60;
		break;
		case 9: //F
		radius *= 0.50;
		break;
		case 10: //Ne
		radius *= 0.50;
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

