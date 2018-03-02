/*****************************************************************

Basic closed-shell spin-restricted DFT-solver for simple molecules using STO-NG

Authors: B. Klumpers (bartkl@live.nl)
		 I.A.W. Filot

Published under GNU General Public License 3.0

Allows for SCF-computation of molecular energies for simple molecules.
Testcases for H, He, H2, HeH+, He2 and CO are included.

*****************************************************************/

#include "quadrature.h"

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
	grid Grid(pos_list.size(), radial_points, angular_points, Basis);

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
		//cout << "atom: " << atom << endl;
		//****************************************************************
		//generate atomic grid:


		Bragg_factor = 0.5*SlaterBragg(atnum_list[atom]); //reset Slater-Bragg radius
		if(atnum_list[atom]==1){
			Bragg_factor = 2.0*Bragg_factor; //do not use the 0.5 prefactor for hydrogen
		}

		//loop across all radial points
		for(int j=1; j<radial_points+1; j++){ //start index at 1 to avoid singularities in the radial coordinate
			//cout << "rpnt: " << j << endl;
			//z_std = j/(radial_points+1); //get z-coordinate of the gridpoint in the transformed Gauss-Chebyshev domain
			x_std = cos(pi*j/(radial_points+1)); //get x-coordinate of the gridpoint in the standard Gauss-Chebyshev domain
			r_point = (1.0 + x_std)/(1.0 - x_std); //transform x-coordinate back to radial coordinate Bragg_factor*
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
	    	//cout << "rpnt: " << rpnt << endl;
	    	for(int apnt=0; apnt<angular_points; apnt++){ //loop across angular points
	    		norm = 0.0; //reset norm
	    		for(int i=0; i<atnum_list.size(); i++){ //loop across all atoms
	    			//cout << "atnum: " << i << endl;
	    			//cout << "set: " << rpnt << "   " << apnt << "    " << atom << endl;
					Gauss = Gauss_weight(i, pos_list, atnum_list, Grid.get_gridpoint(atom, rpnt, apnt)); //compute term of the Gaussian weight -> P_A
					norm += Gauss; //sum over all weights to get normalisation factor
					//cout << "atnum end" << endl;
					if(i == atom){
						nucleic_weight = Gauss; //retrieve Gaussian weight associated with the current atom -> P_n
					}
	    		}
	    		
	    		//map this weight to the grid
	    		if(abs(norm) > 1e-5){ //avoid singularities
	    			Grid.Gauss_weight(atom, rpnt, apnt, nucleic_weight/norm);
	    			//cout << "Pn: " << nucleic_weight/norm << endl;
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
		radius = 0.25; //Becke recommends 0.35 ipv 0.25, evaluate performance before updating case 1 and default
		break;
		case 2: //He
		radius = 0.25;
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
		radius = 0.25;
		break;
	}
	return radius;
}

//End of file
//****************************************************************

