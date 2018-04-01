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

// Transition state: use line search to find maximum along a specific coordinates
//					 search for f^2 = 0 -> f=0 and check that this gives 1 negative eigenvalue for H
// 				  	 NEB
// Also:
//			implement quasi-NR
//			analytical gradient (not 2nd)
//			rotate molecule to apply DOF constraints, perform optim., then rotate back

#include "geo.h"

GEO GEO_HF_nograd(vector<CGF> AO_list, const vector<vec3>& pos_init, const vector<double>& charge_list, const vector<int>& nelec_list){

	GEO geo_optim; //set output structure

	cout << "Commencing geometry optimisation: " << endl;
    cout << "Iter./Energy;" << endl;

	//SCF_E scf_data;
	SCF_E scf_update;
	scf_update = SCF_HF_energy(AO_list, pos_init, charge_list, nelec_list); //try and store matrices in scf_data as well -> assumed parsed

	static const double update_d = 0.1; //go 0.1 Bohr in any direction for updating geo positions

	vector<vec3> pos_temp = pos_init;
	//vec3 pos_new;
	//vec3 pos_opti = pos_init;
	vector<vec3> geoCon = pos_init;

	static const double ETol = -1e-5;
	double Eopti = 0.0;
	double Ediff;
	double Econ = scf_update.energy;

	//do not need this one yet, first gradient-less approach
	//Eigen::MatrixXd geo_derivs = Eigen::MatrixXd:Zero(pos_init.size()-1,pos_init.size()-1); //keep the 1st atom centred on the origin

	bool Convergence = false;
	bool loop_end = false;
	int count = 0;
	int CountC = 0;

	int CountC_max = 5;
	int count_max = 10;

	cout << "Reporting on iteration: " << endl;
	cout << CountC << ":\t" << Econ << endl; //report on first iteration

	//perform below loop until convergence is achieved
	while(!Convergence && CountC < CountC_max){
		CountC++;

		for(int i=1; i<pos_init.size(); i++){ //loop across all atoms that are not the first one
			for(int j=0; j<3; j++){ //go across all 3 coordinates
				//
				if(i==1 && (j==0 || j==1)){
					continue; //freedom of 2nd atom is restricted to stretching vibrations relative to atom1 (e.g. -2 degrees of freedom)
				}
				if(i==2 && j==0){				//see if it is necessary to rotate the molecule to adhere to this criterion
					continue; //if the number of atoms > 2, degrees of freedom are further limited by 1 -> formally assumes atom3 is in the xz-plane as we skip y
				}
				//

				//pos_temp = pos_opti;
				Ediff = 0.0;
				loop_end = false;
				count = 0;
				while(!loop_end && count < count_max){ //move in - direction
					count++;
					pos_temp[i][j] -= update_d; //take a small step in negative j-direction
					scf_update = SCF_HF_energy(AO_list, pos_temp, charge_list, nelec_list); //compute energy at this point
					Ediff = scf_update.energy - Eopti; //see how the energy has changed
					//now if the energy has improved: enew = -1.1, eopti = -1.0, ediff = -0.1 should be at least -0.01 Hartree
					if(Ediff < ETol){
						Eopti = scf_update.energy;
					}
					else{
						loop_end = true; //call for loop to end
						pos_temp[i][j] += update_d;//undo final step which was not favourable
					}
				}
				Ediff = 0.0; //reset params
				loop_end = false;
				count = 0;
				while(!loop_end && count < count_max){ //move in + direction
					count++;
					pos_temp[i][j] += update_d;
					scf_update = SCF_HF_energy(AO_list, pos_temp, charge_list, nelec_list);
					Ediff = scf_update.energy - Eopti;
					if(Ediff < ETol){
						Eopti = scf_update.energy;
					}
					else{
						loop_end = true; //call for loop to end
						pos_temp[i][j] -= update_d;//undo final step which was not favourable
					}
				}
			}
		}

		Ediff = Eopti - Econ;
		if(Ediff < ETol){
			Econ = Eopti;
			geoCon = pos_temp;
		}
		else{
			Convergence = true;
		}
		//end of convergence
		cout << "Reporting on iteration: " << endl;
		cout << CountC << ":\t" << Econ << endl << endl;
	}

	if(CountC == CountC_max){
		cout << "Stopping because maximum number of iterations was reached." << endl;
        cout << "Geometry may not have converged." << endl;
    }
    else{
        cout << "Stopping because convergence was achieved. (geo)" << endl;
    }

	cout << "Optimisation yields: energy = " << Econ << endl;
	cout << "Geometry: " << endl;
	for(int i=0; i<pos_init.size(); i++){
		cout << "Atom[" << i+1 << "]: " << geoCon[i][0] << ", " << geoCon[i][1] << ", " << geoCon[i][2] << endl;
	}
	cout << endl;

	//end of function
	geo_optim.pos = geoCon;
	geo_optim.energy = Econ;
	return geo_optim;
}

GEO GEO_HF_numerical(vector<CGF> AO_list, const vector<vec3>& pos_init, const vector<double>& charge_list, const vector<int>& nelec_list){

	GEO geo_optim; //set output structure

	cout << "Commencing geometry optimisation: " << endl;
    cout << "Iter./Energy;" << endl;

	//SCF_E scf_data;
	SCF_E scf_update;
	scf_update = SCF_HF_energy(AO_list, pos_init, charge_list, nelec_list); //try and store matrices in scf_data as well -> assumed parsed

	static const double update_d = 0.01; //go 0.1 Bohr in any direction for updating geo positions

	vector<vec3> pos_temp = pos_init;
	vector<vec3> pos_store = pos_init;
	//vec3 pos_opti = pos_init;
	vector<vec3> geoCon = pos_init;

	static const double ETol = -1e-5;
	double Eopti = 0.0;
	double Ediff;
	double Econ = scf_update.energy;

	//do not need this one yet, first gradient-less approach
	//Eigen::MatrixXd geo_derivs = Eigen::MatrixXd:Zero(pos_init.size()-1,pos_init.size()-1); //keep the 1st atom centred on the origin

	//Eigen::MatrixXd Hessian; 

	//incorporate exception for constructing B in case of linear molecules
	//catch is pos_list.size < 2 -> do not optimise
	//is size == 2 -> use different m -> get B of 3N x 3N-5 instead -> automatically adjusts Y

	//f(x) vector of 3N, {X,Y,Z} for each atom
	//f(y) [phi] vector of 3N-6 -> f(y) = f(x) * B-1
	//B is 3N x 3N-6 -> solve generating equation -> HOW?
	//get H(Y) by numerically approximating the derivatives of f(y)
	//update Ynew = Yold + H(Yold)-1 * f(Yold)

	//skip atom1, for atom2: skip x and y, for atom3: skip x or y
	//->to limit evaluation to only degrees of freedom instead of nuclear coordinates
	int DOF;
	int Natom = pos_init.size();
	if(Natom < 2){
		cout << "Warning: system only contains 1 (or fewer) atoms." << endl;
		cout << "Not possible to optimise geometry." << endl;
		return geo_optim;
	}
	else if(Natom == 2){
		DOF = 1; //3N - 5 for diatomic is 6 - 5 = 1
	}
	else{
		DOF = 3 * Natom - 6; //3N-6 for general polyatomics
	}
	Eigen::VectorXd f = Eigen::VectorXd::Zero(DOF);
	//f = construct_f(AO_list, pos_temp, charge_list, nelec_list);


	Eigen::VectorXd pos_upd = f; //get all DOF coordinates in one long vector for use in update rules

	bool Convergence = false;
	bool loop_end = false;
	int count = 0;
	int CountC = 0;

	int CountC_max = 100;
	int count_max = 10;

	cout << "Reporting on iteration: " << endl;
	cout << CountC << ":\t" << Econ << endl; //report on first iteration

	//start loop here for updating geometry
	//perform below loop until convergence is achieved
	while(!Convergence && CountC < CountC_max){
		CountC++;

		//construct f
		int index = 0;
		for(int i=1; i<Natom; i++){ //loop across all atoms that are not the first one
			for(int j=0; j<3; j++){ //go across all 3 coordinates
				//
				if(i==1 && (j==0 || j==1)){
					continue; //freedom of 2nd atom is restricted to stretching vibrations relative to atom1 (e.g. -2 degrees of freedom)
				}
				if(i==2 && j==1){				//see if it is necessary to rotate the molecule to adhere to this criterion
					continue; //if the number of atoms > 2, degrees of freedom are further limited by 1 -> formally assumes atom3 is in the xz-plane as we skip y
				}
				pos_temp[i][j] += update_d; //apply perturbation
				scf_update = SCF_HF_energy(AO_list, pos_temp, charge_list, nelec_list); //compute energy for a minor displacement of the structure to positive R
				pos_temp[i][j] -= update_d; //reset
				f[index] = (scf_update.energy - Econ) / update_d; //order of elements: atom1{none}, atom2{z}, atom3{x,z}, atom4{x,y,z}, atom5{x,y,z}, etc.
				index++;
			}
		}
		//f done

		Eigen::MatrixXd H = Eigen::MatrixXd::Zero(DOF, DOF);
		//construct H -> first 2 loops go over all the DOF which need to be differentiated, every time, fill 1 column of H
		Eigen::VectorXd f_new = Eigen::VectorXd::Zero(DOF);
		int Hindex = 0;
		for(int i=1; i<Natom; i++){
			//Hij for fj displaced along DOF i
			for(int j=0; j<3; j++){
				//
				if(i==1 && (j==0 || j==1)){
					continue; //freedom of 2nd atom is restricted to stretching vibrations relative to atom1 (e.g. -2 degrees of freedom)
				}
				if(i==2 && j==1){				//see if it is necessary to rotate the molecule to adhere to this criterion
					continue; //if the number of atoms > 2, degrees of freedom are further limited by 1 -> formally assumes atom3 is in the xz-plane as we skip y
				}
				//
				pos_temp[i][j] += update_d; //update for H
				//get new f
				//construct f
				index = 0;
				for(int iF=1; iF<Natom; iF++){ //loop across all atoms that are not the first one
					for(int jF=0; jF<3; jF++){ //go across all 3 coordinates
						//
						if(iF==1 && (jF==0 || jF==1)){
							continue; //freedom of 2nd atom is restricted to stretching vibrations relative to atom1 (e.g. -2 degrees of freedom)
						}
						if(iF==2 && jF==1){				//see if it is necessary to rotate the molecule to adhere to this criterion
							continue;//if the number of atoms >2, degrees of freedom are further limited by 1 -> formally assumes atom3 is in the xz-plane as we skip y
						}
						pos_temp[iF][jF] += update_d; //apply perturbation
						scf_update = SCF_HF_energy(AO_list, pos_temp, charge_list, nelec_list); //compute energy for a minor displacement of the structure to positive R
						pos_temp[iF][jF] -= update_d; //reset
						f_new[index] = 0.5*(scf_update.energy - Econ) / update_d; //order of elements: atom1{none}, atom2{z}, atom3{x,z}, atom4{x,y,z}, atom5{x,y,z}, etc.
						index++;
					}
				}
				//f done
				pos_temp[i][j] -= update_d; //reset H-update

				//fill one column of H
				for(int row=0; row<DOF; row++){
					H(row,Hindex) = (f_new[row] - f[row]) / update_d;
				}
				//one column of H is done
				Hindex++;
			}
		}
		//H done
		Eigen::MatrixXd H_inverse = H.inverse(); //get inverse of Hessian
		cout << "f" << f << " f_new " << f_new << " H " << H << " Hi " << H_inverse << " C " << H_inverse*f << endl;

		//base pos_upd on pos_temp:
		int idx = 0;
		for(int i=1; i<Natom; i++){
			for(int j=0; j<3; j++){
				//
				if(i==1 && (j==0 || j==1)){
					continue; //freedom of 2nd atom is restricted to stretching vibrations relative to atom1 (e.g. -2 degrees of freedom)
				}
				if(i==2 && j==1){				//see if it is necessary to rotate the molecule to adhere to this criterion
					continue; //if the number of atoms > 2, degrees of freedom are further limited by 1 -> formally assumes atom3 is in the xz-plane as we skip y
				}
				//
				pos_upd[idx] = pos_temp[i][j];
				idx++;
			}
		}
		cout << "upd: " << pos_upd << endl;
		//update geometry according to update rule
		pos_upd -= H_inverse*f*0.1;
		cout << "upd: " << pos_upd << endl;
		//convert upd back to pos_store:
		idx = 0;
		for(int i=1; i<Natom; i++){
			for(int j=0; j<3; j++){
				//
				if(i==1 && (j==0 || j==1)){
					continue; //freedom of 2nd atom is restricted to stretching vibrations relative to atom1 (e.g. -2 degrees of freedom)
				}
				if(i==2 && j==1){				//see if it is necessary to rotate the molecule to adhere to this criterion
					continue; //if the number of atoms > 2, degrees of freedom are further limited by 1 -> formally assumes atom3 is in the xz-plane as we skip y
				}
				//
				pos_temp[i][j] = pos_upd[idx];
				idx++;
			}
		}
		//compute new energy for this geometry
		scf_update = SCF_HF_energy(AO_list, pos_temp, charge_list, nelec_list);
		cout << "iter: " << CountC << ":\t" << scf_update.energy << endl;
		for(int i=0; i<pos_init.size(); i++){
			cout << "Atom[" << i+1 << "]: " << pos_temp[i][0] << ", " << pos_temp[i][1] << ", " << pos_temp[i][2] << endl;
		}

		//proceed to next iteration

		Ediff = scf_update.energy - Econ;
		if(Ediff < ETol){
			Econ = scf_update.energy; //progress energy
			geoCon = pos_temp; //store new geometry data
		}
		else{
			Convergence = true;
		}
		//end of convergence
		cout << "Reporting on iteration: " << endl;
		cout << CountC << ":\t" << Econ << endl << endl;
	}

	if(CountC == CountC_max){
		cout << "Stopping because maximum number of iterations was reached." << endl;
        cout << "Geometry may not have converged." << endl;
    }
    else{
        cout << "Stopping because convergence was achieved. (geo)" << endl;
    }

	cout << "Optimisation yields: energy = " << Econ << endl;
	cout << "Geometry: " << endl;
	for(int i=0; i<pos_init.size(); i++){
		cout << "Atom[" << i+1 << "]: " << geoCon[i][0] << ", " << geoCon[i][1] << ", " << geoCon[i][2] << endl;
	}
	cout << endl;

	//end of function
	geo_optim.pos = geoCon;
	geo_optim.energy = Econ;
	return geo_optim;
}

// VectorXd construct_f(int DOF, vector<CGF> AO_list, const vector<vec3>& pos_temp, const vector<double>& charge_list, const vector<int>& nelec_list){
// 	Eigen::VectorXd f = Eigen::VectorXd::Zero(DOF);
// 	//construct f
// 	int index = 0;
// 	for(int i=1; i<pos_init.size(); i++){ //loop across all atoms that are not the first one
// 		for(int j=0; j<3; j++){ //go across all 3 coordinates
// 			//
// 			if(i==1 && (j==0 || j==1)){
// 				continue; //freedom of 2nd atom is restricted to stretching vibrations relative to atom1 (e.g. -2 degrees of freedom)
// 			}
// 			if(i==2 && j==1){				//see if it is necessary to rotate the molecule to adhere to this criterion
// 				continue; //if the number of atoms > 2, degrees of freedom are further limited by 1 -> formally assumes atom3 is in the xz-plane as we skip y
// 			}
// 			pos_temp[i][j] += update_d; //apply perturbation
// 			scf_update = SCF_HF_energy(AO_list, pos_temp, charge_list, nelec_list); //compute energy for a minor displacement of the structure to positive R
// 			pos_temp[i][j] -= update_d; //reset
// 			f[index] = (scf_update.energy - Econ); //order of elements: atom1{none}, atom2{z}, atom3{x,z}, atom4{x,y,z}, atom5{x,y,z}, etc.
// 			index++;
// 		}
// 	}
// 	//f done
// 	return f;
// }

//call E_HF to compute energy
//call E_HF again, but now with a displacement of 1 of the atoms along x, y or z
//-> repeat for all atoms and every {x,y,z}

//--> will be expensive, but do not care, can be compared with analytical later



//perform " " " " through analytical expressions of the derivatives

//numerical DFT

//analytical DFT