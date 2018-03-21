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

#include "hf.h"

//****************************************************************
//Hartree-Fock energy minimisation:

/*SCF-loop for minimisation of the total energy
//
//{AO_list} is a vectors of CGFs comprising the basis set
//{pos_list} is a vector of the coordinates of the nuclei
//{charge_list} is a vector of the nuclear charges
//{nelec_list} is a vector of the number of electrons associated with each atom
*/
SCF_E SCF_UHF_energy(vector<CGF> AO_list, const vector<vec3>& pos_list, const vector<double>& charge_list, const vector<int>& nelec_list)
{

    SCF_E out; //initialise output struct

    if(pos_list.size()!=charge_list.size()){
        cout << endl << "Error: number of nuclei not equal to supplied list of charges" << endl << endl;
        return out;
    }

    //get total number of electrons in the system
    int nelec = 0;
    for(int i=0; i<nelec_list.size(); i++) {
        nelec += nelec_list[i];
    }
    //assign total number of alpha and beta electrons:
    int N_alpha, N_beta;
    bool spin;
    if(nelec % 2 != 0){
        N_beta = (nelec - 1) / 2;
        N_alpha = nelec - N_beta;
        spin = true;
    }
    else{
        N_alpha = nelec/2;
        N_beta = N_alpha;
        spin = false;
    }

    //initialise matrices for overlap, kinetic and nuclear Fock-matrix components
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //overlap_matrix
    Eigen::MatrixXd Kin = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //kinetic_matrix
    Eigen::MatrixXd Nucl = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //nuclear_matrix

    //calculate respective 1e-integrals: overlap, kinetic and nuclear
    for(int i=0; i<AO_list.size(); i++) {
        for(int j=0; j<AO_list.size(); j++) {
            S(i,j) = overlapCGF(AO_list[i], AO_list[j]);
            Kin(i,j) = kineticCGF(AO_list[i], AO_list[j]);
            for(int n=0; n<charge_list.size(); n++){
                Nucl(i,j) += nuclearCGF(AO_list[i], AO_list[j], charge_list[n], pos_list[n]);
            }
        }
    }
    //get 1e-Hamiltonian
    Eigen::MatrixXd Hamil = Kin + Nucl;

    //calculate all unique 2e-integral using index sorting:
    //
    //vector of unique 2e-integrals, -1 as a reference value for non-calculated integral
    vector<double> Rep(two_electronSort(AO_list.size(),AO_list.size(),AO_list.size(),AO_list.size()), -1.0);
    int ij,kl,current_index;
    for(int i=0; i<AO_list.size(); i++) {
        for(int j=0; j<AO_list.size(); j++) {
            ij = i*(i+1)/2 + j;
            for(int k=0; k<AO_list.size(); k++) {
                for(int l=0; l<AO_list.size(); l++) {
                    kl = k * (k+1)/2 + l;
                    if(ij <= kl) {
                        current_index = two_electronSort(i,j,k,l);
                        //do not compute the integral if it has previously been evaluated (repulsion cannot be negative, so this is a safe default value)
                        if(Rep[current_index] != -1.0) {
                            continue;
                        }
                        Rep[current_index] = two_electronCGF(AO_list[i], AO_list[j], AO_list[k], AO_list[l]);
                    }
                }
            }
        }
    }

    //perform canonical diagonalisation:
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
    Eigen::MatrixXd Diagonal = es.eigenvalues().real().asDiagonal(); //diagonal matrix
    Eigen::MatrixXd Unitary = es.eigenvectors().real(); //unitary matrix

    //take 1/sqrt(D); element-wise since the matrix is diagonal
    for(int i=0; i<AO_list.size(); i++) {
        Diagonal(i,i) = 1.0 / sqrt(Diagonal(i,i));
    }

    //Calculate the transformation matrix
    Eigen::MatrixXd X_transform = Unitary * Diagonal;

    //optimisation parameters:
    static const double alpha = 0.5; //mixing parameter for updating density
    static const double TolEnergy = 1e-5; //tolerance for energy convergence
    static const int loop_max = 100; //maximum number of iterations

    //matrix initialisation for energy-minimisation:
    Eigen::MatrixXd P_density_alpha = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //initialise density-matrix
    Eigen::MatrixXd P_density_beta; //initial guess depends on whether we should look for a unique unrestricted solution
    if(!spin){
        P_density_beta = Eigen::MatrixXd::Constant(AO_list.size(), AO_list.size(), alpha); //initialise density-matrix
    }
    else if(spin){
        P_density_beta = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //initialise density-matrix
    }
    Eigen::MatrixXd P_density_new_alpha = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //initialise update-matrix for density
    Eigen::MatrixXd P_density_new_beta = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //initialise update-matrix for density
    Eigen::MatrixXd G_two_electron_alpha = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //Initialise two-electron Hamiltonian matrix
    Eigen::MatrixXd G_two_electron_beta = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //Initialise two-electron Hamiltonian matrix
    Eigen::MatrixXd C_vector_alpha; //initialise eigenvectors of the Fock-matrix
    Eigen::MatrixXd C_vector_beta; //initialise eigenvectors of the Fock-matrix
    Eigen::VectorXd orbital_energies_alpha; //initialise eigenvalues of the Fock-matrix
    Eigen::VectorXd orbital_energies_beta; //initialise eigenvalues of the Fock-matrix

    Eigen::MatrixXd Fock_alpha = Hamil;
    Eigen::MatrixXd Fock_beta = Hamil;
    Eigen::MatrixXd F_transform_alpha = Hamil;
    Eigen::MatrixXd F_transform_beta = Hamil;

    //initialise update parameters:
    double energy_old = 0.0; //energy obtained from previous iteration for comparison
    if(pos_list.size() == 1){ //correct update for monoatomics
        energy_old = 1.0;
    }
    double energy_difference = 1.0; //difference between new and old energy
    int loop_counter = 0; //initialise counter
    double energy = 0.0; //initialise energy

    cout << "Beginning electronic optimisation: (UHF) " << endl;
    cout << "Iter./Energy;" << endl;

    //SCF-loop for energy optimisation:
    while(energy_difference > TolEnergy && loop_counter < loop_max) {
        loop_counter++; //increment loop counter

        //Calculate two-electron hamiltonian matrix:
        static int indexJ, indexK;
        for(int i=0; i<AO_list.size(); i++) {
            for(int j=0; j<AO_list.size(); j++) {
                G_two_electron_alpha(i,j) = 0; //reset matrix
                G_two_electron_beta(i,j) = 0; //reset matrix
                for(int k=0; k<AO_list.size(); k++) {
                    for(int l=0; l<AO_list.size(); l++) {
                        indexJ = two_electronSort(i,j,l,k);
                        indexK = two_electronSort(i,k,l,j);
                        G_two_electron_alpha(i,j) += (P_density_alpha(k,l) + P_density_beta(k,l)) * Rep[indexJ] - P_density_alpha(k,l) * Rep[indexK];
                        G_two_electron_beta(i,j) += (P_density_alpha(k,l) + P_density_beta(k,l)) * Rep[indexJ] - P_density_beta(k,l) * Rep[indexK];                        
                    }
                }
            }
        }

        //Calculate Fock Matrix
        Fock_alpha = Hamil + G_two_electron_alpha;
        Fock_beta = Hamil + G_two_electron_beta;

        //Transform Fock Matrix
        F_transform_alpha = X_transform.transpose() * Fock_alpha * X_transform;
        F_transform_beta = X_transform.transpose() * Fock_beta * X_transform;

        //Calculate eigenvalues within the transformed basis
        es.compute(F_transform_alpha);
        Eigen::MatrixXd C_vector_transform = es.eigenvectors().real();
        Eigen::MatrixXd orbital_energies_store = es.eigenvalues().real().asDiagonal(); //store values of the orbital energies
        //Transform coefficient matrix back to original basis
        C_vector_alpha = X_transform * C_vector_transform;
        orbital_energies_alpha = orbital_energies_store.diagonal(); //transform matrix back to vector
        //Calculate eigenvalues within the transformed basis
        es.compute(F_transform_beta);
        C_vector_transform = es.eigenvectors().real();
        orbital_energies_store = es.eigenvalues().real().asDiagonal(); //store values of the orbital energies
        //Transform coefficient matrix back to original basis
        C_vector_beta = X_transform * C_vector_transform;
        orbital_energies_beta = orbital_energies_store.diagonal(); //transform matrix back to vector


        //Calculate total electronic energy of the system
        energy = 0.0;
        Eigen::MatrixXd E_alpha = Hamil + Fock_alpha; //non-weighted components of total electronic energy
        Eigen::MatrixXd E_beta = Hamil + Fock_beta;
        for(int i=0; i<AO_list.size(); i++) {
            for(int j=0; j<AO_list.size(); j++) {
                energy += 0.5 * P_density_alpha(j,i) * E_alpha(i,j); //weigh matrix elements
                energy += 0.5 * P_density_beta(j,i) * E_beta(i,j);
            }
        }

        //Add nuclear repulsion to the orbital energies
        for(int i=0; i<charge_list.size(); i++){
            for(int j=0; j<charge_list.size(); j++){
                if(j>i){
                    energy += charge_list[i]*charge_list[j]/(pos_list[i]-pos_list[j]).norm();
                }
            }
        }

        //generate new density matrix
        P_density_new_alpha = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size());
        for(int i=0; i<AO_list.size(); i++) {
            for(int j=0; j<AO_list.size(); j++) {
                for(int k=0; k<N_alpha; k++) {
                    P_density_new_alpha(i,j) += C_vector_alpha(i,k) * C_vector_alpha(j,k);
                }
            }
        }

        //update density matrix by mixing the old and new matrices
        for(int i=0; i<AO_list.size(); i++) {
            for(int j=0; j<AO_list.size(); j++) {
                P_density_alpha(i,j) = (1.0-alpha) * P_density_new_alpha(i,j) + alpha * P_density_alpha(i,j);
            }
        }

        //generate new density matrix
        P_density_new_beta = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size());
        for(int i=0; i<AO_list.size(); i++) {
            for(int j=0; j<AO_list.size(); j++) {
                for(int k=0; k<N_beta; k++) {
                    P_density_new_beta(i,j) += C_vector_beta(i,k) * C_vector_beta(j,k);
                }
            }
        }

        //update density matrix by mixing the old and new matrices
        for(int i=0; i<AO_list.size(); i++) {
            for(int j=0; j<AO_list.size(); j++) {
                P_density_beta(i,j) = (1.0-alpha) * P_density_new_beta(i,j) + alpha * P_density_beta(i,j);
            }
        }

        //report energy and calculate difference with previous energy value
        energy_difference = abs(energy - energy_old);
        energy_old = energy; //store current value for the energy for comparison in the next iteration
        cout << loop_counter << ":\t" << energy << endl;
    }
    //end SCF-loop

    //output SCF results
    cout << "Stopping because energy convergence was achieved." << endl;
    cout << "System Energy = " << energy << endl << endl;

    //sort molecular orbitals
    Eigen::VectorXd e_joined = Eigen::VectorXd::Zero(orbital_energies_alpha.size()+orbital_energies_beta.size());
    e_joined << orbital_energies_alpha, orbital_energies_beta;
    Eigen::VectorXd e_fin = Eigen::VectorXd::Zero(e_joined.size());
    Eigen::MatrixXd c_joined = Eigen::MatrixXd::Zero(C_vector_alpha.rows(), C_vector_alpha.cols() + C_vector_beta.cols());
    for(int i=0; i<C_vector_alpha.cols(); i++){
        c_joined.col(i) = C_vector_alpha.col(i);
    }
    int csize = C_vector_alpha.cols();
    for(int j=csize; j<csize+C_vector_beta.cols(); j++){
        c_joined.col(j) = C_vector_beta.col(j-csize);
    }
    Eigen::MatrixXd c_fin = Eigen::MatrixXd::Zero(c_joined.rows(), c_joined.cols());
    double min, max, temp, idxt;
    min = e_joined[0];
    idxt = 0;
    max = min;
    //get maximum as boundary, get min as first entry
    for(int j=1; j<e_joined.size(); j++){
        temp = e_joined[j];
        if(temp < min){
            min = temp;
            idxt = j;
        }
        if(temp > max){
            max = temp;
        }
    }
    e_fin(0) = e_joined[idxt];
    max += 1.0;
    e_joined[idxt] = max;
    c_fin.col(0) = c_joined.col(idxt);
    //now sort all other entries:
    min = max;
    for(int i=1; i<e_joined.size(); i++){
        for(int j=0; j<e_joined.size(); j++){
            temp = e_joined[j];
            if(temp < min){
                min = temp;
                idxt = j;
            }
        }
        e_fin[i] = e_joined[idxt];
        e_joined[idxt] = max;
        c_fin.col(i) = c_joined.col(idxt);
        min = max;
    }

    //outputting orbital energies and orbital coefficients
    cout << "Orbital Energies:" << endl;
    for(int i=0; i<e_joined.size(); i++){
        cout << "e" << i+1 << " = " << e_fin[i] << endl;
        cout << "coeff: " << endl << c_fin.col(i) << endl << endl;
    }
    // for(int i=0; i<e_joined.size(); i++){
    //     cout << "e" << i+1 << " = " << e_joined[i] << endl;
    //     cout << "coeff: " << endl << c_joined.col(i) << endl << endl;
    // }

    //parse results
    out.SCF_result(energy,e_joined,c_joined);
    return out;
}

//End of file
//****************************************************************
