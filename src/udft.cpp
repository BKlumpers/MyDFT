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

#include "dft.h"

//****************************************************************
//DFT energy minimisation:

/*SCF-loop for minimisation of the total energy
//
//{AO_list} is a vectors of CGFs comprising the basis set
//{pos_list} is a vector of the coordinates of the nuclei
//{charge_list} is a vector of the nuclear charges
//{nelec_list} is a vector of the number of electrons associated with each atom
//{atnum_list} is a vector of the atomic number for each atom
*/
SCF_E SCF_UDFT_energy(vector<CGF> AO_list, const vector<vec3>& pos_list, const vector<double>& charge_list, const vector<int>& nelec_list, const vector<int>& atnum_list)
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
    vector<double> Rep(two_electronSort(AO_list.size(), AO_list.size(), AO_list.size(), AO_list.size()), -1.0);
    if(Poisson == false || Quadrature == false){ //only calculate Rep for non-quadrature, or in case Poisson is not used
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
    Eigen::MatrixXd G_two_electron = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //Initialise two-electron Hamiltonian matrix
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

    //really inefficient grid (Hardcoded)
    static const double Nstep = 30.0; //have integral run in 3d{x,y,z} with Nstep+1 discrete steps along each axis
    static const int steps = int(Nstep) + 1; //turn loop-counter into int to support openmp
    static const double SingularTol = 1e-4; //tolerance for identifying singularities
    static const double DenTol = 1e-10;
    double Density = 0.0; //hold value of the local electron density
    double Density_alpha = 0.0;
    double Density_beta = 0.0;
    double Local = 0.0; //hold value of local AO-density
    double weight = 0.0; //hold value of the weight of the quadratic grid
    Eigen::MatrixXd Exchange_alpha = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //initialise value of the exchange energy integral
    Eigen::MatrixXd Correlation_alpha = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //initialise value of the correlation energy integral
    Eigen::MatrixXd Exchange_beta = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //initialise value of the exchange energy integral
    Eigen::MatrixXd Correlation_beta = Eigen::MatrixXd::Zero(AO_list.size(), AO_list.size()); //initialise value of the correlation energy integral
    double Ex = 0.0; //initialise exchange energy
    double Ec = 0.0; //initialise correlation energy
    double NELEC = 0.0; //initialise counter for electrons
    //define lower and upper bounds for x-axis
    double limset = 3.0;
    double xmin = -limset, xmax = limset;
    //define lower and upper bounds for y-axis
    double ymin = -limset, ymax = limset;
    //define lower and upper bounds for z-axis
    double zmin = -limset, zmax = limset;

    double dV = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)/pow(Nstep,3.0); //mesh-unit volume

    double xpos,ypos,zpos;
    vec3 origin; //mesh-point centre

    //Gauss-Chebyshev Lebedev grid:
    cout << "Generate GCL-grid: " << endl;
    grid GCL_grid = make_grid(pos_list, atnum_list, AO_list.size());
    cout << "Writing wavefunctions: " << AO_list.size() << endl;
    //add line for setting wavefunction amplitudes
    GCL_grid.write_amps(AO_list);
    cout << "Finished generating grid;" << endl;
    
    if(Poisson && Quadrature){
        cout << "Utilising Poisson-solver." << endl;
    }

    cout << "Beginning electronic optimisation: (UDFT) " << endl;
    cout << "Iter./Energy;" << endl;

    double tempos1, tempos2;

    //SCF-loop for energy optimisation:
    while(energy_difference > TolEnergy && loop_counter < loop_max) {
        loop_counter++; //increment loop counter

        Ex = 0.0; //reset energies
        Ec = 0.0;
        NELEC = 0.0;
        //Calculate two-electron hamiltonian matrix:
        if(!Quadrature || !Poisson){
            static int indexJ;
            for(int i=0; i<AO_list.size(); i++) {
                for(int j=0; j<AO_list.size(); j++) {
                    G_two_electron(i,j) = 0; //reset matrix
                    //compute exchange-correlation according to a Cartesian grid
                    if(!Quadrature){
                        Exchange_alpha(i,j) = 0.0;
                        Correlation_alpha(i,j) = 0.0;
                        Exchange_beta(i,j) = 0.0;
                        Correlation_beta(i,j) = 0.0;
                        for(int xx=0; xx<steps; xx++){
                            xpos = xmin + xx*(xmax-xmin)/Nstep;
                            for(int yy=0; yy<steps; yy++){
                                ypos = ymin + yy*(ymax-ymin)/Nstep;
                                for(int zz=0; zz<steps; zz++){
                                    zpos = zmin + zz*(zmax-zmin)/Nstep;
                                    origin << xpos,ypos,zpos;
                                    Density_alpha = density(origin, P_density_alpha, AO_list);
                                    Density_beta = density(origin, P_density_beta, AO_list);
                                    Density = Density_alpha + Density_beta;
                                    if(Density > DenTol){
                                        Local = AO_list[i].getvalue(origin)*AO_list[j].getvalue(origin);
                                        Exchange_alpha(i,j) += dV*exchange_potential(Density_alpha, Density_beta, true)*Local; //Dirac exchange functional
                                        Correlation_alpha(i,j) += dV*VWN_potential(Density_alpha, Density_beta, true)*Local; //Vosko-Wilk-Nusair correlation functional
                                        Exchange_beta(i,j) += dV*exchange_potential(Density_alpha, Density_beta, false)*Local; //Dirac exchange functional
                                        Correlation_beta(i,j) += dV*VWN_potential(Density_alpha, Density_beta, false)*Local; //Vosko-Wilk-Nusair correlation functional
                                        if(i == 0 && j == 0){
                                            //the additional nelec/2 is a scaling factor in restricted systems due to present looping scheme
                                            Ex += dV*exchange_Dirac(Density_alpha, Density_beta)*Density*(nelec/2); //add factor 3/4 to convert potential to energy
                                            Ec += dV*correlation_VWN(Density_alpha, Density_beta)*Density*(nelec/2);
                                            NELEC += dV*Density; //check if nelec-factor is needed here as well
                                        }
                                    }
                                }
                            }
                        }
                    }
                    for(int k=0; k<AO_list.size(); k++) {
                        for(int l=0; l<AO_list.size(); l++) {
                            indexJ = two_electronSort(i,j,l,k);
                            G_two_electron(i,j) += (P_density_alpha(k,l) + P_density_beta(k,l)) * Rep[indexJ];
                        }
                    }
                }
            }
        }

        //****************************************************************
        //
        //Only required change: express XC per spin density
        //                      rewrite loops to contain 2 Fock matrices
        //
        //  F_alpha = Kin + Nucl + P_ttl * J + XC_alpha
        //  F_beta  = Kin + Nucl + P_ttl * J + XC_beta
        //
        //  For   J = integral( rho_ttl * V_poisson )
        //
        //  Where V_poisson is still the same potential as calculated previously
        //
        //****************************************************************

        //compute exchange-correlation terms utilising numerical quadrature
        if(Quadrature){
            Exchange_alpha = Eigen::MatrixXd::Zero(Exchange_alpha.rows(), Exchange_alpha.cols()); //reset XC-matrices
            Correlation_alpha = Eigen::MatrixXd::Zero(Correlation_alpha.rows(), Correlation_alpha.cols());
            Exchange_beta = Exchange_alpha;
            Correlation_beta = Correlation_alpha;
            GCL_grid.set_density(P_density_alpha, P_density_beta); //update electron density map
            if(Poisson){
                G_two_electron = Exchange_alpha; //reset matrix
                GCL_grid.set_Poisson(); //update the Poisson potential
            }
            for(int rpnt=0; rpnt<GCL_grid.get_rpnt(); rpnt++){ //loop across all radial points
                for(int apnt=0; apnt<GCL_grid.get_apnt(); apnt++){ //loop across all angular points
                    for(int atom=0; atom<GCL_grid.get_atomnr(); atom++){ //loop across all atomic grids
                        //get electron density at this gridpoint
                        Density = GCL_grid.get_density(atom, rpnt, apnt);//density(atom, rpnt, apnt, P_density, GCL_grid);
                        Density_alpha = GCL_grid.get_density(atom, rpnt, apnt, true);
                        Density_beta = GCL_grid.get_density(atom, rpnt, apnt, false);
                        if(Density > DenTol){
                            weight = GCL_grid.get_weight(atom,rpnt,apnt);
                            Ex += weight*exchange_Dirac(Density_alpha, Density_beta)*Density; //add factor 3/4 to convert potential to energy
                            Ec += weight*correlation_VWN(Density_alpha, Density_beta)*Density;
                            NELEC += weight*Density;
                            for(int i=0; i<AO_list.size(); i++){ //change nesting order, which limits the number of times the density must be computed
                                for(int j=0; j<AO_list.size(); j++){
                                    Local = GCL_grid.get_amps(atom,rpnt,apnt,i)*GCL_grid.get_amps(atom,rpnt,apnt,j);
                                    Exchange_alpha(i,j) += weight*exchange_potential(Density_alpha,Density_beta,true)*Local; //Dirac exchange functional
                                    Correlation_alpha(i,j) += weight*VWN_potential(Density_alpha,Density_beta,true)*Local; //Vosko-Wilk-Nusair correlation functional
                                    Exchange_beta(i,j) += weight*exchange_potential(Density_alpha,Density_beta,false)*Local; //Dirac exchange functional
                                    Correlation_beta(i,j) += weight*VWN_potential(Density_alpha,Density_beta,false)*Local; //Vosko-Wilk-Nusair correlation functional
                                    if(Poisson){
                                        G_two_electron(i,j) += weight*Local*GCL_grid.get_Poisson(atom,rpnt,apnt);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        //Calculate Fock Matrix
        Fock_alpha = Hamil + G_two_electron + Exchange_alpha + Correlation_alpha;
        Fock_beta = Hamil + G_two_electron + Exchange_beta + Correlation_beta;
        
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

        Eigen::MatrixXd E_total = Hamil + 0.5 * G_two_electron; //non-weighted components of total electronic energy
        for(int i=0; i<AO_list.size(); i++) {
            for(int j=0; j<AO_list.size(); j++) {
                energy += (P_density_alpha(j,i) + P_density_beta(j,i)) * E_total(i,j); //weigh matrix elements
            }
        }
        energy += (Ex + Ec); //include XC-energies

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

    if(Poisson){
        cout << "PC check: " << endl;
        for(int i=0; i<pos_list.size(); i++){
            GCL_grid.PC_check(i);
        }
    }

    //sort molecular orbitals:
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
    double min, max, temp;
    int idxt;
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
    // cout << "Orbital Energies:" << endl;
    // for(int i=0; i<e_joined.size(); i++){
    //     cout << "e" << i+1 << " = " << e_fin[i] << endl;
    //     cout << "coeff: " << endl << c_fin.col(i) << endl << endl;
    // }
    
    //output SCF results
    if(loop_counter == loop_max){
        cout << "Stopping because maximum number of iterations was reached." << endl;
        cout << "Energy may not have converged." << endl;
    }
    else{
        cout << "Stopping because energy convergence was achieved." << endl;
    }
    cout << "System Energy = " << energy << endl;
    cout << "#Electrons = " << NELEC << endl << endl;

    //parse results
    out.SCF_result(energy,e_joined,c_joined);
    return out;
}


//****************************************************************
//Density Functional Theory:


//------------------>//merge below functions:


//calculate localised exchange-part of the DFT exchange-correlation potential
//according to the Dirac expression for a Homogeneous Electron Gas (LDA)
double exchange_Dirac(double Density_alpha, double Density_beta)
{
    static const double pi = 3.14159265359;
    static const double prefactor = -0.75*pow(3.0/pi,1.0/3.0);
    static const double root = sqrt(2.0);
    static const double form = pow(2.0,1.0/3.0) - 1.0;
    static double Density_ttl;
    Density_ttl = Density_alpha + Density_beta;
    static double rho;
    rho = pow(Density_ttl,1.0/3.0);
    static double chi;
    chi = (Density_alpha - Density_beta) / Density_ttl;
    static double f;
    f = (pow(1.0+chi,4.0/3.0) + pow(1.0-chi,4.0/3.0) - 2.0) / (2.0*form);
    return prefactor * rho * (f*form + 1.0);
}
//bool ab is a switch for alpha and beta spin, true for alpha, false for beta
double exchange_potential(double Density_alpha, double Density_beta, bool ab)
{
    static const double pi = 3.14159265359;
    static const double powr = 1.0/3.0;
    static const double prefactor = -0.75*pow(3.0/pi,powr);
    static const double root = sqrt(2.0);
    static const double form = pow(2.0,powr) - 1.0;
    static double Density_ttl;
    Density_ttl = Density_alpha + Density_beta;
    static double rho;
    rho = pow(Density_ttl,powr);
    static double nonspin;
    nonspin = prefactor * rho; // e_x^0
    static double chi;
    chi = (Density_alpha - Density_beta) / Density_ttl;
    static double f;
    f = (pow(1.0+chi,4.0*powr) + pow(1.0-chi,4.0*powr) - 2.0) / (2.0*form);
    static double total;
    total = (f*form + 1.0) * nonspin; //e_x
    static double formderiv;
    formderiv = ( pow(1.0+chi,powr) - pow(1.0-chi,powr) ) * 2.0 * powr/form;
    static double term3;
    if(ab){ //alpha case
        term3 = formderiv * 2.0 * Density_beta / Density_ttl;
    }
    if(!ab){ //beta case
        term3 = -formderiv * 2.0 * Density_alpha / Density_ttl;
    }
    return 4.0*powr*total + form * nonspin * term3;
}

//calculate correlation-part of the DFT exchange-correlation term
//according to the Vosko-Wilk-Nusair correlation functional
double correlation_VWN(double Density_alpha, double Density_beta)
{
    static const double pi = 3.14159265359;
    //Vosko-Wilk-Nusair V parameters (paramagnetic)
    static const double A_p = 0.0621814;
    static const double x0_p = -0.10498;
    static const double b_p = 3.72744;
    static const double c_p = 12.9352;
    //ferromagnetic:
    static const double A_f = A_p/2.0;
    static const double x0_f = -0.32500;
    static const double b_f = 7.06042;
    static const double c_f = 18.0578;   
    //spin-stiffness:
    static const double A_a = -1.0/(3.0*pi*pi);
    static const double x0_a = -0.0047584;
    static const double b_a = 1.13107;
    static const double c_a = 13.0045;
    //end
    static const double form = pow(2.0,1.0/3.0) - 1.0;
    static const double fpp0 = 2.25 * form; // this is 1.0 / fpp(0) | fpp(0) = 4.0/(9.0*form)
    
    double Density_ttl = Density_alpha + Density_beta;
    double chi = (Density_alpha - Density_beta) / Density_ttl;
    double f = (pow(1.0+chi,4.0/3.0) + pow(1.0-chi,4.0/3.0) - 2.0) / (2.0*form);

    double pre = 3.0/(4.0*pi*Density_ttl);
    double x_small = pow(pre,1.0/6.0); //x = sqrt(r0); r0 = (3/(4*pi*Density))^1/3 = pre^1/3
    double X_large_p = x_small*x_small + x_small*b_p + c_p;
    double X_large_f = x_small*x_small + x_small*b_f + c_f;
    double X_large_a = x_small*x_small + x_small*b_a + c_a;
    double Q_p = sqrt(4.0*c_p - b_p*b_p);
    double Q_f = sqrt(4.0*c_f - b_f*b_f);
    double Q_a = sqrt(4.0*c_a - b_a*b_a);
    double Fx_p = log(pow(x_small-x0_p,2.0)/X_large_p) + 2.0*(b_p+2.0*x0_p)*atan(Q_p/(b_p+2.0*x_small))/Q_p;
    double Fx_f = log(pow(x_small-x0_f,2.0)/X_large_f) + 2.0*(b_f+2.0*x0_f)*atan(Q_f/(b_f+2.0*x_small))/Q_f;
    double Fx_a = log(pow(x_small-x0_a,2.0)/X_large_a) + 2.0*(b_a+2.0*x0_a)*atan(Q_a/(b_a+2.0*x_small))/Q_a;
    double ec_p = 0.5*A_p*(    log(x_small*x_small/X_large_p) + 2.0*b_p*atan(Q_p/(2.0*x_small+b_p))/Q_p - b_p*x0_p*Fx_p/(x0_p*x0_p + x0_p*b_p + c_p)    );
    double ec_f = 0.5*A_f*(    log(x_small*x_small/X_large_f) + 2.0*b_f*atan(Q_f/(2.0*x_small+b_f))/Q_f - b_f*x0_f*Fx_f/(x0_f*x0_f + x0_f*b_f + c_f)    );
    double ec_a = 0.5*A_a*(    log(x_small*x_small/X_large_a) + 2.0*b_a*atan(Q_a/(2.0*x_small+b_a))/Q_a - b_a*x0_a*Fx_a/(x0_a*x0_a + x0_a*b_a + c_a)    );
    //multiply exchange potential with Density to get local value of exchange energy, prefactor 0.5 to convert from Rydberg to Hartree
    return ec_p + f * (   ec_a * fpp0 * (1.0 - pow(chi,4.0)) + pow(chi,4.0) * (ec_f - ec_p)   );
}

//calculate correlation-potential as the functional derivative
//of the Vosko-Wilk-Nusair correlation functional
double VWN_potential(double Density_alpha, double Density_beta, bool ab)
{
    static const double pi = 3.14159265359;
    //Vosko-Wilk-Nusair V parameters
    static const double A_p = 0.0621814;
    static const double x0_p = -0.10498;
    static const double b_p = 3.72744;
    static const double c_p = 12.9352;
    //ferromagnetic:
    static const double A_f = A_p/2.0;
    static const double x0_f = -0.32500;
    static const double b_f = 7.06042;
    static const double c_f = 18.0578;   
    //spin-stiffness:
    static const double A_a = -1.0/(3.0*pi*pi);
    static const double x0_a = -0.0047584;
    static const double b_a = 1.13107;
    static const double c_a = 13.0045;
    //end
    static const double form = pow(2.0,1.0/3.0) - 1.0;
    static const double fpp0 = 2.25 * form; // this is 1.0 / fpp(0) | fpp(0) = 4.0/(9.0*form)

    double Density_ttl = Density_alpha + Density_beta;
    double chi = (Density_alpha - Density_beta) / Density_ttl;
    double f = (pow(1.0+chi,4.0/3.0) + pow(1.0-chi,4.0/3.0) - 2.0) / (2.0*form);
    double formderiv = ( pow(1.0+chi,1.0/3.0) - pow(1.0-chi,1.0/3.0) ) * 2.0 / (3.0*form);

    double pre = 3.0/(4.0*pi*Density_ttl);
    double x_small = pow(pre,1.0/6.0); //x = sqrt(r0); r0 = (3/(4*pi*Density))^1/3
    double X_large_p = x_small*x_small + x_small*b_p + c_p;
    double X_large_f = x_small*x_small + x_small*b_f + c_f;
    double X_large_a = x_small*x_small + x_small*b_a + c_a;
    double Q_p = sqrt(4.0*c_p - b_p*b_p);
    double Q_f = sqrt(4.0*c_f - b_f*b_f);
    double Q_a = sqrt(4.0*c_a - b_a*b_a);
    double Fx_p = log(pow(x_small-x0_p,2.0)/X_large_p) + 2.0*(b_p+2.0*x0_p)*atan(Q_p/(b_p+2.0*x_small))/Q_p;
    double Fx_f = log(pow(x_small-x0_f,2.0)/X_large_f) + 2.0*(b_f+2.0*x0_f)*atan(Q_f/(b_f+2.0*x_small))/Q_f;
    double Fx_a = log(pow(x_small-x0_a,2.0)/X_large_a) + 2.0*(b_a+2.0*x0_a)*atan(Q_a/(b_a+2.0*x_small))/Q_a;
    double ec_p = 0.5*A_p*(    log(x_small*x_small/X_large_p) + 2.0*b_p*atan(Q_p/(2.0*x_small+b_p))/Q_p - b_p*x0_p*Fx_p/(x0_p*x0_p + x0_p*b_p + c_p)    );
    double ec_f = 0.5*A_f*(    log(x_small*x_small/X_large_f) + 2.0*b_f*atan(Q_f/(2.0*x_small+b_f))/Q_f - b_f*x0_f*Fx_f/(x0_f*x0_f + x0_f*b_f + c_f)    );
    double ec_a = 0.5*A_a*(    log(x_small*x_small/X_large_a) + 2.0*b_a*atan(Q_a/(2.0*x_small+b_a))/Q_a - b_a*x0_a*Fx_a/(x0_a*x0_a + x0_a*b_a + c_a)    );
    double dF_p = 2.0/(x_small-x0_p)-(2.0*x_small+b_p)/X_large_p-4.0*(2.0*x0_p+b_p)/(pow(2.0*x_small+b_p,2.0)+Q_p*Q_p);
    double dF_f = 2.0/(x_small-x0_f)-(2.0*x_small+b_f)/X_large_f-4.0*(2.0*x0_f+b_f)/(pow(2.0*x_small+b_f,2.0)+Q_f*Q_f);
    double dF_a = 2.0/(x_small-x0_a)-(2.0*x_small+b_a)/X_large_a-4.0*(2.0*x0_a+b_a)/(pow(2.0*x_small+b_a,2.0)+Q_a*Q_a);
    double deriv_p = 0.25*A_p*(2.0/x_small - (2.0*x_small+b_p)/X_large_p - 4.0*b_p/(pow(2.0*x_small+b_p,2.0)+Q_p*Q_p) - (b_p*x0_p/(x0_p*x0_p + x0_p*b_p + c_p)) * dF_p)/x_small;
    double deriv_f = 0.25*A_f*(2.0/x_small - (2.0*x_small+b_f)/X_large_f - 4.0*b_f/(pow(2.0*x_small+b_f,2.0)+Q_f*Q_f) - (b_f*x0_f/(x0_f*x0_f + x0_f*b_f + c_f)) * dF_f)/x_small;
    double deriv_a = 0.25*A_a*(2.0/x_small - (2.0*x_small+b_a)/X_large_a - 4.0*b_a/(pow(2.0*x_small+b_a,2.0)+Q_a*Q_a) - (b_a*x0_a/(x0_a*x0_a + x0_a*b_a + c_a)) * dF_a)/x_small;

    double deriv = deriv_p + f * (   deriv_a * fpp0 * (1.0 - pow(chi,4.0)) + pow(chi,4.0) * (deriv_f - deriv_p)   ); // d/dr_s of e_c
    double mix = 4.0 * chi * chi * chi * f * (ec_f - ec_p - fpp0 * ec_a) + formderiv * (ec_a * fpp0 * (1.0 - pow(chi,4.0)) + pow(chi,4.0) * (ec_f - ec_p)); // d/d chi of e_c

    static double term3;
    if(ab){ //alpha case
        term3 = 2.0 * Density_beta / Density_ttl;
    }
    if(!ab){ //beta case
        term3 = -2.0 * Density_alpha / Density_ttl;
    }

    return ec_p + f * (   ec_a * fpp0 * (1.0 - pow(chi,4.0)) + pow(chi,4.0) * (ec_f - ec_p)   ) - (x_small*x_small/3.0) * deriv + term3 * mix;
}

//End of file
//****************************************************************
