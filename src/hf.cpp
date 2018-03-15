/*****************************************************************

Basic closed-shell spin-restricted DFT-solver for simple molecules using STO-NG

Authors: B. Klumpers
         I.A.W. Filot

Published under GNU General Public License 3.0
Source code available at: https://github.com/BKlumpers/dft

Allows for SCF-computation of molecular energies for simple molecules.
Includes testcases for: H, He, H2, HeH+, He2, CO, and H2O.

*****************************************************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "integrals.h"
#include "orbitals.h"
#include "solvers.h"
#include "dft.h"

typedef Eigen::Vector3d vec3; //define vec3 as the Eigen::Vector3d-object

using namespace std;

const bool visual = false; //process visualisation of orbitals

//switches for setting which testcases to run
const bool H_switch = false;                     
const bool He_switch = false;
const bool mono = true;                    
const bool H2_switch = true;                    
const bool HeHp_switch = false;                 
const bool He2_switch = false;                  
const bool CO_switch = true;             
const bool H2O_switch = true;

//define functional call for visualisation procedure
void plot2D(char const title[], double height, double width, double xmin, double xmax, double zmin, double zmax, vector<CGF> AO_list, const Eigen::MatrixXd& Coeff);

//define the domain which will be visualised
double width = 500.0; //number of gridpoints in the x-direction
double height = 500.0; //number of gridpoints in the z-direction

//x runs from -3,+3 (sigma axis)
double xmin = -3.0, xmax = 3.0;
//z runs -3/+3 apart from each atom (internuclear axis)
double zmin = -3.0, zmax = 3.0;

//****************************************************************
//begin main codeblock

int main() {
    cout << "Initialising" << endl << endl;

    //input parameters for single atoms and simple diatomics:
    double H_Z = 1.0; //core charge for Hydrogen
    double He_Z = 2.0; //core charge for Helium
    double C_Z = 6.0; //core charge for Carbon
    double O_Z = 8.0; //core charge for Oxygen
    int H_nelec = 1; //1 electron on Hydrogen
    int He_nelec = 2; //2 electrons on Helium
    int C_nelec = 6; //6 electrons on Carbon
    int O_nelec = 8; //8 electons on Oxygen
    vec3 pos;
    pos << 0.0,0.0,0.0; //position nucleus at origin
    vec3 pos2;
    pos2 << 0.0,0.0,1.4; //position second atom 1.4A away from the 1st atom (H-H bond length)
    vec3 pos3;
    pos3 << 0.0,0.0,2.116; //position third atom 2.116A away from the 1st atom (C-O bond length in C=O)

    if(H_switch){
        //create contracted gaussian function STO-3G for Hydrogen
        CGF cgfH;
        cgfH.add_gto(0.15432897000000001,3.4252509099999999,0.0,0.0,0.0,pos);
        cgfH.add_gto(0.53532813999999995,0.62391373000000006,0.0,0.0,0.0,pos);
        cgfH.add_gto(0.44463454000000002,0.16885539999999999,0.0,0.0,0.0,pos);
        //Compute integrals for Hydrogen (only 1 electron, so no repulsion)
        double H_overlap = overlapCGF(cgfH,cgfH);
        double H_kinetic = kineticCGF(cgfH,cgfH);
        double H_nuclear = nuclearCGF(cgfH,cgfH,H_Z,pos);
        cout << "Testcase for Hydrogen: " << endl;
        cout << "H_overlap: " << H_overlap << endl;
        cout << "H_kinetic: " << H_kinetic << endl;
        cout << "H_nuclear: " << H_nuclear << endl;
        cout << "H_Energy: " << H_kinetic + H_nuclear << endl << endl;;
    }

    if(He_switch){ 
        //create CGF STO-3G for Helium
        CGF cgfHe;
        cgfHe.add_gto(0.15432897000000001,6.3624213899999997,0.0,0.0,0.0,pos);
        cgfHe.add_gto(0.53532813999999995,1.1589229999999999,0.0,0.0,0.0,pos);
        cgfHe.add_gto(0.44463454000000002,0.31364978999999998,0.0,0.0,0.0,pos);
        //compute integrals for Helium:
        double He_overlap = overlapCGF(cgfHe,cgfHe);
        double He_kinetic = kineticCGF(cgfHe,cgfHe);
        double He_nuclear = nuclearCGF(cgfHe,cgfHe,He_Z,pos);
        double He_repulsion = two_electronCGF(cgfHe,cgfHe,cgfHe,cgfHe);
        cout << "Testcase for Helium: " << endl;
        cout << "He_overlap: " << He_overlap << endl;
        cout << "He_kinetic: " << He_kinetic << endl;
        cout << "He_nuclear: " << He_nuclear << endl;
        cout << "He_repulsion: " << He_repulsion << endl;
        cout << "He_Energy: " << 2.0*He_kinetic + 2.0*He_nuclear + He_repulsion << endl << endl;
    }

    if(mono == true){
        //check SCF for single atom
        vector<CGF> AO_list1; //create object containing all atomic orbitals
        CGF cgfHe;
        cgfHe.add_gto(0.15432897000000001,6.3624213899999997,0.0,0.0,0.0,pos);
        cgfHe.add_gto(0.53532813999999995,1.1589229999999999,0.0,0.0,0.0,pos);
        cgfHe.add_gto(0.44463454000000002,0.31364978999999998,0.0,0.0,0.0,pos);
        AO_list1.push_back(cgfHe);
        vector<vec3> pos_list1; //create list of nucleic positions
        pos_list1.push_back(pos);
        vector<double> charge_list1; //create list of nucleic charges
        charge_list1.push_back(He_Z);
        vector<int> nelec_list1; //create list of electrons for each atom
        nelec_list1.push_back(He_nelec);
        vector<int> atnum_list1 = nelec_list1; //create list of atomic numbers

        cout << "SCF test: He " << endl << endl;

        SCF_E He1_results; //create output struct for H2 SCF energy minimisation
        He1_results = SCF_HF_energy(AO_list1,pos_list1,charge_list1, nelec_list1); //run SCF energy minimisation and get results
        SCF_E He1_DFT; //create output struct for H2 DFT energy minimisation
        He1_DFT = SCF_DFT_energy(AO_list1,pos_list1,charge_list1,nelec_list1,atnum_list1); //run SCF energy minimisation and get results

        CGF cgfBe;
        cgfBe.add_gto(0.15432897000000001,6.3624213899999997,0.0,0.0,0.0,pos);
        cgfBe.add_gto(0.53532813999999995,1.1589229999999999,0.0,0.0,0.0,pos);
        cgfBe.add_gto(0.44463454000000002,0.31364978999999998,0.0,0.0,0.0,pos);

        CGF cgfBe_1S;
        cgfBe_1S.add_gto(0.154329,30.167871,0.0,0.0,0.0,pos);
        cgfBe_1S.add_gto(0.535328,5.495115,0.0,0.0,0.0,pos);
        cgfBe_1S.add_gto(0.444635,1.487193,0.0,0.0,0.0,pos);

        CGF cgfBe_2S;
        cgfBe_2S.add_gto(-0.099967,1.314833,0.0,0.0,0.0,pos);
        cgfBe_2S.add_gto(0.399513,0.305539,0.0,0.0,0.0,pos);
        cgfBe_2S.add_gto(0.700115,0.099371,0.0,0.0,0.0,pos);
        
        CGF cgfBe_2PX;
        cgfBe_2PX.add_gto(0.155916,1.314833,1.0,0.0,0.0,pos);
        cgfBe_2PX.add_gto(0.607684,0.305539,1.0,0.0,0.0,pos);
        cgfBe_2PX.add_gto(0.391957,0.099371,1.0,0.0,0.0,pos);

        CGF cgfBe_2PY;
        cgfBe_2PY.add_gto(0.155916,1.314833,0.0,1.0,0.0,pos);
        cgfBe_2PY.add_gto(0.607684,0.305539,0.0,1.0,0.0,pos);
        cgfBe_2PY.add_gto(0.391957,0.099371,0.0,1.0,0.0,pos);

        CGF cgfBe_2PZ;
        cgfBe_2PZ.add_gto(0.155916,1.314833,0.0,0.0,1.0,pos);
        cgfBe_2PZ.add_gto(0.607684,0.305539,0.0,0.0,1.0,pos);
        cgfBe_2PZ.add_gto(0.391957,0.099371,0.0,0.0,1.0,pos);

        vector<CGF> AO_list_Be; //create object containing all atomic orbitals
        AO_list_Be.push_back(cgfBe_1S);
        AO_list_Be.push_back(cgfBe_2S);
        AO_list_Be.push_back(cgfBe_2PX);
        AO_list_Be.push_back(cgfBe_2PY);
        AO_list_Be.push_back(cgfBe_2PZ);

        vector<double> charge_list_Be; //create list of nucleic charges
        charge_list_Be.push_back(4);

        vector<int> nelec_list_Be; //create list of electrons for each atom
        nelec_list_Be.push_back(4);

        vector<int> atnum_list_Be = nelec_list_Be; //create list of atomic numbers

        cout << "SCF test: Be " << endl << endl;

        SCF_E Be_results; //create output struct for H2 SCF energy minimisation
        Be_results = SCF_HF_energy(AO_list_Be,pos_list1,charge_list_Be, nelec_list_Be); //run SCF energy minimisation and get results
        SCF_E Be_DFT; //create output struct for H2 DFT energy minimisation
        Be_DFT = SCF_DFT_energy(AO_list_Be,pos_list1,charge_list_Be,nelec_list_Be,atnum_list_Be); //run SCF energy minimisation and get results
    }

    if(H2_switch){
        //perform energy minimisation for H2:
        cout << "Testcase for H2: " << endl << endl;
        //create contracted gaussian function STO-3G for Hydrogen
        CGF cgfH_1;
        cgfH_1.add_gto(0.15432897000000001,3.4252509099999999,0.0,0.0,0.0,pos);
        cgfH_1.add_gto(0.53532813999999995,0.62391373000000006,0.0,0.0,0.0,pos);
        cgfH_1.add_gto(0.44463454000000002,0.16885539999999999,0.0,0.0,0.0,pos);
        //H2-molecule
        CGF cgfH_2; //create orbitals for second hydrogen atom at pos2
        cgfH_2.add_gto(0.15432897000000001,3.4252509099999999,0.0,0.0,0.0,pos2);
        cgfH_2.add_gto(0.53532813999999995,0.62391373000000006,0.0,0.0,0.0,pos2);
        cgfH_2.add_gto(0.44463454000000002,0.16885539999999999,0.0,0.0,0.0,pos2);

        vector<CGF> AO_list; //create object containing all atomic orbitals
        AO_list.push_back(cgfH_1);
        AO_list.push_back(cgfH_2);

        vector<vec3> pos_list; //create list of nucleic positions
        pos_list.push_back(pos);
        pos_list.push_back(pos2);

        vector<double> charge_list; //create list of nucleic charges
        charge_list.push_back(H_Z);
        charge_list.push_back(H_Z);

        vector<int> nelec_list; //create list of electrons for each atom
        nelec_list.push_back(H_nelec);
        nelec_list.push_back(H_nelec);

        vector<int> atnum_list = nelec_list; //create list of atomic numbers

        SCF_E H2_results; //create output struct for H2 SCF energy minimisation
        H2_results = SCF_HF_energy(AO_list,pos_list,charge_list, nelec_list); //run SCF energy minimisation and get results

        SCF_E H2_DFT; //create output struct for H2 DFT energy minimisation
        H2_DFT = SCF_DFT_energy(AO_list,pos_list,charge_list,nelec_list,atnum_list); //run SCF energy minimisation and get results

        //visualisation of H2 ground-state orbital
        if(visual){
            //write image of groundstate MO for H2 to file
            plot2D("H2.txt",height,width,xmin,xmax,zmin,zmax+1.4,AO_list,H2_results.C_vector);
        }
    }

    if(HeHp_switch){

        //HeH+ molecule
        cout << "Testcase for HeH+: " << endl << endl;
        //create CGF for Helium
        CGF cgfHe_1p;
        cgfHe_1p.add_gto(0.15432897000000001,6.3624213899999997,0.0,0.0,0.0,pos);
        cgfHe_1p.add_gto(0.53532813999999995,1.1589229999999999,0.0,0.0,0.0,pos);
        cgfHe_1p.add_gto(0.44463454000000002,0.31364978999999998,0.0,0.0,0.0,pos);
        //position hydrogen appropriate distance from helium:
        vec3 pos_HeHp;
        pos_HeHp << 0.0,0.0,0.772;
        CGF cgfH_He; //create orbitals for second hydrogen atom at pos2
        cgfH_He.add_gto(0.15432897000000001,3.4252509099999999,0.0,0.0,0.0,pos_HeHp);
        cgfH_He.add_gto(0.53532813999999995,0.62391373000000006,0.0,0.0,0.0,pos_HeHp);
        cgfH_He.add_gto(0.44463454000000002,0.16885539999999999,0.0,0.0,0.0,pos_HeHp);

        vector<CGF> AO_list_HeHp; //create object containing all atomic orbitals
        AO_list_HeHp.push_back(cgfHe_1p);
        AO_list_HeHp.push_back(cgfH_He);

        vector<vec3> pos_list_HeHp; //create list of nucleic positions
        pos_list_HeHp.push_back(pos);
        pos_list_HeHp.push_back(pos_HeHp);

        vector<double> charge_list_HeHp; //create list of nucleic charges
        charge_list_HeHp.push_back(He_Z);
        charge_list_HeHp.push_back(H_Z);

        vector<int> nelec_list_HeHp; //create list of electrons for each atom
        nelec_list_HeHp.push_back(He_nelec-1); //He+ has one electron less
        nelec_list_HeHp.push_back(H_nelec);

        vector<int> atnum_list_HeHp; //create list of atomic numbers
        atnum_list_HeHp.push_back(He_nelec);
        atnum_list_HeHp.push_back(H_nelec);

        //perform energy minimisation for HeH+:
        SCF_E HeHp_results;
        HeHp_results = SCF_HF_energy(AO_list_HeHp,pos_list_HeHp,charge_list_HeHp,nelec_list_HeHp);

        //perform energy minimisation through DFT:
        SCF_E HeHp_DFT;
        HeHp_DFT = SCF_DFT_energy(AO_list_HeHp,pos_list_HeHp,charge_list_HeHp,nelec_list_HeHp,atnum_list_HeHp);

        if(visual){
            //rescale z-axis to HeH+ bond length, then write image of ground-state MO
            plot2D("HeHp.txt",height,width,xmin,xmax,zmin,zmax+0.772,AO_list_HeHp,HeHp_results.C_vector);
        }
    }

    if(He2_switch){
        //He2 molecule
        cout << "Testcase for He2: " << endl << endl;

        //create CGF for Helium
        CGF cgfHe_1;
        cgfHe_1.add_gto(0.15432897000000001,6.3624213899999997,0.0,0.0,0.0,pos);
        cgfHe_1.add_gto(0.53532813999999995,1.1589229999999999,0.0,0.0,0.0,pos);
        cgfHe_1.add_gto(0.44463454000000002,0.31364978999999998,0.0,0.0,0.0,pos);
        //position second helium appropriate distance from first helium:
        vec3 pos_He2;
        pos_He2 << 0.0,0.0,5.2;
        CGF cgfHe2; //create orbitals for second helium atom at pos2
        cgfHe2.add_gto(0.15432897000000001,6.3624213899999997,0.0,0.0,0.0,pos_He2);
        cgfHe2.add_gto(0.53532813999999995,1.1589229999999999,0.0,0.0,0.0,pos_He2);
        cgfHe2.add_gto(0.44463454000000002,0.31364978999999998,0.0,0.0,0.0,pos_He2);

        vector<CGF> AO_list_He2; //create object containing all atomic orbitals
        AO_list_He2.push_back(cgfHe_1);
        AO_list_He2.push_back(cgfHe2);

        vector<vec3> pos_list_He2; //create list of nucleic positions
        pos_list_He2.push_back(pos);
        pos_list_He2.push_back(pos_He2);

        vector<double> charge_list_He2; //create list of nucleic charges
        charge_list_He2.push_back(He_Z);
        charge_list_He2.push_back(He_Z);

        vector<int> nelec_list_He2; //create list of electrons for each atom
        nelec_list_He2.push_back(He_nelec);
        nelec_list_He2.push_back(He_nelec);

        vector<int> atnum_list_He2 = nelec_list_He2; //create list of atomic numbers

        //perform energy minimisation for He2:
        SCF_E He2_results;
        He2_results = SCF_HF_energy(AO_list_He2,pos_list_He2,charge_list_He2,nelec_list_He2);

        //perform energy minimisation through DFT
        SCF_E He2_DFT;
        He2_DFT = SCF_DFT_energy(AO_list_He2,pos_list_He2,charge_list_He2,nelec_list_He2,atnum_list_He2);

        if(visual){
            //rescale to He2 bond length, then write image of ground-state MO to file
            plot2D("He2.txt",height,width,xmin,xmax,zmin,zmax+5.2,AO_list_He2,He2_results.C_vector);
        }
    }

    if(CO_switch){

        //CO molecule
        cout << "Testcase for CO: " << endl << endl;

        //create CGF for carbon
        CGF cgfC_1S;
        cgfC_1S.add_gto(0.154329,71.616837,0.0,0.0,0.0,pos);
        cgfC_1S.add_gto(0.535328,13.045096,0.0,0.0,0.0,pos);
        cgfC_1S.add_gto(0.444635,3.530512,0.0,0.0,0.0,pos);

        CGF cgfC_2S;
        cgfC_2S.add_gto(-0.099967,2.941249,0.0,0.0,0.0,pos);
        cgfC_2S.add_gto(0.399513,0.683483,0.0,0.0,0.0,pos);
        cgfC_2S.add_gto(0.700115,0.222290,0.0,0.0,0.0,pos);
        
        CGF cgfC_2PX;
        cgfC_2PX.add_gto(0.155916,2.941249,1.0,0.0,0.0,pos);
        cgfC_2PX.add_gto(0.607684,0.683483,1.0,0.0,0.0,pos);
        cgfC_2PX.add_gto(0.391957,0.222290,1.0,0.0,0.0,pos);

        CGF cgfC_2PY;
        cgfC_2PY.add_gto(0.155916,2.941249,0.0,1.0,0.0,pos);
        cgfC_2PY.add_gto(0.607684,0.683483,0.0,1.0,0.0,pos);
        cgfC_2PY.add_gto(0.391957,0.222290,0.0,1.0,0.0,pos);

        CGF cgfC_2PZ;
        cgfC_2PZ.add_gto(0.155916,2.941249,0.0,0.0,1.0,pos);
        cgfC_2PZ.add_gto(0.607684,0.683483,0.0,0.0,1.0,pos);
        cgfC_2PZ.add_gto(0.391957,0.222290,0.0,0.0,1.0,pos);
        
        //position oxygen appropriate distance from carbon:
        CGF cgfO_1S;
        cgfO_1S.add_gto(0.154329,130.709320,0.0,0.0,0.0,pos3);
        cgfO_1S.add_gto(0.535328,23.808861,0.0,0.0,0.0,pos3);
        cgfO_1S.add_gto(0.444635,6.443608,0.0,0.0,0.0,pos3);
        
        CGF cgfO_2S;
        cgfO_2S.add_gto(-0.099967,5.033151,0.0,0.0,0.0,pos3);
        cgfO_2S.add_gto(0.399513,1.169596,0.0,0.0,0.0,pos3);
        cgfO_2S.add_gto(0.700115,0.380389,0.0,0.0,0.0,pos3);
        
        CGF cgfO_2PX;
        cgfO_2PX.add_gto(0.155916,5.033151,1.0,0.0,0.0,pos3);
        cgfO_2PX.add_gto(0.607684,1.169596,1.0,0.0,0.0,pos3);
        cgfO_2PX.add_gto(0.391957,0.380389,1.0,0.0,0.0,pos3);
        
        CGF cgfO_2PY;
        cgfO_2PY.add_gto(0.155916,5.033151,0.0,1.0,0.0,pos3);
        cgfO_2PY.add_gto(0.607684,1.169596,0.0,1.0,0.0,pos3);
        cgfO_2PY.add_gto(0.391957,0.380389,0.0,1.0,0.0,pos3);
        
        CGF cgfO_2PZ;
        cgfO_2PZ.add_gto(0.155916,5.033151,0.0,0.0,1.0,pos3);
        cgfO_2PZ.add_gto(0.607684,1.169596,0.0,0.0,1.0,pos3);
        cgfO_2PZ.add_gto(0.391957,0.380389,0.0,0.0,1.0,pos3);

        vector<CGF> AO_list_CO; //create object containing all atomic orbitals
        AO_list_CO.push_back(cgfC_1S);
        AO_list_CO.push_back(cgfC_2S);
        AO_list_CO.push_back(cgfC_2PX);
        AO_list_CO.push_back(cgfC_2PY);
        AO_list_CO.push_back(cgfC_2PZ);
        AO_list_CO.push_back(cgfO_1S);
        AO_list_CO.push_back(cgfO_2S);
        AO_list_CO.push_back(cgfO_2PX);
        AO_list_CO.push_back(cgfO_2PY);
        AO_list_CO.push_back(cgfO_2PZ);

        vector<vec3> pos_list_CO; //create list of nucleic positions
        pos_list_CO.push_back(pos);
        pos_list_CO.push_back(pos3);

        vector<double> charge_list_CO; //create list of nucleic charges
        charge_list_CO.push_back(C_Z);
        charge_list_CO.push_back(O_Z);

        vector<int> nelec_list_CO; //create list of electrons for each atom
        nelec_list_CO.push_back(C_nelec);
        nelec_list_CO.push_back(O_nelec);

        vector<int> atnum_list_CO = nelec_list_CO; //create list of atomic numbers

        //perform energy minimisation for CO:
        SCF_E CO_results;
        cout << "CO_start (HF): " << endl;
        CO_results = SCF_HF_energy(AO_list_CO,pos_list_CO,charge_list_CO,nelec_list_CO);
        
        cout << "CO_start (DFT): " << endl;
        SCF_E CO_DFT;
        CO_DFT = SCF_DFT_energy(AO_list_CO,pos_list_CO,charge_list_CO,nelec_list_CO,atnum_list_CO);

        if(visual){
            //rescale to CO bond length, then write image of ground-state MO to file
            plot2D("CO.txt",height,width,xmin,xmax,zmin,zmax+2.116,AO_list_CO,CO_results.C_vector);
        }
    }

    if(H2O_switch){

        //H2O molecule
        cout << "Testcase for H2O: " << endl << endl;

        //create CGF for oxygen        
        CGF cgfOx_1S;
        cgfOx_1S.add_gto(0.154329,130.709320,0.0,0.0,0.0,pos);
        cgfOx_1S.add_gto(0.535328,23.808861,0.0,0.0,0.0,pos);
        cgfOx_1S.add_gto(0.444635,6.443608,0.0,0.0,0.0,pos);
        
        CGF cgfOx_2S;
        cgfOx_2S.add_gto(-0.099967,5.033151,0.0,0.0,0.0,pos);
        cgfOx_2S.add_gto(0.399513,1.169596,0.0,0.0,0.0,pos);
        cgfOx_2S.add_gto(0.700115,0.380389,0.0,0.0,0.0,pos);
        
        CGF cgfOx_2PX;
        cgfOx_2PX.add_gto(0.155916,5.033151,1.0,0.0,0.0,pos);
        cgfOx_2PX.add_gto(0.607684,1.169596,1.0,0.0,0.0,pos);
        cgfOx_2PX.add_gto(0.391957,0.380389,1.0,0.0,0.0,pos);
        
        CGF cgfOx_2PY;
        cgfOx_2PY.add_gto(0.155916,5.033151,0.0,1.0,0.0,pos);
        cgfOx_2PY.add_gto(0.607684,1.169596,0.0,1.0,0.0,pos);
        cgfOx_2PY.add_gto(0.391957,0.380389,0.0,1.0,0.0,pos);
        
        CGF cgfOx_2PZ;
        cgfOx_2PZ.add_gto(0.155916,5.033151,0.0,0.0,1.0,pos);
        cgfOx_2PZ.add_gto(0.607684,1.169596,0.0,0.0,1.0,pos);
        cgfOx_2PZ.add_gto(0.391957,0.380389,0.0,0.0,1.0,pos);

        //first hydrogen atom
        vec3 pos_H2O_1;
        pos_H2O_1 << 0.554,0.0,0.709;

        CGF cgfH_1x;
        cgfH_1x.add_gto(0.15432897000000001,3.4252509099999999,0.0,0.0,0.0,pos_H2O_1);
        cgfH_1x.add_gto(0.53532813999999995,0.62391373000000006,0.0,0.0,0.0,pos_H2O_1);
        cgfH_1x.add_gto(0.44463454000000002,0.16885539999999999,0.0,0.0,0.0,pos_H2O_1);

        //second hydrogen atom
        vec3 pos_H2O_2;
        pos_H2O_2 << 0.554,0.0,-0.709;

        CGF cgfH_2x;
        cgfH_2x.add_gto(0.15432897000000001,3.4252509099999999,0.0,0.0,0.0,pos_H2O_2);
        cgfH_2x.add_gto(0.53532813999999995,0.62391373000000006,0.0,0.0,0.0,pos_H2O_2);
        cgfH_2x.add_gto(0.44463454000000002,0.16885539999999999,0.0,0.0,0.0,pos_H2O_2);

        vector<CGF> AO_list_H2O; //create object containing all atomic orbitals
        AO_list_H2O.push_back(cgfOx_1S);
        AO_list_H2O.push_back(cgfOx_2S);
        AO_list_H2O.push_back(cgfOx_2PX);
        AO_list_H2O.push_back(cgfOx_2PY);
        AO_list_H2O.push_back(cgfOx_2PZ);
        AO_list_H2O.push_back(cgfH_1x);
        AO_list_H2O.push_back(cgfH_2x);

        vector<vec3> pos_list_H2O; //create list of nucleic positions
        pos_list_H2O.push_back(pos);
        pos_list_H2O.push_back(pos_H2O_1);
        pos_list_H2O.push_back(pos_H2O_2);

        vector<double> charge_list_H2O; //create list of nucleic charges
        charge_list_H2O.push_back(O_Z);
        charge_list_H2O.push_back(H_Z);
        charge_list_H2O.push_back(H_Z);

        vector<int> nelec_list_H2O; //create list of electrons for each atom
        nelec_list_H2O.push_back(O_nelec);
        nelec_list_H2O.push_back(H_nelec);
        nelec_list_H2O.push_back(H_nelec);

        vector<int> atnum_list_H2O = nelec_list_H2O; //create list of atomic numbers

        //perform energy minimisation for He2:
        SCF_E H2O_results;
        cout << "H2O_start (HF): " << endl;
        H2O_results = SCF_HF_energy(AO_list_H2O,pos_list_H2O,charge_list_H2O,nelec_list_H2O);
        
        cout << "H2O_start (DFT): " << endl;
        SCF_E H2O_DFT;
        H2O_DFT = SCF_DFT_energy(AO_list_H2O,pos_list_H2O,charge_list_H2O,nelec_list_H2O,atnum_list_H2O);

        if(visual){
            //rescale to CO bond length, then write image of ground-state MO to file
            plot2D("H2O.txt",2*height,2*width,xmin-0.554,xmax+0.554,zmin-0.709,zmax+0.709,AO_list_H2O,H2O_results.C_vector);
        }
    }

    cout << "Program has ended" << endl;
    return 0;
}

//end main codeblock
//****************************************************************

//generate an image containing the amplitudes of the ground-state molcular orbital
//this image can be used to generate a colourmap of the MO in an appropriate visualisation tool
//{title} will be the filename
//{AO_list} is a vector of the atomic orbitals
//{Coeff} are the optimised orbital coefficients
//remaining parameters define the domain in which the MO is visualised
void plot2D(char const title[], double height, double width, double xmin, double xmax, double zmin, double zmax, vector<CGF> AO_list, const Eigen::MatrixXd& Coeff)
{
    //visualise the xz-plane
    //
    vec3 origin; //mesh-point centre
    double xpos,ypos,zpos; //initialise position vector components
    ypos = 0.0; //visualise the plane corresponding to y=0
    int orbital_vis = 0; //visualise ground-state orbital
    Eigen::MatrixXd image = Eigen::MatrixXd::Zero(width+1,height+1); //initialise matrix storing pixels

    cout << endl << "Commencing writing to file: " << title << endl;

    //compute amplitude at each point
    for(int i=0; i<(width+1); i++){
        xpos = xmin + i*(xmax-xmin)/width;
        for(int j=0; j<(height+1); j++){
            zpos = zmin + j*(zmax-zmin)/height;
            //j loops across z, i loops across x
            origin << xpos,ypos,zpos;
            for(int orbs=0; orbs<AO_list.size(); orbs++){ //loop across all atomic orbitals
                image(i,j) += AO_list[orbs].getvalue(origin)*Coeff(orbs,orbital_vis);
            }
        }
    }
    // if(string::empty(title)){ //requires #include <string>
    //     title = "2D_plot.txt"; //default filename
    // }
    std::ofstream file(title); //set file as output stream
    file << image << '\n'; //write image to file
    cout << "Finished writing." << endl << endl;
    return;
}

//End of file
//****************************************************************
