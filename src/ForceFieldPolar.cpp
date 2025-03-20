/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: Juen. 10, 2022                                                *
 * -------------------------------------------------------------------------- */

#include "ForceFieldPolar.h"
#include "fftw3.h"

void ForceFieldPolar::init() {
    // * 1 Initialize data members and Hamiltonian.
    ForceFieldBase::init();
    std::cout<<"this is the force field polar part"<<std::endl;
    // * 2 choose whitch TTM0, TTM1, TTM2, TTM3, TTM4
    damping_type = param.getStr("damping_type");
    forcefield_type = param.getStr("forcefield_type");
    iter_max = param.getInt("iter_max");
    TOLERANCE = param.getDouble("TOLERANCE");
    setForceFieldPolarParameters();
    Mu_all.resize(DOFn * 3);
    Pi_all.resize(DOFn * 9);
    Dipole_mol.resize(3);
    cosleft = 0;
    cosright = 0;
    Rdipoleleft = 0;
    numRight = 0;
    numRight = 0;
    numLeft = 0;
}

void ForceFieldPolar::setForceFieldPolarParameters() {
    // set the No palar force field parameters from topology
    // Now construct the polaribility parameters for the polar force field
    int currentNumAtoms = 0;
    setFalseMdPerturb();
    for (int i = 0; i < topologies.molecules.size(); i++) {
        for (int j = 0; j < topologies.molecules[i].moleculeNumber; j++) {
            const Topology::MoleculeTypes& moleculeType = topologies.moleculeTypes[topologies.molecules[i].moleculeTypeIndex]; 
                         
            for (int k = 0; k < moleculeType.polars.size(); k++) {
               int indexI = k + currentNumAtoms;
               double alpha = moleculeType.polars[k].polarizability; //nm3
               ffPolars.push_back({indexI, alpha});     
            } 
            // * Add Bonds to .
            for (int k = 0; k < moleculeType.bonds.size(); k++) {
               int indexI = (moleculeType.bonds[k].atomIndex[0] - 1) + currentNumAtoms;
               int indexJ = (moleculeType.bonds[k].atomIndex[1] - 1) + currentNumAtoms;
               exclusions12.push_back({indexI, indexJ});
            }
            // * Add angles to .
            for (int k = 0; k < moleculeType.angles.size(); k++) {
                int indexI = (moleculeType.angles[k].atomIndex[0] - 1) + currentNumAtoms;
                int indexK = (moleculeType.angles[k].atomIndex[2] - 1) + currentNumAtoms;
                exclusions13.push_back({indexI, indexK});
            }
            // * In order to calculte the currentNumAtoms
            for (int k = 0; k < moleculeType.polars.size(); k++) {
                currentNumAtoms = currentNumAtoms + 1;   
            }
        }
    } 
}

void ForceFieldPolar::calculateDipoleRleftAright() {
    // calculate the dipole size and direction.
    // This is only for dipole size and direction.
    // intial: MolI    :  2236
    // molecule: MolN  : Urea = 200
    // molecule: MolS  : Urea = 8
    // DirecA 
    // DirecB
    // std::cout<<"Test for the dipole calculation"<<std::endl;
    cosleft = 0;
    cosright = 0;
    Rdipoleleft = 0;
    Rdipoleright = 0;
    Vec3 dipole_right;
    Vec3 dipole_left;
    numLeft = 0;
    numRight = 0;
    dirLx = 0;
    dirRx = 0;
    int MolI = param.getInt("MolI");
    int MolN = param.getInt("MolN");
    int MolS = param.getInt("MolS");
    const int DirecA = param.getInt("DirecA");
    const int DirecB = param.getInt("DirecB");
    // This is to calculate the Dipole with R between Ca/Cb and molecule (COM): cos and size.
    // There are three parts: only in the C153 first shell.
    // C153: calculate COM, Ca/Cb with Ch+(COM) and HBD(COM) distance
    // Ch+: calcualte COM, dipole,
    // HBD: calculate COM, dipole
    // C153: molecule=1, atom=36, Ca is and Cb is
    // identify Ca is left, Cb is right.
    Vec3 COM0;
    Vec3 Ca; // 11
    Vec3 Cb; // 14
    Vec3 COM1;
    Vec3 direcX;
    Vec3 direcY;

    for (int a = 0; a < 3; a++) {
        Ca[a] = R[DirecA-1][a];
        Cb[a] = R[DirecB-1][a];
        dipole_right[a] = 0;
        dipole_left[a] = 0;
    }

    for (int a = 0; a < 3; a++) {
        direcX[a] = Cb[a] - Ca[a];
        COM0[a] = (Ca[a] + Cb[a])/2;
    }
    
    bool direction = false; // true: left, false: right
    for (int i = 0; i < MolN; i++) {
        // input Ra, Rb, COM1
        for (int a = 0; a < 3; a++) {
            COM1[a] = 0; 
        }
        double mass1 = 0;
        for (int a = 0; a < MolS; a++) {
            for (int b = 0; b < 3; b++) {
                COM1[b] += R[MolI + (i-1) * MolS + a][b];  
            } 
            //mass1 += ffElec[MolI + (i-1) * MolS + a].mass;
        }
        
        for (int a = 0; a < 3; a++) {
            COM1[a] = COM1[a] / MolS; 
        }

        double dis = 0;
        dis = sqrt((COM0[0]-COM1[0]) * (COM0[0]-COM1[0]) + (COM0[1]-COM1[1]) * (COM0[1]-COM1[1]) + (COM0[2]-COM1[2]) * (COM0[2]-COM1[2]));
        if (dis < 1) {
            double disleft = (Ca[0] - COM1[0]) * (Ca[0] - COM1[0]) + (Ca[1] - COM1[1]) * (Ca[1] - COM1[1]) + (Ca[2] - COM1[2]) * (Ca[2] - COM1[2]);
            double disright = (Cb[0] - COM1[0]) * (Cb[0] - COM1[0]) + (Cb[1] - COM1[1]) * (Cb[1] - COM1[1]) + (Cb[2] - COM1[2]) * (Cb[2] - COM1[2]);
            if (disleft < disright) direction = true;
            // calculate the molecule dipole
            Vec3 dipoleMol;
            for (int a = 0; a < 3; a++) {
                dipoleMol[a] = 0; 
            }
            for (int a = 0; a < MolS; a++) {
                for (int b = 0; b < 3; b++) {
                    if (MolI == 36) dipoleMol[b] += R[MolI + (i-1) * MolS + a][b] * (ffElec[MolI + (i-1) * MolS + a].charge - 1/MolS);
                    else {
                        dipoleMol[b] += R[MolI + (i-1) * MolS + a][b] * ffElec[MolI + (i-1) * MolS + a].charge;
                    }
                }
            }
            double dipoleValue = 0;
            dipoleValue = dipoleMol[0] * dipoleMol[0] + dipoleMol[1] * dipoleMol[1] + dipoleMol[2] * dipoleMol[2];
            if (direction) { // molecule on the left
                numLeft++;
                Rdipoleleft += (COM1[0] - Ca[0]) * dipoleMol[0] + (COM1[1] - Ca[1]) * dipoleMol[1] + (COM1[2] - Ca[2]) * dipoleMol[2];
                cosleft += ((COM1[0] - Ca[0]) * dipoleMol[0] + (COM1[1] - Ca[1]) * dipoleMol[1] + (COM1[2] - Ca[2]) * dipoleMol[2])/sqrt(dipoleValue * disleft);
                for (int a = 0; a < 3; a++) {
                    dipole_left[a] += dipoleMol[a];
                }
            }
            else { // molecule on the right
                numRight++;
                Rdipoleright += (COM1[0] - Cb[0]) * dipoleMol[0] + (COM1[1] - Cb[1]) * dipoleMol[1] + (COM1[2] - Cb[2]) * dipoleMol[2];
                cosright += ((COM1[0] - Cb[0]) * dipoleMol[0] + (COM1[1] - Cb[1]) * dipoleMol[1] + (COM1[2] - Cb[2]) * dipoleMol[2])/sqrt(dipoleValue * disright);
                for (int a = 0; a < 3; a++) {
                    dipole_right[a] += dipoleMol[a];
                }
            }           
        }
    } 
    cosleft = cosleft / numLeft;
    cosright = cosright / numRight;
    for (int a = 0; a < 3; a++) {
        dipole_left[a] = dipole_left[a] / numLeft;
        dipole_right[a] = dipole_right[a] / numRight;
    }
    double xx = 0, yy = 0, zz = 0;
    xx = dipole_left[0] * direcX[0] + dipole_left[1] * direcX[1] + dipole_left[2] * direcX[2];
    yy = sqrt(dipole_left[0] * dipole_left[0] + dipole_left[1]*dipole_left[1]+dipole_left[2]*dipole_left[2]);
    zz = yy * sqrt(direcX[0]*direcX[0]+direcX[1]*direcX[1]+direcX[2]*direcX[2]);
    dirLx = yy * xx/zz; 
    xx = dipole_right[0] * direcX[0] + dipole_right[1] * direcX[1] + dipole_right[2] * direcX[2];
    yy = sqrt(dipole_right[0] * dipole_right[0] + dipole_right[1]*dipole_right[1]+dipole_right[2]*dipole_right[2]);
    zz = yy * sqrt(direcX[0]*direcX[0]+direcX[1]*direcX[1]+direcX[2]*direcX[2]);
    dirRx = yy * xx/zz;
}

void ForceFieldPolar::calculatePiTensor() { //gongxu
    // calculate the many-body polarizability as a function of configuration using iteration method
    // to iter_max order or already converged. Return actual iteration times.
    double pi1[9 * DOFn];
    double pi2[9 * DOFn];
    std::cout<<"begin to calculatePiTensor"<<std::endl;
    double rij[3]={0,0,0};
    double Pi_prev[3][3]={0,0,0,0,0,0,0,0,0};
    double Tij[3][3]={0,0,0,0,0,0,0,0,0};
    double m1[3][3]={0,0,0,0,0,0,0,0,0};  //for swap matrix
    double m2[3][3]={0,0,0,0,0,0,0,0,0};  //for swap matrix
    double tol(1); //tolerance
    int n(0);  //iteration count
    // For Thole model
    double r, r2, rm1, rm2, rm3, rm5;
    double alphaij, Ralphaij, Ralphaij2, Ralphaij4;
    double S0r, S1r, S2r, S3r;
    std::vector<double> S1rP;
    std::vector<double> S2rP;
    std::cout.precision(12);
    double alpha[DOFn];
    S1rP.resize(DOFn*DOFn, 0);
    S2rP.resize(DOFn*DOFn, 0);
    // Set the atom polarizability
    for (int i = 0; i < DOFn; i++) {
       alpha[i] = ffPolars[i].alpha;
    }
    // Initialize Pi_tensor to be 0
    for (int a = 0; a < 3; a++) Pi[a][0] = Pi[a][1] = Pi[a][2] = 0; 
    
    // Initial guess as zero-order alpha
    for ( int i = 0; i < DOFn; i++) { 
       pi2[9 * i] = pi2[9 * i + 4] = pi2[9 * i + 8] = alpha[i];
       pi2[9 * i + 1] = pi2[9 * i + 2] = pi2[9 * i +  + 3] = pi2[9 * i + 5] = pi2[9 * i + 6] = pi2[9 * i + 7] = 0;
    }
    
    std::vector<std::set<int> > bonded12(DOFn); 
    std::vector<std::set<int> > bonded13(DOFn);
    if (damping_type == "TTM3"||damping_type == "TTM4") {
        addTholeExclusionsToSet(bonded12, bonded13);
    }
    for (int i = 0; i < DOFn; i++) {
        for (int j = 0; j < DOFn; j++) {
            if (i == j) {
                continue;
            }
            //calculate Tij
            double deltaR[5];
            const Vec3& atomCoordinatesI = R[i];
            const Vec3& atomCoordinatesJ = R[j];
            // get deltaR 2 atoms, deltaR = Ri - Rj
            GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
            rm1 = 1.0 / deltaR[4];
            rm2 = rm1 * rm1;
            rm3 = rm1 * rm2;
            rm5 = rm2 * rm3;
            // Thole model daming parameter for chosen pair od atom i and j
            // Judge the pair i, j belong to the 12, 13, 14, or noboned part.
            double dampingParam;
            dampingParam = param.getDouble("dampingParam3");
            if (bonded12[i].find(j) == bonded12[i].end()||bonded12[j].find(i) == bonded12[j].end()) 
                dampingParam = param.getDouble("dampingParam1");
            if (bonded13[i].find(j) == bonded13[i].end()||bonded13[j].find(i) == bonded13[j].end()) 
                dampingParam = param.getDouble("dampingParam2");
            // choose whitch thole model damping type.
            // std::cout<<"i  "<<i<<"  j:  "<<j<<"  dampingParam: "<<dampingParam<<std::endl;
            if (damping_type == "TTM0") 
                TTM0Model(S0r, S1r, S2r, S3r);
            else if (damping_type == "TTM1") 
                TTM1Model(S0r, S1r, S2r, S3r);
            //else if (damping_type == "TTM2") 
            //    double alphaEwald, double deltaR, double& S0r, double& S1r, double& S2r, double& S3r) {
            else {
                alphaij = alpha[i] * alpha[j];
                alphaij = pow((double)alphaij,double(1.0/6));   
                alphaij = deltaR[4] * 1.0 / alphaij;
                if (damping_type == "TTM3")
                    TTM3Model(dampingParam, alphaij, S0r, S1r, S2r, S3r);
                else if (damping_type == "TTM4")
                    TTM4Model(dampingParam, alphaij, S0r, S1r, S2r, S3r);
            }
            S1rP[i * DOFn + j] = S1r;
            S2rP[i * DOFn + j] = S2r;
        }
    }
    while ( tol > TOLERANCE && n< iter_max ) {
        n++;  
        // Save Pi_prev<=Pi
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3;b++) {
                Pi_prev[a][b] = Pi[a][b];
            }
        } 
        // Save (old) pi1 <= (new) pi2, and set pi2 zero
        for (int i = 0; i < DOFn; i++) {
            for (int k = 0; k < 9; k++) {
                pi1[9 * i + k] = pi2[9 * i + k];
                pi2[9 * i + k] = 0;
            }
        }
        // Start new iteration for pi2[i=1...N][9]
        for (int i = 0; i < DOFn; i++) {
            for (int j = 0; j < DOFn; j++) {
                if (i == j) continue;
                //calculate Tij
                double deltaR[5];
                const Vec3& atomCoordinatesI = R[i];
                const Vec3& atomCoordinatesJ = R[j];
                // get deltaR 2 atoms, deltaR = Ri - Rj
                GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
                rm1 = 1.0 / deltaR[4];
                rm2 = rm1 * rm1;
                rm3 = rm1 * rm2;
                rm5 = rm2 * rm3;
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {                   
                        Tij[a][b] = 3 * deltaR[a] * deltaR[b] * rm5 * S2rP[i * DOFn + j];
                        if (a == b){
                        Tij[a][b] -= Delta(a,b) * rm3 * S1rP[i * DOFn + j]; 
                        }
                    }
                }                                                               
                // Transfrom pi1[j][9] to a matrix m1[3][3]
                int k = 0;
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {
                        m1[a][b] = pi1[j * 9 + k];
                        k++;
                    }
                }
                MatrixMultiplyPi(Tij, m1, m2);
                // Add matrix m2[3][3] to pi2[i][9]
                k = 0;
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {
                        pi2[i * 9 + k] += m2[a][b];
                        k++;
                    }
                } 
            }
            pi2[i* 9 + 0] += 1;
            pi2[i* 9 + 4] += 1;
            pi2[i* 9 + 8] += 1; 
            for (int k = 0; k < 9; k++) {
                pi2[i * 9 + k] *= alpha[i];
            }   
        } 
        // Calculate many-body polarizability Pi[3][3] (sum up i=1...N)
        // Initialize result tensor to be 0
        for (int a = 0; a < 3; a++) Pi[a][0] = Pi[a][1] = Pi[a][2] = 0; 
        for (int i = 0; i < DOFn; i++) {
            int k = 0;
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++) {
                    Pi[a][b] += pi2[i * 9 + k];
                    k++;
                }
        }
        // Calculate tolerance
        tol = 0;
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                tol += abs(Pi[a][b] - Pi_prev[a][b]);
            }
        }
    }
    // Save pi_i to Pi_all
    for (int i = 0; i < DOFn; i++)
       for (int j = 0; j < 9; j++) {
          Pi_all[i * 9 + j] = pi2[i * 9 + j];
       }                               
}

void ForceFieldPolar::calculateTTM4PolarForceandEnergy(std::vector<Vec3>& forces, 
                                                std::vector<Vec3>& polarforce,
                                                double& totalenergy, 
                                                double& polarenergy, bool includeForces) {
    // We will compare the energy and force with the MBX software.
    // The electrostatic energy can then be written as a sum of four terms, the real, the reciprocal, 
    // the adjusted and the self terms, and the reciprocal part energy and force use the PME methods 
    // to calculate. The system has charges and dipoles.
    polarforce.resize(DOFn, Vec3());
    std::fill(polarforce.begin(), polarforce.end(), Vec3());
    polarenergy = 0.0;
    double totalSelfEwaldEnergy = 0.0;
    double realSpaceEwaldEnergy = 0.0;
    double recipEnergy = 0.0;
    double totalRecipEnergy = 0.0;
    // set the Ewald method some important parameters
    // For Ewald Sum method parameters: alphaEwald                                                                                     
    double Ewald_tolerance = param.getDouble("Ewald_tolerance");
    double alphaEwald;
    double cutlength = param.getDouble("cutoff");;
    alphaEwald = sqrt(-log(2*Ewald_tolerance))/cutlength; 
    double recipCoeff = ONE_4PI_EPS0 * 4 * pi/(periodicBoxVectors[0][0] * periodicBoxVectors[1][1] * periodicBoxVectors[2][2]);
    /**
     * The self part energy.
     */ 
    for (int ii = 0; ii < DOFn; ii++){  
        double charge = ffElec[ii].charge;
        double mu2 = Mu_all[ii*3+0]*Mu_all[ii*3+0] + Mu_all[ii*3+1]*Mu_all[ii*3+1] + Mu_all[ii*3+2]*Mu_all[ii*3+2];
        totalSelfEwaldEnergy -= alphaEwald*ONE_4PI_EPS0/sqrpi*(charge*charge + 2*alphaEwald*alphaEwald*mu2/3);
    } 
    polarenergy += totalSelfEwaldEnergy;    
    /**
     * Real SPACE EWALD ENERGY AND FORCES
     * T = s[0] * 1/r;
     * T[a] = - s[1] * rij[a]/r^3;
     * T[a][b] = s[2] * 3*rij[a]*rij[b]/r^5 -s[1] * delta[ab]/r^3;
     * T[a][b][c] = -s[3]*15*rij[a]*rij[b]*rij[c]/r^7 + s[2]*3*(rij[a]*delta[bc]+rij[b]*delta[ac]+rij[c]*delta[ab])/r^5
     * Using the TTM-4 or TTM3 damping model method.
     */
    double rm1, rm2, rm3, rm5, rm7;
    /** 
     * There is an important thing we should focus on the atom i and j, 
     * the pairs should included in 12, 13, 14, and nonbonded
     * the unit problems should also inportant for the energy and force
     * energy: KJ/mol 
     * force: KJ/mol/nm
     * This is set the pair 12,13, and nonbonded part.
    */
    std::vector<std::set<int> > bonded12(DOFn); 
    std::vector<std::set<int> > bonded13(DOFn);
    addTholeExclusionsToSet(bonded12, bonded13);
    // Set the atom polarizability
    double alpha[DOFn];
    for (int i = 0; i < DOFn; i++) {
       alpha[i] = ffPolars[i].alpha;
    }

    NeighborList neighborList;
    const int nofE = ffExclusions.size();
    std::vector<std::set<int> > Exclusions(nofE);
    std::vector<std::set<int> > bonded(nofE);
    addExclusionsToSet(nofE, bonded, Exclusions);
    bool usePeriodic = true;
    double maxDistance = cutlength;
    bool reportSymmetricPairs = false;
    double minDistance = 0.0;  
    computeNeighborListVoxelHash(neighborList, DOFn, R, Exclusions, periodicBoxVectors, 
                                usePeriodic, maxDistance, minDistance , reportSymmetricPairs);  

    for (auto& pair : neighborList) {
        int ii = pair.first;
        int jj = pair.second;  
        double chargeI = ffElec[ii].charge;
        double chargeJ = ffElec[jj].charge;
        double deltaR[5];
        const Vec3& atomCoordinatesI = R[ii];
        const Vec3& atomCoordinatesJ = R[jj];
        // Get deltaR 2 atoms, deltaR = Ri - Rj
        GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
        rm1 = 1.0 / deltaR[4];
        rm2 = rm1 * rm1;
        rm3 = rm1 * rm2;
        rm5 = rm2 * rm3;
        rm7 = rm5 * rm2;
        double alphaij, Ralphaij2, Ralphaij4;
        alphaij = alpha[ii] * alpha[jj];
        alphaij = pow((double)alphaij,double(1.0/6));            
        alphaij = deltaR[4] * 1.0 / alphaij;
        double S0r, S1r, S2r, S3r;
        // Thole model daming parameter for chosen pair od atom i and j
        double dampingParam;
        // Judge the pair i, j belong to the 12, 13, 14, or noboned part.
        dampingParam = param.getDouble("dampingParam3");
        if (bonded12[ii].find(jj) == bonded12[ii].end()||bonded12[jj].find(ii) == bonded12[jj].end()) 
            dampingParam = param.getDouble("dampingParam1");
        if (bonded13[ii].find(jj) == bonded13[ii].end()||bonded13[jj].find(ii) == bonded13[jj].end()) 
            dampingParam = param.getDouble("dampingParam2");
        // This is the TTM-4 model Exponential damping compare with the MBX software
        TTM4Model(dampingParam, alphaEwald, S0r, S1r, S2r, S3r);
        // This is to struct the T, Ta[a], Tab[a][b], Tabc[a][b][c]
        double T, Ta[3], Tab[3][3], Tabc[3][3][3];
        T = S0r * rm1;
        for (int a = 0; a < 3; a++) {
            Ta[a] = - S1r * deltaR[a] * rm3;
        }
        for (int a = 0; a < 3; a++) 
            for (int b = 0; b < 3; b++) {                   
                Tab[a][b] = 3 * deltaR[a] * deltaR[b] * rm5 * S2r - S1r * rm3 * Delta(a, b);
           } 
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                for (int c = 0; c < 3; c++) {
                    Tabc[a][b][c] = - 15 * deltaR[a] * deltaR[b] * deltaR[c] * S3r * rm7;
                    Tabc[a][b][c] += S2r * 3 * deltaR[a] * rm5 * Delta(b, c); 
                    Tabc[a][b][c] += S2r * 3 * deltaR[b] * rm5 * Delta(a, c) ;  
                    Tabc[a][b][c] += S2r * 3 * deltaR[c] * rm5 * Delta(a, b); 
                }
        // Now begin to calculate the charge-charge, charge-dipole, dipole-dipole interations energy and force.
        // The first part: charge-charge interations energy
        realSpaceEwaldEnergy += ONE_4PI_EPS0 * chargeI * chargeJ * T;

        // The second part: charge-dipole interations energy
        for (int a = 0; a < 3; a++) {
            realSpaceEwaldEnergy += ONE_4PI_EPS0 * Mu_all[ii * 3 + a] * Ta[a] * chargeJ - chargeI * Ta[a] * Mu_all[jj * 3 + a];
        }

        // The third part: dipole-dipole interations energy
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++) {
                realSpaceEwaldEnergy -= ONE_4PI_EPS0 * Mu_all[ii * 3 + a] * Tab[a][b] * Mu_all[jj * 3 + a];
            }
        // This is to calcuate the force at atom i and j.
        double eforce[3];
        if (includeForces) {
            for (int a = 0; a < 3; a++) {
                eforce[a] = ONE_4PI_EPS0 * chargeI * Ta[a] * chargeJ;
                polarforce[ii][a] -= eforce[a];
                polarforce[jj][a] += eforce[a];
                for (int b = 0; b < 3; b++) {
                    eforce[a] = ONE_4PI_EPS0 * chargeI * Tab[a][b] * Mu_all[jj * 3 + b] - chargeJ * Tab[a][b] * Mu_all[ii * 3 + b];
                    polarforce[ii][a] += eforce[a];
                    polarforce[jj][a] -= eforce[a];
                    for (int c = 0; c < 3; c++) {
                        eforce[a] = ONE_4PI_EPS0 * Mu_all[ii * 3 + b] * Tabc[a][b][c] * Mu_all[jj * 3 + c];
                        polarforce[ii][a] += eforce[a];
                        polarforce[jj][a] -= eforce[a];
                    }
                }
            }
        }
    }                 
    polarenergy += realSpaceEwaldEnergy ;
    // This is to use the PME methods to calculate the Rec part energy, 
    // refernece MBX code  helpme::PMEInstance<double> pme_solver_;  
    // Compute the reciprocal space terms, using PME
    // On debug
    // Reciprocal space electric potential and field, in sys order at each site
    // PFF: 这一部分的计算还没有跟fftw合着。需要修改。
    double pme_grid_density = 1.2;
    double pme_spline_order = 6;
    helpme::PMEInstance<double> pme_solver;
    double A = periodicBoxVectors[0][0], B = periodicBoxVectors[1][1], C = periodicBoxVectors[2][2];
    int grid_A = pme_grid_density * A;
    int grid_B = pme_grid_density * B;
    int grid_C = pme_grid_density * C;
    pme_solver.setup(1, alphaEwald, pme_spline_order, grid_A, grid_B, grid_C, 1, 0);
    pme_solver.setLatticeVectors(A, B, C, 90, 90, 90, PMEInstanceD::LatticeType::XAligned);
    // N.B. these do not make copies; they just wrap the memory with some metadata
    // System xyz, not ordered XYZ. xyzxyz...(mon1)xyzxyz...(mon2) ...
    std::vector<double> xyz;
    xyz.resize(3*DOFn);
    for (int i = 0; i < DOFn; i++) 
        for (int a = 0; a < 3; a++) {
            xyz[i * 3 + a] = R[i][a];
        }
    // Charges of each sites Order has to follow mon_type_count.
    std::vector<double> chg;
    chg.resize(DOFn);
    for (int i = 0; i < DOFn; i++) 
        chg[i] = ffElec[i].charge;
    // Induced dipole of each sites Order has to follow mon_type_count.    
    std::vector<double> mu;
    mu.resize(3*DOFn);
    for (int i = 0; i < DOFn; i++) 
        for (int a = 0; a < 3; a++) {
            mu[i * 3 + a] = Mu_all[i * 3 + a];
        }
    // Reciprocal space electric potential and field, in sys order at each site
    std::vector<double> rec_phiAndfield1;
    std::vector<double> rec_phiAndfield2;
    rec_phiAndfield1 = std::vector<double>(DOFn * 4, 0.0);
    rec_phiAndfield2 = std::vector<double>(DOFn * 4, 0.0);
    auto coords = helpme::Matrix<double>(xyz.data(), DOFn, 3);
    auto charges = helpme::Matrix<double>(chg.data(), DOFn, 1);
    auto dipoles = helpme::Matrix<double>(mu.data(), DOFn, 3);
    auto result1 = helpme::Matrix<double>(rec_phiAndfield1.data(), DOFn, 4);
    auto result2 = helpme::Matrix<double>(rec_phiAndfield2.data(), DOFn, 4);
    std::fill(rec_phiAndfield1.begin(), rec_phiAndfield1.end(), 0);
    std::fill(rec_phiAndfield2.begin(), rec_phiAndfield2.end(), 0);
    // This is to calculate the charge part Rec energy and force
    pme_solver.computePRec(0, charges, coords, coords, 1, result1);
    // This is to calculate the dipole part Rec energy and force
    pme_solver.computePRec(-1, dipoles, coords, coords, 2, result2);
    for (int i = 0; i < DOFn; i++) {
        const double *result_ptr1 = result1[i];
        const double *result_ptr2 = result2[i];
        polarenergy += recipCoeff * chg[i] * (result_ptr1[0] + result_ptr2[0]);
        polarforce[i][0] += recipCoeff * chg[i] * (result_ptr1[1] + result_ptr2[1]);
        polarforce[i][1] += recipCoeff * chg[i] * (result_ptr1[2] + result_ptr2[2]);
        polarforce[i][2] += recipCoeff * chg[i] * (result_ptr1[3] + result_ptr2[3]);
    }
    if (includeForces)
        for (int i = 0; i < DOFn; i++) {
            forces[i][0] += polarforce[i][0];
            forces[i][1] += polarforce[i][1];
            forces[i][2] += polarforce[i][2];
        }
    totalenergy += polarenergy;
}

void ForceFieldPolar::calculateTTM3PolarForceandEnergy(std::vector<Vec3>& forces, 
                                                std::vector<Vec3>& polarforce, 
                                                double& totalenergy,
                                                double& polarenergy, bool includeForces) {   
    // The electrostatic energy can then be written as a sum of four terms, the real, the reciprocal, 
    // the adjusted and the self terms, and the reciprocal part energy and force use the PME methods 
    // to calculate. The system has charges and dipoles.
    polarforce.resize(DOFn, Vec3());
    std::fill(polarforce.begin(), polarforce.end(), Vec3());
    polarenergy = 0.0;
    double totalSelfEwaldEnergy = 0.0;
    double realSpaceEwaldEnergy = 0.0;
    double recipEnergy = 0.0;
    double totalRecipEnergy = 0.0;
    // set the Ewald method some important parameters
    // For Ewald Sum method parameters: alphaEwald                                                                                     
    double Ewald_tolerance = param.getDouble("Ewald_tolerance");
    double alphaEwald;
    double cutlength = param.getDouble("cutoff");;
    alphaEwald = sqrt(-log(2*Ewald_tolerance))/cutlength; 
    double recipCoeff = ONE_4PI_EPS0 * 4 * pi/(periodicBoxVectors[0][0] * periodicBoxVectors[1][1] * periodicBoxVectors[2][2]);
    /**
     * The self part energy.
     */    
    for (int ii = 0; ii < DOFn; ii++){  
        double charge = ffElec[ii].charge;
        double mu2 = Mu_all[ii * 3 + 0] * Mu_all[ii * 3 + 0] + Mu_all[ii * 3 + 1] * Mu_all[ii * 3 + 1] 
                   + Mu_all[ii * 3 + 2] * Mu_all[ii * 3 + 2];
        totalSelfEwaldEnergy -= alphaEwald * ONE_4PI_EPS0/sqrpi * (charge * charge + 2 * alphaEwald * alphaEwald * mu2/3);
    } 
    polarenergy += totalSelfEwaldEnergy;
    /**
     * Real SPACE EWALD ENERGY AND FORCES
     * T = s[0] * 1/r;
     * T[a] = - s[1] * rij[a]/r^3;
     * T[a][b] = s[2] * 3*rij[a]*rij[b]/r^5 -s[1] * delta[ab]/r^3;
     * T[a][b][c] = -s[3]*15*rij[a]*rij[b]*rij[c]/r^7 + s[2]*3*(rij[a]*delta[bc]+rij[b]*delta[ac]+rij[c]*delta[ab])/r^5
     * Using the TTM3 damping model method.
     */
    double rm1, rm2, rm3, rm5, rm7;
    // ***** There is an important thing we should focus on the atom i and j, 
    //       the pairs should included in 12, 13, 14, and nonbonded
    // ***** The unit problems should also inportant for the energy and force
    //       energy: KJ/mol 
    //       force: KJ/mol/nm
    // This is set the pair 12,13, and nonbonded part.
    std::vector<std::set<int> > bonded12(DOFn); 
    std::vector<std::set<int> > bonded13(DOFn);
    addTholeExclusionsToSet(bonded12, bonded13);
    // Set the atom polarizability
    double alpha[DOFn];
    for (int i = 0; i < DOFn; i++) {
       alpha[i] = ffPolars[i].alpha/1000;
    }
    for (int ii = 0; ii < DOFn; ii++) 
        for (int jj = 0; jj < DOFn; jj++) {
            double chargeI = ffElec[ii].charge;
            double chargeJ = ffElec[jj].charge;
            double deltaR[5];
            const Vec3& atomCoordinatesI = R[ii];
            const Vec3& atomCoordinatesJ = R[jj];
            // Get deltaR 2 atoms, deltaR = Ri - Rj
            GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
            
            rm1 = 1.0 / deltaR[4];
            rm2 = rm1 * rm1;
            rm3 = rm1 * rm2;
            rm5 = rm2 * rm3;
            rm7 = rm5 * rm2;
            double alphaij, Ralphaij2, Ralphaij4;
            alphaij = alpha[ii] * alpha[jj];
            alphaij = pow((double)alphaij,double(1.0/6));            
            alphaij = deltaR[4] * 1.0 / alphaij;
            double S0r, S1r, S2r, S3r;
            // Thole model daming parameter for chosen pair od atom i and j
            double dampingParam;
            dampingParam = param.getDouble("dampingParam3");
            if (bonded12[ii].find(jj) == bonded12[ii].end()||bonded12[jj].find(ii) == bonded12[jj].end()) 
                dampingParam = param.getDouble("dampingParam1");
            if (bonded13[ii].find(jj) == bonded13[ii].end()||bonded13[jj].find(ii) == bonded13[jj].end()) 
                dampingParam = param.getDouble("dampingParam2");
            // This is the TTM-3 model Exponential damping compare with the Amber and Tinker software
            TTM3Model(dampingParam, alphaEwald, S0r, S1r, S2r, S3r);
            // This is to struct the T, Ta[a], Tab[a][b], Tabc[a][b][c]
            double T, Ta[3], Tab[3][3], Tabc[3][3][3];
            T = S0r * rm1;
            for (int a = 0; a < 3; a++) {
                Ta[a] = - S1r * deltaR[a] * rm3;
            }
            for (int a = 0; a < 3; a++) 
                for (int b = 0; b < 3; b++) {                   
                    Tab[a][b] = 3 * deltaR[a] * deltaR[b] * rm5 * S2r;
                    if (a == b) {
                        Tab[a][b] -= S1r * rm3; 
                    }
               } 
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++)
                    for (int c = 0; c < 3; c++) {
                        Tabc[a][b][c] = - 15 * deltaR[a] * deltaR[b] * deltaR[c] * S3r * rm7;
                        if (b == c) {
                            Tabc[a][b][c] += S2r * 3 * deltaR[a] * rm5; 
                        } 
                        if (a == c) {
                            Tabc[a][b][c] += S2r * 3 * deltaR[b] * rm5; 
                        } 
                        if (a == b) {
                            Tabc[a][b][c] += S2r * 3 * deltaR[c] * rm5; 
                        }  
                    }
        // Now begin to calculate the charge-charge, charge-dipole, dipole-dipole interations energy and force.
        // The first part: charge-charge interations energy
        realSpaceEwaldEnergy += ONE_4PI_EPS0 * chargeI * chargeJ * T;
        // The second part: charge-dipole interations energy
        for (int a = 0; a < 3; a++) {
            realSpaceEwaldEnergy += ONE_4PI_EPS0 * Mu_all[ii * 3 + a] * Ta[a] * chargeJ - chargeI * Ta[a] * Mu_all[jj * 3 + a];
        }
        // The third part: dipole-dipole interations energy
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++) {
                realSpaceEwaldEnergy -= ONE_4PI_EPS0 * Mu_all[ii * 3 + a] * Tab[a][b] * Mu_all[jj * 3 + a];
            }
        // This is to calcuate the force at atom i and j.
        double eforce[3];
        if (includeForces) {
            for (int a = 0; a < 3; a++) {
                eforce[a] = 0.5 * ONE_4PI_EPS0 * chargeI * Ta[a] * chargeJ;
                polarforce[ii][a] -= eforce[a];
                polarforce[jj][a] += eforce[a];
            }
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++) {
                    eforce[a] = 0.5 * ONE_4PI_EPS0 * chargeI * Tab[a][b] * Mu_all[jj * 3 + b] - chargeJ * Tab[a][b] * Mu_all[ii * 3 + b];
                    polarforce[ii][a] += eforce[a];
                    polarforce[jj][a] -= eforce[a];
                }
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++) 
                    for (int c = 0; c < 3; c++) {
                        eforce[a] = 0.5 * ONE_4PI_EPS0 * Mu_all[ii * 3 + b] * Tabc[a][b][c] * Mu_all[jj * 3 + c];
                        polarforce[ii][a] += eforce[a];
                        polarforce[jj][a] -= eforce[a];
                    }
        }                 
    }
    polarenergy += 0.5 * realSpaceEwaldEnergy * ONE_4PI_EPS0;
    // This is to use the PME methods to calculate the Rec part energy, 
    // refernece MBX code  helpme::PMEInstance<double> pme_solver_;  
    // Compute the reciprocal space terms, using PME
    // Reciprocal space electric potential and field, in sys order at each site
    double pme_grid_density = 1.2;
    double pme_spline_order = 6;
    helpme::PMEInstance<double> pme_solver;
    double A = periodicBoxVectors[0][0], B = periodicBoxVectors[1][1], C = periodicBoxVectors[2][2];
    int grid_A = pme_grid_density * A;
    int grid_B = pme_grid_density * B;
    int grid_C = pme_grid_density * C;
    pme_solver.setup(1, alphaEwald, pme_spline_order, grid_A, grid_B, grid_C, 1, 0);
    pme_solver.setLatticeVectors(A, B, C, 90, 90, 90, PMEInstanceD::LatticeType::XAligned);
    // N.B. these do not make copies; they just wrap the memory with some metadata
    // System xyz, not ordered XYZ. xyzxyz...(mon1)xyzxyz...(mon2) ...
    std::vector<double> xyz;
    xyz.resize(3*DOFn);
    for (int i = 0; i < DOFn; i++) 
        for (int a = 0; a < 3; a++) {
            xyz[i * 3 + a] = R[i][a];
        }
    // Charges of each sites Order has to follow mon_type_count.
    std::vector<double> chg;
    chg.resize(DOFn);
    for (int i = 0; i < DOFn; i++) 
        chg[i] = ffElec[i].charge;
    // Induced dipole of each sites Order has to follow mon_type_count.    
    std::vector<double> mu;
    mu.resize(3*DOFn);
    for (int i = 0; i < DOFn; i++) 
        for (int a = 0; a < 3; a++) {
            mu[i * 3 + a] = Mu_all[i * 3 + a];
        }
    std::cout.precision(12);
    // Reciprocal space electric potential and field, in sys order at each site
    std::vector<double> rec_phiAndfield1;
    std::vector<double> rec_phiAndfield2;
    rec_phiAndfield1 = std::vector<double>(DOFn * 4, 0.0);
    rec_phiAndfield2 = std::vector<double>(DOFn * 4, 0.0);
    auto coords = helpme::Matrix<double>(xyz.data(), DOFn, 3);
    auto charges = helpme::Matrix<double>(chg.data(), DOFn, 1);
    auto dipoles = helpme::Matrix<double>(mu.data(), DOFn, 3);
    
    auto result1 = helpme::Matrix<double>(rec_phiAndfield1.data(), DOFn, 4);
    auto result2 = helpme::Matrix<double>(rec_phiAndfield2.data(), DOFn, 4);
    std::fill(rec_phiAndfield1.begin(), rec_phiAndfield1.end(), 0);
    std::fill(rec_phiAndfield2.begin(), rec_phiAndfield2.end(), 0);
    // This is to calculate the charge part Rec energy and force
    pme_solver.computePRec(0, charges, coords, coords, 1, result1);
    // This is to calculate the dipole part Rec energy and force
    //pme_solver.computePRec(-1, dipoles, coords, coords, 2, result2);
    for (int i = 0; i < DOFn; i++) {
        const double *result_ptr1 = result1[i];
        //const double *result_ptr2 = result2[i];
        //polarenergy += recipCoeff * chg[i] * (result_ptr1[0] + result_ptr2[0]);
        //polarforce[i][0] += recipCoeff * chg[i] * (result_ptr1[1] + result_ptr2[1]);
        //polarforce[i][1] += recipCoeff * chg[i] * (result_ptr1[2] + result_ptr2[2]);
        //polarforce[i][2] += recipCoeff * chg[i] * (result_ptr1[3] + result_ptr2[3]);
        polarenergy += recipCoeff * chg[i] * result_ptr1[0];
        polarforce[i][0] += recipCoeff * chg[i] * result_ptr1[1];
        polarforce[i][1] += recipCoeff * chg[i] * result_ptr1[2];
        polarforce[i][2] += recipCoeff * chg[i] * result_ptr1[3];
    }
    if (includeForces)
        for (int i = 0; i < DOFn; i++) {
            forces[i][0] += polarforce[i][0];
            forces[i][1] += polarforce[i][1];
            forces[i][2] += polarforce[i][2];
        }
    totalenergy += polarenergy;
}

// TO DO
//void ForceFieldPolar::calculateESPolarForceAndEnergy(std::vector<Vec3>& forces, 
//                                                    std::vector<Vec3>& polarforce, 
//                                                    double& totalenergy,
//                                                double& polarenergy, bool includeForces) {
//}

void ForceFieldPolar::calculateInducedDipole() {
    // Calculate the many-body Induced dipole as a function of configuration using iteration method
    // to iter_max order or already converged. Return actual iteration times.
    double rij[3] = {0, 0, 0};
    double Mu_prev[3] = {0, 0, 0};
    double Tij[3][3] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    // for swap matrix
    double m1[3] = {0,0,0};  
    double m2[3] = {0,0,0};  
    std::vector<Vec3> mu1;
    std::vector<Vec3> mu2; 
    mu1.resize(DOFn, Vec3());
    mu2.resize(DOFn, Vec3());
    int a, b, i, j, k, l;
    // tolerance
    double tol(1); 
    double r, r2, rm1, rm2, rm3, rm5;
    double alphaij, Ralphaij, Ralphaij2, Ralphaij4;
    double S0r, S1r, S2r, S3r;
    double alpha[DOFn];
    // Initialize the permanent electric field and calculate the E_perm through using E = F / q
    std::vector<Vec3> E_perm;
    E_perm.resize(DOFn, Vec3());
    getElecPerm(E_perm);
    std::vector<std::set<int> > bonded12(DOFn); 
    std::vector<std::set<int> > bonded13(DOFn);
    if (damping_type == "TTM3"||damping_type == "TTM4") {
        addTholeExclusionsToSet(bonded12, bonded13);
    }
    
    // Iteration count
    int n(0);  
    // Set the atom polarizability
    for (i = 0; i < DOFn; i++) {
        alpha[i] = ffPolars[i].alpha /1000;
    }
    for (a = 0; a < 3; a++) Mu[a] = 0; 
    // Initial guess as zero-order μ = α * E
    for (i = 0; i < DOFn; i++) {
        mu2[i][0] = alpha[i] * E_perm[i][0];
        mu2[i][1] = alpha[i] * E_perm[i][1];
        mu2[i][2] = alpha[i] * E_perm[i][2];
    }    

    while ( n< iter_max ) {
    // while ( tol > TOLERANCE && n< iter_max ) {   
        n++;
        // Save Mu_prev<=Mu
        for (a = 0; a < 3; a++) {
           Mu_prev[a] = Mu[a];         
        }
        // Save (old) mu1 <= (new) mu2, and set mu2 zero
        for (i = 0; i < DOFn; i++) 
            for (k = 0; k < 3; k++) {
                mu1[i][k] = mu2[i][k];
                mu2[i][k] = 0;
            }
        
        // Start new iteration for mu2[i=1...DOFn][3]
        for (i = 0; i < DOFn; i++) {
            for (j = 0; j < DOFn; j++) {
                if (i == j) continue;
                double deltaR[5];
                const Vec3& atomCoordinatesI = R[i];
                const Vec3& atomCoordinatesJ = R[j];
                // Get deltaR 2 atoms, deltaR = Ri - Rj
                GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
                rm1 = 1.0 / deltaR[4];
                rm2 = rm1 * rm1;
                rm3 = rm1 * rm2;
                rm5 = rm2 * rm3;
                
                // choose whitch thole model damping type.
                if (damping_type == "TTM0") 
                    TTM0Model(S0r, S1r, S2r, S3r);
                else if (damping_type == "TTM1") 
                    TTM1Model(S0r, S1r, S2r, S3r);
                //else if (damping_type == "TTM2") 
                //    TTM2Model(alphaEwald, deltaR, S0r, S1r, S2r, S3r);
                else {
                    alphaij = alpha[i] * alpha[j];
                    alphaij = pow((double)alphaij,double(1.0/6));            
                    alphaij = deltaR[4] * 1.0 / alphaij;
                    Ralphaij2 = alphaij * alphaij;
                    Ralphaij4 = Ralphaij2 * Ralphaij2;
                    // Thole model daming parameter for chosen pair od atom i and j
                    double dampingParam;
                    // Judge the pair i, j belong to the 12, 13, 14, or noboned part.
                    dampingParam = param.getDouble("dampingParam3");
                    if (bonded12[i].find(j) == bonded12[i].end()||bonded12[j].find(i) == bonded12[j].end()) 
                        dampingParam = param.getDouble("dampingParam1");
                    if (bonded13[i].find(j) == bonded13[i].end()||bonded13[j].find(i) == bonded13[j].end()) 
                        dampingParam = param.getDouble("dampingParam2");

                    if (damping_type == "TTM3")
                        TTM3Model(dampingParam, alphaij, S0r, S1r, S2r, S3r);
                    else if (damping_type == "TTM4")
                        TTM4Model(dampingParam, alphaij, S0r, S1r, S2r, S3r);
                }
                for (int a = 0; a < 3; a++) 
                    for (int b = 0; b < 3; b++) {                   
                        Tij[a][b] = 3 * deltaR[a] * deltaR[b] * rm5 * S2r;
                        if (a == b) {
                            Tij[a][b] -= S1r * Delta(a,b) * rm3; 
                        }
                   }            
                // Transfrom mu1[i][3] to a matrix m1[3][3]
                k = 0;
                for (a = 0; a < 3; a++) {
                    m1[a] = mu1[j][k];
                    k++;
                }          
                
                MatrixMultiplyMu(Tij, m1, m2);              
               // Add matrix m2[3][3] to pi2[i][9]
                for (a = 0; a < 3; a++) {
                    mu2[i][a] += m2[a]; 
                }
            } 
            mu2[i][0] += E_perm[i][0];
            mu2[i][1] += E_perm[i][1];
            mu2[i][2] += E_perm[i][2];
            for (k = 0; k < 3; k++) {
                mu2[i][k] *= alpha[i];
            }   
        } 
        //calculate many-body induced dipole Mu[3] (sum up i=1...N)
        //initialize result tensor to be 0
        for (a = 0; a < 3; a++) Mu[a] = 0; 
        for (l = 0; l < DOFn; l++) {
            k = 0;
            for (a = 0; a < 3; a++) {
                Mu[a] += mu2[l][k];
                k++;
            }
        }
        //calculate tolerance
        tol = 0;
        for (a = 0; a < 3; a++) {
            tol += abs(Mu[a] - Mu_prev[a]);
        }
    }
    //save mu_i to Mu_all
    for (i = 0; i < DOFn; i++)
        for (j = 0; j < 3; j++) {
            Mu_all[i * 3 + j] = mu2[i][j];
        }                   
}

double ForceFieldPolar::getPotentialEnergy(bool includeForces) {
    double PE = 0;
    std::fill(F.begin(), F.end(), Vec3());
    if (force_contribution[1]) 
        calculateBondForceandEnergy(F, forceParts[1], PE, energyParts[1], includeForces);
    if (force_contribution[2]) 
        calculateAngleForceandEnergy(F, forceParts[2], PE, energyParts[2], includeForces);
    if (force_contribution[3]) 
        calculateDihedralForceandEnergy(F, forceParts[3], PE, energyParts[3], includeForces);
    if (force_contribution[4]) 
        calculateLJForceandEnergy(F, forceParts[4], PE, energyParts[4], includeForces);  
    if (forcefield_type == "Polar") {
        // This is to calculate the induced dipole result.
        calculateInducedDipole();
        // use the PME method to calcuate the polar energy and force.
        if (force_contribution[6])
            if (damping_type == "TTM3") {
                calculateTTM3PolarForceandEnergy(F, forceParts[6], PE, energyParts[6], includeForces);
            }
            else if (damping_type == "TTM4") { 
                calculateTTM4PolarForceandEnergy(F, forceParts[6], PE, energyParts[6], includeForces);
            }
            else 
                throw std::runtime_error("ERROR: damping_type " + damping_type); 
    }
    else if (forcefield_type == "nonPolar") {
        if (force_contribution[5]) 
            calculateElecForceandEnergy(F, forceParts[5], PE, energyParts[5], includeForces);   
    }
    else 
        throw std::runtime_error("ERROR: forcefield_type " + forcefield_type); 
    if (md_perturb) {
        std::vector<int> perturb_parameters;
        SplitString(perturb_parameters, param.getStr("perturb_parameters"));
        if (perturb_parameters[1] == 0) {
            std::cout<<"begin to calculate the perturb polar force"<<std::endl;
            calcualtePerturbPolarForce(F); 
        } 
        else 
            calcualtePerturbMuForce(F);      
    } 
    return PE;
}

void ForceFieldPolar::getPiTensor(std::vector<double>& Pi_tensor) {
    for (int a = 0; a < 3; a++) 
        for (int b = 0; b < 3; b++) {
            Pi_tensor[a * 3 + b] = Pi[a][b];
        }
}

void ForceFieldPolar::getDipoleRleftAright(double& DEScosleft, double& DEScosright, int& DESnumLeft, 
                                          int& DESnumRight, double& DESdirLx, double& DESdirRx) {
    DEScosleft = cosleft;
    DEScosright = cosright;
    DESnumLeft = numLeft;
    DESnumRight = numRight;
    DESdirLx = dirLx;
    DESdirRx = dirRx;
}

void ForceFieldPolar::getMuTensor(Vec3& Mu_tensor) {
    for (int a = 0; a < 3; a++) {
        Mu_tensor[a] = Mu[a]; 
    }
}

void ForceFieldPolar::TTM0Model(double& S0r, double& S1r, double& S2r, double& S3r) {
    S0r = 1;
    S1r = 1;
    S2r = 1;
    S3r = 1;
}

void ForceFieldPolar::TTM1Model(double& S0r, double& S1r, double& S2r, double& S3r) {
    S0r = 1;
    S1r = 1;
    S2r = 1;
    S3r = 1;
}

void ForceFieldPolar::TTM2Model(double alphaEwald, double deltaR, double& S0r, double& S1r, double& S2r, double& S3r) {
    // This is the Gaussian ES Damping model
    double Ralpha, Ralpha2, Ralpha3, Ralpha5;
    Ralpha = alphaEwald * deltaR;
    Ralpha2 = Ralpha * Ralpha;
    Ralpha3 = Ralpha2 * Ralpha;
    Ralpha5 = Ralpha3 * Ralpha2;
    
    S0r = erfc(Ralpha); 
    S1r = S0r + 2 * Ralpha * exp(-Ralpha2)/sqrpi;
    S2r = S1r + 4 * Ralpha3 * exp(-Ralpha2)/(3 * sqrpi);
    S3r = S2r + 8 * Ralpha5 * exp(-Ralpha2)/(15 * sqrpi);
}

void ForceFieldPolar::TTM3Model(double dampingParam,double alphaij, double& S0r, double& S1r, double& S2r, double& S3r) {
    // This is the Exponential damping compare with the Amber and OpenMM Amoeba
    double Ralphaij2, Ralphaij3, Ralphaij6, Ralphaij9;
    Ralphaij2 = alphaij * alphaij;
    Ralphaij3 = Ralphaij2 * alphaij;
    Ralphaij6 = Ralphaij3 * Ralphaij3;
    Ralphaij9 = Ralphaij3 * Ralphaij6;
    // Thole model daming parameter for chosen pair od atom i and j, means a.
    //Gamma gamma(2.0/3, dampingParam * Ralphaij3);
    //gamma.SetPrecision(0.001).SetMaxNum(dampingParam * Ralphaij3+10000);// need check
    //double gamma_value = gamma.Exec();
    //S0r = 1 - exp( -dampingParam * Ralphaij3) + alphaij * pow((double)dampingParam,double(1.0/4)) * gamma_value;
    S0r = 1;
    S1r = 1 - exp(- dampingParam * Ralphaij3); //need check 
    S1r = 1 - (1 + dampingParam * Ralphaij3) * exp(-dampingParam * Ralphaij3);
    S3r = 1 - (1 + dampingParam * Ralphaij3 + 3.0 * dampingParam * dampingParam * Ralphaij6/5) * exp(-dampingParam * Ralphaij3);
}

void ForceFieldPolar::TTM4Model(double dampingParam, double alphaij, double& S0r, double& S1r, double& S2r, double& S3r) {
    // This is the TTM-4 model Exponential damping compare with the MBX software
    double Ralphaij2, Ralphaij4;
    Ralphaij2 = alphaij * alphaij;
    Ralphaij4 = Ralphaij2 * Ralphaij2;
    //Gamma gamma(3.0/4, dampingParam * Ralphaij4);
    //gamma.SetPrecision(0.001).SetMaxNum(dampingParam * Ralphaij4+10000);
    //double gamma_value = gamma.Exec();
    double gamma_value = 1;
    S0r = 1 - exp(- dampingParam * Ralphaij4) + pow((double)dampingParam,double(1.0/4)) * alphaij * gamma_value;
    S1r = 1 - exp(- dampingParam * Ralphaij4);
    S2r = S1r - (4.0 * dampingParam * Ralphaij4 / 3.0) * exp(- dampingParam * Ralphaij4);
    S3r = S2r - (4.0 * dampingParam * Ralphaij4 / 15.0) * exp(- dampingParam * Ralphaij4) * (4 * dampingParam * Ralphaij4 - 1);
}

void ForceFieldPolar::addTholeExclusionsToSet(std::vector<std::set<int> >& bonded12, std::vector<std::set<int> >& bonded13) {
    
    int nof12 = exclusions12.size();
    for (int i = 0; i < nof12; ++i) {   
        int atomAIndex = exclusions12[i].atomIndex[0];
        int atomBIndex = exclusions12[i].atomIndex[1];
        bonded12[atomAIndex].insert(atomBIndex);
        bonded12[atomBIndex].insert(atomAIndex);
    }  
    int nof13 = exclusions13.size();
    for (int i = 0; i < nof13; ++i) {   
        int atomAIndex = exclusions13[i].atomIndex[0];
        int atomBIndex = exclusions13[i].atomIndex[1];
        bonded13[atomAIndex].insert(atomBIndex);
        bonded13[atomBIndex].insert(atomAIndex);
    }  
}

void ForceFieldPolar::setTrueMdPerturb() {
    md_perturb = true;
}

void ForceFieldPolar::setFalseMdPerturb() {
    md_perturb = false;
}

void ForceFieldPolar::setPositiveDirection() {
    direction = true;
}

void ForceFieldPolar::setNegativeDirection() {
    direction = false;
}

void ForceFieldPolar::calcualtePerturbPolarForce(std::vector<Vec3>& forces) {
    std::vector<Vec3> perturbforce;
    perturbforce.resize(DOFn, Vec3());
    int conf_diff;
    conf_diff = 1;
    std::vector<double> iPi, iPi_all, Tij_all, DkPi, DPi;
    iPi.resize(9);
    iPi_all.resize(DOFn * 9);
    Tij_all.resize(DOFn * DOFn * 9);
    DkPi.resize(27);
    DPi.resize(DOFn * 27);
    
    Pi_tensor(iPi_all, iPi);
   
    conf_diff = 1;
   
    for (int k = 0; k < DOFn; k++) {
        Dk_Pi_tensor(k, iPi_all, Tij_all, conf_diff, DkPi); //Pi_all or Pi_alpha_all
        
        if (conf_diff == 1) conf_diff = 0;

        for (int u = 0; u < 27; u++) DPi[k * 27 + u] = DkPi[u]; //update Deriv of Pi for all atoms
    }
    
    Perturb_Polar(DPi, perturbforce);
    
    for (int i = 0; i < DOFn; i++) {
        if (direction) {
            forces[i][0] += perturbforce[i][0];
            forces[i][1] += perturbforce[i][1];
            forces[i][2] += perturbforce[i][2];
            }
        else {
            forces[i][0] -= perturbforce[i][0];
            forces[i][1] -= perturbforce[i][1];
            forces[i][2] -= perturbforce[i][2];
        }
    }
}

void ForceFieldPolar::calcualtePerturbMuForce(std::vector<Vec3>& forces) {
    std::vector<Vec3> perturbforce;
    perturbforce.resize(DOFn, Vec3());
    int conf_diff;
    conf_diff = 1;
    std::vector<double> iMu_all, iMu, Tij_all, DkMu, DMu, DkE_perm, DkE;
    iMu_all.resize(DOFn * 3);
    iMu.resize(3);
    Tij_all.resize(DOFn * DOFn * 9);
    DkMu.resize(9);
    DMu.resize(DOFn * 9);
    DkE.resize(9);
    DkE_perm.resize(DOFn * 9);

    Mu_tensor(iMu_all, iMu);

    conf_diff = 1;
    for (int k = 0; k < DOFn; k++) {
        Dk_Mu_tensor(k, iMu_all, Tij_all, conf_diff, DkMu); //Pi_all or Pi_alpha_all
        if (conf_diff == 1) conf_diff = 0;
        for (int u = 0; u < 9; u++) DMu[k * 9 + u] = DkMu[u]; //update Deriv of Pi for all atoms
    }

    // calculate Deriv of the Elec_perm 
    for (int k = 0; k < DOFn; k++) {  
        Dk_E_perm(k, DkE_perm, DkE);
        for (int u = 0; u < 9; u++) DMu[k * 9 + u] += DkE[u]; //update Deriv of E_perm for all atoms
    }
    
}

void ForceFieldPolar::Pi_tensor(std::vector<double>& iPi_all, std::vector<double>& iPi) {
    // calculate the many-body polarizability as a function of configuration using iteration method
    // to iter_max order or already converged. Return actual iteration times.
    double Pi_prev[3][3]={0,0,0,0,0,0,0,0,0};
    double Tij[3][3]={0,0,0,0,0,0,0,0,0};
    double m1[3][3]={0,0,0,0,0,0,0,0,0};  //for swap matrix
    double m2[3][3]={0,0,0,0,0,0,0,0,0};  //for swap matrix
    std::vector<double> pi1, pi2;
    pi1.resize(DOFn * 9); //old
    pi2.resize(DOFn * 9); //new
    double tol(1); //tolerance
    
    double alpha[DOFn]; //assign single molecular polarizability
    for (int i = 0; i < DOFn; i++) {
        alpha[i] = ffPolars[i].alpha/1000;
    }
    
    for (int i = 0; i < DOFn; i++) { // initial guess as zero-order alpha
        pi2[i * 9 + 0] = pi2[i * 9 + 4] = pi2[i * 9 + 8] = alpha[i];
        pi2[i * 9 + 1] = pi2[i * 9 + 2] = pi2[i * 9 + 3] = pi2[i * 9 + 5] = pi2[i * 9 + 6] = pi2[i * 9 + 7] = 0;
    }
    
    std::vector<std::set<int> > bonded12(DOFn); 
    std::vector<std::set<int> > bonded13(DOFn);
    if (damping_type == "TTM3"||damping_type == "TTM4") {
        addTholeExclusionsToSet(bonded12, bonded13);
    }

    int n(0);  //iteration count

    while ( tol > TOLERANCE && n < iter_max ) {
        n++;
        //save Pi_prev<=Pi
        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                Pi_prev[a][b] = iPi[3 * a + b];
            }
        }
        
        //save (old) pi1 <= (new) pi2, and set pi2 zero
        for (int i = 0; i < DOFn; i++) {
            for (int k = 0; k < 9; k++) {
                pi1[i * 9 + k] = pi2[i * 9 + k];
                pi2[i * 9 + k] = 0;
            }
        }
        
        //start new iteration for pi2[i=1...N][9]
        for (int i = 0; i < DOFn; i++) {
            for (int j = 0; j < DOFn; j++) {
                if (i == j) continue;
                //calculate Tij
                double deltaR[5];
                double rm1, rm2, rm3, rm5;
                const Vec3& atomCoordinatesI = R[i];
                const Vec3& atomCoordinatesJ = R[j];
                // get deltaR 2 atoms, deltaR = Ri - Rj
                GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
                rm1 = 1.0 / deltaR[4];
                rm2 = rm1 * rm1;
                rm3 = rm1 * rm2;
                rm5 = rm2 * rm3;
                double S0r, S1r, S2r, S3r;
                // choose whitch thole model damping type.
                if (damping_type == "TTM0") 
                    TTM0Model(S0r, S1r, S2r, S3r);
                else if (damping_type == "TTM1") 
                    TTM1Model(S0r, S1r, S2r, S3r);
                //else if (damping_type == "TTM2") 
                //    double alphaEwald, double deltaR, double& S0r, double& S1r, double& S2r, double& S3r) {
                else {
                    double alphaij;
                    alphaij = alpha[i] * alpha[j];
                    alphaij = pow((double)alphaij,double(1.0/6));            
                    alphaij = deltaR[4] * 1.0 / alphaij;
                    // Thole model daming parameter for chosen pair od atom i and j
                    double dampingParam;
                    // Judge the pair i, j belong to the 12, 13, 14, or noboned part.
                    dampingParam = param.getDouble("dampingParam3");
                    if (bonded12[i].find(j) == bonded12[i].end()||bonded12[j].find(i) == bonded12[j].end()) 
                        dampingParam = param.getDouble("dampingParam1");;
                    if (bonded13[i].find(j) == bonded13[i].end()||bonded13[j].find(i) == bonded13[j].end()) 
                        dampingParam = param.getDouble("dampingParam2");;

                    if (damping_type == "TTM3")
                        TTM3Model(dampingParam, alphaij, S0r, S1r, S2r, S3r);
                    else if (damping_type == "TTM4")
                        TTM4Model(dampingParam, alphaij, S0r, S1r, S2r, S3r);
                }

                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {                   
                        Tij[a][b] = 3 * deltaR[a] * deltaR[b] * rm5 * S2r;
                        if (a == b) {
                            Tij[a][b] -= rm3 * S1r; 
                        }
                    }
                }
                
                //transfrom pi1[i][9] to a matrix m1[3][3]
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {
                        m1[a][b] = pi1[j * 9 + 3 * a + b];
                    }
                }
                
                MatrixMultiplyPi(Tij, m1, m2);           

                //add matrix m2[3][3] to pi2[i][9]
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {
                        pi2[i * 9 + 3 * a + b] += m2[a][b];
                    }
                }
                
            }
            pi2[i * 9 + 0] += 1;
            pi2[i * 9 + 4] += 1;
            pi2[i * 9 + 8] += 1;
            
            for (int k = 0; k < 9; k++) {
                pi2[i * 9 + k] *= alpha[i];
            }
        }
        
        //calculate many-body polarizability Pi[3][3] (sum up i=1...N)
        for (int a = 0; a < 3; a++) 
            iPi[a * 3 + 0] = iPi[a * 3 + 1] = iPi[a * 3 + 2] = 0; //initialize result tensor to be 0

        for (int i = 0; i < DOFn; i++) {
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++) {
                    iPi[a * 3 + b] += pi2[i * 9 + 3 * a + b];
                }
        }
        //calculate tolerance
        tol=0;

        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                tol += abs(iPi[a * 3 + b] - Pi_prev[a][b]);
            }
        }
    }
    //save pi_i to Pi_all
    for (int i = 0; i < DOFn; i++)
        for (int j = 0; j < 9 ;j++) {
            iPi_all[i * 9 + j] = pi2[i * 9 + j];
        }
}

void ForceFieldPolar::Mu_tensor(std::vector<double>& iMu_all, std::vector<double>& iMu) {
    // calculate the many-body polarizability as a function of configuration using iteration method
    // to iter_max order or already converged. Return actual iteration times.
    double Mu_prev[3]={0,0,0};
    double Tij[3][3]={0,0,0,0,0,0,0,0,0};
    double m1[3]={0,0,0};  //for swap matrix
    double m2[3]={0,0,0};  //for swap matrix
    std::vector<double> mu1, mu2;
    mu1.resize(DOFn * 3); //old
    mu2.resize(DOFn * 3); //new
    double tol(1); //tolerance
    
    double alpha[DOFn]; //assign single molecular polarizability
    for (int i = 0; i < DOFn; i++) {
        alpha[i] = ffPolars[i].alpha/1000;
    }
    
    std::vector<Vec3> E_perm;
    E_perm.resize(DOFn, Vec3());
    getElecPerm(E_perm);

    for (int i = 0; i < DOFn; i++) { // initial guess as zero-order alpha
        mu2[i * 3 + 0] = alpha[i] * E_perm[i][0]; 
        mu2[i * 3 + 1] = alpha[i] * E_perm[i][1]; 
        mu2[i * 3 + 2] = alpha[i] * E_perm[i][2];
    }
    
    std::vector<std::set<int> > bonded12(DOFn); 
    std::vector<std::set<int> > bonded13(DOFn);
    if (damping_type == "TTM3"||damping_type == "TTM4") {
        addTholeExclusionsToSet(bonded12, bonded13);
    }

    int n(0);  //iteration count

    while ( tol > TOLERANCE && n < iter_max ) {
        n++;
        //save Mu_prev<=Mu
        for (int a = 0; a < 3; a++) {
            Mu_prev[a] = iMu[a];
        }
        
        //save (old) mu1 <= (new) mu2, and set mu2 zero
        for (int i = 0; i < DOFn; i++) {
            for (int k = 0; k < 3; k++) {
                mu1[i * 3 + k] = mu2[i * 3 + k];
                mu2[i * 3 + k] = 0;
            }
        }
        
        //start new iteration for mu2[]
        for (int i = 0; i < DOFn; i++) {
            for (int j = 0; j < DOFn; j++) {
                if (i == j) continue;
                //calculate Tij
                double deltaR[5];
                double rm1, rm2, rm3, rm5;
                const Vec3& atomCoordinatesI = R[i];
                const Vec3& atomCoordinatesJ = R[j];
                // get deltaR 2 atoms, deltaR = Ri - Rj
                GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
                rm1 = 1.0 / deltaR[4];
                rm2 = rm1 * rm1;
                rm3 = rm1 * rm2;
                rm5 = rm2 * rm3;
                double S0r, S1r, S2r, S3r;
                // choose whitch thole model damping type.
                if (damping_type == "TTM0") 
                    TTM0Model(S0r, S1r, S2r, S3r);
                else if (damping_type == "TTM1") 
                    TTM1Model(S0r, S1r, S2r, S3r);
                //else if (damping_type == "TTM2") 
                //    double alphaEwald, double deltaR, double& S0r, double& S1r, double& S2r, double& S3r) {
                else {
                    double alphaij;
                    alphaij = alpha[i] * alpha[j];
                    alphaij = pow((double)alphaij,double(1.0/6));            
                    alphaij = deltaR[4] * 1.0 / alphaij;
                    // Thole model daming parameter for chosen pair od atom i and j
                    double dampingParam;
                    // Judge the pair i, j belong to the 12, 13, 14, or noboned part.
                    dampingParam = param.getDouble("dampingParam3");
                    if (bonded12[i].find(j) == bonded12[i].end()||bonded12[j].find(i) == bonded12[j].end()) 
                        dampingParam = param.getDouble("dampingParam1");;
                    if (bonded13[i].find(j) == bonded13[i].end()||bonded13[j].find(i) == bonded13[j].end()) 
                        dampingParam = param.getDouble("dampingParam2");;

                    if (damping_type == "TTM3")
                        TTM3Model(dampingParam, alphaij, S0r, S1r, S2r, S3r);
                    else if (damping_type == "TTM4")
                        TTM4Model(dampingParam, alphaij, S0r, S1r, S2r, S3r);
                }

                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {                   
                        Tij[a][b] = 3 * deltaR[a] * deltaR[b] * rm5 * S2r;
                        if (a == b) {
                            Tij[a][b] -= rm3 * S1r; 
                        }
                    }
                }
                
                //transfrom mu1[i][3] to a matrix m1[3]
                for (int a = 0; a < 3; a++) {
                    m1[a] = mu1[j * 3 + a];
                }
                
                MatrixMultiplyMu(Tij, m1, m2);           

                //add matrix m2[3] to mu2[i][3]
                for (int a = 0; a < 3; a++) {
                    mu2[i * 3 + a] += m2[a];
                }
                
            }
            mu2[i * 3 + 0] += E_perm[i][0];
            mu2[i * 3 + 4] += E_perm[i][1];
            mu2[i * 3 + 8] += E_perm[i][2];
            
            for (int k = 0; k < 3; k++) {
                mu2[i * 3 + k] *= alpha[i];
            }
        }
        
        //calculate many-body induced dipole Mu[3] (sum up i=1...N)
        for (int a = 0; a < 3; a++) 
            iMu[a] = 0; //initialize result tensor to be 0

        for (int i = 0; i < DOFn; i++) {
            for (int a = 0; a < 3; a++)
                iMu[a] += mu2[i * 3 + a];
        }
        //calculate tolerance
        tol=0;

        for (int a = 0; a < 3; a++) {
            tol += abs(iMu[a] - Mu_prev[a]);
        }
    }
    //save pi_i to Pi_all
    for (int i = 0; i < DOFn; i++)
        for (int j = 0; j < 3 ;j++) {
            iMu_all[i * 3 + j] = mu2[i * 3 + j];
        }
}

void ForceFieldPolar::Dk_Pi_tensor(int k, std::vector<double>& iPi_all,
                 std::vector<double>& Tij_all, int conf_diff, std::vector<double>& DkPi) {
    //subroutine to calculate partial derivative of the total polairzibility tensor_ab (3*3)
    //with respect to particular rk_u(3), output is (3*3*3) dimensional DkPi[27].
    //e.g.
    //for (a=0; a<3; a++)
    //    for (b=0; b<3; b++)
    //        for (u=0; u<3; u++) Pi_dot[3*a+b] = DkPi[9*a+3*b+u] * V[k][u];
    
    int n_iter(0);
    
    double DkPi_prev[27]; //old DkPi

    std::vector<double> Dkpi1, Dkpi2, DkTkj, DkTik;
    Dkpi1.resize(DOFn * 27);
    Dkpi2.resize(DOFn * 27);
    DkTkj.resize(DOFn * 27);
    DkTik.resize(DOFn * 27);

    double tol(1); //tolerance
    double Tij[3][3]={0,0,0,0,0,0,0,0,0};
    double DkTij[27];
    
    double alpha[DOFn]; //assign single molecular polarizability

    for (int i = 0; i < DOFn; i++) {
        alpha[i] = ffPolars[i].alpha/1000;
    }

    std::vector<std::set<int> > bonded12(DOFn); 
    std::vector<std::set<int> > bonded13(DOFn);
    if (damping_type == "TTM3"||damping_type == "TTM4") {
        addTholeExclusionsToSet(bonded12, bonded13);
    }

    //set Dk_Pi =0
    for (int i = 0; i < 27; i++) DkPi[i] = DkPi_prev[i] = 0;
    
    //initial guess of the Dkpi2
    for (int i = 0; i < DOFn; i++)
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                for (int u = 0; u < 3; u++) {
                    if (a==b && b==u) {
                        Dkpi2[27 * i + 9 * a + 3 * b + u] = 0; //1
                    } 
                    else { 
                        Dkpi2[27 * i + 9 * a + 3 * b + u] = 0;
                    }
                }
    
    //calculate Tij_all
    if (conf_diff == 1) {
        for (int i = 0; i < DOFn ; i++)
            for (int j = 0; j < DOFn ; j++) {
                if (i == j) { //Tii=0
                    for (int m = 0; m < 9; m++) 
                        Tij_all[9 * DOFn * i + 9 * j + m] = 0;
                    continue;
                }
                //calculate Tij
                double deltaR[5];
                double rm1, rm2, rm3, rm5;
                const Vec3& atomCoordinatesI = R[i];
                const Vec3& atomCoordinatesJ = R[j];
                // get deltaR 2 atoms, deltaR = Ri - Rj
                GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
                rm1 = 1.0 / deltaR[4];
                rm2 = rm1 * rm1;
                rm3 = rm1 * rm2;
                rm5 = rm2 * rm3;
                double S0r, S1r, S2r, S3r;
                // choose whitch thole model damping type.
                if (damping_type == "TTM0") 
                    TTM0Model(S0r, S1r, S2r, S3r);
                else if (damping_type == "TTM1") 
                    TTM1Model(S0r, S1r, S2r, S3r);
                //else if (damping_type == "TTM2") 
                //    double alphaEwald, double deltaR, double& S0r, double& S1r, double& S2r, double& S3r) {
                else {
                    double alphaij;
                    alphaij = alpha[i] * alpha[j];
                    alphaij = pow((double)alphaij,double(1.0/6));            
                    alphaij = deltaR[4] * 1.0 / alphaij;
                    // Thole model daming parameter for chosen pair od atom i and j
                    double dampingParam;
                    // Judge the pair i, j belong to the 12, 13, 14, or noboned part.
                    dampingParam = param.getDouble("dampingParam3");
                    if (bonded12[i].find(j) == bonded12[i].end()||bonded12[j].find(i) == bonded12[j].end()) 
                        dampingParam = param.getDouble("dampingParam1");;
                    if (bonded13[i].find(j) == bonded13[i].end()||bonded13[j].find(i) == bonded13[j].end()) 
                        dampingParam = param.getDouble("dampingParam2");;

                    if (damping_type == "TTM3")
                        TTM3Model(dampingParam, alphaij, S0r, S1r, S2r, S3r);
                    else if (damping_type == "TTM4")
                        TTM4Model(dampingParam, alphaij, S0r, S1r, S2r, S3r);
                }
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {                   
                        Tij[a][b] = 3 * deltaR[a] * deltaR[b] * rm5 * S2r;
                        if (a == b) {
                            Tij[a][b] -= Delta(a,b) * rm3 * S1r; 
                        }
                    }
                }
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++) {
                        Tij_all[9 * DOFn * i + 9 * j + 3 * a + b] = Tij[a][b];
                    }
            }
    }
    
    //calculate DkTkj (fixed k, j=0,...N-1)
    for (int j = 0; j < DOFn; j++) {
        Dk_Tij_tensor(k, k, j, DkTij);
        for (int a = 0; a < 27 ; a++) DkTkj[27 * j + a] = DkTij[a];
    }
    
    //calculate DkTik (fixed k, i=0,...N-1)
    for (int i = 0; i < DOFn ;i++) {
        Dk_Tij_tensor(k, i, k, DkTij);
        for (int a = 0; a < 27; a++) DkTik[27 * i + a] = DkTij[a];
    }
    
    //************ iterative calculation of Dkpi1 or 2 ****************
    while ( tol > TOLERANCE && n_iter < iter_max ) {
        //save DkPi_prev <= Dkpi
        for (int m = 0; m < 27; m++) DkPi_prev[m] = DkPi[m];
        
        //Dkpi1 << Dkpi2, and set Dkpi2=0
        for (int i = 0; i < DOFn; i++)
            for (int m = 0; m < 27; m++) {
                Dkpi1[27 * i + m] = Dkpi2[27 * i + m];
                Dkpi2[27 * i + m] = 0;
            }
        
        //start new iteration
        n_iter++;
        for (int i = 0; i < DOFn; i++ ){
            if (i == k) {
                for (int j = 0; j < DOFn; j++) {
                    if (j == k) continue;
                    for (int a = 0; a < 3; a++)
                        for (int b = 0; b < 3; b++)
                            for (int u = 0; u < 3; u++)
                                for (int c = 0; c < 3; c++)
                                    Dkpi2[27 * i + 9 * a + 3 * b + u] +=
                                    DkTkj[27 * j + 9 * a + 3 * c + u] * iPi_all[9 * j + 3 * c + b]
                                    + Tij_all[9 * DOFn * i + 9 * j + 3 * a + c] * Dkpi1 [27 * j + 9 * c + 3 * b + u];
                }
                for (int m = 0; m < 27; m++) 
                    Dkpi2[27 * i + m] *= alpha[i];
            }
            else { // i!=k
                for (int u = 0; u < 3; u++)
                    for (int a = 0; a < 3; a++)
                        for (int b = 0; b < 3; b++)
                            for (int c = 0; c < 3; c++) {
                                Dkpi2[27 * i + 9 * a + 3 * b + u] +=
                                DkTik[27 * i + 9 * a + 3 * c + u] * iPi_all[9 * k + 3 * c + b];
                                for (int j = 0; j < DOFn; j++) {
                                    if (j == i) continue;
                                    Dkpi2[27 * i + 9 * a + 3 * b + u] += Tij_all[9 * DOFn * i + 9 * j + 3 * a + c]
                                    * Dkpi1[27 * j + 9 * c + 3 * b + u];
                                }
                            }
                for (int m = 0; m < 27; m++) 
                    Dkpi2[27 * i + m] *= alpha[i];
            }
        }
        
        //set DkPi[27]=0
        for (int m = 0; m < 27; m++) DkPi[m] = 0;
        
        //sum of i to get new DkPi
        for (int i = 0; i < DOFn; i++)
            for (int m = 0; m < 27; m++) DkPi[m] += Dkpi2[27 * i + m];
        
        //calculate tolerance
        tol = 0;
        for (int m = 0; m < 27 * DOFn; m++) 
            tol += abs(Dkpi2[m] - Dkpi1[m]);
    }
}

void ForceFieldPolar::Dk_Mu_tensor(int k, std::vector<double>& iMu_all,
                 std::vector<double>& Tij_all, int conf_diff, std::vector<double>& DkMu) {
    //subroutine to calculate partial derivative of the total induced dipole tensor_a (3)
    //with respect to particular rk_u(3), output is (3*3) dimensional DkMu[9].
    
    int n_iter(0);
    
    double DkMu_prev[9]; //old DkPi

    std::vector<double> Dkmu1, Dkmu2, DkTkj, DkTik;
    Dkmu1.resize(DOFn * 9);
    Dkmu2.resize(DOFn * 9);
    DkTkj.resize(DOFn * 27);
    DkTik.resize(DOFn * 27);

    double tol(1); //tolerance
    double Tij[3][3]={0,0,0,0,0,0,0,0,0};
    double DkTij[27];
    
    double alpha[DOFn]; //assign single molecular polarizability

    for (int i = 0; i < DOFn; i++) {
        alpha[i] = ffPolars[i].alpha/1000;
    }

    std::vector<std::set<int> > bonded12(DOFn); 
    std::vector<std::set<int> > bonded13(DOFn);
    if (damping_type == "TTM3"||damping_type == "TTM4") {
        addTholeExclusionsToSet(bonded12, bonded13);
    }

    //set Dk_Mu =0
    for (int i = 0; i < 9; i++) DkMu[i] = DkMu_prev[i] = 0;
    
    //initial guess of the Dkmu2
    for (int i = 0; i < DOFn; i++)
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3; b++)
                    if (a==b) {
                        Dkmu2[9 * i + 3 * a + b] = 0; //1
                    } 
                    else { 
                        Dkmu2[9 * i + 3 * a + b] = 0;
                    }
    
    //calculate Tij_all
    if (conf_diff == 1) {
        for (int i = 0; i < DOFn ; i++)
            for (int j = 0; j < DOFn ; j++) {
                if (i == j) { //Tii=0
                    for (int m = 0; m < 9; m++) 
                        Tij_all[9 * DOFn * i + 9 * j + m] = 0;
                    continue;
                }
                //calculate Tij
                double deltaR[5];
                double rm1, rm2, rm3, rm5;
                const Vec3& atomCoordinatesI = R[i];
                const Vec3& atomCoordinatesJ = R[j];
                // get deltaR 2 atoms, deltaR = Ri - Rj
                GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
                rm1 = 1.0 / deltaR[4];
                rm2 = rm1 * rm1;
                rm3 = rm1 * rm2;
                rm5 = rm2 * rm3;
                double S0r, S1r, S2r, S3r;
                // choose whitch thole model damping type.
                if (damping_type == "TTM0") 
                    TTM0Model(S0r, S1r, S2r, S3r);
                else if (damping_type == "TTM1") 
                    TTM1Model(S0r, S1r, S2r, S3r);
                //else if (damping_type == "TTM2") 
                //    double alphaEwald, double deltaR, double& S0r, double& S1r, double& S2r, double& S3r) {
                else {
                    double alphaij;
                    alphaij = alpha[i] * alpha[j];
                    alphaij = pow((double)alphaij,double(1.0/6));            
                    alphaij = deltaR[4] * 1.0 / alphaij;
                    // Thole model daming parameter for chosen pair od atom i and j
                    double dampingParam;
                    // Judge the pair i, j belong to the 12, 13, 14, or noboned part.
                    dampingParam = param.getDouble("dampingParam3");
                    if (bonded12[i].find(j) == bonded12[i].end()||bonded12[j].find(i) == bonded12[j].end()) 
                        dampingParam = param.getDouble("dampingParam1");;
                    if (bonded13[i].find(j) == bonded13[i].end()||bonded13[j].find(i) == bonded13[j].end()) 
                        dampingParam = param.getDouble("dampingParam2");;

                    if (damping_type == "TTM3")
                        TTM3Model(dampingParam, alphaij, S0r, S1r, S2r, S3r);
                    else if (damping_type == "TTM4")
                        TTM4Model(dampingParam, alphaij, S0r, S1r, S2r, S3r);
                }
                for (int a = 0; a < 3; a++) {
                    for (int b = 0; b < 3; b++) {                   
                        Tij[a][b] = 3 * deltaR[a] * deltaR[b] * rm5 * S2r;
                        if (a == b) {
                            Tij[a][b] -= Delta(a,b) * rm3 * S1r; 
                        }
                    }
                }
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++) {
                        Tij_all[9 * DOFn * i + 9 * j + 3 * a + b] = Tij[a][b];
                    }
            }
    }
    
    //calculate DkTkj (fixed k, j=0,...N-1)
    for (int j = 0; j < DOFn; j++) {
        Dk_Tij_tensor(k, k, j, DkTij);
        for (int a = 0; a < 27 ; a++) DkTkj[27 * j + a] = DkTij[a];
    }
    
    //calculate DkTik (fixed k, i=0,...N-1)
    for (int i = 0; i < DOFn ;i++) {
        Dk_Tij_tensor(k, i, k, DkTij);
        for (int a = 0; a < 27; a++) DkTik[27 * i + a] = DkTij[a];
    }
    
    //************ iterative calculation of Dkmu1 or 2 ****************
    while ( tol > TOLERANCE && n_iter < iter_max ) {
        //save DkMu_prev <= Dkmu
        for (int m = 0; m < 9; m++) DkMu_prev[m] = DkMu[m];
        
        //Dkmu1 << Dkmu2, and set Dkmu2=0
        for (int i = 0; i < DOFn; i++)
            for (int m = 0; m < 9; m++) {
                Dkmu1[9 * i + m] = Dkmu2[9 * i + m];
                Dkmu2[9 * i + m] = 0;
            }
        
        //start new iteration
        n_iter++;
        for (int i = 0; i < DOFn; i++ ){
            if (i == k) {
                for (int j = 0; j < DOFn; j++) {
                    if (j == k) continue;
                    for (int a = 0; a < 3; a++)
                        for (int u = 0; u < 3; u++)
                            for (int c = 0; c < 3; c++)
                                Dkmu2[9 * i + 3 * a + u] +=
                                DkTkj[27 * j + 9 * a + 3 * c + u] * iMu_all[3 * j + c]
                                + Tij_all[9 * DOFn * i + 9 * j + 3 * a + c] * Dkmu1 [9 * j + 3 * c + u];
                }
                for (int m = 0; m < 9; m++) 
                    Dkmu2[9 * i + m] *= alpha[i];
            }
            else { // i!=k
                for (int u = 0; u < 3; u++)
                    for (int a = 0; a < 3; a++)
                        for (int c = 0; c < 3; c++) {
                            Dkmu2[9 * i + 3 * a + u] +=
                            DkTik[27 * i + 9 * a + 3 * c + u] * iMu_all[3 * k + c];
                            for (int j = 0; j < DOFn; j++) {
                                if (j == i) continue;
                                Dkmu2[9 * i + 3 * a + u] += Tij_all[9 * DOFn * i + 9 * j + 3 * a + c]
                                * Dkmu1[9 * j + 3 * c + u];
                            }
                        }
                for (int m = 0; m < 9; m++) 
                    Dkmu2[9 * i + m] *= alpha[i];
            }
        }
        
        //set DkMu[9]=0
        for (int m = 0; m < 9; m++) DkMu[m] = 0;
        
        //sum of i to get new DkMu
        for (int i = 0; i < DOFn; i++)
            for (int m = 0; m < 9; m++) DkMu[m] += Dkmu2[9 * i + m];
        
        //calculate tolerance
        tol = 0;
        for (int m = 0; m < 9 * DOFn; m++) 
            tol += abs(Dkmu2[m] - Dkmu1[m]);
    }
}

void ForceFieldPolar::Dk_Tij_tensor(int k, int i, int j, double DkTij[27]) {

    std::vector<std::set<int> > bonded12(DOFn); 
    std::vector<std::set<int> > bonded13(DOFn);
    if (damping_type == "TTM3"||damping_type == "TTM4") {
        addTholeExclusionsToSet(bonded12, bonded13);
    }

    double alpha[DOFn]; //assign single molecular polarizability

    for (int i = 0; i < DOFn; i++) {
        alpha[i] = ffPolars[i].alpha/1000;
    }

    if ((k != i && k != j) || i == j) {
        for (int m = 0; m < 27 ;m++) DkTij[m] = 0;
        return;
    }
    
    if (i != j) {
        double deltaR[5];
        const Vec3& atomCoordinatesI = R[i];
        const Vec3& atomCoordinatesJ = R[j];
        // get deltaR 2 atoms, deltaR = Ri - Rj
        GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
        double rm1, rm2, rm3, rm5, rm7;    
        rm1 = 1.0 / deltaR[4];
        rm2 = rm1 * rm1;
        rm3 = rm1 * rm2;
        rm5 = rm2 * rm3;
        rm7 = rm5 * rm2;
        double alphaij;
        double S0r, S1r, S2r, S3r;
        // choose whitch thole model damping type.
        if (damping_type == "TTM0") 
            TTM0Model(S0r, S1r, S2r, S3r);
        else if (damping_type == "TTM1") 
            TTM1Model(S0r, S1r, S2r, S3r);
        //else if (damping_type == "TTM2") 
        //    double alphaEwald, double deltaR, double& S0r, double& S1r, double& S2r, double& S3r) {
        else {
        // Thole model daming parameter for chosen pair od atom i
            double dampingParam;
            alphaij = alpha[i] * alpha[j];
            alphaij = pow((double)alphaij,double(1.0/6));            
            alphaij = deltaR[4] * 1.0 / alphaij;
            // Judge the pair i, j belong to the 12, 13, 14, or noboned part.
            dampingParam = param.getDouble("dampingParam3");
            if (bonded12[i].find(j) == bonded12[i].end()||bonded12[j].find(i) == bonded12[j].end()) 
                dampingParam = param.getDouble("dampingParam1");;
                if (bonded13[i].find(j) == bonded13[i].end()||bonded13[j].find(i) == bonded13[j].end()) 
                    dampingParam = param.getDouble("dampingParam2");;
                if (damping_type == "TTM3")
                    TTM3Model(dampingParam, alphaij, S0r, S1r, S2r, S3r);
                else if (damping_type == "TTM4")
                    TTM4Model(dampingParam, alphaij, S0r, S1r, S2r, S3r);
        }
            
        for (int a = 0; a < 3; a++)
            for (int b = 0; b < 3;b++)
                for (int c = 0; c < 3; c++) {
                DkTij[9 * a + 3 * b + c] = (Delta(k, i) - Delta(k, j)) *
                    ( 3 * (deltaR[a] * Delta(b, c) + deltaR[b] * Delta(a, c) + deltaR[c] * Delta(a, b)) * rm5 * S2r
                    - 15 * deltaR[a] * deltaR[b] * deltaR[c] * rm7 * S3r);
                }
    } 
    else {
        std::cout << "Dk_Tij_tensor error (i==j)." << std::endl;
    }
}


void ForceFieldPolar::Perturb_Polar(std::vector<double>& DPi, std::vector<Vec3>& perturbforce){
    //Perturb the system with Raman interaction, VR = - 1/2 E1_mu * Pi_munu * E2_nu (* direction)
    //result in change of force of atom k, DeltaF(k)_u = - d VR/ d r(k)_u
    //               =  1/2 E1_mu * d Pi_munu/ dr(k)_u * E2_nu
    double perturb_fields; //perturb energy (epsilon) = perturb_fields(epsilon/sigma^3)* Pi(sigma^3)
    // perturb elec size
    std::vector<int> perturb_parameters;
    // 0 1 2 3 4 5 6 7 8
    SplitString(perturb_parameters, param.getStr("perturb_parameters"));
    int mu, nu;
    int E1, E2;
    if (perturb_parameters[5] == 0) {
        mu = 1;
        nu = 2;
        E1 = perturb_parameters[6];
        E2 = perturb_parameters[7];
    }
    else if (perturb_parameters[6] == 0) {
        mu = 0;
        nu = 2;
        E1 = perturb_parameters[5];
        E2 = perturb_parameters[7];
    }
    else if (perturb_parameters[7] == 0) {
        mu = 0;
        nu = 1;
        E1 = perturb_parameters[5];
        E2 = perturb_parameters[6];
    }

    perturb_fields = - 0.5 * E1 * E2 * 1e20 * Fm2 / mol;// unit KJ/(mol*nm)
    
    int munu;
    munu = 9 * mu + 3 * nu;       //refer Pi_dot[3*a+b] += DkPi[9*a+3*b+u] * V[k][u] (for k=0...N-1)
    for (int k = 0; k < DOFn; k++) {     //atom k
        for (int u = 0; u < 3; u++) { //x,y,z component for the r(k)_u
            perturbforce[k][u] -= DPi[27 * k + munu + u] * perturb_fields;  //F = - dV/dr         
        }
    }
}

void ForceFieldPolar::Perturb_Dipole(std::vector<double>& DMu, std::vector<Vec3>& perturbforce){
    //Perturb the system with Raman interaction, VR = - 1/2 E1_mu * Pi_munu * E2_nu (* direction)
    //result in change of force of atom k, DeltaF(k)_u = - d VR/ d r(k)_u
    //               =  1/2 E1_mu * d Pi_munu/ dr(k)_u * E2_nu
    double perturb_fields; //perturb energy (epsilon) = perturb_fields(epsilon/sigma^3)* Pi(sigma^3)
    // perturb elec size
    std::vector<int> perturb_parameters;
    SplitString(perturb_parameters, param.getStr("perturb_parameters"));
    int mu, nu;
    double E_per[3];
    E_per[0] = perturb_parameters[2];
    E_per[1] = perturb_parameters[3];
    E_per[2] = perturb_parameters[4];
    
    //perturb_fields = - 0.5 * E1 * E2 * 1e20 * Fm2 / mol;// unit KJ/(mol*nm)

    for (int k = 0; k < DOFn; k++) {     //atom k
        for (int u = 0; u < 3; u++) { //x,y,z component for the r(k)_u
            for (int a = 0; a < 3; a++)
                perturbforce[k][u] -= DMu[9 * k + 3 * u + a] * E_per[u];  //F = - dV/dr         
        }
    }
}

void ForceFieldPolar::Dk_E_perm(int k, std::vector<double>& DkE_perm, std::vector<double>& DkE) {
    // This part is to calculate the Deriv of E_perm, and divided into two parts: rec part, and real part.
    static const double epsilon = 1.0;
    // set the Ewald method some important parameters
    // For Ewald Sum method parameters: alphaEwald                                                                                     
    double Ewald_tolerance = param.getDouble("Ewald_tolerance");
    double alphaEwald;
    int kmax[3]={0, 0, 0};
    double cutlength = param.getDouble("cutoff");;
    alphaEwald = sqrt(-log(2*Ewald_tolerance))/cutlength; 
    // get the Kmax on recipical part
    double error;
    error = Ewald_tolerance;
    while (error >= Ewald_tolerance){
        kmax[0]++;
        error = kmax[0]*sqrt(periodicBoxVectors[0][0]*alphaEwald)*exp(-pow((pi*kmax[0]/(periodicBoxVectors[0][0]*alphaEwald)), 2))/20;
    }
    error = Ewald_tolerance;
    while (error >= Ewald_tolerance){
        kmax[1]++;
        error = kmax[1]*sqrt(periodicBoxVectors[1][1]*alphaEwald)*exp(-pow((pi*kmax[1]/(periodicBoxVectors[1][1]*alphaEwald)), 2))/20;
    }
    error = Ewald_tolerance;
    while (error >= Ewald_tolerance){
        kmax[2]++;
        error = kmax[2]*sqrt(periodicBoxVectors[2][2]*alphaEwald)*exp(-pow((pi*kmax[2]/(periodicBoxVectors[2][2]*alphaEwald)), 2))/20;
    }
    double factorEwald = -1 / (4 * alphaEwald * alphaEwald);
    double recipCoeff = ONE_4PI_EPS0 * 4 * pi/(periodicBoxVectors[0][0] * periodicBoxVectors[1][1] * periodicBoxVectors[2][2]) /epsilon;
    double SQRT_PI = sqrt(pi);
    double TWO_PI = 2.0 * pi;

    double Rk[2 * kmax[0] + 1][3];// 2 * pi * k /L
    for (int i = 0; i < 2 * kmax[0] + 1; i++){
        Rk[i][0] = TWO_PI * double(i - kmax[0]) / periodicBoxVectors[0][0];
    }  
    for (int i = 0; i < 2 * kmax[1] + 1; i++){
        Rk[i][1] = TWO_PI * double(i - kmax[1]) / periodicBoxVectors[1][1];
    }    
    for (int i = 0; i < 2 * kmax[2] + 1; i++){
        Rk[i][2] = TWO_PI * double(i - kmax[2]) / periodicBoxVectors[2][2];
    }     
    double KSQ, AK;
    double alpha[DOFn]; //assign single molecular polarizability
    
    for (int i = 0; i < DOFn; i++) {
        alpha[i] = ffPolars[i].alpha/1000;
    }

    for (int i = 0; i < DOFn; i++) {
        for (int j = 0; j < DOFn; j++) {
            if (i == j) continue;
            double delta_ijk = Delta(k,i) - Delta(k,j);
            double charge = ffElec[j].charge;
            double deltaR[5];
            const Vec3& atomCoordinatesI = R[i];
            const Vec3& atomCoordinatesJ = R[j];
            // get deltaR 2 atoms, deltaR = Ri - Rj
            GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
            for (int a = 0; a < 2 * kmax[0] + 1; a++)
                for (int b = 0; a < 2 * kmax[0] + 1; b++)
                    for (int c = 0; a < 2 * kmax[0] + 1; c++) {
                        double abc = a*a + b*b + c*c;
                        if (abc == 0) continue;
                        KSQ = Rk[a][0] * Rk[a][0] + Rk[b][1] * Rk[b][1] + Rk[c][2] * Rk[c][2];                
                        AK = exp(KSQ * factorEwald) / KSQ; 
                        double Rkk[3];
                        Rkk[0] = Rk[a][0];
                        Rkk[1] = Rk[b][1];
                        Rkk[2] = Rk[c][2];        
                        for (int m = 0; m <3; m++)
                            for (int n = 0; n < 3; n++) {
                                DkE_perm[i * 9 + 3 * m + n] += charge * Rkk[m] * Rkk[n] * cos(Rkk[0] * deltaR[0] + Rkk[1] * deltaR[1] + Rkk[2] * deltaR[2]) * delta_ijk;
                                DkE_perm[i * 9 + 3 * m + n] *= recipCoeff * AK;
                            }
                    }
        }
        for (int m = 0; m <3; m++) 
            for (int n = 0; n < 3; n++) {
                DkE_perm[i * 9 + 3 * m + n] *= alpha[i];
            }
    } 
    // This is for the real part 
    double alpha1 = - 4 * alphaEwald * alphaEwald * alphaEwald / sqrpi;
    double alpha2 =  2 * alphaEwald / sqrpi;
    for (int i = 0; i < DOFn; i++) {
        for (int j = 0; j < DOFn; j++) {
            if (i == j) continue;
            double delta_ijk = Delta(k,i) - Delta(k,j);
            double charge = ffElec[j].charge;
            double deltaR[5];
            const Vec3& atomCoordinatesI = R[i];
            const Vec3& atomCoordinatesJ = R[j];
            // get deltaR 2 atoms, deltaR = Ri - Rj
            GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
            double Exp = exp(- alphaEwald * alphaEwald * deltaR[4]); 
            for (int a = 0; a < 3; a++) 
                for (int b = 0; b < 3; b++)
                    DkE_perm[i * 9 + 3 * a + b] += charge * (alpha1 * Exp * deltaR[a] * deltaR[b]/deltaR[4] 
                                                + alpha2 * Exp * (-2 * deltaR[a] * deltaR[b]/(deltaR[4] * deltaR[4]) + Delta(a,b)/deltaR[4])) * delta_ijk
                                                + charge * (alpha2 * Exp * deltaR[a] * deltaR[b]/(deltaR[4]* deltaR[4])
                                                + erfc(alphaEwald * deltaR[3]) * (-3 * deltaR[a] * deltaR[b]/(deltaR[4] * deltaR[4] * deltaR[3]) + Delta(a,b)/(deltaR[4] * deltaR[3]))) * delta_ijk;
        }
        for (int a = 0; a< 3; a++) 
            for (int b = 0; b < 3; b++) 
                DkE_perm[i * 9 + 3 * a + b] *= ONE_4PI_EPS0  * alpha[i];
    }   
    for (int m = 0; m < 9; m++) {
        DkE[m] = 0;
        for (int j=0; j<DOFn; j++) DkE[m] += DkE_perm[j*9 + m];
    }    

}
