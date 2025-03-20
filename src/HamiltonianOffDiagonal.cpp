/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: Jan. 15, 2022                                                *
 * -------------------------------------------------------------------------- */

#include "HamiltonianOffDiagonal.h"

void HamiltonianOffDiagonal::init() {
    // * 1 get the total atom of the system
    DOFn    = param.getInt("DOFn");
    DOFe    = param.getInt("DOFe");
    setNewVdwParameters();
    R.resize(DOFn, Vec3());
}

void HamiltonianOffDiagonal::setNewVdwParameters() {
    for (int i = 0; i < topologies.molecules.size(); i++) {
        for (int j = 0; j < topologies.molecules[i].moleculeNumber; j++) {
            const Topology::MoleculeTypes& moleculeType = topologies.moleculeTypes[topologies.molecules[i].moleculeTypeIndex];
            // * Add particles with mass to System from Atoms.
            for (int k = 0; k < moleculeType.newvdw.size(); k++) {
                int indexI = moleculeType.newvdw[k].atomTypeIndex;
                double vdwparam0 = moleculeType.newvdw[k].vdwparam0; 
                double vdwparam1 = moleculeType.newvdw[k].vdwparam1; 
                double vdwparam2 = moleculeType.newvdw[k].vdwparam2; 
                double vdwparam3 = moleculeType.newvdw[k].vdwparam3; 
                double vdwparam4 = moleculeType.newvdw[k].vdwparam4; 
                double vdwparam5 = moleculeType.newvdw[k].vdwparam5; 
                double vdwparam6 = moleculeType.newvdw[k].vdwparam6; 
                double vdwparam7 = moleculeType.newvdw[k].vdwparam7; 
                double vdwparam8 = moleculeType.newvdw[k].vdwparam8; 
                ffnewvdw.push_back({indexI, vdwparam0, vdwparam1, vdwparam2, vdwparam3, vdwparam4, 
                                    vdwparam5, vdwparam6, vdwparam7, vdwparam8});     
            } 
        }
    }
} 

double HamiltonianOffDiagonal::getDiabaticCoupling(int i, int j) {
    double Diabaticenergy = 0;
    double deltaR[5];
    std::cout.precision(12);
    if (param.getStr("LJ_type") == "LJ8_8") {
        const Vec3& atomCoordinatesI = R[i];
        const Vec3& atomCoordinatesJ = R[j];      
        GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
        double deltaR3 = deltaR[3] * deltaR[4];
        double deltaR4 = deltaR[3] * deltaR[3];
        double deltaR5 = deltaR3 * deltaR[3];
        double deltaR6 = deltaR3 * deltaR3;
        double deltaR7 = deltaR3 * deltaR4;
        double deltaR8 = deltaR4 * deltaR4; 
        Diabaticenergy += ffnewvdw[2].vdwparam8 * deltaR8 + ffnewvdw[2].vdwparam7 * deltaR7 + 
                ffnewvdw[2].vdwparam6 * deltaR6 + ffnewvdw[2].vdwparam5 * deltaR5 + 
                ffnewvdw[2].vdwparam4 * deltaR4 + ffnewvdw[2].vdwparam3 * deltaR3 + 
                ffnewvdw[2].vdwparam2 * deltaR[3] + ffnewvdw[2].vdwparam1 * deltaR[4] + ffnewvdw[2].vdwparam0;                
    //std::cout<<"i j   "<<i<<"   "<<R[i][0]<<" "<<R[i][1]<<"  "<< R[i][2]<<std::endl;
    //std::cout<<"i j   "<<j<<"   "<<R[j][0]<<" "<<R[j][1]<<"  "<< R[j][2]<<std::endl;
    //std::cout<<"Diabaticenergy:  i: "<<i<<" j: "<<j<<" "<<Diabaticenergy<<std::endl;
    
    }
    else 
           throw std::runtime_error("ERROR: Unsupported LJ_type: " + param.getStr("LJ_type")); 
    return Diabaticenergy;
}

std::vector<Vec3> HamiltonianOffDiagonal::getDiabaticCouplingforce(int i, int j) {
    double Diabaticenergy = 0;
    double deltaR[5];
    std::vector<Vec3> offdiagonalforce;
    offdiagonalforce.resize(DOFn, Vec3());
    std::fill(offdiagonalforce.begin(), offdiagonalforce.end(), Vec3());
    std::cout.precision(12);
    if (param.getStr("LJ_type") == "LJ8_8") {
        const Vec3& atomCoordinatesI = R[i];
        const Vec3& atomCoordinatesJ = R[j];      
        GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
        double deltaR3 = deltaR[3] * deltaR[4];
        double deltaR4 = deltaR[3] * deltaR[3];
        double deltaR5 = deltaR3 * deltaR[3];
        double deltaR6 = deltaR3 * deltaR3;
        double deltaR7 = deltaR3 * deltaR4;
        double deltaR8 = deltaR4 * deltaR4; 

        for (int kk = 0; kk < 3; kk++) {
            double force = (8 * ffnewvdw[0].vdwparam8 * deltaR6 + 7 * ffnewvdw[0].vdwparam7 * deltaR5 + 
                            6 * ffnewvdw[0].vdwparam6 * deltaR4 + 5 * ffnewvdw[0].vdwparam5 * deltaR3 + 
                            4 * ffnewvdw[0].vdwparam4 * deltaR[3] + 3 * ffnewvdw[0].vdwparam3 * deltaR[4] + 
                            2 * ffnewvdw[0].vdwparam2 ) * deltaR[kk] + ffnewvdw[0].vdwparam1;     
            offdiagonalforce[i][kk] -= force;
            offdiagonalforce[j][kk] += force;
        }
    }
    else 
           throw std::runtime_error("ERROR: Unsupported LJ_type: " + param.getStr("LJ_type")); 
    return offdiagonalforce;
}

void HamiltonianOffDiagonal::setPositions(const std::vector<Vec3>& positions) {
    R = positions;
}

void HamiltonianOffDiagonal::GetDeltaRPeriodic(const Vec3& atomCoordinatesI, const Vec3& atomCoordinatesJ,
                                       const Vec3 (&periodicperiodicBoxVectors)[3], double deltaR[5]) {                                  
   Vec3 diff = atomCoordinatesJ-atomCoordinatesI;
   diff[2] = diff[2]-periodicBoxVectors[2][2]*floor(diff[2]/periodicBoxVectors[2][2]+0.5);
   diff[1] = diff[1]-periodicBoxVectors[1][1]*floor(diff[1]/periodicBoxVectors[1][1]+0.5);
   diff[0] = diff[0]-periodicBoxVectors[0][0]*floor(diff[0]/periodicBoxVectors[0][0]+0.5);

   deltaR[0] = diff[0];
   deltaR[1] = diff[1];
   deltaR[2] = diff[2];
   deltaR[3] = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
   deltaR[4] = sqrt(deltaR[3]);    
}

void HamiltonianOffDiagonal::setPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c) {  
    periodicBoxVectors[0] = a;
    periodicBoxVectors[1] = b;
    periodicBoxVectors[2] = c;
}

void HamiltonianOffDiagonal::setTopologies(const Topology& fftopologies) {
    topologies = fftopologies;
}