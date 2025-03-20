/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: Jan. 14, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "ForceFieldBase.h"

void ForceFieldBase::init() {
    // * 1 get the total atom of the system
    DOFn    = param.getInt("DOFn");
    std::cout.precision(12);
    // * 2 begin to ge the control of the force contribution from the file
    force_contribution.resize(num_force_parts, true);
    std::vector<std::string> contri;
    std::string turn_off_force_contribution = param.getStr("turn_off_force_contribution");
    SplitString(contri, turn_off_force_contribution);
    // * 3 Determine whether to calculate the Bond, Angle, Dihedral, LJ, Elec.
    for (int i = 0; i < contri.size(); i++) {
        if (contri[i] == "Bond") 
            force_contribution[1] = false;
        else if (contri[i] == "Angle") 
            force_contribution[2] = false;   
        else if (contri[i] == "Dihedral") 
            force_contribution[3] = false; 
        else if (contri[i] == "LJ") 
            force_contribution[4] = false; 
        else if (contri[i] == "Elec") 
            force_contribution[5] = false; 
        else if (contri[i] == "Polar") 
            force_contribution[6] = false; 
        else 
           throw std::runtime_error("ERROR: Unsupported turn_off_force_contribution " + contri[i]);     
    }
    // * 4 Resize all vectors in this class and set their elements to 0.
    setForceFieldParameters();
    R.resize(DOFn, Vec3());
    V.resize(DOFn, Vec3());
    F.resize(DOFn, Vec3()); 
    forceParts.resize(num_force_parts, F);
    energyParts.resize(num_force_parts, 0);  
}

double ForceFieldBase::getPotentialEnergy(bool includeForces) {
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
    if (force_contribution[5]) 
        calculateElecForceandEnergy(F, forceParts[5], PE, energyParts[5], includeForces);     
    return PE;
}

const std::vector<Vec3>& ForceFieldBase::getForces() {
    return F;
}

const std::vector<std::vector<Vec3>>& ForceFieldBase::getForceParts() {
    return forceParts;
}

const std::vector<double>& ForceFieldBase::getEnergyParts() {
    return energyParts;
}

double ForceFieldBase::getPeriodicBoxVolume() const{
    return periodicBoxVectors[0].dot(periodicBoxVectors[1].cross(periodicBoxVectors[2]));
}

void ForceFieldBase::setPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c) {  
    periodicBoxVectors[0] = a;
    periodicBoxVectors[1] = b;
    periodicBoxVectors[2] = c;
}

void ForceFieldBase::setPositions(const std::vector<Vec3>& positions) {
    R = positions;
}

void ForceFieldBase::setVelocities(const std::vector<Vec3>& velocities) {
    V = velocities;
}

void ForceFieldBase::setTopologies(const Topology& fftopologies) {
    topologies = fftopologies;
}

void ForceFieldBase::setForceFieldParameters() {
    int currentNumAtoms = 0;
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
            // * Add particles with charge and LJ parameters to NonbondedForce from AtomTypes.
            for (int k = 0; k < moleculeType.atoms.size(); k++) {
                int indexI = k + currentNumAtoms;
                double charge = moleculeType.atoms[k].charge;
                double mass = topologies.atomTypes[moleculeType.atoms[k].atomTypeIndex].mass;
                double sigma = topologies.atomTypes[moleculeType.atoms[k].atomTypeIndex].sigma;
                double epsilon = topologies.atomTypes[moleculeType.atoms[k].atomTypeIndex].epsilon;
                ffLJ.push_back({indexI, sigma, epsilon});    
                ffElec.push_back({indexI, charge, mass});       
            }   
            // * Add Bonds pair for NonbondedForce calculation 1-4 Interations
            for (int k = 0; k < moleculeType.pairs.size(); k++) {
                int indexI = (moleculeType.pairs[k].atomIndex[0] - 1) + currentNumAtoms;
                int indexJ = (moleculeType.pairs[k].atomIndex[1] - 1) + currentNumAtoms;
                ffExclusions.push_back({indexI, indexJ});
                double sigmaI = ffLJ[indexI].sigma;
                double sigmaJ = ffLJ[indexJ].sigma;
                double epsilonI = ffLJ[indexI].epsilon;
                double epsilonJ = ffLJ[indexJ].epsilon;
                double chargeI = ffElec[indexI].charge;
                double chargeJ = ffElec[indexJ].charge;
                double sigmaLJ, epsilonLJ;
                double chargeProd = chargeI * chargeJ * topologies.forceFiled.Coulomb14Scale;
                // Lorentz-Berthelot rule used by AMBER
                if (topologies.forceFiled.combinationRule == 2) {
                    sigmaLJ = 0.5 * (sigmaI + sigmaJ);
                    epsilonLJ = sqrt(epsilonI * epsilonJ) * topologies.forceFiled.LennardJones14Scale;
                }        
                ff14Intra.push_back({indexI, indexJ, sigmaLJ, epsilonLJ, chargeProd});
            }
            // * Add Bonds to calculateBondsEnergy from Bonds.
            for (int k = 0; k < moleculeType.bonds.size(); k++) {
               int indexI = (moleculeType.bonds[k].atomIndex[0] - 1) + currentNumAtoms;
               int indexJ = (moleculeType.bonds[k].atomIndex[1] - 1) + currentNumAtoms;
               double distance = moleculeType.bonds[k].distance;
               double forceConstant = moleculeType.bonds[k].forceConstant;
               int HydrogenTag = moleculeType.bonds[k].HydrogenTag;
               ffBonds.push_back({indexI, indexJ, HydrogenTag, distance, forceConstant});
               ffExclusions.push_back({indexI, indexJ});
            }
            // * Add angles to calculateAngleForce from Angles.
            for (int k = 0; k < moleculeType.angles.size(); k++) {
                int indexI = (moleculeType.angles[k].atomIndex[0] - 1) + currentNumAtoms;
                int indexJ = (moleculeType.angles[k].atomIndex[1] - 1) + currentNumAtoms;
                int indexK = (moleculeType.angles[k].atomIndex[2] - 1) + currentNumAtoms;
                double angle = moleculeType.angles[k].angle;
                double forceConstant = moleculeType.angles[k].forceConstant;
                bool isHAngle = moleculeType.angles[k].isHAngle;
                ffAngles.push_back({indexI, indexJ, indexK, angle, forceConstant, isHAngle});
                ffExclusions.push_back({indexI, indexK});
            }
            // * Add dihedrals to PeriodicTorsionForce from Dihedrals.
            for (int k = 0; k < moleculeType.dihedrals.size(); k++) {
               int indexI = (moleculeType.dihedrals[k].atomIndex[0] - 1) + currentNumAtoms;
               int indexJ = (moleculeType.dihedrals[k].atomIndex[1] - 1) + currentNumAtoms;
               int indexK = (moleculeType.dihedrals[k].atomIndex[2] - 1) + currentNumAtoms;
               int indexL = (moleculeType.dihedrals[k].atomIndex[3] - 1) + currentNumAtoms;
               int periodicity = moleculeType.dihedrals[k].periodicity;
               double dihedral = moleculeType.dihedrals[k].dihedral;
               double forceConstant = moleculeType.dihedrals[k].forceConstant;
               ffDihedrals.push_back({indexI, indexJ, indexK, indexL, periodicity, dihedral, forceConstant});
               ffExclusions.push_back({indexI, indexL});
            }

            // * In order to calculte the currentNumAtoms
            for (int k = 0; k < moleculeType.atoms.size(); k++) {
                currentNumAtoms = currentNumAtoms + 1;   
            }
        }
    } 
}

void ForceFieldBase::calculateBondForceandEnergy(std::vector<Vec3>& forces, std::vector<Vec3>& bondforce, 
                                                double& totalenergy, double& bondenergy, bool includeForces) { 
    bondforce.resize(DOFn, Vec3());
    std::fill(bondforce.begin(), bondforce.end(), Vec3());
    bondenergy = 0;
    double deltaR[5];                                                  
    const int numberOfBonds = ffBonds.size();  

    for (int ii = 0; ii < numberOfBonds; ii++) { 
        int atomAIndex = ffBonds[ii].atomIndex[0];
        int atomBIndex = ffBonds[ii].atomIndex[1];
        const Vec3& atomCoordinatesI = R[atomAIndex];
        const Vec3& atomCoordinatesJ = R[atomBIndex];
        double bonds0 = ffBonds[ii].distance;
        // get deltaR, R2, and R between 2 atoms
        GetDeltaRPeriodic(atomCoordinatesI, atomCoordinatesJ, periodicBoxVectors, deltaR);

        bondenergy += 0.5 * ffBonds[ii].forceConstant * (deltaR[4]- bonds0) * (deltaR[4]- bonds0);  
        //construct the force part
        if (includeForces) {
            double dEdR = ffBonds[ii].forceConstant *(deltaR[4]- bonds0) ;
            dEdR = deltaR[4] > 0 ? (dEdR/deltaR[4]) : 0;
            bondforce[atomAIndex][0] += dEdR*deltaR[0];
            bondforce[atomAIndex][1] += dEdR*deltaR[1];
            bondforce[atomAIndex][2] += dEdR*deltaR[2];
            bondforce[atomBIndex][0] -= dEdR*deltaR[0];
            bondforce[atomBIndex][1] -= dEdR*deltaR[1];
            bondforce[atomBIndex][2] -= dEdR*deltaR[2];
        }
    }   
    if (includeForces) 
        for (int i = 0; i < DOFn; i++) {
            forces[i][0] += bondforce[i][0];
            forces[i][1] += bondforce[i][1];
            forces[i][2] += bondforce[i][2];
        }
    totalenergy += bondenergy;
}

void ForceFieldBase::calculateAngleForceandEnergy(std::vector<Vec3>& forces, std::vector<Vec3>& angleforce,
                                                    double& totalenergy, double& angleenergy, bool includeForces) {  
    angleforce.resize(DOFn, Vec3());
    std::fill(angleforce.begin(), angleforce.end(), Vec3());
    angleenergy = 0;
    double deltaR1[5];   
    double deltaR2[5];                                                    
    const int numberOfAngles = ffAngles.size();   

    for (int ii = 0; ii < numberOfAngles; ii++) {
        int atomAIndex = ffAngles[ii].atomIndex[0];
        int atomBIndex = ffAngles[ii].atomIndex[1];
        int atomCIndex = ffAngles[ii].atomIndex[2];
        const Vec3& atomCoordinatesI = R[atomAIndex];
        const Vec3& atomCoordinatesJ = R[atomBIndex];
        const Vec3& atomCoordinatesK = R[atomCIndex];

        GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR1);
        GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesK, periodicBoxVectors, deltaR2);
        
        double dot = DOT3(deltaR1, deltaR2);;
        double cosine = dot/(deltaR1[4] * deltaR2[4]);
        double angle; 
        if (cosine >= 1.0)
           angle = 0.0;
        else if (cosine <= -1.0)
            angle = pi;
        else
            angle = acos(cosine);
        
        double angle0 = ffAngles[ii].angle * pi/pi_rad;
        // Compute the force and energy, and apply them to the atoms.
        if (includeForces) {
            double Bterm0 = -(deltaR2[0] + deltaR1[0]) /(deltaR1[4] * deltaR2[4]) + cos(angle) * (deltaR1[0]/deltaR1[3] + deltaR2[0]/deltaR2[3]);
            double Bterm1 = -(deltaR2[1] + deltaR1[1]) /(deltaR1[4] * deltaR2[4]) + cos(angle) * (deltaR1[1]/deltaR1[3] + deltaR2[1]/deltaR2[3]);
            double Bterm2 = -(deltaR2[2] + deltaR1[2]) /(deltaR1[4] * deltaR2[4]) + cos(angle) * (deltaR1[2]/deltaR1[3] + deltaR2[2]/deltaR2[3]);
            
            double Aterm0 = deltaR2[0] /(deltaR1[4] * deltaR2[4]) - cos(angle) * (deltaR1[0]/deltaR1[3]);
            double Aterm1 = deltaR2[1] /(deltaR1[4] * deltaR2[4]) - cos(angle) * (deltaR1[1]/deltaR1[3]);
            double Aterm2 = deltaR2[2] /(deltaR1[4] * deltaR2[4]) - cos(angle) * (deltaR1[2]/deltaR1[3]);
            
            double Cterm0 = deltaR1[0] /(deltaR1[4] * deltaR2[4]) - cos(angle) * (deltaR2[0]/deltaR2[3]);
            double Cterm1 = deltaR1[1] /(deltaR1[4] * deltaR2[4]) - cos(angle) * (deltaR2[1]/deltaR2[3]);
            double Cterm2 = deltaR1[2] /(deltaR1[4] * deltaR2[4]) - cos(angle) * (deltaR2[2]/deltaR2[3]);
            
            double dEdR =  ffAngles[ii].forceConstant * (angle - angle0) /sin(angle);

            angleforce[atomBIndex][0] += dEdR * Bterm0;
            angleforce[atomBIndex][1] += dEdR * Bterm1;
            angleforce[atomBIndex][2] += dEdR * Bterm2;
            
            angleforce[atomAIndex][0] += dEdR * Aterm0;
            angleforce[atomAIndex][1] += dEdR * Aterm1;
            angleforce[atomAIndex][2] += dEdR * Aterm2;

            angleforce[atomCIndex][0] += dEdR * Cterm0;
            angleforce[atomCIndex][1] += dEdR * Cterm1;
            angleforce[atomCIndex][2] += dEdR * Cterm2; 
        }
        angleenergy += 0.5 * ffAngles[ii].forceConstant * (angle - angle0) * (angle - angle0);
    }   
    if (includeForces) 
        for (int i = 0; i < DOFn; i++) {
            forces[i][0] += angleforce[i][0];
            forces[i][1] += angleforce[i][1];
            forces[i][2] += angleforce[i][2];
        }
    totalenergy += angleenergy;
}

void ForceFieldBase::calculateDihedralForceandEnergy(std::vector<Vec3>& forces, std::vector<Vec3>& dihedralforce, 
                                                    double& totalenergy, double& dihedralenergy, bool includeForces) {  
    dihedralforce.resize(DOFn, Vec3());
    std::fill(dihedralforce.begin(), dihedralforce.end(), Vec3());
    dihedralenergy = 0;
    double deltaR1[5];  
    double deltaR2[5];   
    double deltaR3[5];                  
    const int numberOfDihedrals = ffDihedrals.size();   

    for (int ii = 0; ii < numberOfDihedrals; ii++) {
        // get deltaR, R2, and R between three pairs of atoms: [i,j], [j,k], [k,l]
        int atomAIndex = ffDihedrals[ii].atomIndex[0];
        int atomBIndex = ffDihedrals[ii].atomIndex[1];
        int atomCIndex = ffDihedrals[ii].atomIndex[2];
        int atomDIndex = ffDihedrals[ii].atomIndex[3];   
        const Vec3& atomCoordinatesI = R[atomAIndex];
        const Vec3& atomCoordinatesJ = R[atomBIndex];
        const Vec3& atomCoordinatesK = R[atomCIndex];
        const Vec3& atomCoordinatesL = R[atomDIndex];

        GetDeltaRPeriodic(atomCoordinatesI, atomCoordinatesJ, periodicBoxVectors, deltaR1);
        GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesK, periodicBoxVectors, deltaR2);
        GetDeltaRPeriodic(atomCoordinatesK, atomCoordinatesL, periodicBoxVectors, deltaR3); 

        double dihedral0 = ffDihedrals[ii].dihedral * pi/pi_rad; 
        double dihedral;

        // compute force
        if (includeForces) {                   
            double crossjkl[4];
            double crossijk[4];
            double crossjk[3] = {deltaR2[0], deltaR2[1], deltaR2[2]};
            double crosskl[3] = {deltaR3[0], deltaR3[1], deltaR3[2]};
            double crossij[3] = {deltaR1[0], deltaR1[1], deltaR1[2]};

            CrossProductVector4(crossjk, crosskl, crossjkl);
            CrossProductVector4(crossij, crossjk, crossijk);

            double crossvec = (crossijk[0] * crossjkl[0] + crossijk[1] * crossjkl[1] + crossijk[2] * crossjkl[2])/(crossijk[3] * crossjkl[3]);
            if (crossvec >= 1.0)
               dihedral = 0.0;
            else if (crossvec <= -1.0)
               dihedral = pi;
            else
               dihedral = acos(crossvec);
            
            double dEdAngle = -ffDihedrals[ii].forceConstant * ffDihedrals[ii].periodicity * sin(ffDihedrals[ii].periodicity * dihedral - dihedral0)/sin(dihedral);       
            double imple1[3];
            double imple2[3];   
            for (int a = 0; a < 3; a++) {
                imple1[a] = crossjkl[a] /crossjkl[3] - cos(dihedral) * crossijk[a]/crossijk[3];
                imple1[a] /= crossijk[3];
                imple2[a] =  crossijk[a] /crossijk[3] - cos(dihedral) * crossjkl[a]/crossjkl[3];
                imple2[a] /= crossjkl[3];     
            } 
            double forceA[3]; 
            double forceB[3];
            double forceC[3];
            double forceD[3];
            double forceE[3];
            double forceF[3];
            double crossijk1[3];
            double crossjkl2[3];

            for (int a = 0; a <3; a++){
                crossijk1[a] = crossij[a] + crossjk[a];
                crossjkl2[a] = crossjk[a] + crosskl[a];
            }

            CrossProductVector3(imple1, crossjk, forceA); 
            CrossProductVector3(imple1, crossijk1, forceB);
            CrossProductVector3(imple2, crosskl, forceC);
            CrossProductVector3(imple2, crossjkl2, forceD);
            CrossProductVector3(imple1, crossij, forceE);
            CrossProductVector3(imple2, crossjk, forceF);        

            // accumulate forces
            for (int i = 0; i < 3; i++) {
                dihedralforce[atomAIndex][i] += forceA[i] * dEdAngle;
                dihedralforce[atomBIndex][i] -= (forceB[i] -forceC[i]) * dEdAngle;
                dihedralforce[atomCIndex][i] -= (forceD[i] -forceE[i]) * dEdAngle;
                dihedralforce[atomDIndex][i] += forceF[i] * dEdAngle;
            }    
        }
        dihedralenergy += ffDihedrals[ii].forceConstant * (1 + cos(ffDihedrals[ii].periodicity * dihedral - dihedral0));
    }   
    if (includeForces) 
        for (int i = 0; i < DOFn; i++) {
            forces[i][0] += dihedralforce[i][0];
            forces[i][1] += dihedralforce[i][1];
            forces[i][2] += dihedralforce[i][2];
        }
    totalenergy += dihedralenergy;
}

void ForceFieldBase::calculateLJForceandEnergy(std::vector<Vec3>& forces, std::vector<Vec3>& LJforce, 
                                                double& totalenergy, double& LJenergy, bool includeForces) {  
    // this part have two parts: 1-4 LJ Interations, long range dispersion correction and LJ Interations
    LJforce.resize(DOFn, Vec3());
    std::fill(LJforce.begin(), LJforce.end(), Vec3());
    LJenergy = 0;
    double deltaR[5];
    double LJ14Energy = 0.0;
    double LJrealEnergy = 0.0;
    double dispersionEnergy = 0.0;
    double newVDWenergy = 0.0;
    double cutlength = param.getDouble("cutoff");
    bool dispersion_correction = param.getBool("dispersion_correction");
    int i, j;

    if (param.getStr("LJ_type") == "LJ12_6") {
        // This part is the LJ Interations
        NeighborList neighborList;
        const int nofE = ffExclusions.size();
        std::vector<std::set<int> > Exclusions(nofE);
        std::vector<std::set<int> > bonded12(nofE); 
        if (nofE != 0)
            addExclusionsToSet(nofE, bonded12, Exclusions);
        const int Numofatom = DOFn;
        //judge to choose the method of neighborlist
        if (DOFn > 500) {// add EnableNeighborList = true
            bool usePeriodic = true;
            double maxDistance = cutlength;
            bool reportSymmetricPairs = false;
            double minDistance = 0.0;  
    
            computeNeighborListVoxelHash(neighborList, Numofatom, R, Exclusions, periodicBoxVectors, 
                        usePeriodic, maxDistance, minDistance , reportSymmetricPairs);
        }
        else // add EnableNeighborList = false
            computeNoNeighborList(neighborList, Numofatom);
    
        for (auto& pair : neighborList) {   
            int ii = pair.first;
            int jj = pair.second;   
            const Vec3& atomCoordinatesI = R[ii];
            const Vec3& atomCoordinatesJ = R[jj];      
            GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
                
            double r         = deltaR[4];
            double inverseR  = 1.0/(deltaR[4]);
            double sig = 0.5 * (ffLJ[ii].sigma +  ffLJ[jj].sigma);
            double sig2 = inverseR * sig;
            sig2 *= sig2; 
            double sig6 = sig2 * sig2 *sig2;
            double eps = sqrt(ffLJ[ii].epsilon * ffLJ[jj].epsilon);      
            if (r < cutlength) {
                if (includeForces) {      
                    double dEdR = 4 * eps * (12 * sig6 - 6) * sig6 * inverseR * inverseR;
                    for (int kk = 0; kk < 3; kk++) {
                        double force = dEdR * deltaR[kk];
                        LJforce[ii][kk] += force;
                        LJforce[jj][kk] -= force;
                    }
                }
            }
            // This part is the long range dispersion correction
            else if (dispersion_correction) {
                double dinverseR  = 1.0/(deltaR[4]);
                double dsig = 0.5 * (ffLJ[ii].sigma +  ffLJ[jj].sigma);
                double V = periodicBoxVectors[0][0] * periodicBoxVectors[1][1] * periodicBoxVectors[2][2]; 
                double dsig2 = dsig * dsig;
                double dsig6 = dsig2 * dsig2 *dsig2;
                double dsig12 = dsig6 * dsig6;
                double deps = sqrt(ffLJ[ii].epsilon * ffLJ[jj].epsilon);
                double dinverseR3 = dinverseR * dinverseR * dinverseR;
                double dinverseR5 = dinverseR3 * dinverseR * dinverseR;
                double dinverseR9 = dinverseR3 * dinverseR3 * dinverseR3;
                double dinverseR11 = dinverseR9 * dinverseR * dinverseR;
            
                if (includeForces) {
                    double dEdR = 8 * pi * DOFn * DOFn *(9 * deps * dsig12 * dinverseR11 - 3 * deps * dsig2 * dinverseR5) / V;
                    for (int kk = 0; kk < 3; kk++) {
                        double force = dEdR * deltaR[kk];
                        LJforce[ii][kk] += force;
                        LJforce[jj][kk] -= force;
                    }
                }
                dispersionEnergy += 8 * pi * DOFn * DOFn *(deps * dsig12 * dinverseR9 - deps * dsig2 * dinverseR3) / V;
            }       
            LJrealEnergy +=   4 * eps * (sig6 - 1) * sig6;      
        }
        // This part is the 1-4 Interations
        for (int i = 0; i < ff14Intra.size(); i++) {    
            int atomAIndex = ff14Intra[i].atomIndex[0];
            int atomBIndex = ff14Intra[i].atomIndex[1];
            const Vec3& atomCoordinatesI = R[atomAIndex];
            const Vec3& atomCoordinatesJ = R[atomBIndex];      
            GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);       
            double inverseR  = 1.0/(deltaR[4]);
            double sig = ff14Intra[i].sigmaLJ;
            double sig2 = inverseR * sig;
            sig2 *= sig2; 
            double sig6 = sig2 * sig2 *sig2;
            double eps = ff14Intra[i].epsilonLJ * 4;    
            if (includeForces) {
                double dEdR =  eps * (12 * sig6 - 6) * sig6 * inverseR * inverseR;
                for (int kk = 0; kk < 3; kk++) {
                    double force = dEdR * deltaR[kk];
                    LJforce[atomAIndex][kk] += force;
                    LJforce[atomBIndex][kk] -= force;
                }
            }       
            // accumulate energies
            LJ14Energy += eps * (sig6 - 1) * sig6;
        }
    }
    else if (param.getStr("LJ_type") == "LJ8_8") {
        // This is for C60 system new vdw model.
        // 2. The model fitted by Domi is :
        // U(r) = a8 * r^8 + a7 * r^7 + a6 * r^6 + a5 * r^5 + a4 * r^4 + a3 * r^3 + a2 * r^2 + a1 * r^1 + a0;
        int atomElec1 = ffnewvdw[1].atomElec;
        std::cout.precision(8);
        for (int i = 0; i < DOFn - 1; i++) 
           for (int j = i + 1; j < DOFn; j++) {
                const Vec3& atomCoordinatesI = R[i];
                const Vec3& atomCoordinatesJ = R[j];      
                GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
                //GetDeltaR(atomCoordinatesJ, atomCoordinatesI, deltaR);
                
                double deltaR3 = deltaR[3] * deltaR[4];
                double deltaR4 = deltaR[3] * deltaR[3];
                double deltaR5 = deltaR3 * deltaR[3];
                double deltaR6 = deltaR3 * deltaR3;
                double deltaR7 = deltaR3 * deltaR4;
                double deltaR8 = deltaR4 * deltaR4; 
                if (i == atomElec1 || j == atomElec1) {
                    newVDWenergy += ffnewvdw[1].vdwparam8 * deltaR8 + ffnewvdw[1].vdwparam7 * deltaR7 + 
                                    ffnewvdw[1].vdwparam6 * deltaR6 + ffnewvdw[1].vdwparam5 * deltaR5 + 
                                    ffnewvdw[1].vdwparam4 * deltaR4 + ffnewvdw[1].vdwparam3 * deltaR3 + 
                                    ffnewvdw[1].vdwparam2 * deltaR[3] + ffnewvdw[1].vdwparam1 * deltaR[4] + ffnewvdw[1].vdwparam0;
                }
                else {
                    newVDWenergy += ffnewvdw[0].vdwparam8 * deltaR8 + ffnewvdw[0].vdwparam7 * deltaR7 + 
                                    ffnewvdw[0].vdwparam6 * deltaR6 + ffnewvdw[0].vdwparam5 * deltaR5 + 
                                    ffnewvdw[0].vdwparam4 * deltaR4 + ffnewvdw[0].vdwparam3 * deltaR3 + 
                                    ffnewvdw[0].vdwparam2 * deltaR[3] + ffnewvdw[0].vdwparam1 * deltaR[4] + ffnewvdw[0].vdwparam0;
                }
                
                if (includeForces) {
                    for (int kk = 0; kk < 3; kk++) {
                        double force = (8 * ffnewvdw[0].vdwparam8 * deltaR6 + 7 * ffnewvdw[0].vdwparam7 * deltaR5 + 
                                        6 * ffnewvdw[0].vdwparam6 * deltaR4 + 5 * ffnewvdw[0].vdwparam5 * deltaR3 + 
                                        4 * ffnewvdw[0].vdwparam4 * deltaR[3] + 3 * ffnewvdw[0].vdwparam3 * deltaR[4] + 
                                        2 * ffnewvdw[0].vdwparam2 ) * deltaR[kk] + ffnewvdw[0].vdwparam1;
                        LJforce[i][kk] -= force;
                        LJforce[j][kk] += force;
                }
            }  
                         
       }
    }
    else 
           throw std::runtime_error("ERROR: Unsupported LJ_type: " + param.getStr("LJ_type")); 
    if (includeForces)
        for (int i = 0; i < DOFn; i++) {
            forces[i][0] += LJforce[i][0];
            forces[i][1] += LJforce[i][1];
            forces[i][2] += LJforce[i][2];
        }
    LJenergy = LJ14Energy + LJrealEnergy + dispersionEnergy + newVDWenergy;
    totalenergy += LJ14Energy + LJrealEnergy + dispersionEnergy + newVDWenergy;
}

void ForceFieldBase::calculateElecForceandEnergy(std::vector<Vec3>& forces, std::vector<Vec3>& elecforce, 
                                                double& totalenergy, double& elecenergy, bool includeForces) {  

    elecforce.resize(DOFn, Vec3());
    std::fill(elecforce.begin(), elecforce.end(), Vec3());
    elecenergy = 0.0;
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
    while (error >= Ewald_tolerance) {
        kmax[0]++;
        error = kmax[0]*sqrt(periodicBoxVectors[0][0]*alphaEwald)*exp(-pow((pi*kmax[0]/(periodicBoxVectors[0][0]*alphaEwald)), 2))/20;
    }
    error = Ewald_tolerance;
    while (error >= Ewald_tolerance) {
        kmax[1]++;
        error = kmax[1]*sqrt(periodicBoxVectors[1][1]*alphaEwald)*exp(-pow((pi*kmax[1]/(periodicBoxVectors[1][1]*alphaEwald)), 2))/20;
    }
    error = Ewald_tolerance;
    while (error >= Ewald_tolerance) {
        kmax[2]++;
        error = kmax[2]*sqrt(periodicBoxVectors[2][2]*alphaEwald)*exp(-pow((pi*kmax[2]/(periodicBoxVectors[2][2]*alphaEwald)), 2))/20;
    }
    double factorEwald = -1 / (4 * alphaEwald * alphaEwald);
    double SQRT_PI = sqrt(pi);
    double TWO_PI = 2.0 * pi;
    double recipCoeff = ONE_4PI_EPS0 * 4 * pi/(periodicBoxVectors[0][0] * periodicBoxVectors[1][1] * periodicBoxVectors[2][2]) /epsilon;
    double totalSelfEwaldEnergy = 0.0;
    double realSpaceEwaldEnergy = 0.0;
    double recipEnergy = 0.0;
    double recipDispersionEnergy = 0.0;
    double totalRecipEnergy = 0.0;

    double V = periodicBoxVectors[0][0] * periodicBoxVectors[1][1] * periodicBoxVectors[2][2];
    int m, n, l, kk;
    double KSQ, AK;
    /**
     * SELF ENERGY
     * totalSelfEwaldEnergy = (alpha_Ew / sqrpi) * Sum_i=0_i=N(q[i] * q[i]);
     */    
    for (int ii = 0; ii < DOFn; ii++) {  
        double charge = ffElec[ii].charge;
        totalSelfEwaldEnergy -=  alphaEwald * charge * charge * ONE_4PI_EPS0 / sqrpi;
    } 
    elecenergy += totalSelfEwaldEnergy;
    /** 
     * RECIPROCAL SPACE EWALD ENERGY AND FORCES
     * totalRecipEnergy = 1 / 2V  * 
     *                   Sum( 4 * pi /k^2 * exp(-k * k / 4 * alpha_Ew * alpha_Ew) * 
     *                   (sum_i=1_i=N(q[i] * exp(i * K * R[i]))));
     *  Erec = Sum{q[j] / V * 
     *        Sum( 4 * pi * K /k * k * exp(-k * k / 4 * alpha_Ew * alpha_Ew) * sin(K * Rij))};
     *  Frec = q[i] * Sum{q[j] / V
     *        * Sum( 4 * pi * K /k * k * exp(-k * k / 4 * alpha_Ew * alpha_Ew) * sin(K * Rij))};
     */
    // to get more quick speed to calculate the energy and force
    typedef std::complex<double> d_complex;
    double recipBoxSize[3] = { TWO_PI / periodicBoxVectors[0][0], TWO_PI / periodicBoxVectors[1][1], TWO_PI / periodicBoxVectors[2][2]};
    double totalRecipEnergy1 = 0;
    #define EIR(x, y, z) eir[(x)*DOFn*3+(y)*3+z]
        std::vector<d_complex> eir(kmax[0]*DOFn*3);
        std::vector<d_complex> tab_xy(DOFn);
        std::vector<d_complex> tab_qxyz(DOFn);
        if (kmax[0] < 1)
            throw std::runtime_error("ERROR: kmax for Ewald summation < 1 ");

        for (int i = 0; (i < DOFn); i++) {
            for (int m = 0; (m < 3); m++)
                EIR(0, i, m) = d_complex(1,0);

            for (int m = 0; (m < 3); m++)
                EIR(1, i, m) = d_complex(cos(R[i][m] * recipBoxSize[m]),
                                         sin(R[i][m] * recipBoxSize[m]));

            for (int j = 2; (j < kmax [0]); j++)
                for (int m = 0; (m < 3); m++)
                    EIR(j, i, m) = EIR(j - 1, i, m) * EIR(1, i, m);
        }

        // calculate reciprocal space energy and forces

        int lowry = 0;
        int lowrz = 1;

        for (int rx = 0; rx < kmax[0]; rx++) {

            double kx = rx * recipBoxSize[0];

            for (int ry = lowry; ry < kmax[1]; ry++) {

                double ky = ry * recipBoxSize[1];

                if (ry >= 0) {
                    for (int n = 0; n < DOFn; n++)
                        tab_xy[n] = EIR(rx, n, 0) * EIR(ry, n, 1);
                }

                else {
                    for (int n = 0; n < DOFn; n++)
                        tab_xy[n]= EIR(rx, n, 0) * std::conj (EIR(-ry, n, 1));
                }

                for (int rz = lowrz; rz < kmax[2]; rz++) {

                    if (rz >= 0) {
                        for (int n = 0; n < DOFn; n++) {
                            double charge = ffElec[n].charge;
                            tab_qxyz[n] = charge * (tab_xy[n] * EIR(rz, n, 2));
                        }
                    }

                    else {
                        for (int n = 0; n < DOFn; n++) {
                            double charge = ffElec[n].charge;
                            tab_qxyz[n] = charge * (tab_xy[n] * std::conj(EIR(-rz, n, 2)));
                    
                        }
                    }

                    double cs = 0.0f;
                    double ss = 0.0f;

                    for (int n = 0; n < DOFn; n++) {
                        cs += tab_qxyz[n].real();
                        ss += tab_qxyz[n].imag();
                    }

                    double kz = rz * recipBoxSize[2];
                    double k2 = kx * kx + ky * ky + kz * kz;
                    double ak = exp(k2*factorEwald) / k2;

                    if (includeForces)
                        for (int n = 0; n < DOFn; n++) {
                            double force = ak * (cs * tab_qxyz[n].imag() - ss * tab_qxyz[n].real());
                            elecforce[n][0] += 2 * recipCoeff * force * kx ;
                            elecforce[n][1] += 2 * recipCoeff * force * ky ;
                            elecforce[n][2] += 2 * recipCoeff * force * kz ;
                        }

                    recipEnergy       = recipCoeff * ak * (cs * cs + ss * ss);
                    totalRecipEnergy += recipEnergy;
                    
                    lowrz = 1 - kmax[2];
                }
                lowry = 1 - kmax[1];
            }
        }   
    elecenergy += totalRecipEnergy;     

    /**
     * Real SPACE EWALD ENERGY AND FORCES
     */
    std::vector<Vec3> realforce;
    realforce.resize(DOFn, Vec3());
    double totalRealSpaceEwaldEnergy = 0.0f;
    NeighborList neighborList;
    const int nofE = ffExclusions.size();
    std::vector<std::set<int> > Exclusions(nofE);
    std::vector<std::set<int> > bonded12(nofE);
    addExclusionsToSet(nofE, bonded12, Exclusions);
   
    bool usePeriodic = true;
    double maxDistance = cutlength;
    bool reportSymmetricPairs = false;
    double minDistance = 0.0;  
    computeNeighborListVoxelHash(neighborList, DOFn, R, Exclusions, periodicBoxVectors, 
                                usePeriodic, maxDistance, minDistance , reportSymmetricPairs);  

    for (auto& pair : neighborList) {
        int ii = pair.first;
        int jj = pair.second;  
        double deltaR1[5];
        const Vec3& atomCoordinatesI = R[ii];
        const Vec3& atomCoordinatesJ = R[jj];
        GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR1);

        double r         = deltaR1[4];
        double inverseR  = 1.0/(deltaR1[4]);
        double alphaR = alphaEwald * r;

        if (includeForces) {
            double dEdR = ONE_4PI_EPS0 * ffElec[ii].charge * ffElec[jj].charge * inverseR * inverseR * inverseR;
            dEdR = dEdR * (erfc(alphaR) + 2 * alphaR * exp (- alphaR * alphaR) / SQRT_PI);
            // accumulate forces
            for (int kk = 0; kk < 3; kk++) {
                double force = dEdR * deltaR1[kk];
                realforce[ii][kk] += force;
                realforce[jj][kk] -= force;
                elecforce[ii][kk] += force;
                elecforce[jj][kk] -= force;
            }
        }
        // accumulate energies
        realSpaceEwaldEnergy = ONE_4PI_EPS0 * ffElec[ii].charge * ffElec[jj].charge * inverseR * erfc(alphaR);     
        totalRealSpaceEwaldEnergy  += realSpaceEwaldEnergy;
    }
    elecenergy +=  totalRealSpaceEwaldEnergy ;
    // Now subtract off the exclusions, since they were implicitly included in the reciprocal space sum.
    double totalExclusionEnergy = 0.0f;
    const double TWO_OVER_SQRT_PI = 2/sqrt(pi);
    for (int i = 0; i < DOFn; i++)
        for (int exclusion : Exclusions[i]) {
            int ii = i;
            int jj = exclusion;
            if (exclusion > i) {
                int ii = i;
                int jj = exclusion;
                const Vec3& atomCoordinatesI = R[ii];
                const Vec3& atomCoordinatesJ = R[jj];
                double deltaR1[5];
                double deltaR2[5];
                GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR1);
                double r = deltaR1[4];
                double inverseR  = 1.0/(deltaR1[4]);
                double alphaR    = alphaEwald * r;
                if (erf(alphaR) > 1e-6) {
                    if (includeForces) {
                        double dEdR = ONE_4PI_EPS0 * ffElec[ii].charge * ffElec[jj].charge * inverseR * inverseR * inverseR;
                        dEdR = dEdR * (erf(alphaR) - 2 * alphaR * exp (- alphaR * alphaR) / SQRT_PI);
                        // accumulate forces
                        for (int kk = 0; kk < 3; kk++) {
                            double force = dEdR * deltaR1[kk];
                            realforce[ii][kk] += force;
                            realforce[jj][kk] -= force;
                            elecforce[ii][kk] -= force;
                            elecforce[jj][kk] += force;
                        }
                    }
                    // accumulate energies
                    realSpaceEwaldEnergy = ONE_4PI_EPS0 * ffElec[ii].charge * ffElec[jj].charge * inverseR * erf(alphaR);
                }
                else {
                    realSpaceEwaldEnergy = alphaEwald*TWO_OVER_SQRT_PI*ONE_4PI_EPS0* ffElec[ii].charge * ffElec[jj].charge;
                }              
                totalExclusionEnergy += realSpaceEwaldEnergy;
            }
        }    
    elecenergy -= totalExclusionEnergy;

    double total14realEnergy;
    double real14SpaceEwaldEnergy;
    for (int i = 0; i < ff14Intra.size(); i++) {
        double deltaR[5];
        int atomAIndex = ff14Intra[i].atomIndex[0];
        int atomBIndex = ff14Intra[i].atomIndex[1];
        const Vec3& atomCoordinatesI = R[atomAIndex];
        const Vec3& atomCoordinatesJ = R[atomBIndex];

        GetDeltaRPeriodic(atomCoordinatesJ, atomCoordinatesI, periodicBoxVectors, deltaR);
        double r         = deltaR[4];
        double inverseR  = 1.0/(deltaR[4]);
        double alphaR    = alphaEwald * r;
        
        if (includeForces) {
            double dEdR = ONE_4PI_EPS0 * ff14Intra[i].chargeProd * inverseR * inverseR*inverseR;
            // accumulate forces
            for (int kk = 0; kk < 3; kk++) {
                double force = dEdR * deltaR[kk];
                realforce[atomAIndex][kk] += force;
                realforce[atomBIndex][kk] -= force;
                elecforce[atomAIndex][kk] += force;
                elecforce[atomBIndex][kk] -= force;
            }
        }
        // accumulate energies
        real14SpaceEwaldEnergy = ONE_4PI_EPS0 * ff14Intra[i].chargeProd  * inverseR;
        total14realEnergy += real14SpaceEwaldEnergy; 
    }  
    elecenergy += total14realEnergy ;   
      
    if (includeForces)
        for (int i = 0; i < DOFn; i++) {
            forces[i][0] += elecforce[i][0];
            forces[i][1] += elecforce[i][1];
            forces[i][2] += elecforce[i][2];
        }
    totalenergy += elecenergy;
}

void ForceFieldBase::addExclusionsToSet(const int nofE, std::vector<std::set<int> >& bonded12, 
                                        std::vector<std::set<int> >& exclusions) {
    for (int i = 0; i < nofE; ++i) {   
        int atomAIndex = ffExclusions[i].atomIndex[0];
        int atomBIndex = ffExclusions[i].atomIndex[1];
        bonded12[atomAIndex].insert(atomBIndex);
        bonded12[atomBIndex].insert(atomAIndex);
    }
    for (int i = 0; i < DOFn; i++)
        for (int bonded : bonded12[i]) {          
            int ii = i;
            int jj = bonded;
            exclusions[ii].insert(jj);
        }    
}

void ForceFieldBase::GetDeltaRPeriodic(const Vec3& atomCoordinatesI, const Vec3& atomCoordinatesJ,
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

void ForceFieldBase::GetDeltaR(const Vec3& atomCoordinatesI, const Vec3& atomCoordinatesJ,
                               double deltaR[5]) {
   deltaR[0]    = atomCoordinatesJ[0] - atomCoordinatesI[0];
   deltaR[1]    = atomCoordinatesJ[1] - atomCoordinatesI[1];
   deltaR[2]    = atomCoordinatesJ[2] - atomCoordinatesI[2];
   deltaR[3]   = DOT3(deltaR, deltaR);
   deltaR[4]    = sqrt(deltaR[3]);
}

void ForceFieldBase::CrossProductVector3(double* vectorX, double* vectorY, double* vectorZ) {
   vectorZ[0]  = vectorX[1]*vectorY[2] - vectorX[2]*vectorY[1];
   vectorZ[1]  = vectorX[2]*vectorY[0] - vectorX[0]*vectorY[2];
   vectorZ[2]  = vectorX[0]*vectorY[1] - vectorX[1]*vectorY[0];
   return;
}

void ForceFieldBase::CrossProductVector4(double* vectorX, double* vectorY, double* vectorZ) {
   vectorZ[0]  = vectorX[1]*vectorY[2] - vectorX[2]*vectorY[1];
   vectorZ[1]  = vectorX[2]*vectorY[0] - vectorX[0]*vectorY[2];
   vectorZ[2]  = vectorX[0]*vectorY[1] - vectorX[1]*vectorY[0];
   vectorZ[3] = sqrt(vectorZ[0] * vectorZ[0] + vectorZ[1]*vectorZ[1]+ vectorZ[2]*vectorZ[2]);
   return;
}

double ForceFieldBase::GetDihedralAngleBetweenThreeVectors(double*  vector1,
                                                            double*  vector2, 
                                                            double*  vector3, 
                                                            double** outputCrossProduct , 
                                                            double*  cosineOfAngle     , 
                                                            double*  signVector        , 
                                                            double*  signOfAngle       ,       
                                                            int          hasREntry = 0) {
   double   tempVectors[6]         = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   // get cross products between vectors and then angle between cross product vectors
   double* crossProduct[2];
   if (outputCrossProduct) {
      crossProduct[0] = outputCrossProduct[0];
      crossProduct[1] = outputCrossProduct[1];
   } else {
      crossProduct[0] = tempVectors;
      crossProduct[1] = tempVectors + 3;
   }
   
   CrossProductVector3(vector1, vector2, crossProduct[0]);
   CrossProductVector3(vector2, vector3, crossProduct[1]);

   double angle = GetAngleBetweenTwoVectors(crossProduct[0], crossProduct[1], cosineOfAngle, 0);
   // take care of sign of angle
   if (signVector) {
      double dotProduct = DOT3(signVector, crossProduct[1]);
      double sign       = dotProduct < 0.0 ? -1.0 : 1.0; 
      if (signOfAngle) {
         *signOfAngle = sign;
      }
      angle *= sign;
   }
   return angle;
}

double ForceFieldBase::GetAngleBetweenTwoVectors(double* vector1, double* vector2, 
                                                double* outputDotProduct = NULL,
                                                int hasREntry = 0) {
    // get dot product betweenn vectors and then angle
   double dotProduct = GetNormedDotProduct(vector1, vector2, hasREntry);
   double angle;
   if (dotProduct > 0.99 || dotProduct < -0.99) {
       // We're close to the singularity in acos(), so take the cross product and use asin() instead.
       double cross[3];
       CrossProductVector3(vector1, vector2, cross);
       double scale = DOT3(vector1, vector1)*DOT3(vector2, vector2);
       angle = asin(sqrt(DOT3(cross, cross)/scale));
       if (dotProduct < 0.0)
           angle = M_PI-angle;
    } 
    else {
      angle = acos(dotProduct);
    }
    if (outputDotProduct) {
        *outputDotProduct = dotProduct;
    }
    return angle;
}

double ForceFieldBase::GetNormedDotProduct(double* vector1, double* vector2,int hasREntry = 0) {
    double dotProduct = DOT3(vector1, vector2);
    if (dotProduct != 0.0) {
        if (hasREntry) {
            dotProduct   /= (vector1[4]*vector2[4]);
        } 
        else {
            double norm1  = DOT3(vector1, vector1);
            double norm2  = DOT3(vector2, vector2);
            dotProduct   /= sqrt(norm1*norm2);
        }
    }   
    // clamp dot product to [-1,1]
    if (dotProduct > 1.0) {
        dotProduct = 1.0;
    } 
    else if (dotProduct < -1.0) {
        dotProduct = -1.0;
    }
    return dotProduct;
}

double ForceFieldBase::getmass(int index) {
    return ffElec[index].mass;         
}

void ForceFieldBase::getElecPerm(std::vector<Vec3>& E_perm) {
    std::vector<Vec3> F_perm;
    double energy_perm;
    calculateElecForceandEnergy(F_perm, F_perm, energy_perm, energy_perm, true);
    for (int i = 0; i < DOFn; i++) {
        for (int j = 0; j < 3; j ++) {
            E_perm[i][j] = F_perm[i][j] / ffElec[i].charge ; 
        }
    }
}

