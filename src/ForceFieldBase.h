/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: Jan. 14, 2021                                                *
 * -------------------------------------------------------------------------- */
#pragma once
#include "Tools.h"
#include "Parameters.h"
#include "Topology.h"
#include "Vec3.h"
#include "ForceFieldNeighborList.h"

class ForceFieldBase {

public: 
    /**
     * Construct a ForceFieldBase object.
     *
     * @param param   the global paramters
     */    
    ForceFieldBase(Parameters& param) : param(param) {}
    virtual ~ForceFieldBase() {}
    /**
     * Make preparations for calculating the forces and energy.
     * 1. Initialize force_contribution to control the calculation of force part 
     * 2. Initialize data members such as R, V, F, forceParts and energyParts.
     */
    virtual void init();
    /**
     * Compute and return potential energy of given state. If includeForces = true,
     * the forces will be updated, and you can use getForces() to retrieve the 
     * forces that were calculated after calling this.
     * @param includeForces  true if forces should be calculated
     * @return               the potential energy
     */
    virtual double getPotentialEnergy(bool includeForces); 
    /**
     * Calculete the induced dipole tensor
     */
    virtual void calculateInducedDipole() {};
    /**
     * Get the Induced dipole tensor from the ForceFieldPolar
     */
    virtual void getMuTensor(Vec3& Mu_tensor) {}; 
    /**
     * Calculete the polarizability tensor
     */ 
    virtual void calculatePiTensor() {}; 
    /**
     * Calculete the dipole size and direction
     */
    virtual void calculateDipoleRleftAright() {}; 
    /**
     * Get the polarizable tensor from the ForceFieldPolar
     */
    virtual void getPiTensor(std::vector<double>& Pi_tensor) {}; 
    virtual void getDipoleRleftAright(double& DEScosleft, double& DEScosright, int& DESnumLeft, 
                                          int& DESnumRight,double& DESdirLx, double& DESdirRx) {}; 
    virtual void setTrueMdPerturb() {};
    virtual void setFalseMdPerturb() {};
    virtual void setPositiveDirection() {};
    virtual void setNegativeDirection() {};
    /**
     * Get the forces F (in kj/mol/nm) from ForceFieldBase and return
     * this is for HamiliationForceFieldBase to call out the forces
     * @return reference to the F in Hamiltonian object [read only]
     */
    const std::vector<Vec3>& getForces();
    /**
     * Get the forceParts (in kj/mol/nm) from ForceFieldBase and return
     * this is for HamiliationForceFieldBase to call out the forceforceParts
     * @return reference to the F in Hamiltonian object [read only]
     */
    const std::vector<std::vector<Vec3>>& getForceParts();
    /**
     * Get the EnergyParts (in kj/mol) from ForceFieldBase and return
     * this is for HamiliationForceFieldBase to call out the EnergyParts
     * @return reference to the F in Hamiltonian object [read only]
     */
    const std::vector<double>& getEnergyParts();    
    /**
     * Get the box Volume of the system.
     */
    double getPeriodicBoxVolume() const;
    /**
     * Set the default values of the vectors defining the axes of the periodic 
     * box (measured in nm).  Any newly created system will have its box vectors 
     * set to these.  They will affect any Force added to the System that uses 
     * periodic boundary conditions.
     *
     * Triclinic boxes are supported, but the vectors must satisfy certain 
     * requirements.  In particular, a must point in the x direction, b must 
     * point "mostly" in the y direction, and c must point "mostly"
     * in the z direction.  See the documentation for details.
     *
     * @param a      the vector defining the first edge of the periodic box
     * @param b      the vector defining the second edge of the periodic box
     * @param c      the vector defining the third edge of the periodic box
     */
    void setPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c);
    /**
     * set the posiation from ForceFieldBase 
     */    
    void setPositions(const std::vector<Vec3>& positions);
    /**
     * set the velocity from ForceFieldBase
     */
    void setVelocities(const std::vector<Vec3>& velocities);
    /**
     * set the Topologies from ForceFieldBase
     */
    void setTopologies(const Topology& fftopologies);
    /**
     * Fill the each part with the force field parameters from topology.
     * @param topology     the Topology object loaded from topology file
     * @param param        the global simulation parameters
     */
    void setForceFieldParameters();
    /**
     * Compute the Bond force and energy components
     * @param forces           total force that to compute
     * @param bondforce        bondforce that to compute
     * @param totalenergy      totalenergy that to compute
     * @param bondenergy       bondenergy that to compute
     * @param includeForces    true if forces should be calculated
     */
    void calculateBondForceandEnergy(std::vector<Vec3>& forces, std::vector<Vec3>& bondforce, 
                                    double& totalenergy, double& bondenergy, bool includeForces);
    /**
     * Compute the Angle force and energy cosmponents
     * @param forces           total force that to compute
     * @param agleforce        angleforce that to compute
     * @param totalenergy      totalenergy that to compute
     * @param angleenergy      angleenergy that to compute
     * @param includeForces    true if forces should be calculated
     */
    void calculateAngleForceandEnergy(std::vector<Vec3>& forces, std::vector<Vec3>& angleforce, 
                                        double& totalenergy, double& angleenergy, bool includeForces); 
    /**
     * Compute the Dihedral force and energy components
     * @param forces           total force that to compute
     * @param dihedralforce    dihedralforce that to compute
     * @param totalenergy      totalenergy that to compute
     * @param dihedralenergy   dihedralenergy that to compute
     * @param includeForces    true if forces should be calculated
     */                            
    void calculateDihedralForceandEnergy(std::vector<Vec3>& forces, std::vector<Vec3>& dihedralforce,
                                        double& totalenergy, double& dihedralenergy, bool includeForces);
    /**
     * Compute the LJ force and energy components
     * @param forces           total force that to compute
     * @param LJforce          LJforce that to compute
     * @param totalenergy      totalenergy that to compute
     * @param LJenergy         LJenergy that to compute
     * @param includeForces    true if forces should be calculated
     */                               
    void calculateLJForceandEnergy(std::vector<Vec3>& forces, std::vector<Vec3>& LJforce, 
                                    double& totalenergy, double& LJenergy, bool includeForces);     
    /**
     * Compute the Elec force and energy components
     * @param forces           total force that to compute
     * @param elecforce        elecforce that to compute
     * @param totalenergy      totalenergy that to compute
     * @param elecenergy       elecenergy that to compute
     * @param includeForces    true if forces should be calculated
     */
    void calculateElecForceandEnergy(std::vector<Vec3>& forces, std::vector<Vec3>& elecforce,
                                     double& totalenergy, double& elecenergy, bool includeForces);

    /**
     * Set the exclusions to calculate the nonbonded the force and energy components
     */
    void addExclusionsToSet(const int NofE, 
                            std::vector<std::set<int>>& bonded12, 
                            std::vector<std::set<int>>& exclusions);
    /**
     * Get the dR in PBC system
     */ 
    void GetDeltaRPeriodic(const Vec3& atomCoordinatesI, const Vec3& atomCoordinatesJ,
                        const Vec3 (&periodicBoxVectors)[3], double deltaR[5]);
    /**
     * Get the dR in NoPBC system
     */ 
    void GetDeltaR(const Vec3& atomCoordinatesI, const Vec3& atomCoordinatesJ, double deltaR[5]);
    /**
     * Get the vectorZ[x] = vectorX[y]*vectorY[z] - vectorX[z]*vectorY[y];
     */ 
    void CrossProductVector3(double* vectorX, double* vectorY, double* vectorZ);
    /**
     * Get the vectorZ[x] = vectorX[y]*vectorY[z] - vectorX[z]*vectorY[y] and
     * the HI
     */ 
    void CrossProductVector4(double* vectorX, double* vectorY, double* vectorZ);    
    /**
     * Get dihedral angle between three vectors
     * @param  vector1            first vector
     * @param  vector2            second vector
     * @param  vector3            third vector
     * @param  outputCrossProduct output cross product vectors
     * @param  cosineOfAngle      cosine of angle (output)
     * @param  signVector         vector to test sign (optional)
     * @param  signOfAngle        sign of angle (output) (optional) 
     * @param  hasREntry          if set, then vector1[ReferenceForce::RIndex] = norm of vector 
     *                            defaults to 0
     * @return cosine of dihedral angle in radians 
     */                                   
    double GetDihedralAngleBetweenThreeVectors(double* vector1, double* vector2, 
                                            double* vector3, double** outputCrossProduct, 
                                            double* cosineOfAngle, double* signVector, 
                                            double* signOfAngle, int hasREntry);                                                                                                          
    /**
     * Get angle between two vectors
     * @param  vector1            first vector
     * @param  vector2            second vector
     * @param  outputDotProduct   output cosine of angle between two vectors (optional)
     * @param  hasREntry          if set, then vector1[ReferenceForce::RIndex] = norm of vector
     *                            defaults to 0 -> R unavailable
     * @return cosine of angles in radians
     */    
    double GetAngleBetweenTwoVectors(double* vector1, double* vector2, 
                                    double* outputDotProduct, int hasREntry);
    /**
     * Get normed dot product between two vectors
     * Do computation in double?
     * @param  vector1            first vector
     * @param  vector2            second vector
     * @param  hasREntry          if set, then vector1[ReferenceForce::RIndex] = norm of vector
     * defaults to 0 (i.e., R unavailable)
     * @return dot product
     */       
    double GetNormedDotProduct(double* vector1, double* vector2, int hasREntry);
    /**
     * Get the mass of atom
     */
    double getmass(int index);   
    /**
     * Get Elec only charge system field
     */
    void getElecPerm(std::vector<Vec3>& E_perm);                                                                                                                                                                                                      
    /**
    * Internal class of MoleculeTypes used to records the parameters of the Harmonic
    * bond stretching (2-body) interaction. The form for Harmonic bond potential
    * is: E = 1/2*k(x-x0)^2, in which k is force constant (in kj/mol/nm^2), and x0
    * is equilibrium distance (in nm). Note that the atom index counts from 1.
    * Here, HydrogenTag records which atom is Hydrogen: 0: no-H; 1: atomI; 2: atomJ.
    */
    struct FFBonds { int atomIndex[2], HydrogenTag; double distance, forceConstant; }; 
    /**
     * Internal class of MoleculeTypes used to records the parameters of the Harmonic
     * bond angle (3-body) interaction. The form for Harmonic angle potential is:
     * E = 1/2*k(c-c0)^2, in which k is force constant (in kj/mol/rad^2), and c0 is
     * equilibrium angle (in degree). Here, isHAngle means H-X-H or H-O-X (where X
     * is an arbitrary atom), which can be constrainted.
     */
    struct FFAngles { int atomIndex[3]; double angle, forceConstant; bool isHAngle; };
    /**
     * Internal class of MoleculeTypes used to records the parameters of periodic type
     * dihedral (4-body) interaction. The form for periodic type dihedral potential is:
     * E = k(1+cos(n*c-c0)), where c is the dihedral angle (in degree) formed by the
     * four particles, c0 is the phase offset, n is the periodicity, and k is the
     * force constant (in kj/mol/rad^2).
     * Note that in AMBER/GAFF proper and improper dihedral angles are defined by
     * the same functional form (periodic type).
     */
    struct FFDihedrals { int atomIndex[4], periodicity; double dihedral, forceConstant; };
    /**
     * Internal class of MoleculeTypes used to records the parameters of 1-4 (2-body)
     * interaction which should be calculated by multiplying the scale factor for Coulomb
     * and Lennard Jones interactions, for example, the scale factors used for Coulomb
     * and Lennard Jones interactions by AMBER/GAFF force filed are 0.5 and 0.833333,
     * respectively. And they should be excluded from the normal nonbonded interactions.
     */ 
    struct FF14Intra{int atomIndex[2]; double sigmaLJ, epsilonLJ, chargeProd;};
    /**
     * Internal class of Topology used to records atomIndex and Lennard Jones interactions 
     * parameters (sigma/nm and epsilon/kJ·mol-1) of atom in a system. sigma is the 
     * distance at which the energy equals zero and epsilon sets the strength of the interaction.
     */
    struct FFLJ{ int atomIndex; double sigma, epsilon; }; 
    /**
     * Internal class of Topology used to records atomIndex and Coulomb interactions parameters 
     * of atom in a system. charge and mass. 
     */ 
    struct FFElec{ int atomIndex; double charge, mass; };
    /**
     * It starts with one atom index, followed by onen or more atom indices. All non-bonded 
     * interactions between the first atom and the other atoms will be excluded for 
     * generating exclusions.
     */
    struct FFExclusions{ int atomIndex[2];};
    /** 
     * This is to get the atom or molecule polarizability
     */
    struct FFNewVDW{ int atomElec; double vdwparam0; double vdwparam1; double vdwparam2; double vdwparam3; 
                    double vdwparam4; double vdwparam5; double vdwparam6; double vdwparam7; double vdwparam8;};    

protected:
    // Parameters object controls the simulation
    Parameters& param;
    // force field parameters of the system
    Topology topologies; 
    std::vector<FFBonds>           ffBonds;
    std::vector<FFAngles>         ffAngles;
    std::vector<FFDihedrals>   ffDihedrals;
    std::vector<FFLJ>                 ffLJ;
    std::vector<FFElec>             ffElec;
    std::vector<FF14Intra>       ff14Intra;
    std::vector<FFExclusions> ffExclusions;
    std::vector<FFNewVDW>         ffnewvdw;
    // box size
    Vec3 periodicBoxVectors[3];
    // nuclear positions (in nm) that to be propagated
    std::vector<Vec3> R;
    // nuclear velocities (in nm/ps) that to be propagated
    std::vector<Vec3> V;
    // nuclear forces (in kj/mol/nm) used for nuclear propagation.
    // Note that the meaning of it may be different, for example, in mapping
    // dynamics, it can be called effective forces which needs F_all, F_avg, and
    // electronic mapping variable q and p to be updated (this updating of F is
    // implemented in DynamicsMapping object).
    // Anyway, it should be the forces which are used to do nuclear propagation.
    std::vector<Vec3> F;
    // [0]: all  [1]: Bond [2]: Angle [3]: Dihedral [4]: LJ [5]: Elec [6]: Polar
    std::vector<bool> force_contribution;
    // This is the number of the all contributions of energy or force
    int num_force_parts = 7;
    // These are every part of the force result: DOFn * 3 * num_force_parts
    std::vector<std::vector<Vec3>> forceParts;
    // These are every part of the energy result: num_force_parts
    std::vector<double> energyParts;
    // this is the total atom in the sysytem
    int DOFn;
};