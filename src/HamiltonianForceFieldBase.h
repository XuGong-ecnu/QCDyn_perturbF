/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: Jan. 15, 2022                                                *
 * -------------------------------------------------------------------------- */
#pragma once

#include "HamiltonianBase.h"
#include "StructureVec3.h"
#include "Topology.h"
#include "ForceFieldBase.h"
#include "HamiltonianOffDiagonal.h"
#include "ForceFieldPolar.h"
#include "Vec3.h" 
#include <vector>
#include "Tools.h"

class HamiltonianForceFieldBase : public HamiltonianBase{
    friend class DynamicsBase;
    friend class DynamicsMD;
    friend class DynamicsReadTraj;
    // This is for test DynamicsLSC for mapping dynamics
    //friend class DynamicsLSC;
    //friend class DynamicsMQCBase;
    //friend class DynamicsMF;
    //friend class DynamicsSQC;
    //friend class DynamicsSPM;
    //friend class DynamicsCMM;
    //friend class DynamicsECMM;
    //friend class DynamicsFSSH;
    //friend class DynamicsTBSH;
    //friend class DynamicsMFRDM;
    //friend class DynamicsECMMCV;
    //friend class DynamicsMixPES;
public:
/**
     * Construct a HamiltonianForceFieldBase object.
     *
     * @param param   the global paramters
     */
    HamiltonianForceFieldBase(Parameters& param) : HamiltonianBase(param) {}
    ~HamiltonianForceFieldBase() {}
    /**
     * Make preparations for an all-atom simulation.
     * 1. Create Stucture object from structure file.
     * 2. Create Topology object from topology file.
     * 3. Create and initialize ForceFieldBase objects.
     * 4. Initialize data members such as masses, R, V, and so on.
     */
    void init(); 
    /**
     * Update the index Pi tensor.
     */
    void updatePiTensor(int index, std::vector<double>& Pi_tensor);
    /**
     * Get the index Pi tensor.
     */
    void getPiTensor(int index, std::vector<double>& Pi_tensor);
    /**
     * Get the index Dipole.
     */
    void getDipoleRleftAright(int index, double& cosleft, double& cosright, 
                            int& numLeft, int& numRight, double& dirLx, double& dirRx);
    /**
     * Update the index Dipole.
     */
    void updateDipoleRleftAright(int index, double& cosleft, double& cosright, 
                                int& numLeft, int& numRight, double& dirLx, double& dirRx);
    /**
     * Update the index Mu tensor.
     */
    void updateMuTensor(int index, Vec3& Mu_tensor);
    /**
     * Get the index Mu tensor.
     */
    void getMuTensor(int index, Vec3& Mu_tensor);
    /**
     * Compute and return potential energy (in kj/mol) of given state.
     *
     * @param index          the index of state, default is 0.
     * @return               the potential energy (in kj/mol)
     */
    void calculatePotentialEnergyandForces(int index, bool includeForces);                                                                 
    /**
     * Return potential energy (in kj/mol) of given state.
     *
     * @param index          the index of state, default is 0.
     * @return               the potential energy (in kj/mol)
     */    
    double getPotentialEnergy(int index);
    /**
     * To get the each part energy of the every state.
     */
    const std::vector<double>& getEnergyGroups(int index);
    /**
     * To get the each part force of the every state.
     */    
    const std::vector<std::vector<Vec3>>& getForceGroups(int index);
    /**
     * Update all states potential energy (in kj/mol) and forces.
     */
    void updateAllForces();
    /**
     * if DOFe >1, update all states EffectiveForces (in kj/mol) and forces.
     */
    void updateEffectiveForces();
    /**
     * TO DO: This is not supported
     * Update Hamiltonian in diabatic representation (potential energies, forces
     * of all states, and averaged forces) when running nonadiabatic dynamics
     * with CustomForceField.
     * It will update data: H, H_old, Heff, Heff_old, F_all, F_avg
     * Note, Hamiltonian (H) denotes the standard one. If the averege potential
     * energy (H_avg) is removed from the diagonal elements of Hamiltonian (H),
     * then it is called effective Hamiltonian (Heff = H - H_avg), which is
     * commonly used in diabatic propagation of some dynamics methods.
     */
    void updateDiabaticHamiltonian();
    /**
     * TO DO: This is not supported
     */
    void updateAdiabaticHamiltonian();
    /**
     * TO DO: This is not supported
     * Update Hamiltonian in quasi-diabatic representation, which is not
     * supported by OpenMM interface.
     */
    void updateQuasiDiabaticHamiltonian();
    /**
     * Get kinetic energy (in kj/mol). All states share same velocities.
     *
     * Note: For leapfrog-like integrator, the velocitiy at each step is delayed
     * haf step, that is V(t-0.5*DT), in this case the kinetic energy should be
     * compueted by shifting: 0.5*m(V(t-0.5*DT)+0.5DT*F/m)^2, not 0.5*mV(t)^2.
     * Thus, the forces based on current positions should be calculated.
     * @return  the kinetic energy (in kj/mol)
     */
    double getKineticEnergy();
    /**
     * Get the index state of the Volue from the ForceFieldBase
     */
    double getPeriodicBoxVolume(int index) const;
    /**
     * Get the position (in nm) of each particle (read only).
     *
     * @return reference to the R in Hamiltonian object
     */
    void getPositions(std::vector<Vec3>& positions);
    /**
     * Set the position (in nm) of each particle (read only).
     *
     * @return reference to the R in Hamiltonian object
     */
    void setPositions(std::vector<Vec3>& positions);
    /**
     * Get the velocity (in nm/ps) of each particle (read only).
     *
     * @return reference to the V in Hamiltonian object
     */
    void getVelocities(std::vector<Vec3>& velocities);
    /**
     * Set the velocity (in nm/ps) of each particle (read only).
     *
     * @return reference to the V in Hamiltonian object
     */
    void setVelocities(std::vector<Vec3>& velocities);
    /**
     * Get the force (in kj/mol/nm) of each particle (read only).
     *
     * @return reference to the F in Hamiltonian object
     */
    void getForces(std::vector<Vec3>& forces);
    /**
     * Get the number of degrees of freedom (DOF) of system.
     *
     * @return  the number of degrees of freedom (DOF) of system
     */
    int getSystemDOF() const;
    /**
     * Compute and return the number of degrees of freedom of system, which is used
     * to compute the instantaneous temperature according to kinetic energy.
     *
     * systemDOF = 3*numAtoms - numConstraints - 3(if removing COM motion)
     * Note that if running a multi-state simulation, the DOF of system is
     * computed from the first one.
     *
     * @return  the number of degrees of freedom of system
     */
    int computeSystemDOF() const;
    /**
     * Set the periodicBoxVectors (measured in nm) of all particles to ForceFieldBase.
     *
     * @param index          the index of state, default is 0.
     */
    void setPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c);
    /**
     * Get the periodicBoxVectors (measured in nm) of all particles to output.
     *
     */  
    void getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c);
    /**
     * Get the mass of system (in atomic mass unit).
     * 1 atomic mass unit = 1.660539040(20)e-27 kg (CODATA2014).
     *
     * @return  the mass of system (in atomic mass unit)
     */
    double getSystemMass() const;
    /**
     * Get the atom information (a string including resid, resname, atomname
     * and atomindex in Gromacs gro format) of each particle (read only),
     * which can be used to output trajectory or structure file, such as
     * Gromacs gro file.
     *
     * @return reference to the atominfo in Hamiltonian object [read only]
     */
    const std::vector<std::string>& getAtomInfo(); 
    /**
     * If the gromacs input file doesn't have the velocities, we can use 
      initial_velocity == "Boltzmann" to set the velocities of each atom.
     */
    void setBoltzmannVelocities();
    /**
     * set the perturb MD simulation
     */  
    void setTrueMdPerturb();
    /**
     * set the perturb MD simulation
     */  
    void setFalseMdPerturb();
    /**
     * set the perturb MD simulation
     */  
    void setPositiveDirection();
    /**
     * set the perturb MD simulation
     */  
    void setNegativeDirection();
    /**
     * save the eq CONFIG
     */ 
    void saveTheEqConfig();
    /**
     * set the eq CONFIG
     */ 
    void setTheEqConfig();
    /**
     * get the direction + or -
     */ 
    bool getDirection() {
        return direction;
    }
    /**
     * The step is the Obs and Ha
     */
    void addTheOneStep();
    int getcurrentStep() {
        return currentStep;
    }
    const std::vector<double>& getPE() const {
        return PE;
    }
private:   
    /**
     * Initialize the positions and velocities of atoms by loading them from
     * structure for the ForceFieldBase in this class.
     *
     * @param structrue  the Structure object loaded from structure file
     * @param param      the global simulation parameters
     */
    void loadInitialStructure(const StructureVec3& structure, const Parameters& param);
  
public:
    //** 
    // The smarter pointer to HamiltonianElec object, which stores the electronic
    // mapping variables of Hamiltonian.
    std::shared_ptr<HamiltonianOffDiagonal> HaOffDiagonal;
    // Each string in this vector stores the information of atoms, includes resid,
    // resname, atomname and atomindex which is loaded form structure file and
    // will be used when saving trajectory.
    std::vector<std::string> atominfo;
    // Potential energy (in kj/mol) of each state, the size of vector is DOFe.
    // The potential energy of each state is computed by function of base class:
    std::vector<double> PE;
    // Get kinetic energy (in kj/mol).
    double KE;
    // unit is ps.
    double DT;
    // Integrator can only do a general MD simulation , such as leap-frog Verlet, velocity
    // Verlet, Langevin integrators.
    std::string integrator;
    // Control the force field polar and nonpolar
    std::string forcefield_type;
    // Control the polarizability calculation
    // masses of each particle (in atomic mass unit)
    // Note that the atomic mass in System will be modified sometimes, e.g.,
    // heavyHydrogenMass or set to 0 to freeze a atom.
    // 1 atomic mass unit = 1.660538921e-27 kg. (Reference: Gromacs manual)
    // Here, inverseMasses is 1.0 / mass of each particle.
    std::vector<double> masses, inverseMasses;
    // total physical mass of system (in atomic mass unit) based on original
    // topology, which is used to compute density of system. Note that if running
    // a multi-state simulation, the system mass is computed from the first one.
    // Don't compute total mass by getting partical masses from OpenMM::System.
    double systemMass;
    // systemDOF is the number of degrees of freedom of system, which is used
    // to compute the instantaneous temperature according to T = 2 * Ek / DOF / R,
    // , R is gas constant, R = 8.3144621e-3 kj/mol/K from Gromacs manual.
    // systemDOF = 3*numAtoms - numConstraints - 3(if removing COM motion)
    // Note that if running a multi-state simulation, the DOF of system is
    // computed from the first one.
    // Note that the usage of systemDOF is is differnent from DOFn. For an all-atom
    // simulation, DOFn is the number of atoms, which is used to propagate positions,
    // velocities of each particle in system.
    double systemDOF;
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
    // box size.
    Vec3 periodicBoxVectors[3];
    // This is to save the each part energy, such as bond, angle, dihedral, vdw, 
    // Elec and polar.
    std::vector<double> energygroups;
    // FFList stores the force filed parameters loaded from Topology of all state 
    // DOFe. If the topology files are more than one
    std::vector<std::shared_ptr<ForceFieldBase>> FFList;
    // The nuclear forces (in kj/mol/nm) matrix of system (stored as vector).
    // Here, the size of F_all is DOFe * DOFe, since it is the negative gradient
    // of Hamiltonian matrix. The digonal forces are forces of each state.
    // When Condon approximation is used, the off-diganol of Hamiltonian is
    // constant and thus the off-dignoal forces is zero.
    std::vector<std::vector<Vec3>> F_all;
    // The average nuclear forces (in kj/mol/nm) of of each atom of all states
    std::vector<Vec3> F_avg; 
    // Temperarure for the intial velocities
    double temperature;
    // random number generator, which is used in nuclear sampling.
    std::mt19937 nucl_gen;
    // This is the eq config
    // nuclear positions (in nm) that to be propagated
    std::vector<Vec3> R_eq;
    // nuclear velocities (in nm/ps) that to be propagated
    std::vector<Vec3> V_eq;
    // nuclear forces (in kj/mol/nm) used for nuclear propagation.
    // Note that the meaning of it may be different, for example, in mapping
    // dynamics, it can be called effective forces which needs F_all, F_avg, and
    // electronic mapping variable q and p to be updated (this updating of F is
    // implemented in DynamicsMapping object).
    // Anyway, it should be the forces which are used to do nuclear propagation.
    std::vector<Vec3> F_eq;
    double currentStep;
    bool direction;
};
