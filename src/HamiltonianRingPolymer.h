/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: July. 20, 2022                                                *
 * -------------------------------------------------------------------------- */
#pragma once

#include "HamiltonianBase.h"
#include "StructureVec3.h"
#include "Topology.h"
#include "ForceFieldBase.h"
#include "ForceFieldPolar.h"
#include "Vec3.h" 
#include <vector>
#include "Tools.h"

class HamiltonianRingPolymer : public HamiltonianBase{
    friend class DynamicsBase;
    friend class DynamicsRPMD;
   
public:
/**
     * Construct a HamiltonianRingPolymer object.
     *
     * @param param   the global paramters
     */
    HamiltonianRingPolymer(Parameters& param) : HamiltonianBase(param) {}
    ~HamiltonianRingPolymer() {}
    /**
     * Make preparations for F * DOFe .
     */
    void init();   
    /**
     * Update all states potential energy (in kj/mol) and forces.
     */
    void updateAllForces();                                                              
    /**
     * Update Hamiltonian in diabatic representation (potential energies, forces
     * of all states, and averaged forces).
     *
     * It will update data: H, H_old, Heff, Heff_old, F_all, F_avg
     *
     * Note, Hamiltonian (H) denotes the standard one. If the averege potential
     * energy (H_avg) is removed from the diagonal elements of Hamiltonian (H),
     * then it is called effective Hamiltonian (Heff = H - H_avg), which is
     * commonly used in diabatic propagation of some dynamics methods.
     */
    void updateDiabaticHamiltonian();
    /**
     * Update Hamiltonian in adiabatic representation (potential energies, forces
     * of all states, and nonadiabatic coupling (NAC) vectors).
     *
     * It will update data: H, H_old, F_all, NAC
     */
    void updateAdiabaticHamiltonian();
    /**
     * Update Hamiltonian in quasi-diabatic representation.
     *
     * It will update data: H, H_old, F_all
     */
    void updateQuasiDiabaticHamiltonian();
    /**
     * return potential energy.
     *
     * @param index  the index of state
     * @return       the potential energy
     */
    double getPotentialEnergy(int index);
    /**
     * Get kinetic energy. All states share same velocities.
     *
     * @return       the kinetic energy
     */
    double getKineticEnergy();
    /**
     * Compute and return electronic diabatic coupling (in au) i.e., the off-diagonal
     * element of Hamiltonian matrix, H_ij (i != j). Here, assuming the
     * Hamiltonian is a real symmetric matrix, and thus H_ij = H_ji (i != j).
     *
     * If Condon_approximation is true, then the diabatic coupling is a constant.
     * Otherwise, it is R-dependent and thus should be computed based on current
     * positions.
     *
     * @param i  the zero-based state index i
     * @param j  the zero-based state index j, i != j
     * @return   the diabatic coupling in au
     */
    double getDiabaticCoupling(int i, int j);
    void loadInitialStructure(const StructureVec3& structure, const Parameters& param);
    int computeSystemDOF() const;
    void setPeriodicBoxVectors(int index, Vec3& a, Vec3& b, Vec3& c);
    void calculatePotentialEnergyandForces(int index, bool includeForces);
    double getPeriodicBoxVolume(int index) const;
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
     * Get the periodicBoxVectors (measured in nm) of all particles to output.
     *
     */  
    void getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c);
    /**
     * Get the position (in nm) of each particle (read only).
     *
     * @return reference to the R in Hamiltonian object
     */
    void getPositions(std::vector<std::vector<Vec3>>& positions);
    /**
     * Get the velocity (in nm/ps) of each particle (read only).
     *
     * @return reference to the V in Hamiltonian object
     */
    void getVelocities(std::vector<std::vector<Vec3>>& velocities);
    /**
     * Get the force (in kj/mol/nm) of each particle (read only).
     *
     * @return reference to the F in Hamiltonian object
     */
    void getForces(std::vector<std::vector<Vec3>>& forces);
    /**
     * Get the mass of system (in atomic mass unit).
     * 1 atomic mass unit = 1.660539040(20)e-27 kg (CODATA2014).
     *
     * @return  the mass of system (in atomic mass unit)
     */
    double getSystemMass() const;
    /**
     * Get the number of degrees of freedom (DOF) of system.
     *
     * @return  the number of degrees of freedom (DOF) of system
     */
    int getSystemDOF() const;
    /**
     * To get the each part energy of the every state.
     */
    const std::vector<double>& getEnergyGroups(int index);

private:
    int nbeads;
    // FFList stores the force filed parameters loaded from Topology of all state 
    // DOFe. If the topology files are more than one
    std::vector<std::shared_ptr<ForceFieldBase>> FFList;
    std::vector<double> PE;
    double Hspr;
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
    std::vector<std::vector<Vec3>> R_RP, V_RP, F_RP;
    double beta_n, omega_n;
    // Each string in this vector stores the information of atoms, includes resid,
    // resname, atomname and atomindex which is loaded form structure file and
    // will be used when saving trajectory.
    std::vector<std::string> atominfo;
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
    // Integrator can only do a general MD simulation , such as leap-frog Verlet, velocity
    // Verlet, Langevin integrators.
    std::string integrator;
    // box size.
    Vec3 periodicBoxVectors[3];
    double DT;
};
