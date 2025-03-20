/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "HamiltonianFMO.h"
#include "HamiltonianMSH.h"
#include "HamiltonianLVC.h"
#include "HamiltonianQVC.h"
#include "HamiltonianOpenMM.h"
#include "HamiltonianSpinBoson.h"
#include "HamiltonianNCMorse.h"
#include "HamiltonianForceFieldBase.h"
#include "HamiltonianRingPolymer.h"

/**
 * Dynamics class defines a method for simulating a system (defined by Hamiltonian
 * class) by integrating the equations of motion of nucleus and/or electron.
 *
 * Note: the dynamics can be on the diabatic or adiabatic representation,
 * which is depended on the system Hamiltonian or dynamics method.
 *
 * Each Dynamics object should be bound to a particular Hamiltonian object which
 * it integrates. This connection is specified by passing the HamiltonianBase as
 * an argument to the constructor of the Dynamics. The data (positions, velocities,
 * and forces) that being propagated are stored in Hamiltonian object (except
 * running a general MD simulation with OpenMM internal Integrator).
 *
 * The nuclear sampling and movement of nuclear coordinates and momenta are
 * defined here, since they are almost same in most dynamics methods
 *
 * DynamicsBase is an abstract class as an interface.
 * The subclasses of it define particular dynamics methods.
 */
class DynamicsBase {
public:
    virtual ~DynamicsBase() {}
    /**
     * Initialize data members and check the legality of parameters.
     */
    virtual void init();
    /**
     * Make preparations for starting to propagate a trajectory.
     *
     * You must call it before running a trajectory.
     */
    virtual void beforeOneTraj() = 0;
    /**
     * Advance a simulation through time by taking a series of time steps.
     *
     * @param steps the number of time steps to take silently
     */
    virtual void dynamics(int steps) = 0;
    /**
     * Make conclusion for a trajectory. (Maybe need for some methods).
     *
     * You must call it after propagating a trajectory.
     */
    virtual void afterOneTraj() = 0;
    /**
     * Get the current step of nuclear propagation.
     */
    int getStep() const {
        return step;
    }
    /**
     * Get the time step size (DT) used for nuclear propagation.
     *
     * Note that the unit of it depends on the system and dynamcis, which
     * may different from the one in input control file.
     * For model system, unit is dimentionless or atomic unit (a.u.).
     * For OpenMM dynamics simulation, unit is ps.
     */
    double getStepSize() const {
        return DT;
    }
    /**
     * Get the current time of nuclear propagation.
     *
     * The unit of it depends on the unit of time step size (DT).
     */
    double getTime() const {
        return DT * step;
    }
    /*
    * Get current on-the-fly calculatation result （SNAPSHOT） of the reduced density matrix elements.
    *    * In DynamicsBase class, this function is served as an interface only, which
    *    * is required to be defined explicitly in nonadiabatic dynamcis method.
    *    * If there is no RDM results for this dynamics method, it will throw an
    *    * exception.
    *
    * @return the constant reference to the sigma (or sigma_*LSC*) in Dynamics object (if any)
    * 
    * Here we return the sigma matrix which is a 2-dimension matrix，sigma[i][j*DOFe+k]
    * whose first dimension is the order of the initial state (states if the inital_state == all),
    * and the second dimension is the order of the time correlation of all the electronic states. 
    * For example, sigma[i][j*DOFe+k] is the time correlation between the i^{th} initial state and 
    * electronic state |j><k| at the given time snapshot. 
    * Nota Bonne: fur LSC method, we have 5 different kinds of RDM come out of one trajectory together, 
    * which should be notificated. 
    */
    virtual const std::vector<std::vector<Complex>>& getRDM() {
        // Please define it in the subclasses of nonadiabatic dynamcis method.
        throw std::runtime_error("ERROR: This dynamics method doesn't have electronic reduced density matrix.");
    }
    virtual const std::vector<std::vector<Complex>>& getRDMLSC1() {
        // Please define it in the subclasses of nonadiabatic dynamcis method (DynamicsLSC).
        throw std::runtime_error("ERROR: This dynamics method doesn't have electronic reduced density matrix.");
    }
    virtual const std::vector<std::vector<Complex>>& getRDMLSC2() {
        // Please define it in the subclasses of nonadiabatic dynamcis method (DynamicsLSC).
        throw std::runtime_error("ERROR: This dynamics method doesn't have electronic reduced density matrix.");
    }
    virtual const std::vector<std::vector<Complex>>& getRDMRILSC1() {
        // Please define it in the subclasses of nonadiabatic dynamcis method (DynamicsLSC).
        throw std::runtime_error("ERROR: This dynamics method doesn't have electronic reduced density matrix.");
    }
    virtual const std::vector<std::vector<Complex>>& getRDMRILSC2() {
        // Please define it in the subclasses of nonadiabatic dynamcis method (DynamicsLSC).
        throw std::runtime_error("ERROR: This dynamics method doesn't have electronic reduced density matrix.");
    }
    virtual const std::vector<std::vector<Complex>>& getRDMRILSC3() {
        // Please define it in the subclasses of nonadiabatic dynamcis method (DynamicsLSC).
        throw std::runtime_error("ERROR: This dynamics method doesn't have electronic reduced density matrix.");
    }
    virtual const std::vector<std::vector<Complex>>& getAdiabaticRDM() {
        // Please define it in the subclasses of nonadiabatic dynamcis method.
        throw std::runtime_error("ERROR: This dynamics method doesn't have electronic reduced density matrix.");
    }
    virtual const std::vector<std::vector<Complex>>& getAdiabaticRDMLSC1() {
        // Please define it in the subclasses of nonadiabatic dynamcis method (DynamicsLSC).
        throw std::runtime_error("ERROR: This dynamics method doesn't have electronic reduced density matrix.");
    }
    virtual const std::vector<std::vector<Complex>>& getAdiabaticRDMLSC2() {
        // Please define it in the subclasses of nonadiabatic dynamcis method (DynamicsLSC).
        throw std::runtime_error("ERROR: This dynamics method doesn't have electronic reduced density matrix.");
    }
    virtual const std::vector<std::vector<Complex>>& getAdiabaticRDMRILSC1() {
        // Please define it in the subclasses of nonadiabatic dynamcis method (DynamicsLSC).
        throw std::runtime_error("ERROR: This dynamics method doesn't have electronic reduced density matrix.");
    }
    virtual const std::vector<std::vector<Complex>>& getAdiabaticRDMRILSC2() {
        // Please define it in the subclasses of nonadiabatic dynamcis method (DynamicsLSC).
        throw std::runtime_error("ERROR: This dynamics method doesn't have electronic reduced density matrix.");
    }
    virtual const std::vector<std::vector<Complex>>& getAdiabaticRDMRILSC3() {
        // Please define it in the subclasses of nonadiabatic dynamcis method (DynamicsLSC).
        throw std::runtime_error("ERROR: This dynamics method doesn't have electronic reduced density matrix.");
    }

protected:
    /**
     * This function is only used in the constructor of subclass, since you can
     * not construct an object for an abstract class.
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     */
    DynamicsBase(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : param(param), Ha(Ha) {}
    /**
     * Get initial sampling of nuclear DOF (positions and velocities).
     * For a simulation of model, Classical or Wigner method can be used.
     * For an all-atom simulation, the Wigner nuclear sampling for all atom needs
     * Hessian, the double precision of Gromacs can do it [not implemented].
     * And the Classical sampling can be done in internal of program via running
     * a general NVE/NVT/NPT MD simulation firstly. Or you can do it yourself
     * externally (nucl_sample = None), then run mapping dynamics one trajectory
     * by one trajectory with one input configuration.
     * TODO: Classical for all-atom in internal
     * TODO: Wigner for all-atom
     */
    void samplingNucl();
    /**
     * Update nuclear velocities (in nm/ps) with a half step.
     * i.e., V(t+0.5*DT) = V(t) + 0.5*DT*F(t)/m.
     * This function is used for all-atom simulation.
     *
     * @param masses the mass of each atom (in atomic mass unit)
     * @param F      the forces (in kj/mol/nm) used to propagation
     * @param V      the velocities (in nm/ps) to be updated
     */
    void MOVE_Half_V(const std::vector<double>& masses, const std::vector<OpenMM::Vec3>& F, std::vector<OpenMM::Vec3>& V);
        /**
     * Update nuclear velocities (in nm/ps) with a half step.
     * i.e., V(t+0.5*DT) = V(t) + 0.5*DT*F(t)/m.
     * This function is used for all-atom simulation.
     *
     * @param masses the mass of each atom (in atomic mass unit)
     * @param F      the forces (in kj/mol/nm) used to propagation
     * @param V      the velocities (in nm/ps) to be updated
     */
    void MOVE_Half_V(const std::vector<double>& masses, const std::vector<Vec3>& F, std::vector<Vec3>& V);
    /**
     * Update nuclear velocities (in nm/ps) with a full step.
     * i.e., V(t+DT) = V(t) + DT*F(t)/m.
     * This function is used for an all-atom simulation.
     *
     * @param masses the mass of each atom (in atomic mass unit)
     * @param F      the forces (in kj/mol/nm) used to propagation
     * @param V      the velocities (in nm/ps) to be updated
     */
    void MOVE_V(const std::vector<double>& masses, const std::vector<OpenMM::Vec3>& F, std::vector<OpenMM::Vec3>& V);
    /**
     * Update nuclear velocities (in nm/ps) with a full step.
     * i.e., V(t+DT) = V(t) + DT*F(t)/m.
     * This function is used for an all-atom simulation.
     *
     * @param masses the mass of each atom (in atomic mass unit)
     * @param F      the forces (in kj/mol/nm) used to propagation
     * @param V      the velocities (in nm/ps) to be updated
     */
    void MOVE_V(const std::vector<double>& masses, const std::vector<Vec3>& F, std::vector<Vec3>& V);
    /**
     * Update nuclear positions (in nm) with a full step.
     * i.e., R(t+DT) = R(t) + DT*V. (Here, V is V(t+DT) or V(t+0.5*DT))
     * This function is used for an all-atom simulation.
     *
     * @param masses the mass of each atom (in atomic mass unit)
     * @param V      the velocities (in nm/ps)
     * @param R      the positions (in nm) to be updated
     */
    void MOVE_R(const std::vector<double>& masses, const std::vector<OpenMM::Vec3>& V, std::vector<OpenMM::Vec3>& R);
    /**
     * Update nuclear positions (in nm) with a full step.
     * i.e., R(t+DT) = R(t) + DT*V. (Here, V is V(t+DT) or V(t+0.5*DT))
     * This function is used for an all-atom simulation.
     *
     * @param masses the mass of each atom (in atomic mass unit)
     * @param V      the velocities (in nm/ps)
     * @param R      the positions (in nm) to be updated
     */
    void MOVE_R(const std::vector<double>& masses, const std::vector<Vec3>& V, std::vector<Vec3>& R);
    /**
     * Update nuclear velocities with a half step.
     * i.e., V(t+0.5*DT) = V(t) + 0.5*DT*F(t)/m.
     * This function is used for a simulation of model (mass=1).
     *
     * @param F      the forces used to propagation
     * @param V      the velocities to be updated
     */
    void MOVE_Half_V(const std::vector<double>& F, std::vector<double>& V);
    /**
     * Update nuclear velocities with a full step.
     * i.e., V(t+DT) = V(t) + DT*V(t).
     * This function is used for a simulation of model (mass=1).
     *
     * @param F      the forces used to propagation
     * @param V      the velocities to be updated
     */
    void MOVE_V(const std::vector<double>& F, std::vector<double>& V);
    /**
     * Update nuclear positions with a full step.
     * i.e., R(t+DT) = R(t) + DT*V. (Here, V is V(t+DT) or V(t+0.5*DT))
     * This function is used for a simulation of model (mass=1).
     *
     * @param V      the velocities
     * @param R      the positions to be updated
     */
    void MOVE_R(const std::vector<double>& V, std::vector<double>& R);
    // This is for RPMD
    void RP_MOVE_Half_V(const std::vector<double>& masses, std::vector<std::vector<Vec3>>& F, std::vector<std::vector<Vec3>>& V);
    void RP_MOVE_V(const std::vector<double>& masses, std::vector<std::vector<Vec3>>& F, std::vector<std::vector<Vec3>>& V);
    
protected:
    // Parameters object controls the simulation
    Parameters& param;
    // Base class of Hamiltonian. We need to cast to the pointer to specific
    // subclass sometimes.
    std::shared_ptr<HamiltonianBase> Ha;
    // The system type: model, onthefly, allatom.
    // model: anlytical model Hamiltonian, SB, MSH, FMO and so on
    // onthefly: electronic structure interface,
    // allatom: calssical force filed MD, OpenMM
    std::string system_type;
    // dyn_type defines the type of dynamics, such as LSC, SQC, OpenMM.
    // means a general MD simultion which will be carried out by OpenMM
    // directly (No Dynamics object is needed); LSC and SQC is mapping dynamics.
    std::string dyn_type;
    // allatom_type defines the type of the allatom, such as OpenMM, QCDyn.
    std::string allatom_type;
    std::string PIMD_type;
    // Equations of motion (EOM) in diabatic, adiabatic, and
    // quasi-diabatic representation.
    std::string representation;
    // total nuclear DOF, read from input for model,
    // for all atom, it is number of atoms
    int DOFn;
    // electronic DOF = number of states/topologies/surfaces.
    // For nonadiabatic dynamics, DOFe should be greater than 1.
    int DOFe;
    int nbeads;
    // the number of trajectories should be propagated.
    // For a general MD simultion, ntraj always equals 1 regardless of input value.
    // For an all-atom mapping dynamics, ntraj should be equal to the number of
    // nuclear sampling (if nucl_sample=none, then ntraj=1).
    int ntraj;
    // total steps for one trajectory
    int nsteps;
    // the nuclear time step size.
    // Note that the unit of it depends on the system and dynamics, which
    // may different from the one in input control file.
    // For model system, unit is reduced unit or atomic unit (a.u.), when real
    // unit (ps) is used in the input control file, convert it to au to propagate.
    // For dynamics simulation with OpenMM, unit is ps.
    double DT;
    // the current step in dynamics, starting from 0 (if it is not a restart job)
    // the current time will be step*DT.
    int step;
    // The nuclear sampling method, allowed values are Wigner, Classical, None,
    // and File (load from file, which is used to debug)
    std::string nucl_sample;
    // real temperature in Kelvin. 1 kT = 3.168429139e-6 a.u.
    double temperature;
    // beta = 1/kT, which is used in nuclear sampling of model.
    // The unit of it is a.u. except the dimentionless SB/GOA model.
    double beta;
    // random number generator, which is used in nuclear sampling.
    std::mt19937 nucl_gen;
};