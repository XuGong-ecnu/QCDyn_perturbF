/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 14, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "DynamicsMF.h"
#include "DynamicsLSC.h"
#include "DynamicsSQC.h"
#include "DynamicsSPM.h"
#include "DynamicsCMM.h"
#include "DynamicsTTM.h"
#include "DynamicsECMM.h"
#include "DynamicsFSSH.h"
#include "DynamicsTBSH.h"
#include "DynamicsECMMCV.h"
#include "DynamicsMFRDM.h"
#include "DynamicsOpenMM.h"
#include "DynamicsMD.h"
#include "DynamicsRPMD.h"
#include "DynamicsReadTraj.h"
#include "HamiltonianFMO.h"
#include "HamiltonianMSH.h"
#include "HamiltonianLVC.h"
#include "HamiltonianQVC.h"
#include "HamiltonianSpinBoson.h"
#include "HamiltonianForceFieldBase.h"
#include "HamiltonianNCMorse.h"
#include "Observables.h"
#include "Vec3.h"

/**
 * This opaque class is used to init/run/report a simulation, which includes
 * smarter pointers to the interface of Hamiltonian and Dynamics objects.
 *
 * Some I/O is also defined here, such as when and how to output simulation
 * information (such as energies, temperature, and so on) and trajectory.
 */
class Simulation {
public:
    /**
     * Construct a Simulation object.
     *
     * @param param the global simulation parameters
     */
    Simulation(Parameters& param) : param(param) {}
    ~Simulation() {}
    /**
     * Make preparations for a simulation:
     * 1. Load parameters from input control file and command line.
     * 2. Create and initialize Hamiltonian and Dynamics objects.
     *
     * @param argc the number of arguments from command line
     * @param argv the arguments from command line, argc[0] is program itself
     */
    void init(int argc, char *argv[]);
    /**
     * Run a dynamics simulation.
     */
    void run();
    /**
     * Report current simulation information, such as progress, energies,
     * temperature, positions, velocities, and so on.
     *
     * @param traj  the current traj, 1, ..., ntraj
     * @param step  the current step, 0, 1, ..., nsteps
     */
    void report(int traj, int step);

private:
    /**
     * Run classical OpenMM MD simulation.
     *
     * @param traj      the current traj, 1, ..., ntraj
     * @param nsteps    the total number of steps to run
     */
    void runOpenMM(int traj, int nsteps);
    /**
     * Run CustomForceField MD simulation.
     *
     * @param traj      the current traj, 1, ..., ntraj
     * @param nsteps    the total number of steps to run
     */
    void runPerturbOpenMM(int traj, int nsteps);//add by Xu

    void runMD(int traj, int nsteps);
    /**
     * Run CustomForceField MD simulation for perturb MD.
     *
     * @param traj      the current traj, 1, ..., ntraj
     * @param nsteps    the total number of steps to run
     */   
    void runPerturbMD(int traj, int nsteps);
    /**
     * Read traj MD simulation.
     *
     * @param traj      the current traj, 1, ..., ntraj
     * @param nsteps    the total number of steps to run
     */
    void runTrajMD(int traj, int nsteps);
    /**
     * Get the frerquency to call report() to reprot energies, trajectory and
     * so on for simulation with OpenMM base on the global paramters.
     * It will check the validity of the paramters to control the frequency.
     *
     * @return the skip_steps to call report()
     */
    int getOpenMMReportFrequency();
    /**
     * Run simulated annealing when performing classical MD simulation with OpenMM.
     * This is a part of runOpenMM().
     *
     * @param step        the current step [will change]
     * @param skip_steps  the frequency to report simulation status
     */
    void runOpenMMAnnealing(int& step, int skip_steps);
    /**
     * Print current simulation progress to STDOUT.
     *
     * @param traj  the current traj, 1, ..., ntraj
     * @param step  the current step, 0, 1, ..., nsteps
     */
    void reportProgress(int traj, int step);
    /**
     * Save the energies and so on to file when perfroming nonadiabatic
     * simulation for model.
     *
     * @param traj  the current traj, 1, ..., ntraj
     * @param step  the current step, 0, 1, ..., nsteps
     */
    void reportModelData(int traj, int step);
    /**
     * Save energy/structure/trajecory/checkpoint files when perfroming MD
     * simulation with OpenMM.
     *
     * @param traj  the current traj, 1, ..., ntraj
     * @param step  the current step, 0, 1, ..., nsteps
     */
    void reportOpenMMData(int traj, int step);
    /**
     * Save energy/structure/trajecory/checkpoint files when perfroming MD
     * simulation with ForceField 
     *
     * @param traj  the current traj, 1, ..., ntraj
     * @param step  the current step, 0, 1, ..., nsteps
     */
     void reportMDData(int traj, int step);    /**
     * Save energy/structure/trajecory/checkpoint files when perfroming MD
     * simulation with ForceField 
     *
     * @param traj  the current traj, 1, ..., ntraj
     * @param step  the current step, 0, 1, ..., nsteps
     */
    void reportRPMDData(int traj, int step);
    
    /**
     * Write the energy/structure/trajecory/checkpoint files with Custom force field.
     *
     * @param traj  the current traj, 1, ..., ntraj
     * @param step  the current step, 0, 1, ..., nsteps
     */
    void writeMDData(const std::string& file, std::vector<double>& data, int step, double time);
    /**
     * Get the energy/structure/trajecory/checkpoint files with Custom force field.
     *
     * @param traj  the current traj, 1, ..., ntraj
     * @param step  the current step, 0, 1, ..., nsteps
     */
    void getMDData(int traj, std::vector<double>& data, bool includeForces = false, bool energy_decompose = false);
    /**
     * Get the energy/structure/trajecory/checkpoint files with Custom force field.
     *
     * @param traj  the current traj, 1, ..., ntraj
     * @param step  the current step, 0, 1, ..., nsteps
     */
    void getRPMDData(int traj, std::vector<double>& data, bool includeForces = false, bool energy_decompose = false);
    /**
     * Save energy/structure/trajecory/checkpoint files when perfroming ReadTraj
     * simulation with ForceField 
     *
     * @param traj  the current traj, 1, ..., ntraj
     * @param step  the current step, 0, 1, ..., nsteps
     */
    void reportReadTrajData(int traj, int step);
    /**
     * Get the polarizable tensor with Custom force field.
     *
     * @param traj  the current traj, 1, ..., ntraj
     * @param data
     * @param step  the current step, 0, 1, ..., nsteps
     * @param time
     */
    void writePiState(const std::string& file, std::vector<double>& data, int step, double time);
    /**
     * Get the polarizable tensor with Custom force field.
     *
     * @param traj  the current traj, 1, ..., ntraj
     * @param date  
     */
    void getPiState(int traj, std::vector<double>& data);
    /**
     * Get the frerquency to call report() to reprot energies, trajectory and
     * so on for simulation with OpenMM base on the global paramters.
     * It will check the validity of the paramters to control the frequency.
     *
     * @return the skip_steps to call report()
     */
    int getReportFrequency();
    /**
     * Get current state(s) data from HamiltonianOpenMM object.
     *
     * The data in vector should be (energies in kj/mol):
     * data[0]: Temperature (in K), data[1]: Volume (in nm^3),
     * data[2]: Density (in kg/m^3), data[3]: TotalEnergy, data[4]: KineticEnergy,
     * data[5-12]: total potential energy of state 0, and the energy components
     * of potential energy: Nonbond, Bond, Angle, Dihedral, Others, Coulomb, vdW.
     * Here, the Others means from Custom Force such as potsition restraint.
     * data[13-25]: same as above, but for second state (state 1).
     * data[26-38]: same as above, but for third state (state 2).
     * ...
     *
     * If energy_decompose = true, then energy components of potential energy
     * will be computed, otherwise, they are zero.
     *
     * Note that, the total energy is the sum of kinetic energy and potential
     * energy; Nonbond energy is the sum of Coulomb and vdW; potential energy is
     * the sum of Nonbond, Bond, Angle, Dihedral, and Others
     *
     * @param traj             the current traj, 1, ..., ntraj
     * @param data             the vector used to store the state(s) data
     * @param includeForces    whether to compute force, default is false
     * @param energy_decompose whether to do energy decompose, default is false
     */
    void getStateData(int traj, std::vector<double>& data, bool includeForces = false, bool energy_decompose = false);
    /**
     * Save current state(s) data from HamiltonianOpenMM object to data file.
     *
     * This function is used when perfroming OpenMM dynamics simulation.
     *
     * The data to be written are as following (one line):
     * current step, current time (in ps), state index (from 0).
     * Temperature (in K), Volume (in nm^3), Density (in kg/m^3),
     * TotalEnergy, KineticEnergy, PotentialEnergy,. (in kj/mol)
     * Energy components: Nonbonded_E, Bond_E, Angle_E, Dihedral_E, Others_E,
     * Coulomb_E and Lennard-Jones_E (in kj/mol). [if energy_decompose = true]
     *
     * Note:
     * 1. The total energy is the sum of kinetic energy and potential energy.
     * 2. The potential energy is the sum of Nonbonded_E, Bond_E, Angle_E, Dihedral_E, Others_E.
     * 3. The nonbonded energy is the sum of Coulomb and Lennard-Jones energy
     * 4. If it is a multi-state simulation, the data of states will
     *    printed other lines (one line for one state).
     * 5. For a different state, the values of temperature, volume, density, and
     *    kinetic energy are same (since the propagated trajectory is same one).
     *
     * If an old file with the same name already exists, its contents are discarded.
     *
     * @param file  the filename to save data
     * @param data  the vector that stores the state(s) data
     * @param step  the current MD step
     * @param time  the current MD time in ps
     */
    void writeStateData(const std::string& file, std::vector<double>& data, int step, double time);
    /**
     * Save multi-frame average electronic reduced density matrix to file.
     * The density matrix is stroed as DOFe^2-dimentional vector in RDM.
     *
     * This function is used when perfroming nonadiabatic dynamics simulation.
     *
     * If an old file with the same name already exists, its contents are discarded.
     *
     * @param RDM  the time-dependent reduced density matrix to be outputted
     * @param file the filename to save this RDM.
     */
    void writeTTMFile(const Complex_Matrix& RDM, const std::string& file);

private:
    // Parameters object controls the simulation
    Parameters& param;
    // Hamiltonian object defines the system to be simulated
    std::shared_ptr<HamiltonianBase> Ha;
    // Dynamics object defines a method for simulating the system
    std::shared_ptr<DynamicsBase> Dy;
    // A object contains the keys and points to objects of observables
    std::shared_ptr<Observables> Obs; 
    // The system type: model, onthefly, allatom.
    // model: anlytical model Hamiltonian, SB, MSH, FMO and so on
    // onthefly: electronic structure interface,
    // allatom: calssical force filed MD, OpenMM
    // forcefield: polar, nopolar
    std::string system_type, model_type, onthefly_type, allatom_type;
    // The simulation job type: dynamics, preprocess, postprocess.
    // dyn_type: dynamics method, OpenMM, LSC, SQC, ...
    // prep_type: preprocess type, model_gen
    // post_type: postprocess type, TTM
    std::string job_type, dyn_type, prep_type, post_type, forcefield_type, PIMD_type;
    // Control the polarizability calculation
    bool polar_obs;
    // RDM steps: between how many steps we observe the system
    int RDM_steps;
    // maximum tolerance of the failed trajectory friction
    double max_friction_fail_traj;
    // minumum step beginning to judge if the trajectory is available
    int min_judge_traj_avail;
};
