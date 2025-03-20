/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 14, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "Topology.h"
#include "Structure.h"
#include "HamiltonianBase.h"
#include "ForceFieldBase.h"
#include "ForceFieldPolar.h"      
#include "HamiltonianForceFieldBase.h"

class ForceFieldPolar;
class HamiltonianForceFieldBase;
class ForceFieldBase;
/**
 * An HamiltonianOpenMM class defines the Hamiltonian of an all-atom system (
 * nuclear information only), which includes smarter pointers to OpenMM core
 * objects (System, Integrator, Context) required by an all-atom simulation.
 *
 * This class is an OpenMM interface based on the modified OpenMM C++ library.
 * The charming feature of QCDyn-OpenMM is that it can do multi-state energy
 * calculation when propagation on one specific state and nonadiabatic dynamics
 * simulation for complex condesed phase system with mixed quantum-classical methods.
 * To realize some special functionalties such as multi-state energy calculation,
 * the modification of OpenMM source code may be required.
 *
 * It can be used to do a classical OpenMM MD simultion with being bound to a
 * DynamicsOpenMM object using the internal Integrator of OpenMM.
 * A general MD simulation requires four OpenMM objects: System, Integrator,
 * Platfrom and Context. System stores the force filed parameters loaded from
 * Topology; Integrator implements an algorithm for advancing the simulation
 * through time; Platform defines an implementation of all the kernels needed
 * to perform some calculation. Context stores the complete state of a simulation:
 * positions and velocities, whose initial values are loaded from Structure object.
 *
 * Moreover, it also can be used to do a nonadiabatic dynamics simulation within
 * diabatic representation with being bound to an object of mixed quantum-classical
 * dynamics object, including surface hopping, mean-field Ehrenfest, and mapping
 * dynamics. In this case, the nulcear propagation is still performed with OpenMM
 * integrator in the OpenMM platform, but the forces used to propagation are
 * provided externally. Therefore, the data (such as positions, velocities, and
 * forces) are stored as OpenMM::Vec3 format and OpenMM real units are used for
 * them. While, within Hamiltonian matrix, the atomic unit (au) is used to keep
 * consistense with model Hamiltonian.
 *
 * Note that all the data that being propagated is stored in the OpenMM internal
 * Platform when running simulation with OpenMM. If you want to get positions,
 * you can call getPostions() to get a copy into data memebr.
 */
class HamiltonianOpenMM : public HamiltonianBase {
    friend class DynamicsBase;
    friend class DynamicsOpenMM;
    friend class DynamicsMQCBase;
    friend class DynamicsMF;
    friend class DynamicsLSC;
    friend class DynamicsSQC;
    friend class DynamicsSPM;
    friend class DynamicsCMM;
    friend class DynamicsECMM;
    friend class DynamicsFSSH;
    friend class DynamicsTBSH;
    friend class DynamicsMFRDM;
    friend class DynamicsECMMCV;
    friend class DynamicsMixPES;
    friend class DynamicsReadTraj;


public:
    /**
     * Construct a HamiltonianOpenMM object.
     *
     * @param param   the global paramters
     */
    HamiltonianOpenMM(Parameters& param) : HamiltonianBase(param) {}
    ~HamiltonianOpenMM() {}
    /**
     * Make preparations for an all-atom simulation.
     * 1. Create Stucture object from structure file.
     * 2. Create Topology object from topology file.
     * 3. Create and initialize OpenMM objects (System, Integrator, Platfrom, Context).
     * 4. Initialize data members such as masses, R, V, and so on.
     *
     * Note that we don't hold Stucture and Topology objects in this class, since
     * the data of them have been stored in OpenMM objects or data members.
     * But we must hold three OpenMM objects (System, Integrator, Context) to
     * keep them alive until the termination of program. The Platfrom is a static
     * object, we don't need to control it.
     */
    void init();
    /**
     * Update Hamiltonian in diabatic representation (potential energies, forces
     * of all states, and averaged forces) when running nonadiabatic dynamics
     * with OpenMM interface.
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
     * Update Hamiltonian in adiabatic representation, which is not supported
     * by OpenMM interface.
     */
    void updateAdiabaticHamiltonian();
    /**
     * Update Hamiltonian in quasi-diabatic representation, which is not
     * supported by OpenMM interface.
     */
    void updateQuasiDiabaticHamiltonian();
    /**
     * Compute and return potential energy (in kj/mol) of given state.
     *
     * @param index          the index of state, default is 0.
     * @return               the potential energy (in kj/mol)
     */
    double getPotentialEnergy(int index = 0);
    /**
     * Compute and return potential energy (in kj/mol) of given state. If includeForces=true,
     * the forces in Context will be updated, and you can use getForces() to
     * retrieve the forces that were calculated after calling this.
     *
     * It should be noted that the forces in GPU platform always change as long
     * as the energy calculation is required, but if includeForces is not true,
     * the forces got by getForces() is useless (I don't know what it is).
     * However, for CPU Reference platform, doing energy calculation only
     * without including forces, will not change the forces data in CPU memory
     * (the forces of previous step are restored).
     *
     * @param index          the index of state.
     * @param includeForces  true if forces should be calculated
     * @return               the potential energy (in kj/mol)
     */
    double getPotentialEnergy(int index, bool includeForces);
    /**
     * Compute and return potential energy (in kj/mol) of given force groups. If includeForces=true,
     * the forces in Context will be updated, and you can use getForces() to
     * retrieve the forces that were calculated after calling this.
     *
     * It should be noted that the forces in GPU platform always change as long
     * as the energy calculation is required, but if includeForces is not true,
     * the forces got by getForces() is useless (I don't know what it is).
     * However, for CPU Reference platform, doing energy calculation only
     * without including forces, will not change the forces data in CPU memory
     * (the forces of previous step are restored).
     *
     * @param groups        a set of bit flags for which force groups to include.
     *                      Group i in state j will be included if (groups[j]&(1<<(i-32*j)) != 0.
     * @param includeForces true if forces should be calculated
     * @return              the potential energy (in kj/mol)
     */
    double getPotentialEnergy(const std::vector<int>& groups, bool includeForces);
    /**
     * Get kinetic energy (in kj/mol). All states share same velocities.
     *
     * Note: For leapfrog-like integrator, the velocitiy at each step is delayed
     * haf step, that is V(t-0.5*DT), in this case the kinetic energy should be
     * compueted by shifting: 0.5*m(V(t-0.5*DT)+0.5DT*F/m)^2, not 0.5*mV(t)^2.
     * Thus, the forces based on current positions should be calculated.
     * Whether an integrator need pre-computed forces can be known by
     * integrator->kineticEnergyRequiresForce(). As I know, the Langevin/leapfrog
     * needs forces, while velocityVerlet/LangevinMiddle/NoseHooverChain doesn't
     * need.
     *
     * This funtion won't calculated forces before compute kinetic energy, it
     * just compute kinetic energy according to current velocities and forces in
     * Conetxt, sometime the forces in the Context is not the forces of current
     * positions, therefore, the returned kinetic energy may not right, if the
     * calculation of kinetic energy require the froces based on current positions.
     * In this case, make sure call getPotentialEnergy() with includeForces=true
     * firstly to store correct forces in Context, then call this to compute
     * kinetic energy.
     *
     * @return  the kinetic energy (in kj/mol)
     */
    double getKineticEnergy();
    /**
     * Get whether getKineticEnergy() expects forces to have been pre-computed.
     * This is decided from OpenMM integrator. Generally, non-leapfrog integrators
     * such as velocityVerlet, return false, while leapfrog-like integrators will
     * return true, such as leapfrog and Langevin.
     *
     * @return  true if pre-computed forces is required for getKineticEnergy().
     */
    bool kineticEnergyRequiresForce();
    /**
     * Get the atom information (a string including resid, resname, atomname
     * and atomindex in Gromacs gro format) of each particle (read only),
     * which can be used to output trajectory or structure file, such as
     * Gromacs gro file.
     *
     * @return reference to the atominfo in Hamiltonian object [read only]
     */
    const std::vector<std::string>& getAtomInfo() const;
    /**
     * Get the positions R (in nm) from OpenMM Context and return
     * the reference to R in Hamiltonian object (read only).
     *
     * @return reference to the R in Hamiltonian object [read only]
     */
    const std::vector<OpenMM::Vec3>& getPositions();
    /**
     * Get the velocities V (in nm/ps) from OpenMM Context and return
     * the reference to V in Hamiltonian object (read only).
     *
     * @return reference to the V in Hamiltonian object [read only]
     */
    const std::vector<OpenMM::Vec3>& getVelocities();
    /**
     * Get the forces F (in kj/mol/nm) from OpenMM Context and return
     * the reference to F in Hamiltonian object (read only).
     *
     * You should call getPotentialEnergy() firstly with includeForces=True,
     * then the forces got from OpenMM Context is what you want.
     *
     * @return reference to the F in Hamiltonian object [read only]
     */
    const std::vector<OpenMM::Vec3>& getForces();
    /**
     * Get the position (in nm) of each particle (read only).
     *
     * @return reference to the R in Hamiltonian object
     */
    const std::vector<OpenMM::Vec3>& getPositions() const {
        return R;
    }
    /**
     * Get the velocity (in nm/ps) of each particle (read only).
     *
     * @return reference to the V in Hamiltonian object
     */
    const std::vector<OpenMM::Vec3>& getVelocities() const {
        return V;
    }
    /**
     * Get the force (in kj/mol/nm) of each particle (read only).
     *
     * @return reference to the F in Hamiltonian object
     */
    const std::vector<OpenMM::Vec3>& getForces() const {
        return F;
    }
    /**
     * Get the potential energy of each state (in kj/mol/nm), which stored in
     * the data member PE.
     *
     * It is computed in the function of updateDiabaticHamiltonina(),
     * which is only used for nonadiabatic simulation. For classical MD
     * simulation, it is always zero.
     *
     * @return reference to the PE in Hamiltonian object
     */
    const std::vector<double>& getPE() const {
        return PE;
    }
    /**
     * Get the positions (in nm) from OpenMM Context to a given vector.
     *
     * @param positions the vector to store positions
     */
    void getPositions(std::vector<OpenMM::Vec3>& positions) const;
    /**
     * Get the velocities (in nm/ps) from OpenMM Context to a given vector.
     *
     * @param velocities the vector to store velocities
     */
    void getVelocities(std::vector<OpenMM::Vec3>& velocities) const;
    /**
     * Get the forces (in kj/mol/nm) from OpenMM Context to a given vector.
     *
     * You should call getPotentialEnergy() firstly with includeForces=True,
     * then the forces got from OpenMM Context is what you want.
     *
     * @param forces the vector to store forces
     */
    void getForces(std::vector<OpenMM::Vec3>& forces) const;
    /**
     * Set the positions (measured in nm) of all particles to OpenMM Context.
     *
     * @param positions the vector contains positions
     */
    void setPositions(const std::vector<OpenMM::Vec3>& positions) const;
    /**
     * Set the velocities (measured in nm/picosecond) of all particles to OpenMM Context.
     *
     * @param velocities the vector contains velocities
     */
    void setVelocities(const std::vector<OpenMM::Vec3>& velocities) const;
    /**
     * Set the forces (measured in kj/mol/nm) of all particles to OpenMM Context.
     *
     * @param forces the vector contains forces
     */
    void setForces(const std::vector<OpenMM::Vec3>& forces) const;
    /**
     * Upload the positions (data member R in this object) to OpenMM Context.
     */
    void uploadPositions() const;
    /**
     * Upload the velocities (data member V in this object) to OpenMM Context.
     */
    void uploadVelocities() const;
    /**
     * Upload the forces (data member F in this object) to OpenMM Context.
     */
    void uploadForces() const;
    /**
     * Get the vectors defining the axes of the periodic box (measured in nm)
     * from OpenMM Context.
     *
     * Note: If the PBC is not used in this system, the a, b, c will be 0.
     *
     * @param[out] a  the vector defining the first edge of the periodic box
     * @param[out] b  the vector defining the second edge of the periodic box
     * @param[out] c  the vector defining the third edge of the periodic box
     */
    void getPeriodicBoxVectors(OpenMM::Vec3& a, OpenMM::Vec3& b, OpenMM::Vec3& c) const;
    /**
     * Set the vectors defining the axes of the periodic box (measured in nm)
     * in OpenMM Context.
     *
     * @param[out] a  the vector defining the first edge of the periodic box
     * @param[out] b  the vector defining the second edge of the periodic box
     * @param[out] c  the vector defining the third edge of the periodic box
     */
    void setPeriodicBoxVectors(const OpenMM::Vec3& a, const OpenMM::Vec3& b, const OpenMM::Vec3& c) const;
    /**
     * Returns whether or not this System use periodic boundaries.
     *
     * @return  true if System uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const;
    /**
     * Get the volume (in nm^3) of the periodic box from OpenMM Context.
     *
     * Note: If the PBC is not used in this system, it returns 0.
     *
     * @return  the volume (in nm^3) of the periodic box if System uses PBC
     *          and return 0 otherwise
     */
    double getPeriodicBoxVolume() const;
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
     * Get the number of atoms of system from OpenMM System.
     *
     * @return  the number of atoms of system
     */
    int getNumAtoms() const;
    /**
     * Create a checkpoint recording the current state of the Context. This should
     * be treated as an opaque block of binary data. See loadCheckpoint() for more details.
     *
     * When this function attemps to write each successive checkpoint file, it first
     * renames the previous checkpoint with the suffix _prev, so that even if something goes
     * wrong while writing the new checkpoint file, only recent progress can be lost.
     *
     * @param file  the filename of checkpoint file to save
     */
    void createCheckpoint(const std::string& file);
    /**
     * Load a checkpoint that was written by createCheckpoint().
     *
     * A checkpoint contains not only publicly visible data such as the particle
     * positions and velocities, but also internal data such as the states of
     * random number generators. Ideally, loading a checkpoint should restore
     * the Context to an identical state to when it was written, such that
     * continuing the simulation will produce an identical trajectory.
     * This is not strictly guaranteed to be true, however, and should not be
     * relied on.  For most purposes, however, the internal state should be close
     * enough to be reasonably considered equivalent.
     *
     * A checkpoint contains data that is highly specific to the Context from
     * which it was created. It depends on the details of the System,
     * the Platform being used, and the hardware and software of the computer
     * it was created on.  If you try to load it on a computer with different
     * hardware, or for a System that is different in any way, loading is likely
     * to fail. Checkpoints created with different versions of OpenMM are also
     * often incompatible.
     * If a checkpoint cannot be loaded, that is signaled by throwing an exception.
     *
     * @param file  the filename of checkpoint file to load
     */
    void loadCheckpoint(const std::string& file);
    /**
     * Perform energy minimization for current Context by creating a object of
     * OpenMM::LocalEnergyMinimizer.
     *
     * The required paramters for minimization (converge tolerance of focres and
     * maximum number of iterations) will be got from global parameters.
     *
     * Once finished, the Context will have been updated with the new positions.
     *
     * Note, if a multi-state simulation is performed, then the state to do energy
     * minimization is same as the state to do dynamcis (propagate_state).
     *
     * @param file   the structure file name to save minimized struture,
     *               if empty or set to "none", no file will generate.
     * @param force_tolerance this specifies how precisely the energy minimum
     *                        must be located. Minimization will be halted once
     *                        the root-mean-square value of all force components
     *                        reaches this tolerance. The default value is 10.
     * @param max_iterations  the maximum number of iterations to perform. If
     *                        this is 0, minimation is continued until the results
     *                        converge without regard to how many iterations it
     *                        takes. The default value is 10000.
     */
    void energyMinimize(const std::string& file = "", double force_tolerance = 10.0, int max_iterations = 10000);
    /**
     * This should be called at the start of each time step. It will call the
     * updateContextState() on each ForceImpl in the system, allowing them to
     * modify the values of state variables. For example, an AndersenThermostat
     * can randomize velocities or a MonteCarloBarostat can scale particle
     * positions, or a CMMotionRemover will adjusts the individual particle
     * velocities to make it zero.
     *
     * @return true if the state was modified in any way that would cause the
     *         forces on particles to change, false otherwise
     */
    bool updateContextState();
    /**
     * Get the temperature of the heat bath (in Kelvin).
     *
     * @return the temperature (in K) of heat bath
     */
    double getTemperature() const;
    /**
     * Set the temperature of the heat bath (in Kelvin).
     *
     * @param temperature the temperature (in K) of heat bath
     */
    void setTemperature(double temperature);
    /**
     * Get the current time of the simulation (in picoseconds) in Context.
     *
     * @return the current time in ps
     */
    double getTime() const;
    /**
     * Set the current time of the simulation (in picoseconds) in Context.
     *
     * @param time the current time in ps
     */
    void setTime(double time);
    /**
     * Get the current step in the platform data.
     */
    int getStep() const;
    /**
     * Set the current step by modifying the platform data.
     */
    void setStep(int step);
    /**
     * Get whether the system has distance constraints.
     *
     * @return true is system has distance constraints.
     */
    bool hasConstraints();
    /**
     * Update the positions of particles so that all distance constraints are
     * satisfied.  This also recomputes the locations of all virtual sites.
     *
     * The distance tolerance within which constraints must be satisfied will be
     * get from integrator.
     */
    void applyConstraints();
    /**
     * Update the velocities of particles so the net velocity of each constrained
     * distance is zero.
     *
     * The velocity tolerance within which constraints must be satisfied will be
     * get from integrator.
     */
    void applyVelocityConstraints();
    /**
     * Get whether the system has virtual sites.
     *
     * @return true is system has virtual sites.
     */
    bool hasVirtualSites();
    /**
     * Recompute the locations of all virtual sites.  There is rarely a reason
     * to call this, since virtual sites are also updated by applyConstraints().
     * This is only for the rare situations when you want to enforce virtual
     * sites but not constraints.
     */
    void computeVirtualSites();
    /**
     * When a Context is created, it caches information about the System being simulated
     * and the Force objects contained in it.  This means that, if the System or Forces are then
     * modified, the Context does not see the changes.  Call reinitialize() to force
     * the Context to rebuild its internal representation of the System and pick up any changes
     * that have been made.
     *
     * This is an expensive operation, so you should try to avoid calling it too frequently.
     * Most Force classes have an updateParametersInContext() method that provides a less expensive
     * way of updating certain types of information.  However, this method is the only way to
     * make some types of changes, so it is sometimes necessary to call it.
     *
     * By default, reinitializing a Context causes all state information (positions, velocities,
     * etc.) to be discarded.  You can optionally tell it to try to preserve state information.
     * It does this by internally creating a checkpoint, then reinitializing the Context, then
     * loading the checkpoint.  Be aware that if the System has changed in a way that prevents
     * the checkpoint from being loaded (such as changing the number of particles), this will
     * throw an exception and the state information will be lost.
     */
    void reinitialize(bool preserveState=false);

private:
    /**
     * Initialize the OpenMM::System object in this class by adding the force
     * filed parameters from topologies and periodic box vectors from structure.
     * Before calling this function, you should create a smarter pointer to the
     * System object in heap space firstly: std::make_shared<OpenMM::System>().
     *
     * @param topologies the vector of Topology objects loaded from topology files
     * @param structrue  the Structure object loaded from structure file
     * @param param      the global simulation parameters
     */
    void initializeOpenMMSystem(const std::vector<Topology>& topologies, const Structure& structure, const Parameters& param);
    /**
     * Fill the System and Force Objects with the force field parameters from topology.
     * This is an internal function of createOpenMMSystem().
     *
     * @param topology     the Topology object loaded from topology file
     * @param param        the global simulation parameters
     * @param system       the pointer to System object that to be set
     * @param nonbond      the pointer to NonbondedForce object that to be set
     * @param bondStretch  the pointer to HarmonicBondForce object that to be set
     * @param bondBend     the pointer to HarmonicAngleForce object that to be set
     * @param bondTorsion  the pointer to PeriodicTorsionForce object that to be set
     */
    void setForceFieldParameters(const Topology& topology,
                                 const Parameters& param,
                                 OpenMM::System* system,
                                 OpenMM::NonbondedForce* nonbond,
                                 OpenMM::HarmonicBondForce* bondStretch,
                                 OpenMM::HarmonicAngleForce* bondBend,
                                 OpenMM::PeriodicTorsionForce* bondTorsion);
    /**
     * Set nonbonded method and parameters for a NonbondedForce object.
     * This is an internal function of createOpenMMSystem().
     *
     * @param nonbond the pointer to NonbondedForce object that to be set
     * @param param   the global simulation parameters
     */
    void setNonbondedMethod(OpenMM::NonbondedForce* nonbond, const Parameters& param);
    /**
     * Add a restraint force via CustomExternalForce to current system object.
     * The restraint force is used to apply position restraint with harmonic
     * force to given atoms. The restraint depends on four parameters: the
     * spring force constant k and equilibrium coordinates x0, y0, and z0
     * (reference position).
     *
     * Note: When barostat is used, the restrained atoms must belong to one
     * molecule. And the structure is translated so that the centroid of
     * restrained molecule is in the origin to avoid the barostat rejections.
     *
     * Note: You should call this function before creating Context and after
     * calling initializeOpenMMSystem().
     *
     * Reference:
     * [1] http://docs.openmm.org/latest/api-c++/generated/CustomExternalForce.html
     * [2] https://openmmtools.readthedocs.io/en/0.18.1/_modules/openmmtools/forcefactories.html
     *
     * @param structure the Structure object loaded from structure file,
     *                  the positions of it may be translated if barostat is used.
     * @param param     the global simulation parameters
     */
    void addRestraintForce(Structure& structure, const Parameters& param);
    /**
     * Initialize the OpenMM::Integrator object in this class according the
     * settings in global parameters.
     * Before calling this function, you should create a smarter pointer to the
     * Integrator object in heap space firstly: std::make_shared<OpenMM::Integrator>().
     *
     * @param param   the global simulation parameters
     * @return        the smarter pointer to the Integrator object in heap space
     */
    void initializeOpenMMIntegrator(const Parameters& param);
    /**
     * Create an OpenMM::Platform static object.
     * We don't need to create a smarter pointer to hold the Platform object,
     * since it is a global static object after being created. It always exists
     * before the termination of program.
     *
     * @param param   the global simulation parameters
     * @return        the pointer to the Platform object
     */
    OpenMM::Platform* createOpenMMPlatform(const Parameters& param);
    /**
     * Initialize the positions and velocities of atoms by loading them from
     * structure for the Context object in this class.
     *
     * @param structrue  the Structure object loaded from structure file
     * @param param      the global simulation parameters
     */
    void initializeOpenMMContext(const Structure& structure, const Parameters& param);
    /**
     * Print the Platform and all Platform-specific properties that will be used.
     */
    void printPlatformInfo() const;
    /**
     * Print the parameters being used for (LJ)PME in a particular Context.
     *
     * Because some platforms have restrictions on the allowed grid sizes,
     * the values that are actually used may be slightly different from those
     * specified with setPMEParameters(), or the standard values calculated
     * based on the Ewald error tolerance.
     */
    void printPMEParameters() const;
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
     * Get perturbation forces calculated by ForceFieldPolar
     *  
     * @param forceFieldIndex the index of the ForceField to use (typically 0)
     * @param perturbForces vector to store the calculated perturbation forces
     * @return true if forces were successfully calculated, false otherwise
     **/
bool getForceFieldForces(std::vector<OpenMM::Vec3>& perturbForces);

    /**
     * Get access to the ForceField object at the specified index
     *   
     * @param index the index of the ForceField to retrieve
     * @return pointer to the ForceField object, or nullptr if not found
     **/
//     ForceFieldBase* getForceField(int index = 0);


private:
    std::vector<std::shared_ptr<ForceFieldPolar>> polarForceFields;
    // System is an core object in OpenMM. System stores the force filed
    // parameters loaded from Topology. If the topology files are more than one
    // (a multi-surface or mapping dynamics simulation), we don't create diferrent
    // systems, but create other Force objects which are belong to different
    // force groups to represent the other system or state with bond and nonbonded
    // terms only.
    std::shared_ptr<OpenMM::System> system;
    // Integrator is an core object in OpenMM. Integrator implements an algorithm
    // for advancing the simulation through time, such as leap-frog Verlet, velocity
    // Verlet, Langevin integrators.
    // Note that the internal OpenMM::Integrator can only do a general MD simulation
    // (the forces that are used to propagate is only one state.)
    std::shared_ptr<OpenMM::Integrator> integrator;
    // Context is an core object in OpenMM, which is consist of System, Integrator,
    // Platform and stores the complete state of a simulation: positions and
    // velocities. We use it to compute potential energies and forces in mapping
    // dynamics or to run a general MD simulation.
    std::shared_ptr<OpenMM::Context> context;
    // Each string in this vector stores the information of atoms, includes resid,
    // resname, atomname and atomindex which is loaded form structure file and
    // will be used when saving trajectory.
    std::vector<std::string> atominfo;
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
    // force groups of each state which represents the forces should be added in
    // the calculation of the phsical total potential energy or integration.
    // In my code, each force object will belong to one single force group (
    // a force group can have several force objects). And force group is a
    // interger (>= 0). I have modified the original OpenMM code so that limitation
    // of force group (0<=group<31) dissappers. And in my code, group index 0 is
    // the common share group for all states, such as pseudo-forces: CMMotionRemover,
    // thermostat, barostat, they are no contribution to potential energy, and
    // some custom force such as restriant_force is also included group 0, which
    // has non-zero energy, so we should include group 0 in the energy calculation
    // of each state and propagation. Then, group 1~31 is belong to state 0 (
    // from the first topology file); group 32~63 is belong to the state 1, ....
    // Currently, for AMBER force filed, in each state, group 1 is nonbond,
    // 2 is bond, 3 is angel, 4 is dihedral, and 5 is a fake LennardJonesForce
    // (if any), which is used to do energy decompose only (to get isolated vdW
    // energy), and group 5 should be excluded in the propagation and calculation
    // of total potential energy for one state. Therefore, the gounp index for
    // state i (starts from 0), the group of its nonbond is (1 + i * 32), and
    // bond is (2 + i * 32), angel: (3 + i * 32); dihedral: (4 + i * 32);
    // fake LennardJonesForce (if any): (5 + i * 32). group 6~32 is undefined now.
    // Here, we use a vector of interger to represents which force groups should
    // be included in one specific state (the index of vector is state index).
    // Since the forces from group j will be included if (groups&(1<<j)) != 0,
    // -1 means all force groups will be included, then for example, the 1<<0,
    // the the forces of common group 0 will be included and the value like (1<<0)
    // is the value of element of this vector.
    std::vector<std::vector<int>> forceGroups;
    // Potential energy (in kj/mol) of each state, the size of vector is DOFe.
    // The potential energy of each state is computed by function of base class:
    // HamiltonianOpenMM::getPotentialEnergy(), in which, the function of Context
    // context->getPotentialEnergy() will be called. It is computed in the function
    // of updateDiabaticHamiltonina(), which is only used for nonadiabatic simulation.
    // For classical MD simulation, it is always zero.
    std::vector<double> PE;
    // masses of each particle (in atomic mass unit)
    // Note that the atomic mass in OpeMM::System will be modified sometimes, e.g.,
    // heavyHydrogenMass or set to 0 to freeze a atom.
    // 1 atomic mass unit = 1.660538921e-27 kg. (Reference: Gromacs manual)
    // Here, inverseMasses is 1.0 / mass of each particle.
    std::vector<double> masses, inverseMasses;
    // nuclear positions (in nm) that to be propagated
    std::vector<OpenMM::Vec3> R;
    // nuclear velocities (in nm/ps) that to be propagated
    std::vector<OpenMM::Vec3> V;
    // nuclear forces (in kj/mol/nm) used for nuclear propagation.
    // Note that the meaning of it may be different, for example, in mapping
    // dynamics, it can be called effective forces which needs F_all, F_avg, and
    // electronic mapping variable q and p to be updated (this updating of F is
    // implemented in DynamicsMapping object).
    // Anyway, it should be the forces which are used to do nuclear propagation.
    std::vector<OpenMM::Vec3> F;
    // The nuclear forces (in kj/mol/nm) matrix of system (stored as vector).
    // Here, the size of F_all is DOFe * DOFe, since it is the negative gradient
    // of Hamiltonian matrix. The digonal forces are forces of each state.
    // When Condon approximation is used, the off-diganol of Hamiltonian is
    // constant and thus the off-dignoal forces is zero.
    std::vector<std::vector<OpenMM::Vec3>> F_all;
    // The average nuclear forces (in kj/mol/nm) of of each atom of all states
    std::vector<OpenMM::Vec3> F_avg;
};
