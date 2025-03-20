/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 14, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "HamiltonianOpenMM.h"

void HamiltonianOpenMM::init() {
//    polarForceFields.resize(DOFe);
//    for (int i = 0; i < DOFe; ++i) {
//        polarForceFields[i] = std::make_shared<ForceFieldPolar>(param);
//        polarForceFields[i]->init();
//    }
    // * 0 Initialize date members of HamiltonianBase.
    HamiltonianBase::init();
    if (system_type != "allatom")
        throw std::runtime_error("ERROR: Unsupported system_type=" + system_type + " for OpenMM Hamiltonian.");
    if (allatom_type != "OpenMM")
        throw std::runtime_error("ERROR: Unsupported allatom_type=" + model_type + " for OpenMM Hamiltonian.");
    if (representation != "diabatic")
        throw std::runtime_error("ERROR: Only diabatic representation is allowed for OpenMM Hamiltonian.");
    if (!Condon_approximation)
        throw std::runtime_error("ERROR: Non-Condon simulation is not yet supported for OpenMM Hamiltonian.");

    // * 1 Check if it is a restart simulation. The checkpoint file created by
    // createCheckpoint() contains particle positions and velocities, and random
    // number generators. But it doesn't contain the information of System,
    // Integrator, and Platform. Therefore, it is nessesary to rebuild the System
    // Integrator, and Platform objects using same input file as previous.
    // Then, create the Context and loadCheckpoint to restore the state.
    // It will allow you to exactly continue a simulation if the need arises
    // (though whether the simulation is deterministic depends on platform and
    // methods. However, this binary can only be used to restart simulations on
    // machines with the same hardware and the same OpenMM version as the one
    // that saved it.
    const std::string restart_file = param.getStr("restart");
    const bool isRestart = restart_file == "none" ? false : true;

    // * 2 Create Stucture object from structure file and get atominfo.
    std::cout << "Start to initialize OpenMM objects since " << CurrentTime() << ".\n" << std::endl;
    Structure structure;
    structure.loadStructure(param.getStr("structure"));
    atominfo = structure.getAtomInfo();
    std::cout << "Loaded structure from file: " << param.getStr("structure") << "\n" << std::endl;

    // * 3 Create Topology objects from topology files and get total mass.
    std::vector<std::string> files;
    SplitString(files, param.getStr("topology"), ',');
    // Reset DOFe to the number of topology files (input value is ignored),
    // since number of states (DOFe) must equal to the number of topology files.
    DOFe = files.size();
    std::vector<Topology> topologies(DOFe);
    for (int i = 0; i < DOFe; ++i) {
        topologies[i].loadTopologyFile(files[i]);
        std::cout << "Loaded topology from file: " << files[i] << std::endl;
    }
    std::cout << std::endl;
    systemMass = topologies[0].computeSystemMass();

    // * 4 Create OpenMM core objects and initialize them.
    // 4.1 Create smarter pointer to OpenMM System and initialize, and get systemDOF.
    system = std::make_shared<OpenMM::System>();
    initializeOpenMMSystem(topologies, structure, param);
    systemDOF = computeSystemDOF();
    // 4.1.5 (optional) Apply position restraint with harmonic force
    // Note: call this before creating Context and after initializing System.
    if (param.getBool("position_restraint"))
        addRestraintForce(structure, param);
    // 4.1.6 Precompute the force groups should be included in the intergration
    // and calculation of total potential enenrgy of each state.
    // Note: Do this before initializing intergrator
    for (int state = 0; state < DOFe; ++state) {
        // the index of this vector denotes state, and the value represents
        // the force groups of this state should be included.
        std::vector<int> groups((state+1), 0);
        // group 0 is common group for all states such as CMMotionRemover,
        // barostat, positiona restriant, and should be included always.
        // The forces from group i will be included if (groups&(1<<i)) != 0.
        groups[0] = 1 << 0;
        // group 1-4 is physical forces from topology: nonbond, bond, angel,
        // dihedral of each state. and should be included.
        for (int i = 1; i <= 4; ++i)
            groups[state] += 1 << i;
        forceGroups.push_back(groups);
        // group 5 is is a fake LennardJonesForce (if any), which is used to do
        // energy decompose only (to get isolated vdW energy), and not inculded
        // here. group 6~32 is undefined currently.
    }
    // 4.2 Create smarter pointer to OpenMM Integrator and initialize the object.
    // Integrator is an abstract class, we need create a smarter pointer to
    // subclass firstly and assign it to the class data member: integrator.
    initializeOpenMMIntegrator(param);
    // 4.3 Create OpenMM Platform static object.
    OpenMM::Platform* platform = createOpenMMPlatform(param);
    // 4.3 Create OpenMM Context based on the System, Integrator and Platform objects.
    // If platform is NULL, the Context will choose the best available Platform with default properties.
    context = std::make_shared<OpenMM::Context>(*system, *integrator, *platform);
    // 4.4 Initialize the positions and velocities in the Context
    if (!isRestart)
        initializeOpenMMContext(structure, param);
    // if this is a restart job, the simulation state will load from checkpoint file.
    else {
        loadCheckpoint(restart_file);
        std::cout << "Load simulation state from checkpoint file: " << restart_file << "\n" << std::endl;
    }
    // 4.5 Print infomation (optional).
    printPlatformInfo();
    printPMEParameters();
    std::cout << "The initialization of OpenMM objects is finished at " << CurrentTime() << ".\n" << std::endl;

    // * 5 Reset DOFn. For all atom, the DOFn is number of atoms.
    DOFn = system->getNumParticles();

    // * 6 Get mass of each particle from OpenMM System.
    masses.resize(DOFn, 0);
    inverseMasses.resize(DOFn, 0);
    for (int j = 0; j < DOFn; j++) {
        masses[j] = system->getParticleMass(j);
        inverseMasses[j] = masses[j] == 0 ? 0 : (1.0/masses[j]);
    }

    // * 7 Do energy minimization (if any) before dynamics.
    // This will update the positions in Context directly, so do it after Context
    // is created. And if it is a restart simulation, this will not to do.
    if (param.getBool("energy_minimization") && !isRestart) {
        std::cout << "Start to do energy minimization since " << CurrentTime() << "." << std::endl;
        const double force_tolerance = param.getDouble("force_tolerance");
        const int max_iterations = param.getInt("max_iterations");
        const std::string file = param.getStr("minimized_structure");
        energyMinimize(file, force_tolerance, max_iterations);
        // Note, if a multi-state simulation is performed, then the state to do energy
        // minimization is same as the state to do dynamcis (propagate_state).
        if (DOFe > 1) // multi-state topology files are provided
            std::cout << "The zero-based index of state to do energy minimization is " <<
                param.getInt("propagate_state") <<  "." << std::endl;
        std::cout << "The converge tolerance of force components is " << param.getStr("force_tolerance") <<
            " kj/mol/nm." << std::endl;

        if (max_iterations == 0)
            std::cout << "Since max_iterations=0, the minimation will be continued until the results "
                "converge without regard to how many iterations it takes." << std::endl;
        else
            std::cout << "The maximum number of iterations to perform energy minimization is " <<
                param.getStr("max_iterations") << "." << std::endl;
        if (!file.empty() && file != "none")
            std::cout << "The minimized structure is saved to file: " << file << std::endl;
        std::cout << "The energy minimization is finished at " << CurrentTime() << ".\n" << std::endl;
    }

    // * 8 Resize all vectors in this class and set their elements to 0.
    // We will update them only when we need.
    PE.resize(DOFe, 0);
    R.resize(DOFn, OpenMM::Vec3());
    V.resize(DOFn, OpenMM::Vec3());
    F.resize(DOFn, OpenMM::Vec3());
    F_all.resize(DOFe*DOFe, std::vector<OpenMM::Vec3>(DOFn, OpenMM::Vec3()));
    F_avg.resize(DOFn, OpenMM::Vec3());
    

    polarForceFields.resize(DOFe);
    for (int i = 0; i < DOFe; ++i) {
        polarForceFields[i] = std::make_shared<ForceFieldPolar>(param);
        polarForceFields[i]->init();
    }

}

void HamiltonianOpenMM::updateDiabaticHamiltonian() {
    // This function is very simillar as HamiltonianModelBase::updateDiabaticHamiltonian().
    // Save current H as H_old before update
    H_old = H;
    Heff_old = Heff;
    // Remember to reset H_avg to 0.
    // the average potential energy of all states, which is removed in
    // effective energy matrix (Heff)
    double H_avg = 0.0;
    for (int i = 0; i < DOFe; i++) {
        // Compute potential energy (in kj/mol) each state
        PE[i] = getPotentialEnergy(i, true);
        // The total potential energy should add energy correction
        // Note the unit of Hamiltonian matrix is au.
        H[i][i] = PE[i] * kj2au  + epsilon[i];
        // Compute H_avg
        H_avg += H[i][i] / DOFe;
        // Get forces of current state, which is diganol element of force matrix
        getForces(F_all[i * DOFe + i]);
        // Get the off-diagonal element of Hamiltonian matrix and forces
        // When using Condon approximation, coupling is constant and force is 0.
        // So, we only need to compute them in the non-Condon case.
        // Here, the Hamiltonian is a real symmetrix matrix, i.e., H_ij = H_ji
        // and F[i][j] = F[j][i] = - dH_ij/dR, i != j.
        // Note, the unit of force is kj/mol/nm.
        // TODO: undefined for all-atom
        if (!Condon_approximation)
            for (int j = i + 1; j < DOFe; ++j) {
                H[j][i] = H[i][j] = getDiabaticCoupling(i, j);
                //getForces(i, j, F_all[i * DOFe + j]); // undefined
                F_all[j * DOFe + i] = F_all[i * DOFe + j];
            }
    }
    // Get effective Hamiltioan matrix (Heff) with removing H_avg
    Heff = H;
    for (int i = 0; i < DOFe; i++)
        Heff[i][i] -= H_avg;
    // Reset all elements to 0, then compute average forces (F_avg)
    std::fill(F_avg.begin(), F_avg.end(), OpenMM::Vec3());
    for (int j = 0; j < DOFn; j++) {
        for (int i = 0; i < DOFe; i++)
            F_avg[j] += F_all[i * DOFe + i][j];
        F_avg[j] /= DOFe;
    }
}

void HamiltonianOpenMM::updateAdiabaticHamiltonian() {
    throw std::runtime_error("ERROR: Nonadiabatic dynamics simulation within "
        "adiabatic representation is not supported by OpenMM interface.");
}

void HamiltonianOpenMM::updateQuasiDiabaticHamiltonian() {
    throw std::runtime_error("ERROR: Nonadiabatic dynamics simulation within "
        "quasi-diabatic representation is not supported by OpenMM interface.");
}

double HamiltonianOpenMM::getPotentialEnergy(int index) {
    if (index >= 0 && index < DOFe)
        return context->getPotentialEnergy(false, forceGroups[index]);
    else
        throw std::runtime_error("ERROR: Called getPotentialEnergy() with wrong state index.");
}

double HamiltonianOpenMM::getPotentialEnergy(int index, bool includeForces) {
    if (index >= 0 && index < DOFe)
        return context->getPotentialEnergy(includeForces, forceGroups[index]);
    else
        throw std::runtime_error("ERROR: Called getPotentialEnergy() with wrong state index.");
}

double HamiltonianOpenMM::getPotentialEnergy(const std::vector<int>& groups, bool includeForces) {
    if (groups.size() > 0 && groups.size() <= DOFe)
        return context->getPotentialEnergy(includeForces, groups);
    else
        throw std::runtime_error("ERROR: Called getPotentialEnergy() with wrong force groups.");
}

double HamiltonianOpenMM::getKineticEnergy() {
    return context->getKineticEnergy();
}

bool HamiltonianOpenMM::kineticEnergyRequiresForce() {
    return context->kineticEnergyRequiresForce();
}

const std::vector<OpenMM::Vec3>& HamiltonianOpenMM::getPositions() {
    static const bool enforce_box = param.getBool("enforce_box");
    static const bool usePBC = usesPeriodicBoundaryConditions();
    context->getPositions(R, usePBC && enforce_box);
    return R;
}

const std::vector<OpenMM::Vec3>& HamiltonianOpenMM::getVelocities() {
    context->getVelocities(V);
    return V;
}

const std::vector<OpenMM::Vec3>& HamiltonianOpenMM::getForces() {
    context->getForces(F);
    return F;
}

void HamiltonianOpenMM::getPositions(std::vector<OpenMM::Vec3>& positions) const {
    static const bool enforce_box = param.getBool("enforce_box");
    static const bool usePBC = usesPeriodicBoundaryConditions();
    context->getPositions(positions, usePBC && enforce_box);
}

void HamiltonianOpenMM::getVelocities(std::vector<OpenMM::Vec3>& velocities) const {
    context->getVelocities(velocities);
}

void HamiltonianOpenMM::getForces(std::vector<OpenMM::Vec3>& forces) const {
    context->getForces(forces);
}

void HamiltonianOpenMM::setPositions(const std::vector<OpenMM::Vec3>& positions) const {
    context->setPositions(positions);
}

void HamiltonianOpenMM::setVelocities(const std::vector<OpenMM::Vec3>& velocities) const {
    context->setVelocities(velocities);
}

void HamiltonianOpenMM::setForces(const std::vector<OpenMM::Vec3>& forces) const {
    context->setForces(forces);
}

void HamiltonianOpenMM::uploadPositions() const {
    context->setPositions(this->R);
}

void HamiltonianOpenMM::uploadVelocities() const {
    context->setVelocities(this->V);
}

void HamiltonianOpenMM::uploadForces() const {
    context->setForces(this->F);
}

void HamiltonianOpenMM::getPeriodicBoxVectors(OpenMM::Vec3& a, OpenMM::Vec3& b, OpenMM::Vec3& c) const {
    if (usesPeriodicBoundaryConditions())
        context->getPeriodicBoxVectors(a, b, c);
    else
        a = b = c = OpenMM::Vec3(0, 0, 0);
}

void HamiltonianOpenMM::setPeriodicBoxVectors(const OpenMM::Vec3& a, const OpenMM::Vec3& b, const OpenMM::Vec3& c) const {
    if (usesPeriodicBoundaryConditions())
        context->setPeriodicBoxVectors(a, b, c);
    else
        throw std::runtime_error("ERROR: Called setPeriodicBoxVectors() for non-PBC simulation.");
}

bool HamiltonianOpenMM::usesPeriodicBoundaryConditions() const {
    return system->usesPeriodicBoundaryConditions();
}

double HamiltonianOpenMM::getPeriodicBoxVolume() const {
    if (usesPeriodicBoundaryConditions())
        return context->getPeriodicBoxVolume();
    else
        return 0.0;
}

double HamiltonianOpenMM::getSystemMass() const {
    return systemMass;
}

int HamiltonianOpenMM::getSystemDOF() const {
    return systemDOF;
}

int HamiltonianOpenMM::getNumAtoms() const {
    return system->getNumParticles();
}

const std::vector<std::string>& HamiltonianOpenMM::getAtomInfo() const {
    return atominfo;
}

bool HamiltonianOpenMM::updateContextState() {
    return context->updateContextState();
}

double HamiltonianOpenMM::getTemperature() const {
    static const std::string thermostat = param.getStr("thermostat");
    if (thermostat == "none")
        throw std::runtime_error("ERROR: Called getTemperature() without using thermostat.");
    if (thermostat == "Andersen") {
        return context->getParameter("AndersenTemperature");
    }
    else if (thermostat == "Langevin") {
        static std::shared_ptr<OpenMM::LangevinIntegrator> Langevin =
            std::static_pointer_cast<OpenMM::LangevinIntegrator>(integrator);
        return Langevin->getTemperature();
    }
    else if (thermostat == "LangevinMiddle") {
        static std::shared_ptr<OpenMM::LangevinMiddleIntegrator> LangevinMiddle =
            std::static_pointer_cast<OpenMM::LangevinMiddleIntegrator>(integrator);
        return LangevinMiddle->getTemperature();
    }
    else if (thermostat == "NoseHooverChain") {
        static std::shared_ptr<OpenMM::NoseHooverIntegrator> NoseHoover =
            std::static_pointer_cast<OpenMM::NoseHooverIntegrator>(integrator);
        return NoseHoover->getTemperature();
    }
}

void HamiltonianOpenMM::setTemperature(double temperature) {
    static const std::string thermostat = param.getStr("thermostat");
    static const std::string barostat = param.getStr("barostat");
    if (thermostat == "none")
        throw std::runtime_error("ERROR: Called setTemperature() without using thermostat.");
    if (thermostat == "Andersen") {
        // This is the only way to change the temperature of heat bath in Context
        // Here, "AndersenTemperature" is the name of the parameter which stores
        // the current temperature in Context. It is known from source code.
        context->setParameter("AndersenTemperature", temperature);
    }
    else if (thermostat == "Langevin") {
        static std::shared_ptr<OpenMM::LangevinIntegrator> Langevin =
            std::static_pointer_cast<OpenMM::LangevinIntegrator>(integrator);
        Langevin->setTemperature(temperature);
    }
    else if (thermostat == "LangevinMiddle") {
        static std::shared_ptr<OpenMM::LangevinMiddleIntegrator> LangevinMiddle =
            std::static_pointer_cast<OpenMM::LangevinMiddleIntegrator>(integrator);
        LangevinMiddle->setTemperature(temperature);
    }
    else if (thermostat == "NoseHooverChain") {
        static std::shared_ptr<OpenMM::NoseHooverIntegrator> NoseHoover =
            std::static_pointer_cast<OpenMM::NoseHooverIntegrator>(integrator);
        NoseHoover->setTemperature(temperature);
    }
    // If barostat is used, the temperature of it should be modified, too.
    if (barostat == "MonteCarlo")
        context->setParameter("MonteCarloTemperature", temperature);
}

double HamiltonianOpenMM::getTime() const {
    return context->getTime();
}

void HamiltonianOpenMM::setTime(double time) {
    context->setTime(time);
}

int HamiltonianOpenMM::getStep() const {
    return context->getStep();
}

void HamiltonianOpenMM::setStep(int step) {
    context->setStep(step);
}

bool HamiltonianOpenMM::hasConstraints() {
    if (system->getNumConstraints() > 0)
        return true;
    else
        return false;
}

void HamiltonianOpenMM::applyConstraints() {
    context->applyConstraints(integrator->getConstraintTolerance());
}

void HamiltonianOpenMM::applyVelocityConstraints() {
    context->applyVelocityConstraints(integrator->getConstraintTolerance());
}

bool HamiltonianOpenMM::hasVirtualSites() {
    for (int i = 0; i < system->getNumParticles(); ++i)
        if (system->isVirtualSite(i))
            return true;
    return false;
}

void HamiltonianOpenMM::computeVirtualSites() {
    context->computeVirtualSites();
}

void HamiltonianOpenMM::reinitialize(bool preserveState) {
    context->reinitialize(preserveState);
}

void HamiltonianOpenMM::createCheckpoint(const std::string& file) {
    std::string old_file = GetFilePrefix(file) + "_prev." + GetFileSuffix(file);
    // If the file exists, it will rename it to old_file (if old_file exists too,
    // the old_file will be override); if the file doesn't exist, do nothing.
    std::rename(file.c_str(), old_file.c_str());
    std::ofstream chkFile(file);
    if (!chkFile)
        throw std::runtime_error("ERROR: Can't create the checkpoint file: " + file);
    context->createCheckpoint(chkFile);
    chkFile.close();
}

void HamiltonianOpenMM::loadCheckpoint(const std::string& file) {
    std::ifstream chkFile(file);
    if (!chkFile)
        throw std::runtime_error("ERROR: Can't open the checkpoint file: " + file);
    context->loadCheckpoint(chkFile);
    chkFile.close();
}

void HamiltonianOpenMM::energyMinimize(const std::string& file, double force_tolerance, int max_iterations) {
    if (force_tolerance <= 0)
        throw std::runtime_error("ERROR: Called energyMinimize() with illegal force_tolerance=" + std::to_string(force_tolerance));
    if (max_iterations < 0)
        throw std::runtime_error("ERROR: Called energyMinimize() with illegal max_iterations=" + std::to_string(max_iterations));
    // Once finished, the Context will have been updated with the new positions.
    OpenMM::LocalEnergyMinimizer em;
    em.minimize(*context, force_tolerance, max_iterations);
    // Save the optimized structure to a file
    // TODO: using structure object in Hamiltonina directly
    if (!file.empty() && file != "none") {
        // Create a structure object and get data from Hamiltonian obejct
        Structure structure;
        std::vector<std::string>&  atominfo   = structure.getAtomInfo();
        std::vector<OpenMM::Vec3>& positions  = structure.getPositions();
        OpenMM::Vec3         (&boxVectors)[3] = structure.getBoxVectors();
        atominfo = this->atominfo;
        structure.setNumAtoms(atominfo.size());
        getPositions(positions);
        if (usesPeriodicBoundaryConditions())
            getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        structure.saveStructure(file, 0, 0);
    }
}

void HamiltonianOpenMM::initializeOpenMMSystem(const std::vector<Topology>& topologies, const Structure& structure, const Parameters& param) {
    // * Set periodic boundary conditions (PBC)  for this System.
    std::string PBC_type = param.getStr("PBC");
    if (PBC_type == "xyz") { // PBC in all directions
        // Box vectors is loaded from structure file and has been converted to
        // the OpenMM reduced form (if nessesary) in structure object.
        OpenMM::Vec3 a, b, c;
        structure.getBoxVectors(a, b, c);
        // It will affect any Force added to the System that uses periodic boundary conditions.
        system->setDefaultPeriodicBoxVectors(a, b, c);
        std::cout << "The periodic boundary conditions (PBC) in all directions will be used in the simulation.\n";
        std::cout << "And the box vectors (in nm) are loaded from structure file and reduced to the form required by OpenMM: \n";
        std::cout << "a = (" << std::fixed << std::setprecision(2) <<
            a[0] << ", " << a[1] << ", " << a[2] << "), b = (" <<
            b[0] << ", " << b[1] << ", " << b[2] << "), c = (" <<
            c[0] << ", " << c[1] << ", " << c[2] << ").\n" << std::endl;
    }
    else if (PBC_type == "none") // any Force added to the System won't use PBC
        std::cout << "The periodic boundary conditions (PBC) will be not used in the simulation.\n" << std::endl;
    else
        throw std::runtime_error("ERROR: Unsupported periodic boundary conditions type, PBC=" + PBC_type);

    // * Set whether to remove motion of centor of mass (COM) to prevent drifting.
    // if do, set the frequency (in time steps) to be removed.
    // The default frequency used in Gromacs is 100 and in Amber is 1000
    // The frequency used in OpenMM (Python) is 1.
    if (param.getBool("remove_COMMotion")) {
        OpenMM::CMMotionRemover* removeMotion = new OpenMM::CMMotionRemover(param.getInt("remove_steps"));
        system->addForce(removeMotion);
        std::cout << "The motion for center of mass (COM) will be removed at each " <<
            removeMotion->getFrequency() << " time steps to prevent the COM from drifting with time.\n" << std::endl;
    }

    // * Set which thermostat to use.
    // Curently, It supports the Langevin/LangevinMiddle/NoseHooverChain
    // integrator-builtin thermostats and Andersen thermostat.
    std::string thermostat_type = param.getStr("thermostat");
    if (thermostat_type != "none") {
        const std::string integrator = param.getStr("integrator");
        // For the integrator-builtin thermostat, the same name intgrator must
        // be used. So settings of them will be done in initializeOpenMMIntegrator().
        if (thermostat_type == "Langevin" || thermostat_type == "LangevinMiddle" ||
            thermostat_type == "NoseHooverChain") {
            if (integrator != thermostat_type)
                throw std::runtime_error("ERROR: When using the " + thermostat_type +
                    " integrator-builtin thermostat, the same name intgrator must be used.");
        }
        else if (thermostat_type == "Andersen") { // Andersen Thermostat
            if (integrator == "Langevin" && integrator == "LangevinMiddle" && integrator == "NoseHooverChain")
                throw std::runtime_error("ERROR: When using Andersen thermostat, the "
                    "integrator cannot be Langevin/LangevinMiddle/NoseHooverChain integrator.");
            // Note: don't use smarter pointer for Force object, the destructor
            // of System object will delete it at the end of program. If using
            // smarter pointer, a Segmentation fault (core dumped) error occurs.
            OpenMM::AndersenThermostat* thermostat =
                new OpenMM::AndersenThermostat(param.getDouble("temperature"), param.getDouble("collision_frequency"));
            // Note: for AndersenThermostat, the precise meaning of this parameter is undefined,
            // and is left up to each Platform to interpret in an appropriate way.
            // And, this behavior will throw an exception when creating context:
            // "EXCEPTION: IntegrationUtilities::initRandomNumberGenerator():
            // Requested two different values for the random number seed"
            // The function is in platforms/common/src/IntegrationUtilities.cpp
            // Therefore, it is no way to produce same results for two differnent
            // simulations when using the AndersenThermostat
            // Note: the random number seed for it cannot specified by user.
            if (param.getInt("thermostat_seed") != 0)
                std::cout << "WARNING: The random number seed for Andersen thermostat "
                    "can't be specified by user, which will be decided by platfrom." << std::endl;
            // Set a name for it, it is useful to locate it in System when nessary.
            thermostat->setName("AndersenThermostat");
            system->addForce(thermostat); // add it to System
            std::cout << "The Andersen thermostat in " << std::fixed << std::setprecision(2) <<
                thermostat->getDefaultTemperature() << " K with a collision frequency " <<
                thermostat->getDefaultCollisionFrequency() << " ps^-1 will be applied.\n" << std::endl;
        }
        else
            throw std::runtime_error("ERROR: Unsupported thermostat=" + thermostat_type);
    }

    // * Set whether to use barostat
    // Curently, only the isotropic Monte Carlo barostat is implemented.
    std::string barostat_type = param.getStr("barostat");
    if (barostat_type != "none") {
        if (param.getStr("PBC") == "none")
            throw std::runtime_error("ERROR: Barostat can only be used for a PBC simulation.\n");
        if (param.getStr("thermostat") == "none")
            throw std::runtime_error("ERROR: When using barostat, a thermostat should be used, too.\n");
        // Currently, only the isotropic Monte Carlo barostat is implemented.
        // When using Monte Carlo barostat (isotropic), assumes the simulation is
        // also being run at constant temperature, which means a NPT simulation.
        if (barostat_type == "MonteCarlo") {
            OpenMM::MonteCarloBarostat* barostat =
                new OpenMM::MonteCarloBarostat(param.getDouble("pressure"),
                    param.getDouble("temperature"), param.getInt("pressure_steps"));
            if (param.getInt("barostat_seed") != 0)
                barostat->setRandomNumberSeed(param.getInt("barostat_seed"));
            barostat->setName("MonteCarloBarostat"); // set a name for it
            system->addForce(barostat); // add it to System
            std::cout << "The isotropic Monte Carlo barostat with a random number seed of " <<
                barostat->getRandomNumberSeed() << "will be used.\n";
            std::cout << "If the random number seed is 0 (default value), which means a unique seed is chosen.\n";
            std::cout << "The reference pressure for coupling is " << std::fixed << std::setprecision(2) <<
                barostat->getDefaultPressure() << "bar. And the frequency for coupling the pressure is " <<
                barostat->getFrequency() << " (in time steps).\n" << std::endl;
        }
        else
            throw std::runtime_error("ERROR: Unsupported barostat=" + barostat_type);
    }

    // * Set whether to use implicit solvation GBSA-OBC model.
    // Note that Gromacs doesn't support this since 2019. And there is no implicit
    // solvent parameters in topology file now.
    // TODO: The implicit solvation model is not yet supported.
    if (param.getBool("implicit_solvation")) {
        OpenMM::GBSAOBCForce* gbsa = new OpenMM::GBSAOBCForce();
        system->addForce(gbsa);
        // Specify dielectrics for GBSA implicit solvation.
        gbsa->setSolventDielectric(param.getDouble("solvent_dielectric"));
        gbsa->setSoluteDielectric(param.getDouble("solute_dielectric"));
        std::cout << "The GBSA-OBC implicit solvation model with solvent dielectric (" <<
            gbsa->getSolventDielectric() << ") and solute dielectric (" << gbsa->getSoluteDielectric() <<
            ") will be used.\n" << std::endl;
        throw std::runtime_error("ERROR: Sorry, implicit solvation model is not yet supported.\n");
    }

    // * Check which bonds angles should be implemented with constraints
    // Allowed values are none, HBonds, AllBonds, or HAngles.
    // The setting of constraints will be done in setForceFieldParameters().
    // Here, only print the information of constraints.
    std::string constraints = param.getStr("constraints");
    if (constraints == "none")
        std::cout << "No distance constraints will be applied in the simulation.\n" << std::endl;
    else if (constraints == "HBonds")
        std::cout << "The distances of the bonds that involve a hydrogen atom will be constrained in the simulation.\n" << std::endl;
    else if (constraints == "AllBonds")
        std::cout << "The distances of all bonds will be constrained in the simulation.\n" << std::endl;
    else if (constraints == "HAngles")
        std::cout << "The distances of all bonds and the angles of the form H-X-H or H-O-X (X is an arbitrary atom) will be constrained in the simulation.\n" << std::endl;
    else
        throw std::runtime_error("ERROR: Illegal value for constraints: \n" + constraints);

    // * Check the setting of rigid water model or flexible water model
    // The setting of constraints will be done in setForceFieldParameters().
    // Here, only print the information of water model.
    for (int i = 0; i < topologies[0].moleculeTypes.size(); i++)
        if (!topologies[0].moleculeTypes[i].settles.empty()) {
            if (param.getBool("rigid_water"))
               std::cout << "Found water molecules and will be fully rigid in the simulation regardless of the value of constraints.\n";
            else if (topologies[0].moleculeTypes[i].angles.empty())
               throw std::runtime_error("ERROR: Found water molecules and flexible water is requested, but there is no angle parameters.\n");
            else if (constraints == "none")
               std::cout << "WARNING: Found water molecules and flexible water is requested, but the value of constraints is not none, which may not reasonable.\n";
            else if (int(param.getInt("DT") * 10000) >= 10)
               std::cout << "WARNING: Found water molecules and flexible water is requested, but step size is larger than 1 fs, which may not reasonable.\n";
            else
               std::cout << "Found water molecules and flexible water will be used in the simulation, which may require to use a small step size such as 0.5 fs.\n";
            std::cout << std::endl;
            break;
        }

    // * Add the physical forces into System, such as bond, angel, dihedral, nonbond
    // If the topology files are more than one, we don't create a new system, but
    // create other Force objects which are belong to different force groups to
    // represent the other system or state. But, create bond and nonbonded terms only.
    const int numTopologies = topologies.size();
    std::vector<OpenMM::NonbondedForce*> nonbonds(numTopologies);
    std::vector<OpenMM::HarmonicBondForce*> bondStretches(numTopologies);
    std::vector<OpenMM::HarmonicAngleForce*> bondBends(numTopologies);
    std::vector<OpenMM::PeriodicTorsionForce*> bondTorsions(numTopologies);
    for (int i = 0; i < numTopologies; ++i) {
        // Create Force objects within the System in heap space.
        // Note that the System owns the force objects and will take care of deleting them.
        // We don't need to delete Force objects.
        nonbonds[i]      = new OpenMM::NonbondedForce();
        bondStretches[i] = new OpenMM::HarmonicBondForce();
        bondBends[i]     = new OpenMM::HarmonicAngleForce();
        bondTorsions[i]  = new OpenMM::PeriodicTorsionForce();
        // Add all Force objects into the System that we want to return.
        system->addForce(nonbonds[i]);
        system->addForce(bondStretches[i]);
        system->addForce(bondBends[i]);
        system->addForce(bondTorsions[i]);
        // Set force group for each Force object so as to evaluate their energies independently.
        // If we don't do this explicitly, all of them would be 0.
        // The force groups of first topology is 0 ~ 31.
        // The force groups of second topology is 32 ~ 63.
        nonbonds[i]->setForceGroup(1 + i * 32);
        bondStretches[i]->setForceGroup(2 + i * 32);
        bondBends[i]->setForceGroup(3 + i * 32);
        bondTorsions[i]->setForceGroup(4 + i * 32);
        // The atomic mass, constraints are from the first (master) topology.
        if (i == 0)
            setForceFieldParameters(topologies[i], param, system.get(), nonbonds[i], bondStretches[i], bondBends[i], bondTorsions[i]);
        // The other topologies only provide the force field paramters of bond/nonbonded terms.
        else {
            // Create a dummy System without Force which is used to record atom information only.
            OpenMM::System dummySystem = OpenMM::System();
            setForceFieldParameters(topologies[i], param, &dummySystem, nonbonds[i], bondStretches[i], bondBends[i], bondTorsions[i]);
        }
        // Set nonbonded method for nonbonded force, if multi-topology is loaded
        // the method for them are same.
        setNonbondedMethod(nonbonds[i], param);
    }

    // * Add fake Coulomb and/or Lennard-Jones Force for energy docompose
    if (param.getBool("energy_decompose")) { // TODO: Coulomb Force for Ewald electric field
        // We create fake Lennard-Jones Forces by copying the original
        // NonbondedForce object and setting the all charges to zero.
        // So the energy calculated from this Force object only inlcuding
        // van der Waals energy.
        std::vector<OpenMM::NonbondedForce*> LJForces(numTopologies);
        for (int i = 0; i < numTopologies; ++i) {
            LJForces[i] = new OpenMM::NonbondedForce();
            *LJForces[i] = *nonbonds[i]; // this is a copy
            // Get nonbonded parameters and reset them with zero charge
            for (int index = 0; index < LJForces[i]->getNumParticles(); ++index) {
                double charge, sigma, epsilon;
                nonbonds[i]->getParticleParameters(index, charge, sigma, epsilon);
                LJForces[i]->setParticleParameters(index, 0.0, sigma, epsilon);
            }
            // Get exception parameters and reset them with zero charge
            for (int index = 0; index < LJForces[i]->getNumExceptions(); ++index) {
                int particle1, particle2;
                double chargeProd, sigma, epsilon;
                nonbonds[i]->getExceptionParameters(index, particle1, particle2, chargeProd, sigma, epsilon);
                LJForces[i]->setExceptionParameters(index, particle1, particle2, 0.0, sigma, epsilon);
            }
            // Set NonbondedMethod to cutoff for Coulomb, which is faster
            // While in the case of non-PBC calculation, keep the original method.
            if (param.getStr("PBC") != "none")
                LJForces[i]->setNonbondedMethod(OpenMM::NonbondedForce::CutoffPeriodic);
            LJForces[i]->setName("LennardJonesForce"); // zero Coulomb energy
            LJForces[i]->setForceGroup(5 + i * 32); // force group 5 in each state
            system->addForce(LJForces[i]);
        }
    }
}

void HamiltonianOpenMM::setForceFieldParameters(const Topology& topology,
                                                 const Parameters& param,
                                                 OpenMM::System* system,
                                                 OpenMM::NonbondedForce* nonbond,
                                                 OpenMM::HarmonicBondForce* bondStretch,
                                                 OpenMM::HarmonicAngleForce* bondBend,
                                                 OpenMM::PeriodicTorsionForce* bondTorsion) {
    // Fill the System and Forces with the force field parameters from topology
    // with looping over molecule types and the specified number of each type.
    for (int i = 0; i < topology.molecules.size(); i++) {
        for (int j = 0; j < topology.molecules[i].moleculeNumber; j++) {
            // Get molecule type of current molecule
            const Topology::MoleculeTypes& moleculeType = topology.moleculeTypes[topology.molecules[i].moleculeTypeIndex];
            // getNumParticles() return the number of atoms has been added to System,
            // record this number to define the atom index in whole system later.
            int currentNumAtoms = system->getNumParticles();
            // This will be used for generating exclusions later.
            std::vector<std::pair<int,int>> bondPairs;
            // * Add particles with mass to System from Atoms.
            // * Add particles with charge and LJ parameters to NonbondedForce from AtomTypes.
            for (int k = 0; k < moleculeType.atoms.size(); k++) {
                double charge = moleculeType.atoms[k].charge;
                double mass = topology.atomTypes[moleculeType.atoms[k].atomTypeIndex].mass;
                double sigma = topology.atomTypes[moleculeType.atoms[k].atomTypeIndex].sigma;
                double epsilon = topology.atomTypes[moleculeType.atoms[k].atomTypeIndex].epsilon;
                system->addParticle(mass);
                nonbond->addParticle(charge, sigma, epsilon);
                // TODO: implicit solvation model
                // If implicit solvation model is used, the gbsa parameters should be added.
                //if (param.getBool("implicit_solvation"))
                    //gbsa.addParticle(atom.charge, atype.gbsaRadiusInNm, atype.gbsaScaleFactor);
            }
            // * Add bonds that without constraints to HarmonicBondForce from Bonds.
            // * Add distance constraints to System from Bonds if constraints applied.
            for (int k = 0; k < moleculeType.bonds.size(); k++) {
                // Note that the atom index from .top file is in its molecule
                // and start from 1 , while the atom index used to add OpenMM
                // System should be the index in whole system and start form 0.
                // The relation between them is:
                int indexI = (moleculeType.bonds[k].atomIndex[0] - 1) + currentNumAtoms;
                int indexJ = (moleculeType.bonds[k].atomIndex[1] - 1) + currentNumAtoms;
                double distance = moleculeType.bonds[k].distance;
                double forceConstant = moleculeType.bonds[k].forceConstant;
                int HydrogenTag = moleculeType.bonds[k].HydrogenTag;
                // Append bond pairs, will be use to generate nonbond exclusions.
                bondPairs.emplace_back(std::make_pair(indexI, indexJ));
                // In Gromacs topology file, only water has a section called [ settles ].
                // Here, we use it to check if it is a water molecule.
                // If rigid water is used, water molecules will be fully rigid
                // regardless of the value passed for the constraints argument.
                std::string constraints = param.getStr("constraints");
                if (!moleculeType.settles.empty() && param.getBool("rigid_water"))
                    system->addConstraint(indexI, indexJ, distance);
                // Constraint level is AllBonds or HAngles:
                else if (constraints == "AllBonds" || constraints == "HAngles")
                    system->addConstraint(indexI, indexJ, distance);
                // Constraint level is HBonds and this is a bond involving a hydrogen atom.
                else if (constraints == "HBonds" && HydrogenTag != 0)
                    system->addConstraint(indexI, indexJ, distance);
                // No constraint (constraints == none)
                else
                    bondStretch->addBond(indexI, indexJ, distance, forceConstant);
                // Set whether to use heavy hydrogen mass.
                // The mass hydrogens that are bonded to heavy atoms can be increased
                // to slow down the fast motions of hydrogens, while keeps their
                // total mass constant. When combined with constraints (typically
                // AllBonds), this allows a further increase in integration step size.
                // Here, we just use a scale factor of 4, which is the most used.
                // And, this behaviour also can be done with "gmx pdb2gmx -heavyh",
                // then the mass in structure file (.gro) is changed.
                // This also works for water molecule.
                // Reference: J. Chem. Theory Comput. 2015, 11, 1864−1874
                // TODO：HMASS check Hydrogen mass repartitioning (HMR)
                if (param.getBool("repart_HMass") && HydrogenTag != 0) {
                    if (constraints == "HBonds" || constraints == "none")
                        std::cout << "WARNING: You are using heavy hydrogen mass, while the constraint level is not AllBonds or HAngles.\n" << std::endl;
                    // The hydrogen is atom I.
                    if (HydrogenTag == 1) {
                        system->setParticleMass(indexI, 4 * 1.008);
                        system->setParticleMass(indexJ, system->getParticleMass(indexJ) - 3 * 1.008);
                    }
                    // The hydrogen is atom J.
                    else {
                        system->setParticleMass(indexJ, 4 * 1.008);
                        system->setParticleMass(indexI, system->getParticleMass(indexI) - 3 * 1.008);
                    }
                }
            }
            // * Add extra constraints to System from Constraints.
            // This is explicitly defined in the section [ constraints ] in .top file.
            // Constraints that are explicitly specified will always be included
            // regardless of the value passed for the constraints argument.
            for (int k = 0; k < moleculeType.constraints.size(); k++) {
                int indexI = (moleculeType.constraints[k].atomIndex[0] - 1) + currentNumAtoms;
                int indexJ = (moleculeType.constraints[k].atomIndex[1] - 1) + currentNumAtoms;
                double distance = moleculeType.constraints[k].distance;
                system->addConstraint(indexI, indexJ, distance);
                // If constraints atoms connected by a chemical bond, then them
                // also will be used to generate nonbond exclusions.
                if (moleculeType.constraints[k].useGenExclusion)
                    bondPairs.emplace_back(std::make_pair(indexI, indexJ));
            }
            // * Add angles to HarmonicAngleForce from Angles.
            // * Add HAngles constraints by adding distance constraints of I-K to System.
            for (int k = 0; k < moleculeType.angles.size(); k++) {
                int indexI = (moleculeType.angles[k].atomIndex[0] - 1) + currentNumAtoms;
                int indexJ = (moleculeType.angles[k].atomIndex[1] - 1) + currentNumAtoms;
                int indexK = (moleculeType.angles[k].atomIndex[2] - 1) + currentNumAtoms;
                double angle = moleculeType.angles[k].angle;
                double forceConstant = moleculeType.angles[k].forceConstant;
                bool isHAngle = moleculeType.angles[k].isHAngle;
                // H-X-H or H-O-X (where X is an arbitrary atom) angles can be constrained.
                // And if HAngles is used, the lengths of all bonds are constrained, too.
                // And if rigid water model is used, the angle of water will be constrained.
                if ((!moleculeType.settles.empty() && param.getBool("rigid_water"))
                    || (param.getStr("constraints") == "HAngles" && isHAngle)) {
                    // A: length of I-J; B: J-K; C: I-K
                    double A = 0, B = 0, C = 0;
                    bool getA = false, getB = false;
                    for (int l = 0; l < moleculeType.bonds.size(); l++) {
                        if (moleculeType.angles[k].atomIndex[1] == moleculeType.bonds[l].atomIndex[0]) {
                            if (moleculeType.angles[k].atomIndex[0] == moleculeType.bonds[l].atomIndex[1]) {
                                    A = moleculeType.bonds[l].distance; getA = true;
                                    if (getB) break;
                                }
                            if (moleculeType.angles[k].atomIndex[2] == moleculeType.bonds[l].atomIndex[1]) {
                                    B = moleculeType.bonds[l].distance; getB = true;
                                    if (getA) break;
                                }
                        }
                        else if (moleculeType.angles[k].atomIndex[1] == moleculeType.bonds[l].atomIndex[1]) {
                            if (moleculeType.angles[k].atomIndex[0] == moleculeType.bonds[l].atomIndex[0]) {
                                    A = moleculeType.bonds[l].distance; getA = true;
                                    if (getB) break;
                                }
                            if (moleculeType.angles[k].atomIndex[2] == moleculeType.bonds[l].atomIndex[0]) {
                                    B = moleculeType.bonds[l].distance; getB = true;
                                    if (getA) break;
                                }
                        }
                    }
                    // Compute distance of I-K (C) and add a distance constraint.
                    C = sqrt(A * A + B * B - 2 * A * B * cos(angle * OpenMM::RadiansPerDegree));
                    system->addConstraint(indexI, indexK, C);
                }
                else
                    // Convert the degree to radian for angle.
                    bondBend->addAngle(indexI, indexJ, indexK, angle * OpenMM::RadiansPerDegree, forceConstant);
            }
            // * Add dihedrals to PeriodicTorsionForce from Dihedrals.
            // Note that both proper and improper dihedrals are defined as same
            // form (periodic type) in AMBER force filed.
            for (int k = 0; k < moleculeType.dihedrals.size(); k++) {
                int indexI = (moleculeType.dihedrals[k].atomIndex[0] - 1) + currentNumAtoms;
                int indexJ = (moleculeType.dihedrals[k].atomIndex[1] - 1) + currentNumAtoms;
                int indexK = (moleculeType.dihedrals[k].atomIndex[2] - 1) + currentNumAtoms;
                int indexL = (moleculeType.dihedrals[k].atomIndex[3] - 1) + currentNumAtoms;
                int periodicity = moleculeType.dihedrals[k].periodicity;
                double dihedral = moleculeType.dihedrals[k].dihedral;
                double forceConstant = moleculeType.dihedrals[k].forceConstant;
                // Note that in OpenMM (Python), the torsion will not be added if
                // its force constant is zero when processing gromacs topology file.
                // But, we add all dihedrals that loaded from topology file to
                // match the number of dihedrals in topology file.
                bondTorsion->addTorsion(indexI, indexJ, indexK, indexL, periodicity,
                    dihedral * OpenMM::RadiansPerDegree, forceConstant);
            }
            // * Add exceptions with createExceptionsFromBonds() to NonbondedForce from bondPairs.
            // OpenMM provides createExceptionsFromBonds() to create the common exceptions:
            // (1) for atoms connected by one bond (1-2) or connected by just one
            // additional bond (1-3), Coulomb and van der Waals terms do not apply;
            // (2) for atoms connected by three bonds (1-4), Coulomb and van der Waals
            // terms apply but are reduced by a force-field dependent scale factor.
            // If there is no [ pairs ] for this molecule, which means to exclude
            // 1-2, 1-3 bonded atoms from nonbonded forces, and calculate 1-4
            // interactions in common way. So we can use this function directly.
            if (moleculeType.pairs.empty() && !bondPairs.empty())
                nonbond->createExceptionsFromBonds(bondPairs, topology.forceFiled.Coulomb14Scale, topology.forceFiled.LennardJones14Scale);
            // If the vector of Pairs is not empty, which means to exclude 1-2,
            // 1-3, 1-4 bonded atoms from nonbonded forces, and 1-4 interactions
            // defined in [ pairs ] should be calculated.
            // Note that this is the most case in Gromacs topology file.
            // Here, we use this function to exclude 1-2, 1-3, but 1-4 interactions
            // are reduced to zero. Since we will create exception for 1-4 interactions
            // that in [ pairs ] additionally (will overide 1-4 interaction defined here).
            else if (!moleculeType.pairs.empty() && !bondPairs.empty())
                nonbond->createExceptionsFromBonds(bondPairs, 0.0, 0.0);
            // * Add extra exceptions (1-4 interactions) manually to NonbondedForce from Pairs.
            // The [ pairs ] defined in .top file describe extra 1-4 interactions explicitly.
            // If this section is not empty (nrexcl = 3) in .top file,
            // 1-4 interactions should excluded from the non-bonded interactions.
            // For AMBER force filed, pair charge/LJ parameters will get from atoms/
            // atomtypes, and calculate nonbonded force normally multiplied by scale factor.
            for (int k = 0; k < moleculeType.pairs.size(); k++) {
                int indexI = (moleculeType.pairs[k].atomIndex[0] - 1) + currentNumAtoms;
                int indexJ = (moleculeType.pairs[k].atomIndex[1] - 1) + currentNumAtoms;
                // Get atom parameters from OpenMM that has been added.
                double chargeI, sigmaI, epsilonI, chargeJ, sigmaJ, epsilonJ;
                nonbond->getParticleParameters(indexI, chargeI, sigmaI, epsilonI);
                nonbond->getParticleParameters(indexJ, chargeJ, sigmaJ, epsilonJ);
                double chargeProd = chargeI * chargeJ * topology.forceFiled.Coulomb14Scale;
                double sigmaComb, epsilonComb;
                // Lorentz-Berthelot rule used by AMBER
                if (topology.forceFiled.combinationRule == 2) {
                    sigmaComb = 0.5 * (sigmaI + sigmaJ);
                    epsilonComb = sqrt(epsilonI * epsilonJ) * topology.forceFiled.LennardJones14Scale;
                }
                // Jorgensen rule used by OPLS-AA force filed
                // TODO: custom combination rule
                else if (topology.forceFiled.combinationRule == 3) {
                    sigmaComb = sqrt(sigmaI + sigmaJ);
                    epsilonComb = sqrt(epsilonI * epsilonJ) * topology.forceFiled.LennardJones14Scale;
                    throw std::runtime_error("ERROR: Sorry, OPLS-AA force filed are not yet supported.");
                }
                else
                    throw std::runtime_error("ERROR: Sorry, other combination rules are not yet supported.");
                // Here, "Replace = true", which means to overide the exceptions that
                // previous created by createExceptionsFromBonds().
                nonbond->addException(indexI, indexJ, chargeProd, sigmaComb, epsilonComb, true);
            }
            // * Add extra exceptions (ignore nonbonded interactions) manually to NonbondedForce from Exclusions.
            // The extra exclusion explicitly defined in the section [ exclusion ] in .top file.
            // Nonbonded interactions between the first and other atoms are excluded.
            // Note that each water model in Gromacs has this section.
            for (int k = 0; k < moleculeType.exclusions.size(); k++) {
                int indexI = (moleculeType.exclusions[k].atomIndex[0] - 1) + currentNumAtoms;
                // Loop from the second atom (n = 1).
                for (int n = 1; n < moleculeType.exclusions[k].atomIndex.size(); n++) {
                    int indexJ = (moleculeType.exclusions[k].atomIndex[n] - 1) + currentNumAtoms;
                    // Set both charge and epsilon to 0, means interaction to
                    // be completely omitted from force and energy calculations.
                    // If there is already an exception for the same two particles
                    // only the last one is valid.
                    nonbond->addException(indexI, indexJ, 0.0, 1.0, 0.0, true);
                }
            }
        }
    }
}

void HamiltonianOpenMM::setNonbondedMethod(OpenMM::NonbondedForce* nonbond, const Parameters& param) {
    // * Set the method used for handling long range nonbonded interactions.
    // Ref: http://docs.openmm.org/latest/api-c++/generated/NonbondedForce.html
    // Enum: OpenMM::NonbondedForce::NonbondedMethod
    // (0) NoCutoff. Computed exactly for non-PBC.
    // (1) CutoffNonPeriodic. Cutoff method for non-PBC.
    // (2) CutoffPeriodic. Cutoff method for PBC.
    // (3) Ewald. Ewald for PBC (expensive) to compute Coulomb interaction.
    // (4) PME. PME for PBC, but cutoff is used for LJ interaction.
    // (5) LJPME. PME for PBC, both for Coulomb and LJ interaction.
    const bool usePBC = param.getStr("PBC") != "none" ? true : false;
    int nonbondedMethod; // using integer represent NonbondedMethod
    std::string nonbonded_method = param.getStr("nonbonded_method");
    if (nonbonded_method == "NoCutoff")
        nonbondedMethod = 0;
    else if (nonbonded_method == "CutoffNonPeriodic")
        nonbondedMethod = 1;
    else if (nonbonded_method == "CutoffPeriodic")
        nonbondedMethod = 2;
    else if (nonbonded_method == "Ewald")
        nonbondedMethod = 3;
    else if (nonbonded_method == "PME")
        nonbondedMethod = 4;
    else if (nonbonded_method == "LJPME")
        nonbondedMethod = 5;
    else
        throw std::runtime_error("ERROR: Unsupported nonbonded_method=" + nonbonded_method);
    switch (nonbondedMethod) { // set nonbondedMethod
        case OpenMM::NonbondedForce::NoCutoff:
            if (!usePBC) {
                nonbond->setNonbondedMethod(OpenMM::NonbondedForce::NoCutoff);
                std::cout << "The NoCutoff nonbonded method will be used to compute nonbonded interactions exactly.\n" << std::endl;
            }
            else
                throw std::runtime_error("ERROR: NoCutoff method is not supported when PBC is used.\n");
            break;
        case OpenMM::NonbondedForce::CutoffNonPeriodic:
            if (!usePBC) {
                nonbond->setNonbondedMethod(OpenMM::NonbondedForce::CutoffNonPeriodic);
                std::cout << "The CutoffNonPeriodic nonbonded method will be used to compute the nonbonded interactions. Interactions beyond the cutoff distance are ignored.\n";
                std::cout << "Coulomb interactions closer than the cutoff distance are modified using the reaction field method.\n" << std::endl;
            }
            else
                throw std::runtime_error("ERROR: CutoffNonPeriodic method is not supported when PBC is used.\n");
            break;
        case OpenMM::NonbondedForce::CutoffPeriodic:
            if (usePBC) {
                nonbond->setNonbondedMethod(OpenMM::NonbondedForce::CutoffPeriodic);
                std::cout << "The CutoffPeriodic nonbonded method will be used to compute the nonbonded interactions. Interactions beyond the cutoff distance are ignored.\n";
                std::cout << "Coulomb interactions closer than the cutoff distance are modified using the reaction field method.\n" << std::endl;
            }
            else
                throw std::runtime_error("ERROR: CutoffPeriodic nonbonded method is not supported when PBC is not used.\n");
            break;
        case OpenMM::NonbondedForce::Ewald:
            if (usePBC) {
                nonbond->setNonbondedMethod(OpenMM::NonbondedForce::Ewald);
                std::cout << "The Ewald nonbonded method will be used to compute the long range Coulomb interactions.\n" << std::endl;
            }
            else
                throw std::runtime_error("ERROR: Ewald nonbonded method is not supported when PBC is not used.\n");
            break;
        case OpenMM::NonbondedForce::PME:
            if (usePBC) {
                nonbond->setNonbondedMethod(OpenMM::NonbondedForce::PME);
                std::cout << "The PME nonbonded method will be used to compute the long range Coulomb interactions.\n" << std::endl;
            }
            else
                throw std::runtime_error("ERROR: PME nonbonded method is not supported when PBC is not used.\n");
            break;
        case OpenMM::NonbondedForce::LJPME:
            if (usePBC) {
                nonbond->setNonbondedMethod(OpenMM::NonbondedForce::LJPME);
                std::cout << "The LJPME nonbonded method will be used to compute the both long range Coulomb and Lennard-Jones interactions.\n" << std::endl;
            }
            else
                throw std::runtime_error("ERROR: LJPME nonbonded method is not supported when PBC is not used.\n");
            break;
        default:
            throw std::runtime_error("ERROR: Unsupported nonbonded method.\n");
    }
    // Set the cutoff distance (in nm) being used for nonbonded interactions.
    // Valid for all methods except NoCutoff. The default value is 1.0 nm.
    if (nonbondedMethod != 0) {
        nonbond->setCutoffDistance(param.getDouble("cutoff"));
        std::cout << "A cutoff distance of " << std::fixed << std::setprecision(2) <<
            nonbond->getCutoffDistance() << " nm will be used in the calculations for nonbonded interactions.\n" << std::endl;
    }
    // Add a contribution to the energy which approximates the effect of all
    // Lennard-Jones interactions beyond the cutoff in a periodic system.
    // When running a simulation at constant pressure, this can improve the quality
    // of the result. This correction can also be applied to force in Gromacs.
    // Only valid when do PBC calculation.
    // This is enabled by default.
    if (param.getBool("dispersion_correction"))
        if (usePBC) {
            nonbond->setUseDispersionCorrection(true);
            std::cout << "A long range dispersion correction to the energy which approximates the effect of all "
                "Lennard-Jones interactions beyond the cutoff in a periodic system will be applied.\n" << std::endl;
        }
        else
            throw std::runtime_error("ERROR: The dispersion correction only works for PBC=xyz.");
    else
        nonbond->setUseDispersionCorrection(false);
    // Set the error tolerance for Ewald summation. This corresponds to the
    // fractional error in the forces which is acceptable. This value is used to
    // select the reciprocal space cutoff and separation parameter so that the
    // average error level will be less than the tolerance. However, there is not a
    // rigorous guarantee that all forces on all atoms will be less than the tolerance.
    // Only valid for Ewald/PME/LJPME method.
    // The default value is 0.0005.
    if (nonbondedMethod > 2) {
        nonbond->setEwaldErrorTolerance(param.getDouble("Ewald_tolerance"));
        std::cout << "The Ewald error tolerance is " << std::setprecision(6) <<
            nonbond->getEwaldErrorTolerance() << ". The parameters for " << nonbonded_method <<
            " calculations are chosen based on this tolerance.\n" << std::endl;
    }
    // Set the parameters for PME calculations by user-defined.
    // PMEalpha：the separation parameter.
    // PMEnx/y/z: the number of grid points along the X/Y/Z axis.
    // Only valid for (LJ)PME method.
    // By default (PMEalpha = 0), the parameters of (LJ)PME nonbonded method are
    // based on Ewald error tolerance. We don't set them explicitly.
    if (param.getStr("PME_parameters") != "0") {
        if (nonbondedMethod < 4)
            throw std::runtime_error("ERROR: The PME parameters is specified explicitly, "
                "while the nonbonded method is not PME or LJPME.");
        std::vector<std::string> PME_parameters;
        SplitString(PME_parameters, param.getStr("PME_parameters"));
        if (PME_parameters.size() != 4)
            throw std::runtime_error("ERROR: Illeagal PME_parameters=" + param.getStr("PME_parameters") +
                ". It should be 4 values representing (alpha,nx,ny,nz).");
        nonbond->setPMEParameters(std::stod(PME_parameters[0]), std::stod(PME_parameters[1]),
            std::stod(PME_parameters[2]), std::stod(PME_parameters[3]));
        std::cout << "The PME parameters are specified by user, and the values of them are as following:\n";
        std::cout << "The separation parameter is " << PME_parameters[0] << ", " <<
            "and the number of grid points along the X/Y/Z axis are " << PME_parameters[1] <<
            ", " << PME_parameters[2] << ", and" << PME_parameters[3] << ", respectively.\n" << std::endl;
    }
    // Set the parameters for the dispersion term in LJPME calculations by user-defined.
    // Only valid for LJPME method.
    if (param.getStr("LJPME_parameters") != "0") {
        if (nonbondedMethod != 5)
            throw std::runtime_error("ERROR: The LJPME parameters is specified explicitly, "
                "while the nonbonded method is not LJPME.");
        std::vector<std::string> LJPME_parameters;
        SplitString(LJPME_parameters, param.getStr("LJPME_parameters"));
        if (LJPME_parameters.size() != 4)
            throw std::runtime_error("ERROR: Illeagal LJPME_parameters=" + param.getStr("LJPME_parameters") +
                ". It should be 4 values representing (alpha,nx,ny,nz).");
        nonbond->setLJPMEParameters(std::stod(LJPME_parameters[0]), std::stod(LJPME_parameters[1]),
            std::stod(LJPME_parameters[2]), std::stod(LJPME_parameters[3]));
        std::cout << "The LJPME parameters are specified by user, and the values of them are as following:\n";
        std::cout << "The separation parameter is " << LJPME_parameters[0] << ", " <<
            "and the number of grid points along the X/Y/Z axis are " << LJPME_parameters[1] <<
            ", " << LJPME_parameters[2] << ", and" << LJPME_parameters[3] << ", respectively.\n" << std::endl;
    }
    // When using a cutoff, by default Lennard-Jones interactions are sharply
    // truncated at the cutoff distance. Optionally you can instead use a switching
    // function to make the interaction smoothly go to zero over a finite distance range.
    // You must also specify the distance at which the interaction should begin
    // to decrease. The switching distance must be less than the cutoff distance.
    // Only valid for CutoffPeriodic/PME/Ewald method.
    // By default, this is false.
    if (param.getBool("switching_function")) {
        if (nonbondedMethod == 0)
            throw std::runtime_error("ERROR: If using switching function, the "
                "nonbonded method can not be NoCutoff.");
        double switching_distance = param.getDouble("switching_distance");
        if (switching_distance >= nonbond->getCutoffDistance())
            throw std::runtime_error("ERROR: If using switching function, the "
                "switch distance must be less than the cutoff distance.");
        nonbond->setUseSwitchingFunction(true);
        nonbond->setSwitchingDistance(switching_distance);
    }
    // TODO: When GBSA implict solvation is enabled ?
    // When using CutoffNonPeriodic and CutoffPeriodic nonbonded method, Coulomb
    // interactions closer than the cutoff distance are modified using the reaction
    // field method. Therefore, you may want to set the the dielectric constant
    // to use for the solvent in the reaction field approximation.
    // Only valid for CutoffNonPeriodic or CutoffPeriodic method.
    // The default value is 78.3. This is used when using implict solvation model.
    // Note: when using implict solvation GBSA model, you should call
    // setReactionFieldDielectric(1.0) on the NonbondedForce to turn off the
    // reaction field approximation, which does not produce correct results
    // when combined with GBSA.
    if (param.getBool("reaction_field")) {
        if (nonbondedMethod != 1 && nonbondedMethod != 2)
            throw std::runtime_error("ERROR: If using reaction field approximation, the "
                "nonbonded method should be CutoffNonPeriodic or CutoffPeriodic.");
        nonbond->setReactionFieldDielectric(param.getDouble("RF_dielectric"));
        std::cout << "The dielectric constant to use for the solvent in the reaction field approximation is " <<
            std::setprecision(2) << nonbond->getReactionFieldDielectric() << ".\n" << std::endl;
    }
}

void HamiltonianOpenMM::addRestraintForce(Structure& structure, const Parameters& param) {
    // 1. Get indices of restraint_atoms
    std::vector<int> restraint_atoms; // list of indices of restraint atoms
    GetIndexList(restraint_atoms, param.getStr("restraint_atoms"));
    if (restraint_atoms.size() == 0)
        throw std::runtime_error("ERROR: position_restraint is true while restraint_atoms is empty.");
    // 2. Get restraint force constant
    double k = param.getDouble("restraint_force"); // force constant, in kj/mol/mn^2
    if (k <= 0)
        throw std::runtime_error("ERROR: Illegal restraint_force=" + std::to_string(k));
    // 3. Get reference positions of restraint_atoms from gromacs structure file
    const std::string reference_file = param.getStr("reference_structure");
    if (reference_file.empty())
        throw std::runtime_error("ERROR: position_restraint is true while reference_structure is empty.");
    std::vector<OpenMM::Vec3> reference_positions;
    // If the files of structure and reference_structure are same, then the
    // reference position is from the initial structure.
    if (param.getStr("structure") == reference_file)
        reference_positions = structure.getPositions(); // this is a copy
    else { // the reference_structure is different from initial structure
        Structure reference_structure; // a object
        reference_structure.loadStructure(reference_file);
        // Check if they are same system (atominfo).
        if (structure.getAtomInfo() != reference_structure.getAtomInfo())
            throw std::runtime_error("ERROR: The atominfo in reference_structure=" + reference_file +
                " is inconsistent with that in structure=" + param.getStr("structure"));
        reference_positions = reference_structure.getPositions(); // this is a copy
    }
    // If barostat is used, we need to translate the initial structure
    // positions so that the centroid of the restrained molecule is at
    // the origin, which can avoid MonteCarloBarostat rejections (see openmm#1854).
    // However, in this case, the restrained atoms must belong to one molecule.
    if (param.getStr("barostat") != "none") {
        // Create a context to get the molecule infomation.
        // Note: Do this after adding the particals into system.
        OpenMM::Platform& platform = OpenMM::Platform::getPlatformByName("Reference");
        OpenMM::VerletIntegrator verlert_integrator(0.001);
        OpenMM::Integrator* integrator = &(verlert_integrator);
        OpenMM::Context context(*system, *integrator, platform);
        // Get a description of how the particles in the system are grouped
        // into molecules. Two particles are in the same molecule if they
        // are connected by constraints or bonds.
        // It is descriped by atom indices in assending order.
        const std::vector<std::vector<int>>& molecules_atoms = context.getMolecules();
        // Make sure the atoms to restrain belong only to a single molecule.
        std::vector<int> restrained_molecule_atoms; // all atoms in the restrained molecule
        for (auto molecule_atoms : molecules_atoms) {
            // if max_index of restraint_atoms is <= the max_index of this molecule,
            // then min_index of restraint_atoms should be >= min_index of this molecule,
            // otherwise, the restraint_atoms don't belong to same molecule.
            if (restraint_atoms[restraint_atoms.size()-1] <=  molecule_atoms[molecule_atoms.size()-1])
                if (restraint_atoms[0] < molecule_atoms[0])
                    throw std::runtime_error("The restraint_atoms=" + param.getStr("restraint_atoms") +
                        " don't belong to a single molecule, which is not supported when using a barostat.");
                else { // belong to current molecule, then exit the loop
                    restrained_molecule_atoms = molecule_atoms;
                    break;
                }
        }
        // Translate structure so that the centroid of restrained molecule is
        // in the origin to avoid the barostat rejections. In this case,
        // the barostat won't change it, but only works for single molecule.
        // 1. Translate the reference structure
        OpenMM::Vec3 reference_centroid; // centroid of restrained molecule (reference)
        for (int index : restrained_molecule_atoms)
            reference_centroid += reference_positions[index];
        reference_centroid /= restrained_molecule_atoms.size();
        for (int index = 0; index < reference_positions.size(); ++index)
            reference_positions[index] -= reference_centroid;
        // 2. (optional) If the structure and reference_structure are not
        // same, we need to translate the structure by itself, too.
        // Then, update the initial positions in structure object.
        if (param.getStr("structure") != reference_file) {
            std::vector<OpenMM::Vec3> positions = structure.getPositions();
            OpenMM::Vec3 centroid; // centroid of restrained molecule (initial)
            for (int index : restrained_molecule_atoms)
                centroid += positions[index];
            centroid /= restrained_molecule_atoms.size();
            for (int index = 0; index < positions.size(); ++index)
                positions[index] -= centroid;
            structure.setPositions(positions);
        }
        else
            structure.setPositions(reference_positions);
    }
    // 4. Create OpenMM CustomExternalForce with harmonic energy_expression
    std::string energy_expression; // energy_expression for restraint_force
    if (system->usesPeriodicBoundaryConditions()) // periodic distance
        energy_expression = "(k/2)*periodicdistance(x, y, z, x0, y0, z0)^2";
    else // non-periodic distance
        energy_expression = "(k/2)*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)";
    OpenMM::CustomExternalForce* restraint_force = new OpenMM::CustomExternalForce(energy_expression);
    // add the four parameters needed for position restraint
    restraint_force->addGlobalParameter("k", k);
    restraint_force->addPerParticleParameter("x0");
    restraint_force->addPerParticleParameter("y0");
    restraint_force->addPerParticleParameter("z0");
    std::vector<double> parameters(3, 0); // x0, y0, z0
    for (int index : restraint_atoms) {
        for (int i = 0; i < 3; ++i)
            parameters[i] = reference_positions[index][i];
        restraint_force->addParticle(index, parameters);
    }
    // 5. Add this restraint force to system
    // By default, the force group of it is 0, (the group 0 is a common group
    // for all states, which will be included in propagtaion)
    // We can set a name for it, by default, the name is class name.
    restraint_force->setName("HarmonicRestraintForce");
    system->addForce(restraint_force);
    std::cout << "The position restraint will be applied to atoms: " << param.getStr("restraint_atoms") <<
        " with a force constant: " << param.getStr("restraint_force") << " kj/mol/nm^2.\n";
    std::cout << "The reference positions of restraint atoms are loaded from structure file: " <<
        reference_file << ".\n";
    if (param.getStr("barostat") != "none") {
        // Save the translated reference structure which can be used to
        // in the subsequent simulation (following this NPT simulation). If the
        // original reference structure is used in the subsequent simulation,
        // the simulation may be not stable due to the fact that the position of
        // restrained molecule of last frame may be far away from the original
        // reference structure and thus resulting in large restraint forces.
        const std::string translated_file = GetFilePrefix(reference_file) +
            "_translated." + GetFileSuffix(reference_file);
        structure.saveStructure(translated_file, 0, 0);
        std::cout << "Since barostat is used, the structure is translated so that the centroid "
            "of restrained molecule is in the origin to avoid the barostat rejections.\n";
        std::cout << "And the translated structure is saved to the file: " << translated_file << "\n";
        std::cout << "One can use this translated structure as reference structure in the subsequent simulation.\n";
    }
    std::cout << std::endl;
}

void HamiltonianOpenMM::initializeOpenMMIntegrator(const Parameters& param) {
    const std::string integrator_type = param.getStr("integrator");
    const double      DT              = param.getDouble("DT"); // time step in ps
    const int         propagate_state = param.getInt("propagate_state"); // starts from 0
    // 1. NVT or NPT simulation with Langevin Dynamics as well as Langevin thermostat.
    if (integrator_type == "Langevin") {
        if (param.getStr("thermostat") != "Langevin")
            throw std::runtime_error("ERROR: If Langevin integrator is used, the thermostat must be specified as Langevin, too.");
        std::shared_ptr<OpenMM::LangevinIntegrator> Langevin =
            std::make_shared<OpenMM::LangevinIntegrator>(param.getDouble("temperature"),
                param.getDouble("friction_coefficient"), DT);
        // Set the random number seed. If seed is set to 0 (which is the default value assigned),
        // a unique seed is chosen when a Context is created from this Force. This is done to ensure
        // that each Context receives unique random seeds without you needing to set them explicitly.
        if (param.getInt("thermostat_seed") != 0)
            Langevin->setRandomNumberSeed(param.getInt("thermostat_seed"));
        std::cout << "The OpenMM LangevinIntegrator with builtin thermostat will be used for NVT or NPT simulation.\n";
        std::cout << "The temperature of the heat bath is " << std::setprecision(2) <<
            Langevin->getTemperature() << " K.\n";
        std::cout << "The friction coefficient which couples the system to the heat bath is " <<
            Langevin->getFriction() << " ps^-1.\n";
        std::cout << "The random number seed is " << Langevin->getRandomNumberSeed() << ".\n";
        if (Langevin->getRandomNumberSeed() == 0)
            std::cout <<"The random number seed is 0 (default value), which means a unique seed is chosen.\n";
        std::cout << std::endl;
        integrator = Langevin;
    }
    // 2. Langevin dynamics with the LFMiddle discretization, which is tend to
    // produce more accurate configurational sampling than LangevinIntegrator
    // and can use the a bigger step size, such as 4 fs when using HBonds constriants.
    // Note, it returns half step (leapfrog) velocities.
    // Reference: J. Phys. Chem. A 2019, 123, 28, 6056-6079
    else if (integrator_type == "LangevinMiddle") {
        if (param.getStr("thermostat") != "LangevinMiddle")
            throw std::runtime_error("ERROR: If LangevinMiddle integrator is used, "
                "the thermostat must be specified as LangevinMiddle, too.");
        std::shared_ptr<OpenMM::LangevinMiddleIntegrator> LangevinMiddle =
            std::make_shared<OpenMM::LangevinMiddleIntegrator>(param.getDouble("temperature"),
                param.getDouble("friction_coefficient"), DT);
        // Set the random number seed. If seed is set to 0 (which is the default value assigned),
        // a unique seed is chosen when a Context is created from this Force. This is done to ensure
        // that each Context receives unique random seeds without you needing to set them explicitly.
        if (param.getInt("thermostat_seed") != 0)
            LangevinMiddle->setRandomNumberSeed(param.getInt("thermostat_seed"));
        std::cout << "The OpenMM LangevinMiddleIntegrator with builtin thermostat will be used for NVT or NPT simulation.\n";
        std::cout << "The temperature of the heat bath is " << std::setprecision(2) <<
            LangevinMiddle->getTemperature() << " K.\n";
        std::cout << "The friction coefficient which couples the system to the heat bath is " <<
            LangevinMiddle->getFriction() << " ps^-1.\n";
        std::cout << "The random number seed is " << LangevinMiddle->getRandomNumberSeed() << ".\n";
        if (LangevinMiddle->getRandomNumberSeed() == 0)
            std::cout <<"The random number seed is 0 (default value), which means a unique seed is chosen.\n";
        std::cout << std::endl;
        integrator = LangevinMiddle;
    }
    // 3. The "middle" leapfrog propagation (J. Phys. Chem. A 2019, 123, 6056-6079)
    // using one or more Nose Hoover chain thermostats. The advantage of it is
    // that several thermostat can be created to control the different subsystems.
    // TODO: multi thermostate for different subsystems is not supported
    else if (integrator_type == "NoseHooverChain") {
        if (param.getStr("thermostat") != "NoseHooverChain")
            throw std::runtime_error("ERROR: If NoseHoover integrator is used, "
                "the thermostat must be specified as NoseHooverChain, too");
        std::shared_ptr<OpenMM::NoseHooverIntegrator> NoseHooverChain;
        // Set the detail parameters for NoseHooverChain thermostat with three integers.
        // They represent the number of beads in the Nose-Hoover chain, the number of step
        // in the multiple time step chain propagation algorithm, the number of terms in
        // the Yoshida-Suzuki multi time step decomposition used in the chain propagation
        // algorithm (must be 1, 3, 5, or 7), respectively.
        const std::string NHC_parameters = param.getStr("NHC_parameters");
        if (NHC_parameters.empty() || NHC_parameters == "3,3,7") // default parameters
            NoseHooverChain = std::make_shared<OpenMM::NoseHooverIntegrator>(
                param.getDouble("temperature"), param.getDouble("collision_frequency"), DT);
        else {
            std::vector<int> parameters;
            SplitString(parameters, NHC_parameters);
            if (parameters.size() != 3)
                throw std::runtime_error("ERROR: Use NoseHooverChain thermostat "
                    "with wrong number of NHC_parameters=" + NHC_parameters +
                    ", only three non-zero integers are accepted.");
            NoseHooverChain = std::make_shared<OpenMM::NoseHooverIntegrator>(
                param.getDouble("temperature"), param.getDouble("collision_frequency"), DT,
                parameters[0], parameters[1], parameters[2]);
        }
        // Note: the random number seed for it cannot specified by user.
        if (param.getInt("thermostat_seed") != 0)
            std::cout << "WARNING: The random number seed for NoseHooverChain thermostat "
                "can't be specified by user, which will be decided by platfrom.\n";
        std::cout << "The OpenMM NoseHooverIntegrator with a simple Nose-Hoover Chain thermostat "
            "to control the temperature of the full system will be used for simulation.\n";
        std::cout << "The temperature of the heat bath is " << std::setprecision(2) <<
            NoseHooverChain->getTemperature() << " K.\n";
        std::cout << "The friction coefficient which couples the system to the heat bath is " <<
            NoseHooverChain->getCollisionFrequency() << " ps^-1.\n";
        if (NHC_parameters.empty() || NHC_parameters == "3,3,7")
            std::cout << "The NoseHooverChain parameters (number of beads, number of steps "
                "in the multiple time step chain propagation, and number of terms in the "
                "Yoshida-Suzuki multi time step decomposition) are: 3,3,7, respectively.\n";
        else
            std::cout << "The NoseHooverChain parameters (number of beads, number of steps "
                "in the multiple time step chain propagation, and number of terms in the "
                "Yoshida-Suzuki multi time step decomposition) are: " << NHC_parameters << ", respectively.\n";
        if (param.getInt("thermostat_seed") != 0)
            std::cout << "WARNING: The random number seed for Nose-Hoover Chain thermostat "
                "can't be specified by user, which will be decided by platfrom.\n";
        std::cout << std::endl;
        integrator = NoseHooverChain;
    }
    // 4. NVE simulation with leapfrog Verlet Integrator, is the default integrator in Gromacs.
    else if (integrator_type == "leapfrog") {
        integrator = std::make_shared<OpenMM::VerletIntegrator>(DT);
        std::cout << "The OpenMM VerletIntegrator using the leapfrog Verlet algorithm will be used for simulation.\n" << std::endl;
    }
    // 5. NVE simulation with velocity Verlet Integrator.
    else if (integrator_type == "velocityVerlet") {
        // velocity Verlet Integrator for single surface, which is a CustomIntegrator.
        if (param.getStr("dyn_type") == "OpenMM") {
            std::shared_ptr<OpenMM::CustomIntegrator> custom = std::make_shared<OpenMM::CustomIntegrator>(DT);
            custom->addPerDofVariable("x1", 0);
            custom->addUpdateContextState();
            custom->addComputePerDof("v", "v+0.5*dt*f/m");
            custom->addComputePerDof("x", "x+dt*v");
            custom->addComputePerDof("x1", "x");
            custom->addConstrainPositions();
            custom->addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt");
            custom->addConstrainVelocities();
            std::cout << "The velocity Verlet Integrator will be used for simulation.\n" << std::endl;
            integrator = custom;
        }
        // Use CompoundIntegrator with two CustomIntegrator to do velocity Verlet
        // for nonadiabatic dynamics with external forces.
        // Note: in this case, the ConstraintTolerance and IntegrationForceGroups
        // should be set separatedly for each integrator.
        else {
            OpenMM::CustomIntegrator* vv1 = new OpenMM::CustomIntegrator(DT);
            vv1->addPerDofVariable("x1", 0);
            vv1->addComputePerDof("v", "v+0.5*dt*f/m");
            vv1->addComputePerDof("x", "x+dt*v");
            vv1->addComputePerDof("x1", "x");
            vv1->addConstrainPositions();
            vv1->addComputePerDof("v", "v+(x-x1)/dt");
            vv1->setConstraintTolerance(param.getDouble("constraint_tolerance"));
            vv1->setIntegrationForceGroups(forceGroups[propagate_state]);
            OpenMM::CustomIntegrator* vv2 = new OpenMM::CustomIntegrator(DT);
            vv2->addComputePerDof("v", "v+0.5*dt*f/m");
            vv2->addConstrainVelocities();
            vv2->setConstraintTolerance(param.getDouble("constraint_tolerance"));
            vv2->setIntegrationForceGroups(forceGroups[propagate_state]);
            std::shared_ptr<OpenMM::CompoundIntegrator> compound = std::make_shared<OpenMM::CompoundIntegrator>();
            compound->addIntegrator(vv1);
            compound->addIntegrator(vv2);
            std::cout << "The velocity Verlet CompoundIntegrator will be used for simulation.\n" << std::endl;
            integrator = compound;
        }
    }
    else
        throw std::runtime_error("ERROR: Unsupported integrator: " + integrator_type);
    // Set the distance tolerance for constraint (a fraction of the constrained distance).
    // Note that this doesn't apply constraints to the start configuration.
    // Default value is 0.00001.
    if (param.getBool("rigid_water") || param.getStr("constraints") != "none") {
        integrator->setConstraintTolerance(param.getDouble("constraint_tolerance"));
        std::cout << "The distance tolerance within which constraints are maintained is " <<
            std::setprecision(6) << integrator->getConstraintTolerance() << ", as a fraction "
            "of the constrained distance.\n" << std::endl;
    }
    std::cout << "The step size with which to integrate the system is " << std::setprecision(5) <<
        integrator->getStepSize() << " ps.\n";
    // Set which state to do propagation by specifying the force groups which
    // should be included in the propagation.
    if (propagate_state < 0 || propagate_state >= DOFe)
        throw std::runtime_error("ERROR: Illegal value for propagate_state=" + param.getStr("propagate_state"));
    integrator->setIntegrationForceGroups(forceGroups[propagate_state]);
    std::cout << "The nuclear dynamics propagation will be performed on the surface of state " <<
        propagate_state << ".\n" << std::endl;
}

OpenMM::Platform* HamiltonianOpenMM::createOpenMMPlatform(const Parameters& param) {
    std::cout << "The version number of the OpenMM library is " << OpenMM::Platform::getOpenMMVersion() <<
        ".\n" << std::endl;
    // 1. Load all plugins before specifying the Platform.
    if (param.getStr("plugins_directory").empty()) {
        const std::string& plugins_directory = OpenMM::Platform::getDefaultPluginsDirectory();
        OpenMM::Platform::loadPluginsFromDirectory(plugins_directory);
        std::cout << "The default directory from which to load OpenMM plugins is " <<
            plugins_directory << ".\n" << std::endl;
    }
    else {
        OpenMM::Platform::loadPluginsFromDirectory(param.getStr("plugins_directory"));
        std::cout << "The directory defined by user from which to load OpenMM plugins is " <<
            param.getStr("plugins_directory") << ".\n" << std::endl;
    }
    // 2. Specify the Platform and Platform-specific properties, if requested.
    // Note that we need specify properties before creating a Context.
    // Since the properties are read-only after a Context is created.
    OpenMM::Platform* platform = nullptr;
    const std::string platformName = param.getStr("platform");
    if (platformName == "CPU" || platformName == "CUDA" || platformName == "OpenCL" || platformName == "Reference") {
        // Platform is a single static object, we don't need to create it with "new".
        platform = &(OpenMM::Platform::getPlatformByName(platformName));
        // Set Platform-specific properties.
        // Only valid properties which are different from default value for specific
        // platform will be set, and others will be ignored.
        if (platformName == "CPU") {
            if (!param.getStr("CPU_threads").empty())
                platform->setPropertyDefaultValue("Threads", param.getStr("CPU_threads"));
            if (param.getStr("deterministic_forces") == "true")
                platform->setPropertyDefaultValue("DeterministicForces", "true");
            if (param.getStr("precision") == "mixed" || param.getStr("precision") == "double")
                std::cout << "WARNING: The CPU platform supports single precision only!" << std::endl;
        }
        else if (platformName == "CUDA" || platformName == "OpenCL") {
            // Properties used both for CUDA and OpenCL
            if (param.getStr("precision") == "mixed" || param.getStr("precision") == "double")
                platform->setPropertyDefaultValue("Precision", param.getStr("precision"));
            if (param.getStr("CPU_PME") == "true")
                platform->setPropertyDefaultValue("UseCpuPme", "true");
            if (param.getStr("disable_PMEStream") == "true")
                platform->setPropertyDefaultValue("DisablePmeStream", "true");
            if (!param.getStr("device_index").empty())
                platform->setPropertyDefaultValue("DeviceIndex", param.getStr("device_index"));
            // CUDA properties only
            if (platformName == "CUDA") {
                if (param.getStr("deterministic_forces") == "true")
                    platform->setPropertyDefaultValue("DeterministicForces", "true");
                if (param.getStr("blocking_sync") == "false")
                    platform->setPropertyDefaultValue("UseBlockingSync", "false");
                if (!param.getStr("CUDA_compiler").empty())
                    platform->setPropertyDefaultValue("CudaComplier", param.getStr("CUDA_compiler"));
                if (!param.getStr("temp_directory").empty())
                    platform->setPropertyDefaultValue("TempDirectory", param.getStr("temp_directory"));
            }
            // OpenCL properties only
            if (platformName == "OpenCL")
                if (!param.getStr("OpenCL_index").empty())
                    platform->setPropertyDefaultValue("OpenCLPlatformIndex", param.getStr("OpenCL_index"));
        }
    }
    // If not specified, the Context will choose the best available Platform.
    else if (platformName == "auto" || platformName.empty())
        platform = nullptr;
    else
        throw std::runtime_error("ERROR: Unsupported platform: " + platformName);
    return platform;
}

void HamiltonianOpenMM::initializeOpenMMContext(const Structure& structure, const Parameters& param) {
    const std::vector<OpenMM::Vec3>& positions = structure.getPositions();
    const std::vector<OpenMM::Vec3>& velocities = structure.getVelocities();
    const int numAtoms = positions.size();
    const int numVelocities = velocities.size();
    // 1. Set initial positions.
    context->setPositions(positions);
    // 2. Set initial velocities. The allowed value: 0, file, Boltzmann, auto.
    const std::string initial_velocity = param.getStr("initial_velocity");
    // 2.1 Apply a half step correction to initial velocities when using leapfrog integrator.
    if (param.getBool("velocity_correction")) {
        // TODO: How to define the kinetic energy (same)
        // TODO: used for switching between a leapfrog and non-leapfrog integrator
        if (param.getStr("integrator") != "leapfrog")
            throw std::runtime_error("ERROR: The half step correction to initial "
                "velocities can only be used for leapfrog Verlet integrator.");
        if (initial_velocity== "Boltzmann")
            throw std::runtime_error("ERROR: The half step correction to initial "
                "velocities can't be used when random velocities is set.");
        // Note: the forces in this correction is from the propagate state,
        // which is the same as the force groups in the integrator.
        // Here, use the State object (a snapshot of Context) to get forces.
        // The another way to get the forces is to call firstly
        // getPotentialEnergy(integrator->getIntegrationForceGroups(), true)
        // then call getForces(std::vector<OpenMM::Vec3>& forces) to get froces
        // directly from OpenMM::Context.
        // These two ways give you identical forces. The latter is more
        // efficient but the former is more safety. Since we do this only
        // at the initial step, so the former is used.
        const OpenMM::State& state = context->getState(OpenMM::State::Forces,
            false, integrator->getIntegrationForceGroups());
        const std::vector<OpenMM::Vec3>& forces = state.getForces();
        std::vector<OpenMM::Vec3> newVelocities = velocities;
        if (numVelocities == 0)
            newVelocities.resize(numAtoms, OpenMM::Vec3());
        for (int i = 0; i < numAtoms; i++) {
            double mass = context->getSystem().getParticleMass(i);
            if (mass != 0)
                for (int j = 0; j < 3; ++j)
                    newVelocities[i][j] -= 0.5*forces[i][j]*param.getDouble("DT")/mass;
        }
        context->setVelocities(newVelocities);
        std::cout << "A half step correction to initial velocities is applied: V_new = V_old - 0.5*dt*F/m. "
            "V_old is from structure file, if not exist, V_old is zero.\n" << std::endl;
    }
    // 2.2 Set initial velocities automatically: if found velocities in structure
    // file then use them, otherwise the initial velocitie are zero.
    else if (initial_velocity== "auto") { // This is default
        if (numVelocities == numAtoms) {
            context->setVelocities(velocities);
            std::cout << "The initial velocities are loaded from structure file and applied.\n" << std::endl;
        }
        else // don't set the velocitie, then they are zero
            std::cout << "The initial velocities of all atoms are set to zero.\n" << std::endl;
    }
    // 2.2 Set initial velocities from structure file.
    else if (initial_velocity == "file") {
        if (numVelocities !=  numAtoms)
            throw std::runtime_error("ERROR: Set initial_velocity = file while "
                "no velocities found in structure file: " + param.getStr("structure"));
        context->setVelocities(velocities);
        std::cout << "The initial velocities are loaded from structure file and applied.\n" << std::endl;
    }
    // 2.3 Set random initial velocities according to Boltzmann distribution.
    else if (initial_velocity== "Boltzmann") {
        // Get the Boltzmann temperature and random number seed.
        const double Boltzmann_temperature = param.getDouble("Boltzmann_temperature");
        const int Boltzmann_seed = param.getInt("Boltzmann_seed");
        if (Boltzmann_seed == 0) // automatically choosed
            context->setVelocitiesToTemperature(Boltzmann_temperature);
        else // spcified by user
            context->setVelocitiesToTemperature(Boltzmann_temperature, Boltzmann_seed);
        std::cout << "The initial velocities are set to random values chosen from"
            " a Boltzmann distribution at a temperature: " << std::setprecision(2) <<
            param.getDouble("temperature") << " K.\n" << std::endl;
    }
    // 2.4 Set initial velocities to zero with do nothing.
    else if (initial_velocity == "0") {
        std::cout << "The initial velocities of all atoms are set to zero.\n" << std::endl;
    }
    else
        throw std::runtime_error("ERROR: Illegal value of initial_velocity: " + initial_velocity);
    // 3. Set whether to apply constraints to the start configuration and reset shells.
    // If "true", the positions of particles will be updated so that all distance
    // constraints are satisfied, which is default behavior in Gromacs. This also
    // recomputes the locations of all virtual sites.
    // This should be "false" for exact continuation and reruns.
    if (param.getBool("constraint_start")) {
        context->applyConstraints(param.getDouble("constraint_tolerance"));
        std::cout << "The constraints will be applied to the starting coordinates (step 0).\n"
            "(And the locations of all virtual sites will be recomputed, if have.)\n" << std::endl;
    }
}

void HamiltonianOpenMM::printPlatformInfo() const {
    if (param.getStr("platform").empty())
        std::cout << "The platform is not specifed by user and will be chosen by program with default properties.\n";
    const OpenMM::Platform& platform = context->getPlatform();
    std::cout << "The OpenMM Platform that will be used in this simulation is: " <<
        platform.getName() << ".\n";
    std::cout << "And all Platform-specific properties of it that will be used are as following:\n";
    const std::vector<std::string>& properties = platform.getPropertyNames();
    for (int i = 0; i < properties.size(); i++)
        std::cout << properties[i] << ": " <<  platform.getPropertyValue(*context, properties[i]) << "\n";
    std::cout << std::endl;
}

void HamiltonianOpenMM::printPMEParameters() const {
    // Usually, the index of NonbondedForce is 0.
    for (int i = 0; i < system->getNumForces(); i++)
        // if the Force* cannot be converted to NonbondedForce*, the expression is 0
        if (const OpenMM::NonbondedForce* nonbond = dynamic_cast<OpenMM::NonbondedForce*>(const_cast<OpenMM::Force*>(&(system->getForce(i)))))
            // PME or LJPME method
            if (nonbond->getNonbondedMethod() > 3) {
                double alpha;
                int nx, ny, nz;
                nonbond->getPMEParametersInContext(*context, alpha, nx, ny, nz);
                std::cout << "The values of PME parameters that are actually used by Platform are as following:\n";
                std::cout << "The separation parameter is " << alpha << ", and the number of grid points along "
                    "the X/Y/Z axis are " << nx << ", " << ny << ", " << "and " << nz << ", respectively.\n";
                // LJPME method
                if (nonbond->getNonbondedMethod() == 5) {
                    nonbond->getLJPMEParametersInContext(*context, alpha, nx, ny, nz);
                    std::cout << "The values of LJPME parameters that are actually used by Platform are as following:\n";
                    std::cout << "The separation parameter is " << alpha << ", and the number of grid points along "
                        "the X/Y/Z axis are " << nx << ", " << ny << ", " << "and " << nz << ", respectively.\n";
                }
                std::cout << "Note that the above values may be slightly different from standard values calculated based on the Ewald error "
                    "tolerance or specified by user due to restrictions on the allowed grid sizes from some platforms.\n" << std::endl;
                break;
            }
}

int HamiltonianOpenMM::computeSystemDOF() const {
    int systemDOF = 0;
    const int numAtoms = system->getNumParticles();
    for(int i = 0; i < numAtoms; i++)
        if (system->getParticleMass(i) != 0)
            systemDOF += 3;
    systemDOF -= system->getNumConstraints();
    if (param.getBool("remove_COMMotion"))
            systemDOF -= 3;
    return systemDOF;
}


/*
bool HamiltonianOpenMM::getForceFieldForces(int forceFieldIndex, std::vector<OpenMM::Vec3>& perturbForces) {
    
    perturbForces.resize(getNumAtoms());
    
    
    ForceFieldBase* ff = getForceField(forceFieldIndex);
    if (ff == nullptr) {
        std::cerr << "Error: Cannot find ForceField at index " << forceFieldIndex << std::endl;
        return false;
    }
    
    
    ForceFieldPolar* ffPolar = dynamic_cast<ForceFieldPolar*>(ff);
    if (ffPolar == nullptr) {
        std::cerr << "Error: ForceField at index " << forceFieldIndex << " is not of type ForceFieldPolar" << std::endl;
        return false;
    }
   
    std::vector<Vec3> convertedForces;

    convertedForces.reserve(perturbForces.size());
    for (const auto& force : perturbForces) {
        convertedForces.push_back(Vec3(force[0], force[1], force[2]));
    }
    bool success = ffPolar->calcualtePerturbPolarForce(convertedForces);



 
//    bool success = ffPolar->calcualtePerturbPolarForce(perturbForces);
    
    if (!success) {
        std::cerr << "Error: Failed to calculate perturbation forces" << std::endl;
        return false;
    }
    
    return true;
}*/

//version 2
/*
bool HamiltonianOpenMM::getForceFieldForces(int forceFieldIndex, std::vector<OpenMM::Vec3>& perturbForces) {
   
    perturbForces.resize(getNumAtoms());
    
    
    ForceFieldBase* ff = getForceField(forceFieldIndex);
    if (ff == nullptr) {
        std::cerr << "Error: Cannot find ForceField at index " << forceFieldIndex << std::endl;
        return false;
    }
    
    
    ForceFieldPolar* ffPolar = dynamic_cast<ForceFieldPolar*>(ff);
    if (ffPolar == nullptr) {
        std::cerr << "Error: ForceField at index " << forceFieldIndex << " is not of type ForceFieldPolar" << std::endl;
        return false;
    }
    
    
    std::vector<Vec3> convertedForces(getNumAtoms(), Vec3(0, 0, 0));
    
    
    ffPolar->calcualtePerturbPolarForce(convertedForces);
    
    
    for (size_t i = 0; i < convertedForces.size() && i < perturbForces.size(); ++i) {
        perturbForces[i] = OpenMM::Vec3(convertedForces[i][0], 
                                       convertedForces[i][1], 
                                       convertedForces[i][2]);
    }
    
    return true;
}






ForceFieldBase* HamiltonianOpenMM::getForceField(int index) {
   
    HamiltonianForceFieldBase* ffBase = dynamic_cast<HamiltonianForceFieldBase*>(this);
    if (ffBase != nullptr && index >= 0 && index < ffBase->FFList.size()) {
        return ffBase->FFList[index].get();
    }
       
    return nullptr;
}
*/


bool HamiltonianOpenMM::getForceFieldForces(std::vector<OpenMM::Vec3>& perturbForces) {
    int activeState = 0;  
    try {
        activeState = param.getInt("propagate_state");
    } catch (...) {
    }
    
    if (activeState < 0 || activeState >= polarForceFields.size()) {
        std::cerr << "Error: Invalid state index: " << activeState << std::endl;
        return false;
    }
    
    auto ffPolar = polarForceFields[activeState];
    if (!ffPolar) {
        std::cerr << "Error: ForceField at state " << activeState << " is not initialized" << std::endl;
        return false;
    }
    
    std::vector<Vec3> tempForces(DOFn, Vec3(0, 0, 0));
    ffPolar->calcualtePerturbPolarForce(tempForces);
    
    perturbForces.resize(tempForces.size());
    for (size_t i = 0; i < tempForces.size(); ++i) {
        perturbForces[i] = OpenMM::Vec3(tempForces[i][0], 
                                       tempForces[i][1], 
                                       tempForces[i][2]);
    }
    
    return true;
}






