/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: Jan. 15, 2022                                                *
 * -------------------------------------------------------------------------- */

#include "HamiltonianForceFieldBase.h"

void HamiltonianForceFieldBase::init() {
    // * 1 Initialize data members and Hamiltonian.
    HamiltonianBase::init();
    if (system_type != "allatom")
        throw std::runtime_error("ERROR: Unsupported system_type=" + system_type + " for ForceField Hamiltonian.");
    if (allatom_type != "CustomForceField")
        throw std::runtime_error("ERROR: Unsupported allatom_type=" + allatom_type + " for ForceField Hamiltonian.");
    if (representation != "diabatic")
        throw std::runtime_error("ERROR: Only diabatic representation is allowed for ForceField Hamiltonian.");
    // The non Condon_approximation has supported for ForceField Hamiltonian applied to the C60 system.
    // added by Xiaofang
    //if (!Condon_approximation)
    //    throw std::runtime_error("ERROR: Non-Condon simulation is not yet supported for ForceField Hamiltonian.");

    // * 2 Check if it is a restart simulation. The checkpoint file created by
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

    // * 3 Create Stucture object from structure file and get atominfo.
    std::cout << "Start to initialize ForceField objects since " << CurrentTime() << ".\n" << std::endl;
    StructureVec3 structure;
    structure.loadStructure(param.getStr("structure"));
    std::cout << "Loaded structure from file: " << param.getStr("structure") << "\n" << std::endl;

    // * 4 load the force field parameters from topology files and get total masses.
    std::vector<std::string> files; 
    SplitString(files, param.getStr("topology"), ','); 
    // Reset DOFe to the number of topology files (input value is ignored),
    // since number of states (DOFe) must equal to the number of topology files.
    DOFe = param.getInt("DOFe");

    // * 5 Reset DOFe to the number of FFList.
    std::vector<Topology> topologies(DOFe); 
    forcefield_type = param.getStr("forcefield_type");
    static const bool polar_obs = param.getBool("polar_obs");
    static const bool perturb = param.getBool("perturb");
    if (forcefield_type == "nonPolar") {
        if (polar_obs || perturb) {
            for(int i = 0; i < DOFe; ++i) {
                FFList.push_back(std::make_shared<ForceFieldPolar>(param));
            }
        }
        else {
            for(int i = 0; i < DOFe; ++i) {
                FFList.push_back(std::make_shared<ForceFieldBase>(param));
            }
        }
    }
    else if (forcefield_type == "Polar") {
        for(int i = 0; i < DOFe; ++i){
            FFList.push_back(std::make_shared<ForceFieldPolar>(param));
        }
    }
    else 
        throw std::runtime_error("ERROR: Unsupported forcefield_type: " + forcefield_type);
    for (int i = 0; i < DOFe; ++i) {
        topologies[i].loadTopologyFile(files[i]);
        FFList[i]->setTopologies(topologies[i]);
        FFList[i]->init();
        std::cout << "Loaded topology from file: " << files[i] << std::endl;
    } 
    if (!Condon_approximation) {
        HaOffDiagonal = std::make_shared<HamiltonianOffDiagonal>(param);
        HaOffDiagonal->setTopologies(topologies[0]);
        HaOffDiagonal->init();
    }
    // * 6 get the system DOFn, integrater, DT, temperature parameters.
    DOFn = param.getInt("DOFn");
    integrator = param.getStr("integrator");
    DT = param.getDouble("DT");
    temperature = param.getDouble("temperature");
    if (integrator == "leapfrog") {
        std::cout << "The ForceField VerletIntegrator using the leapfrog Verlet algorithm will be used for simulation.\n" << std::endl;
    }
    else if (integrator == "velocityVerlet") {
        std::cout << "The velocity Verlet Integrator will be used for simulation.\n" << std::endl;
    }
    else
        throw std::runtime_error("ERROR: Unsupported integrator: " + integrator);
    // * 7 Set random number seed for initial velocities sampling.
    const int nucl_seed = param.getInt("nucl_seed");
    if (nucl_seed == 0) { // Use system-time seed for real simulation.
        std::random_device rd;
        nucl_gen.seed(rd());
    }
    else // Use deterministic seed for debug.
        nucl_gen.seed(nucl_seed);
    // * 8 Get mass of each particle from topologies.
    systemMass = topologies[0].computeSystemMass();
    systemDOF = computeSystemDOF();
    // * 9 Resize all vectors in this class and set their elements to 0.
    // We will update them only when we need.
    PE.resize(DOFe, 0);
    KE = 0;
    R.resize(DOFn, Vec3());
    V.resize(DOFn, Vec3());
    F.resize(DOFn, Vec3());
    R_eq.resize(DOFn, Vec3());
    V_eq.resize(DOFn, Vec3());
    F_eq.resize(DOFn, Vec3());
    F_all.resize(DOFe*DOFe, std::vector<Vec3>(DOFn, Vec3()));
    F_avg.resize(DOFn, Vec3());
    // * 10 Get mass of each particle from CustomForceField System.
    atominfo = structure.getAtomInfo();
    inverseMasses.resize(DOFn, 0);
    masses.resize(DOFn, 0);
    for (int j = 0; j < DOFn; j++) {
        masses[j] = FFList[0]->getmass(j);
        inverseMasses[j] = masses[j] == 0 ? 0 : (1.0/masses[j]);
    }
    // * 11 Initialize the positions, velocities and Box vectors in the ForceFieldBase.
    if (!isRestart) 
         loadInitialStructure(structure, param);
    // if this is a restart job, the simulation state will load from checkpoint file.
    else {
        //loadCheckpoint(restart_file);
        std::cout << "Load simulation state from checkpoint file: " << restart_file << "\n" << std::endl;
    }
}

void HamiltonianForceFieldBase::saveTheEqConfig() {
    R_eq = R;
    V_eq = V;
    F_eq = F;
}

void HamiltonianForceFieldBase::setTheEqConfig() {
    R = R_eq;
    V = V_eq;
    F = F_eq;
}

void HamiltonianForceFieldBase::updatePiTensor(int index, std::vector<double>& Pi_tensor) {
    FFList[index]->calculatePiTensor();
    FFList[index]->getPiTensor(Pi_tensor); 
}

void HamiltonianForceFieldBase::updateDipoleRleftAright(int index, double& cosleft,
    double& cosright, int& numLeft, int& numRight, double& dirLx, double& dirRx) {
    FFList[index]->calculateDipoleRleftAright();
    FFList[index]->getDipoleRleftAright(cosleft, cosright, numLeft, numRight, dirLx, dirRx); 
}

void HamiltonianForceFieldBase::updateMuTensor(int index, Vec3& Mu_tensor) {
    FFList[index]->calculateInducedDipole();
    FFList[index]->getMuTensor(Mu_tensor);
}

void HamiltonianForceFieldBase::getPiTensor(int index, std::vector<double>& Pi_tensor) {
    FFList[index]->getPiTensor(Pi_tensor); 
}

void HamiltonianForceFieldBase::getDipoleRleftAright(int index, double& cosleft,
    double& cosright, int& numLeft, int& numRight, double& dirLx, double& dirRx) {
    FFList[index]->getDipoleRleftAright(cosleft, cosright, numLeft, numRight, dirLx, dirRx); 
}


void HamiltonianForceFieldBase::getMuTensor(int index, Vec3& Mu_tensor) {
    FFList[index]->getMuTensor(Mu_tensor); 
}

void HamiltonianForceFieldBase::updateAllForces() {
    // Calculate the one state[index] especially for propagate_state enery and force 
    // for propagating.
    int propagate_state = param.getInt("propagate_state");
    // Calculate the propagate_state energy and force.
    calculatePotentialEnergyandForces(propagate_state, true);
    for (int i = 0; i < DOFe; i++) {
        if (i == propagate_state) continue;
        calculatePotentialEnergyandForces(i, true);
    }
}


void HamiltonianForceFieldBase::calculatePotentialEnergyandForces(int index, bool includeForces) {
    // Judging the correctness of the calculated state.
    if (index >= 0 && index < DOFe) {
        // Set the R for each state.
        FFList[index]->setPositions(R);
        // Get and calculate the potential energy.
        PE[index] = FFList[index]->getPotentialEnergy(includeForces);
        // Get the all states forces.
        F_all[index * DOFe + index] = FFList[index]->getForces();
        // Only propagate state force should be in HamiltonForceFieldBase.cpp.
        int propagate_state = param.getInt("propagate_state");
        if (propagate_state == index) F = FFList[index]->getForces();
    }
    else
        throw std::runtime_error("ERROR: Called calculatePotentialEnergyandForces() with wrong state index.");
}

double HamiltonianForceFieldBase::getPotentialEnergy(int index){
    return PE[index];
}

const std::vector<double>& HamiltonianForceFieldBase::getEnergyGroups(int index){
    return FFList[index]->getEnergyParts();
}

const std::vector<std::vector<Vec3>>& HamiltonianForceFieldBase::getForceGroups(int index){
    return FFList[index]->getForceParts();
}

void HamiltonianForceFieldBase::getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) {
    a = periodicBoxVectors[0];
    b = periodicBoxVectors[1];
    c = periodicBoxVectors[2];
}

void HamiltonianForceFieldBase::updateEffectiveForces() {
    throw std::runtime_error("ERROR: Unsupported effective forces [to be constructed].");
}

void HamiltonianForceFieldBase::updateDiabaticHamiltonian() {
    // This function is very simillar as HamiltonianModelBase::updateDiabaticHamiltonian().
    // Save current H as H_old before update.
    //HaOffDiagonal.setPositions(R);
    H_old = H;
    Heff_old = Heff;
    // Remember to reset H_avg to 0.
    // the average potential energy of all states, which is removed in
    // effective energy matrix (Heff).
    double H_avg = 0.0;
    HaOffDiagonal->setPositions(R);
    for (int i = 0; i < DOFe; i++) {
        FFList[i]->setPositions(R);
        // Compute potential energy (in kj/mol) each state
        PE[i] = FFList[i]->getPotentialEnergy(true);
        // The total potential energy should add energy correction
        // Note the unit of Hamiltonian matrix is au.
        // H[i][i] = PE[i] * kj2au  + epsilon[i];
        // This is for the C60 system degree 8 vdw model.
        H[i][i] = PE[i] + epsilon[i];
        // Compute H_avg
        H_avg += H[i][i] / DOFe;
        // Get forces of current state, which is diganol element of force matrix
        F_all[i * DOFe + i] = FFList[i]->getForces();
        // Get the off-diagonal element of Hamiltonian matrix and forces
        // When using Condon approximation, coupling is constant and force is 0.
        // So, we only need to compute them in the non-Condon case.
        // Here, the Hamiltonian is a real symmetrix matrix, i.e., H_ij = H_ji
        // and F[i][j] = F[j][i] = - dH_ij/dR, i != j.
        // Note, the unit of force is kj/mol/nm.
        // TODO: undefined for all-atom
        // Test
        if (!Condon_approximation)
            for (int j = i + 1; j < DOFe; ++j) {
                H[j][i] = H[i][j] = HaOffDiagonal->getDiabaticCoupling(i, j);
                F_all[j * DOFe + i] = F_all[i * DOFe + j]=HaOffDiagonal->getDiabaticCouplingforce(i, j);;
            }
    }
    // Get effective Hamiltioan matrix (Heff) with removing H_avg.
    Heff = H;
    for (int i = 0; i < DOFe; i++)
        Heff[i][i] -= H_avg;
    // Reset all elements to 0, then compute average forces (F_avg).
    std::fill(F_avg.begin(), F_avg.end(), Vec3());
    for (int j = 0; j < DOFn; j++) {
        for (int i = 0; i < DOFe; i++)
            F_avg[j] += F_all[i * DOFe + i][j];
        F_avg[j] /= DOFe;
    }
}

void HamiltonianForceFieldBase::updateAdiabaticHamiltonian() {
    throw std::runtime_error("ERROR: Nonadiabatic dynamics simulation within "
        "adiabatic representation is not supported by ForceField interface.");
}

void HamiltonianForceFieldBase::updateQuasiDiabaticHamiltonian() {
    throw std::runtime_error("ERROR: Nonadiabatic dynamics simulation within "
        "quasi-diabatic representation is not supported by ForceField interface.");
}

double HamiltonianForceFieldBase::getKineticEnergy() {
    // To get the kinetic energy of the system at time t
    KE = 0;
    // from V(t-0.5dt) to the V(t) which means dtV.
    double dtV[DOFn][3];
    // For velocity Verlet integrator
    if (integrator == "velocityVerlet") {
        for (int i = 0; i < DOFn; i++) {
            KE += 0.5 * masses[i] * (V[i][0]*V[i][0]+V[i][1]*V[i][1]+V[i][2]*V[i][2]);
        }
    }
    // For leap frog integrator
    else if (integrator == "leapfrog") {
        for (int i = 0; i < DOFn; i++) {
            for (int j = 0; j < 3; j++) {
                dtV[i][j] = V[i][j] + 0.5 * DT * F[i][j] / masses[i];
            }
            KE += 0.5 * masses[i] * (dtV[i][0]*dtV[i][0]+dtV[i][1]*dtV[i][1]+dtV[i][2]*dtV[i][2]);
        }
    } 
    return  KE;  
}

double HamiltonianForceFieldBase::getPeriodicBoxVolume(int index) const {
    // To judge the PBC type an calculate the box volume. 
    std::string PBC_type = param.getStr("PBC");
    if (PBC_type == "xyz") {
        return FFList[index]->getPeriodicBoxVolume();
    }
    else
        return 0.0;
}

void HamiltonianForceFieldBase::getPositions(std::vector<Vec3>& positions) {
    positions = R;
}

void HamiltonianForceFieldBase::getVelocities(std::vector<Vec3>& velocities) {
    velocities = V;
}

void HamiltonianForceFieldBase::setPositions(std::vector<Vec3>& positions) {
    for (int index = 0; index < DOFe; index++) {
        FFList[index]->setPositions(positions); 
    }
    if (!Condon_approximation) {
        HaOffDiagonal->setPositions(positions);  
    }
    R = positions;
}

void HamiltonianForceFieldBase::setVelocities(std::vector<Vec3>& velocities) {
    V = velocities;
}

void HamiltonianForceFieldBase::getForces(std::vector<Vec3>& forces) {
    forces = F;
}

int HamiltonianForceFieldBase::getSystemDOF() const {
    return systemDOF;
}

int HamiltonianForceFieldBase::computeSystemDOF() const {
    int systemDOF = 3 * DOFn;
    return systemDOF;
}

void HamiltonianForceFieldBase::setPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) {
    for (int i = 0; i < DOFe; i++) {
        FFList[i]->setPeriodicBoxVectors(a, b, c);
    }
    if (!Condon_approximation) {
        HaOffDiagonal->setPeriodicBoxVectors(a,b,c);
    }
}

double HamiltonianForceFieldBase::getSystemMass() const {
    return systemMass;
}

// PFF
const std::vector<std::string>& HamiltonianForceFieldBase::getAtomInfo(){
    return atominfo;
}

void HamiltonianForceFieldBase::setBoltzmannVelocities() { 
    // This function is to construct the intial velocities for every atoms.
    std::vector<double> Sum_V(3, 0);
    std::vector<double> sigma_V(DOFn, 0);
    // Boltzmann velocities distribution
    std::normal_distribution<double> normal_dist(0.0, 1.0);
    temperature = param.getDouble("Boltzmann_temperature");
    int nucl_seed = param.getInt("nucl_seed");
    nucl_gen.seed(nucl_seed);
    // Construct the Boltzmann velocities temperature distributions.
    double Ke = 0;
    double Kt = 0;
    for (int j = 0 ; j < DOFn ; j++) {  
        sigma_V[j] = sqrt(Rbar  * temperature / masses[j]);
        V[j][0] = normal_dist(nucl_gen) * sigma_V[j];
        V[j][1] = normal_dist(nucl_gen) * sigma_V[j];
        V[j][2] = normal_dist(nucl_gen) * sigma_V[j];
        Sum_V[0] += V[j][0];
        Sum_V[1] += V[j][1];
        Sum_V[2] += V[j][2];
    }
    Sum_V[0] /= DOFn;
    Sum_V[1] /= DOFn;
    Sum_V[2] /= DOFn;
    // Make the velocity 0 in all directions. 
    for (int j = 0 ; j < DOFn ; j++) {
        V[j][0] -= Sum_V[0];
        V[j][1] -= Sum_V[1];
        V[j][2] -= Sum_V[2];
    } 
    // Use the Kt parameters to get the right temperature distribution velocities.
    for (int j = 0; j < DOFn; j++) {
        Ke += 0.5 * masses[j] * ((V[j][0])*(V[j][0]) + (V[j][1])*(V[j][1])  +(V[j][2])*(V[j][2]));
    }
    Kt = sqrt(0.5 * Rbar * systemDOF * temperature/Ke);
    for (int j = 0 ; j < DOFn ; j++) {
        V[j][0] *= Kt;
        V[j][1] *= Kt;
        V[j][2] *= Kt;
    } 
}

void HamiltonianForceFieldBase::loadInitialStructure(const StructureVec3& structure, const Parameters& param) {
    const std::vector<Vec3>& positions = structure.getPositions();
    const std::vector<Vec3>& velocities = structure.getVelocities();
    const int numAtoms = positions.size();
    const int numVelocities = velocities.size();
    DOFe = param.getInt("DOFe");
    // 1. Set initial positions.
    R = positions;
    // 2. Set initial velocities. The allowed values: 0, file, Boltzmann, auto.
    const std::string initial_velocity = param.getStr("initial_velocity");
    // 2.1 Set initial velocities automatically: if found velocities in structure
    // file then use them, otherwise the initial velocitie are zero.
    if (initial_velocity== "zero") { // This is default
        std::cout << "The initial velocities of all atoms are set to zero.\n" << std::endl;
    }
    // 2.2 Set initial velocities from structure file.
    else if (initial_velocity == "file") {
        if (velocities.size() != DOFn) 
            throw std::runtime_error("ERROR: The initial velocities are not found in file.");
        else {
            V = velocities;
            std::cout << "The initial velocities are loaded from structure file and applied.\n" << std::endl;
        }
    }
    // 2.3 Set random initial velocities according to Boltzmann distribution.
    else if (initial_velocity == "Boltzmann") {
        setBoltzmannVelocities();
    }
    // 3. Set periodic boundary conditions (PBC)  for this System.
    std::string PBC_type = param.getStr("PBC");
    if (PBC_type == "xyz") { // PBC in all directions
        // Box vectors is loaded from structure file and has been converted to
        // the ForceFieldBase reduced form (if nessesary) in structure object.
        Vec3 a, b, c;
        structure.getBoxVectors(a, b, c);
        periodicBoxVectors[0] = a;
        periodicBoxVectors[1] = b;
        periodicBoxVectors[2] = c;
        // It will affect any Force added to the System that uses periodic boundary conditions.
        setPeriodicBoxVectors(a, b, c);
        std::cout << "The periodic boundary conditions (PBC) in all directions will be used in the simulation.\n";
        std::cout << "And the box vectors (in nm) are loaded from structure file and reduced to the form required by ForceField: \n";
        std::cout << "a = (" << std::fixed << std::setprecision(2) <<
            a[0] << ", " << a[1] << ", " << a[2] << "), b = (" <<
            b[0] << ", " << b[1] << ", " << b[2] << "), c = (" <<
            c[0] << ", " << c[1] << ", " << c[2] << ").\n" << std::endl;
    }
    else if (PBC_type == "none") // any Force added to the System won't use PBC
        std::cout << "The periodic boundary conditions (PBC) will be not used in the simulation.\n" << std::endl;
    else
        throw std::runtime_error("ERROR: Unsupported periodic boundary conditions type, PBC=" + PBC_type);
}

void HamiltonianForceFieldBase::setTrueMdPerturb() {
    for (int i = 0; i < DOFe; i++)
        FFList[i]->setTrueMdPerturb();
}

void HamiltonianForceFieldBase::setFalseMdPerturb() {
    for (int i = 0; i < DOFe; i++)
        FFList[i]->setFalseMdPerturb();
}

void HamiltonianForceFieldBase::setPositiveDirection() {
    currentStep = 0;
    direction = true;
    for (int i = 0; i < DOFe; i++)
        FFList[i]->setPositiveDirection();
}

void HamiltonianForceFieldBase::setNegativeDirection() {
    currentStep = 0;
    direction = false;
    for (int i = 0; i < DOFe; i++)
        FFList[i]->setNegativeDirection();
}

void HamiltonianForceFieldBase::addTheOneStep() {
    currentStep += 1;
}