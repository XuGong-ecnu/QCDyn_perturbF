/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author:Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: July. 20, 2022                                                *
 * -------------------------------------------------------------------------- */

#include "HamiltonianRingPolymer.h"

void HamiltonianRingPolymer::init() {
    // * 1 Initialize data members and Hamiltonian.
    HamiltonianBase::init();
    std::cout.precision(12);
    if (PIMD_type != "RPMD" || PIMD_type != "PIMD" || PIMD_type != "CMD")// 
        throw std::runtime_error("ERROR: Unsupported PIMD_type=" + PIMD_type + " for PIMD Hamiltonian.");
    nbeads = param.getInt("nbeads");
    // * 2 Create Stucture object from structure file and get atominfo.
    std::cout << "Start to initialize ForceField objects since " << CurrentTime() << ".\n" << std::endl;
    StructureVec3 structure;
    structure.loadStructure(param.getStr("structure"));
    std::cout << "Loaded structure from file: " << param.getStr("structure") << "\n" << std::endl;
    // * 3 load the force field parameters from topology files and get total masses.
    std::vector<std::string> files; 
    SplitString(files, param.getStr("topology"), ','); 
    // Reset DOFe to the number of topology files (input value is ignored),
    // since number of states (DOFe) must equal to the number of topology files.
    DOFe = param.getInt("DOFe");
    DT = param.getDouble("DT");
    // * 4 Reset DOFe to the number of FFList.
    std::vector<Topology> topologies(DOFe); 
    // * 5 To set the nbeads ForceFieldBase
    for(int i = 0; i < DOFe; ++i) {
        FFList.push_back(std::make_shared<ForceFieldBase>(param));
    }
    for (int i = 0; i < DOFe; ++i) {
        topologies[i].loadTopologyFile(files[i]);
        FFList[i]->setTopologies(topologies[i]);
        FFList[i]->init();
        std::cout << "Loaded topology from file: " << files[i] << std::endl;
    } 
    // * 6 Get mass of each particle from topologies.
    atominfo = structure.getAtomInfo();
    integrator = param.getStr("integrator");
    systemMass = nbeads * topologies[0].computeSystemMass();
    systemDOF = computeSystemDOF();
    // * 7 Resize all vectors in this class and set their elements to 0.
    // We will update them only when we need.
    PE.resize(DOFe, 0);
    Hspr = 0;
    R.resize(DOFn, Vec3());
    V.resize(DOFn, Vec3());
    F.resize(DOFn, Vec3());
    R_RP.resize(nbeads, std::vector<Vec3>(DOFn, Vec3()));
    V_RP.resize(nbeads, std::vector<Vec3>(DOFn, Vec3()));
    F_RP.resize(nbeads, std::vector<Vec3>(DOFn, Vec3()));
    //F_all.resize(DOFe * DOFe , std::vector<Vec3>(DOFn, Vec3()));
    //F_avg.resize(DOFn, Vec3());  
    // * 8 Initialize the positions, velocities and Box vectors in the ForceFieldBase.
    loadInitialStructure(structure, param);
    double temperature = param.getDouble("temperature");
    double beta = 1.0/(kB * temperature );
    // kB   = 1.380649e-23; 
    // double hbar_SI = 1.054571800e-34; // J.s
    std::cout<<"temperature: "<<temperature<<std::endl;
    beta_n = beta / nbeads;
    omega_n =  1E-12 / (beta_n * hbar_SI); // unit: ps-1
    atominfo = structure.getAtomInfo();
    inverseMasses.resize(DOFn, 0);
    masses.resize(DOFn , 0);
}

void HamiltonianRingPolymer::calculatePotentialEnergyandForces(int index, bool includeForces) {
    // Judging the correctness of the calculated state.
    if (index >= 0 && index < DOFe) {
        // Set the R for each state. 
        for (int n = 0; n < nbeads; n++) {
            FFList[index]->setPositions(R_RP[n]);// Get and calculate the potential energy.
            PE[index] = FFList[index]->getPotentialEnergy(includeForces);
            // Get the all states forces.
            F_RP[n] = FFList[index]->getForces();
            // Only propagate state force should be in HamiltonForceFieldBase.cpp.
            int propagate_state = param.getInt("propagate_state");
        }
    }
    else
        throw std::runtime_error("ERROR: Called calculatePotentialEnergyandForces() with wrong state index.");
}

void HamiltonianRingPolymer::updateAllForces() {
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

void HamiltonianRingPolymer::updateDiabaticHamiltonian() {
}

void HamiltonianRingPolymer::updateAdiabaticHamiltonian() {
    throw std::runtime_error("ERROR: Nonadiabatic dynamics simulation within "
        "adiabatic representation is not supported by OpenMM interface.");
}

void HamiltonianRingPolymer::updateQuasiDiabaticHamiltonian() {
    throw std::runtime_error("ERROR: Nonadiabatic dynamics simulation within "
        "quasi-diabatic representation is not supported by OpenMM interface.");
}

double HamiltonianRingPolymer::getPotentialEnergy(int index) {
    return PE[index];
}

double HamiltonianRingPolymer::getKineticEnergy() {
    // calculate the H0 energy
    Hspr = 0;
    // from V(t-0.5dt) to the V(t) which means dtV.
    double dtV[nbeads][DOFn][3];
    double dq2;
    double dx, dy, dz;
    // For velocity Verlet integrator
    std::cout<<"masses[i]:0,1,10   "<<masses[0]<<"  "<<masses[1]<<"   "<<masses[10]<<std::endl;
    std::cout<<"V_RP[0][10][0]:  "<<V_RP[0][10][0]<<std::endl;
    std::cout<<" nbeads  "<<nbeads<<std::endl;
    std::cout<<"DOFn:   "<<DOFn<<std::endl;
    double xf1 = 0;
    for (int i = 0; i < DOFn; i++) {
        //std::cout<<"vvvvv::::  "<<i<<"   "<<V_RP[0][i][0]<<"  "<<V_RP[0][i][1]<<"  "<<V_RP[0][i][2]<<std::endl;
        //std::cout<<"vvvvv2::::  "<<i<<"   "<<V_RP[0][i][0]* V_RP[0][i][0]<<"  "<<V_RP[0][i][1]* V_RP[0][i][1]<<"  "<<V_RP[0][i][2]* V_RP[0][i][2]<<std::endl;
        //std::cout<<"ddddl::::: "<<i<<"  "<<0.5 * masses[i] * (V_RP[0][i][0] * V_RP[0][i][0] + V_RP[0][i][1] * V_RP[0][i][1] + V_RP[0][i][2] *  V_RP[0][i][2])<<std::endl;
        xf1 += 0.5 * masses[i] * (V_RP[0][i][0] * V_RP[0][i][0] + V_RP[0][i][1] * V_RP[0][i][1] + V_RP[0][i][2] *  V_RP[0][i][2]);         
        //std::cout<<"xf1:::::   "<<xf1<<std::endl;
    }
    std::cout<<"TTTTTTTTTTTTT:   "<<2.0 * xf1 / (3 * DOFn) / Rbar;
    if (integrator == "velocityVerlet") {
        for (int j = 0; j < nbeads; j++)
            for (int i = 0; i < DOFn; i++) {
                Hspr += 0.5 * masses[i] * (V_RP[j][i][0] * V_RP[j][i][0] + V_RP[j][i][1] * V_RP[j][i][1] + V_RP[j][i][2] *  V_RP[j][i][2]);
                //if (j == 0) {
                //    dx = R_RP[j][i][0] - R_RP[nbeads-1][i][0];
                //    dy = R_RP[j][i][1] - R_RP[nbeads-1][i][1];
                //    dz = R_RP[j][i][2] - R_RP[nbeads-1][i][2];
                //}
                //else {
                //    dx = R_RP[j][i][0] - R_RP[j-1][i][0];
                //    dy = R_RP[j][i][1] - R_RP[j-1][i][1];
                //    dz = R_RP[j][i][2] - R_RP[j-1][i][2];
                //}
                //dq2 = dx * dx + dy * dy + dz * dz;
                //Hspr += 0.5 * masses[i] * omega_n * omega_n * dq2;
        }
    }
    // For leap frog integrator
    else if (integrator == "leapfrog") {
        for (int j = 0; j < nbeads; j++) {
            for (int i = 0; i < DOFn; i++) {
                for (int k = 0; k < 3; k++) {
                    dtV[j][i][k] = V_RP[j][i][k] + 0.5 * DT * F_RP[j][i][k] / masses[i];
                }
                Hspr += 0.5 * masses[i] * (dtV[j][i][0] * dtV[j][i][0] + dtV[j][i][1] * dtV[j][i][1] + dtV[j][i][2] * dtV[j][i][2]);
                //if (j == 0) {
                //    dx = R_RP[0][i][0] - R_RP[nbeads-1][i][0];
                //    dy = R_RP[0][i][1] - R_RP[nbeads-1][i][1];
                //    dz = R_RP[0][i][2] - R_RP[nbeads-1][i][2];
                //}
                //else {
                //    dx = R_RP[j][i][0] - R_RP[j-1][i][0];
                //    dy = R_RP[j][i][1] - R_RP[j-1][i][1];
                //    dz = R_RP[j][i][2] - R_RP[j-1][i][2];
                //}
                //dq2 = dx * dx + dy * dy + dz * dz;
                //Hspr += 0.5 * masses[i] * omega_n * omega_n * dq2;
            }
        }
    }
    return  Hspr;  
}

// To Do
double HamiltonianRingPolymer::getDiabaticCoupling(int i, int j) {
    return 0.0;
}

void HamiltonianRingPolymer::loadInitialStructure(const StructureVec3& structure, const Parameters& param) {
    const std::vector<Vec3>& positions = structure.getPositions();
    const std::vector<Vec3>& velocities = structure.getVelocities();
    const int numAtoms = positions.size();
    const int numVelocities = velocities.size();
    DOFe = param.getInt("DOFe");
    // 1. Set initial positions.
    R = positions;
    for (int i = 0; i < nbeads; i++)
        R_RP[i] = R;
    // 2. Set initial velocities. The allowed value: 0, file, Boltzmann, auto.
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
            for (int i = 0; i < nbeads; i++)
                V_RP[i] = velocities;
            std::cout << "The initial velocities are loaded from structure file and applied.\n" << std::endl;
        }
    }
    // 2.3 Set random initial velocities according to Boltzmann distribution.
    else if (initial_velocity == "Boltzmann") {
        //setBoltzmannVelocities();
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
        for (int i = 0; i < DOFe; i++) {
            setPeriodicBoxVectors(i, a, b, c);
        }
        std::cout << "The periodic boundary conditions (PBC) in all directions will be used in the simulation.\n";
        std::cout << "And the box vectors (in nm) are loaded from structure file and reduced to the form required: \n";
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

int HamiltonianRingPolymer::computeSystemDOF() const {
    int systemDOF = 3 * DOFn * nbeads;
    return systemDOF;
}

void HamiltonianRingPolymer::setPeriodicBoxVectors(int index, Vec3& a, Vec3& b, Vec3& c) {
    FFList[index]->setPeriodicBoxVectors(a, b, c);
}

double HamiltonianRingPolymer::getPeriodicBoxVolume(int index) const {
    // To judge the PBC type an calculate the box volume.  
    std::string PBC_type = param.getStr("PBC");
    if (PBC_type == "xyz") {
        return FFList[index]->getPeriodicBoxVolume();
        }
    else 
        return 0.0;
}

const std::vector<std::string>& HamiltonianRingPolymer::getAtomInfo(){
    return atominfo;
}

void HamiltonianRingPolymer::getPeriodicBoxVectors(Vec3& a, Vec3& b, Vec3& c) {
    a = periodicBoxVectors[0];
    b = periodicBoxVectors[1];
    c = periodicBoxVectors[2];
}

void HamiltonianRingPolymer::getPositions(std::vector<std::vector<Vec3>>& positions) {
    positions = R_RP;
}

void HamiltonianRingPolymer::getVelocities(std::vector<std::vector<Vec3>>& velocities) {
    velocities = V_RP;
}

void HamiltonianRingPolymer::getForces(std::vector<std::vector<Vec3>>& forces) {
    forces = F_RP;
}

double HamiltonianRingPolymer::getSystemMass() const {
    return systemMass;
}

int HamiltonianRingPolymer::getSystemDOF() const {
    return systemDOF;
}

const std::vector<double>& HamiltonianRingPolymer::getEnergyGroups(int index){
    return FFList[index]->getEnergyParts();
}
