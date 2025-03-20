/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * First created: Nov. 23, 2021                                               *
 * Last updated: Nov. 30, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "HamiltonianLVC.h"

void HamiltonianLVC::init() {
    // * 0. Initialize data members and Hamiltonian.
    HamiltonianModelBase::init();
    if (model_type != "LVC")
        throw std::runtime_error("ERROR: Undefined model_type=" + model_type + " for LVC model.");
    if (Condon_approximation)
        throw std::runtime_error("ERROR: Please set Condon_approximation to false for LVC model.");
    // For LVC model, number of normal modes, N = DOFn.
    N = DOFn;

    // * 1. Initialize and get model parameters by building or loading from file
    omega.resize(N, 0); // frequecy
    shift.resize(N, 0); // used for initial nuclear sampling
    c.resize(DOFe, std::vector<double>(N, 0)); // system-bath coupling (the linear shifts d_j in literature)
    gamma.resize((DOFe-1)*DOFe/2, std::vector<double>(N, 0)); // linear coupling coefficient
    initializeModelParameters();
    getSamplingShifts();
    // * 2. Compute and report reaction free energy and reorganization energy
    if (DOFe == 2) { // currently, works for two-state LVC model only
        const double DeltaE = getReactionFreeEnergy();
        const double Er = getReorganizationEnergy();
        std::cout << "The donor-to-acceptor reaction free energy (DeltaE) for the LVC model is " <<
            std::setprecision(6) << DeltaE << " au.\n";
        std::cout << "The reorganization energy (Er) for the LVC model is " <<
            std::setprecision(6) << Er << " au.\n";
        std::cout << "Er > |DeltaE| and Er < |DeltaE| correspond to the Marcus normal "
            "and inverted regions, respectively.\n";
        if (Er > abs(DeltaE))
            std::cout << "So, this LVC model corresponds to the Marcus normal region.\n";
        else
            std::cout << "So, this LVC model corresponds to the Marcus inverted region.\n";
        std::cout << std::endl;
    }
}

void HamiltonianLVC::buildModelParameters() {
    throw std::runtime_error("ERROR: Contructing LVC model parameters from spectral density is not supported yet.");
}

void HamiltonianLVC::loadModelParameters(const std::string& loadfile) {
    std::vector<std::vector<double>*> data = {&omega};
    for (int i = 0; i < DOFe; ++i)
        data.emplace_back(&c[i]);
    for (int i = 0; i < DOFe; ++i)
        for (int j = i+1; j < DOFe; ++j)
            data.emplace_back(&gamma[(2*DOFe-i-1)*i/2 + j-i - 1]);
    LoadDataFile(loadfile, data);
}

void HamiltonianLVC::saveModelParameters(const std::string& savefile) {
    std::vector<std::vector<double>*> data = {&omega};
    std::vector<std::string> headers = {"omega"};
    for (int i = 0; i < DOFe; ++i) {
        data.emplace_back(&c[i]);
        headers.emplace_back("c" + std::to_string(i));
    }
    for (int i = 0; i < DOFe; ++i)
        for (int j = i+1; j < DOFe; ++j) {
            data.emplace_back(&gamma[(2*DOFe-i-1)*i/2 + j-i - 1]);
            headers.emplace_back("gamma" + std::to_string(i) + std::to_string(j));
        }
    SaveDataFile(savefile, headers, data);
}

double HamiltonianLVC::getPotentialEnergy(int index) {
    if (index >= DOFe || index < 0)
        throw std::runtime_error("ERROR: State index out of scope.");
    double PE = 0.0;
    for (int j = 0; j < DOFn; ++j) {
        PE += 0.5 * omega[j] * omega[j] * R[j] * R[j];
        PE += c[index][j] * R[j];
    }
    return PE + epsilon[index];
}

void HamiltonianLVC::getForces(int i, int j, std::vector<double>& F) {
    if (i >= DOFe || i < 0 || j >= DOFe || j < 0)
        throw std::runtime_error("ERROR: Called getForces() with wrong state index.");
    F.resize(DOFn, 0);
    if (i == j) // froces of each state from diagonal nuclear Hamiltonian
        for (int k = 0; k < DOFn; ++k)
            F[k] = -omega[k] * omega[k] * R[k] - c[i][k];
    else { // forces from off-diagonal Hamiltonian (diabatic coupling)
        // H_ji = H_ij = γ_ij * R_k, where γ_ij is the linear coupling coefficient.
        // The index of γ_ij in vector gamma is (2*DOFe-i-1)*i/2 + j-i - 1, i < j
        // And only the i < j case is stored since γ_ji = γ_ij.
        // So force is F[k] = -dH_ij/dR = - γ_ij[k]
        int index = i < j ? ((2*DOFe-i-1)*i/2 + j-i - 1) : ((2*DOFe-j-1)*j/2 + i-j - 1);
        for (int k = 0; k < DOFn; ++k)
            F[k] = -gamma[index][k];
    }
}

double HamiltonianLVC::getNonCondonCoupling(int i, int j) {
    if (i >= DOFe || i < 0 || j >= DOFe || j < 0 || i == j)
        throw std::runtime_error("ERROR: Called getNonCondonCoupling() with wrong state index.");
    // H_ji = H_ij = Σγ_ij * R_k, where γ_ij is the linear coupling coefficient.
    // The index of γ_ij in vector gamma is (2*DOFe-i-1)*i/2 + j-i - 1, i < j
    // And only the i < j case is stored since γ_ji = γ_ij.
    int index = i < j ? ((2*DOFe-i-1)*i/2 + j-i - 1) : ((2*DOFe-j-1)*j/2 + i-j - 1);
    double coupling = 0;
    for (int k = 0; k < DOFn; ++k)
        coupling += gamma[index][k] * R[k];
    return coupling;
}

double HamiltonianLVC::getReactionFreeEnergy() {
    // Formula see Eq. 3 in Reference: J. Chem. Theory Comput. 16, 4479−4488 (2020).
    // Note, for donor-to-acceptor ΔE, a negative sign should be added.
    if (DOFe > 2)
        throw std::runtime_error("ERROR: The reaction free energy for multi-state "
            "LVC model is not defined yet.");
    double DeltaE = 0;
    for (int k = 0; k < DOFn; ++k)
        DeltaE += (c[0][k] * c[0][k] - c[1][k] * c[1][k]) / (omega[k] * omega[k]);
    return 0.5 * DeltaE + epsilon[1] - epsilon[0];
}

double HamiltonianLVC::getReorganizationEnergy() {
    // Formula see Eq. 3 in Reference: J. Chem. Theory Comput. 16, 4479−4488 (2020).
    if (DOFe > 2)
        throw std::runtime_error("ERROR: The reorganization energy for multi-state "
            "LVC model is not defined yet.");
    double Er = 0;
    for (int k = 0; k < DOFn; ++k)
        Er += (c[0][k] - c[1][k]) * (c[0][k] - c[1][k]) / (omega[k] * omega[k]);
    return 0.5 * Er;
}

void HamiltonianLVC::getSamplingShifts() {
    // 1. Load the equilbrium shifts of R for nuclear sampling from file.
    const std::string loadfile = param.getStr("shift_load");
    if (!loadfile.empty()) {
        // Requirments of the shifts file (.csv or .dat file):
        // 1. First line is headers: index,shift
        // 2. The number of lines is DOFn+1.
        // 3. The unit is au.
        std::vector<std::vector<double>*> data = {&shift};
        LoadDataFile(loadfile, data);
    }
    // 2. Get shifts from the system-bath coupling c_j of specified sample state.
    // The minimal point of PES of state i is at (Rj+c(i)_j/omega_j^2), then
    // shift is c_j/omega_j^2. If sample_state = -1, the shifts are 0 (default).
    else {
        const int sample_state = param.getInt("sample_state");
        if (sample_state >= 0 || sample_state < DOFe)
            for (int k = 0; k < DOFn; ++k)
                shift[k] = c[sample_state][k] / (omega[k]*omega[k]);
        else if (sample_state != -1)
            throw std::runtime_error("ERROR: Illegal sample_state=" + std::to_string(sample_state));
    }
}