/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Nov. 23, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "HamiltonianFMO.h"

void HamiltonianFMO::init() {
    // * 0. Initialize data members and Hamiltonian.
    HamiltonianModelBase::init();
    if (model_type != "FMO" && model_type != "FMO7" && model_type != "FMO4" && model_type != "FMO5")
        throw std::runtime_error("ERROR: Undefined model_type=" + model_type + " for FMO model.");
    if (DOFn % DOFe != 0)
        throw std::runtime_error("ERROR: For FMO model, total nuclear DOFn should be divided exactly by DOFe.");
    // For FMO model, number of normal modes, N = DOFn/DOFe.
    N = DOFn/DOFe;
    // The parameters for standard 7-state or 4- and 5-state reduced-dimensional
    // FMO model can found in J. Chem. Phys. 150, 104101 (2019) Appendix
    // For these FMO models, the couplings and epsilons can not be modified by
    // user, i.e., the input value of gamma_DA and epsilon will be ignored.
    if (model_type != "FMO") {
        if (std::stoi(model_type.substr(3)) != DOFe)
            throw std::runtime_error("ERROR: For model_type=" + model_type + ", DOFe should be equal to " + model_type.substr(3));
        // Standard 7-state FMO (FMO7) model. (in cm-1)
        const Real_Matrix FMO7 = {{12410, -87.7,   5.5,  -5.9,   6.7, -13.7,  -9.9},
                                  {-87.7, 12530,  30.8,   8.2,   0.7,  11.8,   4.3},
                                  {  5.5,  30.8, 12210, -53.5,  -2.2,  -9.6,   6.0},
                                  { -5.9,   8.2, -53.5, 12320, -70.7, -17.0, -63.3},
                                  {  6.7,   0.7,  -2.2, -70.7, 12480,  81.1,  -1.3},
                                  {-13.7,  11.8,  -9.6, -17.0,  81.1, 12630,  39.7},
                                  { -9.9,   4.3,   6.0, -63.3,  -1.3,  39.7, 12440}};
        const std::string init_state = param.getStr("init_state");
        // A reasonable 4-state FMO model is given by cancelling the state 5-7,
        // sicne states 5–7 are relatively unimportant when the init_state is 1.
        // So, for FMO4 model, initial state must be state 1 (index=0 in FMO4)
        if (model_type == "FMO4")
            if (init_state != "0,0")
                throw std::runtime_error("ERROR: For reduced FMO4 model, a reasonable dynamcis is starting from init_state=0,0.");
        // A reasonable 5-state FMO model is given by cancelling the state 1-2,
        // since the dynamics does not significantly involve state 1-2, when
        // the initial excitation is on the state 6 (the state 4 in FMO5).
        // So, for FMO5 model, initial state must be state 6 (index=3 in FMO5)
        else if (model_type == "FMO5")
            if (init_state != "3,3")
                throw std::runtime_error("ERROR: For reduced FMO5 model, a reasonable dynamcis is starting from init_state=3,3.");
        // shift is the start index of state in FMO7 matrix
        // For FMO4 and FMO7, it is 0, and for FMO5 it is 2.
        const int shift = model_type == "FMO5" ? 2 : 0;
        // Set epsilon and couplings and convert from cm-1 to au.
        int index = 0;
        for (int i = 0; i < DOFe; ++i) {
            H[i][i] = epsilon[i] = FMO7[i+shift][i+shift] * cm2au;
            for (int j = 0; j < DOFe; ++j)
                if (i < j) {
                    H[j][i] = H[i][j] = FMO7[i+shift][j+shift] * cm2au;
                    diabatic_coupling[index++] = H[i][j];
                }
        }
    }

    // * 1. Initialize and get model parameters by building or loading
    omega.resize(N, 0);
    c.resize(N, 0);
    req.resize(N, 0);
    shift.resize(DOFn, 0);
    initializeModelParameters();

    // * 2. Extend parameters for nulcear sampling and energy/forces calculations
    // Note that after ModelBase::init(), the model parameters have been created
    // or loaded, but for convenience and general, we extend the length of omega,
    // req, and set shift accordingly.
    // Let omega size = DOFn, the R of multi-state share same set omega j.
    if (omega.size() == N) {
        omega.resize(DOFn, 0);
        for (int i = 1; i < DOFe; i++)
            for (int j = 0; j < N; j++)
                omega[i*N+j] = omega[j];
    }
    // Extend the size of req with same set req, if only one set req is builded.
    if (req.size() == N) {
        req.resize(DOFn, 0);
        for (int i = 1; i < DOFe; i++)
            for (int j = 0; j < N; j++)
                req[i*N+j] = req[j];
    }
    // Set equilibrium shifts of each R (used for nulcear sampling)
    // For FMO model, shift = 0 (initial values are 0) (ground state).

    // * 3. Save Hamiltonian matrix from file (if required) in input energy unit
    // Do this in subclasses since the subclasses init() may modify it.
    const std::string H_save= param.getStr("H_save");
    if (!H_save.empty())
        saveHamiltonianMatrix(H_save);
}

void HamiltonianFMO::buildModelParameters() {
    // Get the cutoff frequency or characteristic frequency with unit converison
    const double omega_c = param.getDouble("omega_c") * unit2au;
    // eta is Kondo parameter for Ohmic (dimensionless) bath (useless here.)
    const double eta = param.getDouble("eta");
    // lambda is reorganization energy for Debye bath, describe the system-bath coupling
    // strength. Generally, using one value for all states, which is done in paper.
    // But, here, we can load a list of lambda (the order of lambda follow the order of
    // state index) and generate different couplings, but with one set omega.
    // Note that in most/all of publications, same one lambda is used for all states.
    std::vector<double> lambda;
    SplitString(lambda, param.getStr("lambda"));
    for (int i = 0; i < lambda.size(); i++)
        lambda[i] *= unit2au;
    // Get omega, c by discretization of Debye spectral density.
    // This will give you one set omega and c.
    discretizeSpectralDensity("Debye", N, omega_c, eta, lambda[0], omega, c);
    if (lambda.size() == DOFe) { // using different sets c for differnet states
        c.resize(DOFn, 0); // For FMO model, total DOFn = DOFe*N
        // This discretization scheme used in Jian's CMM paper
        // Ref1: X. He, J. Liu, J. Chem. Phys. 2019, 151, 024105. eq. 40, 46, 47
        // Get other sets of c according to given lambda
        for (int i = 1; i < DOFe; i++)
            for (int j = 0; j < N; j++)
                c[i*N+j] = sqrt(2*lambda[i]/(1.0+N)) * omega[j];
    }
    else if (lambda.size() != 1)
        throw std::runtime_error("ERROR: For FMO model, the number of values for lambda should be 1 or DOFe.");
    // Get displacements req according to req = c/omega^2
    // Ref: J. Phys. Chem. A, 124, 11006−11016 (2020).
    // section: Frenkel 7-Exciton Model for FMO
    req.resize(c.size(), 0);
    for (int i = 0; i < lambda.size(); i++)
        for (int j = 0; j < N; j++)
            req[i*N+j] = c[i*N+j] / (omega[j]*omega[j]);
}

void HamiltonianFMO::loadModelParameters(const std::string& loadfile) {
    // We need know the number of lines of file firstly to reszie the vector
    // if DOFe sets c/req data are provided.
    const double lineNumber = GetLineNumber(loadfile);
    if (lineNumber == (DOFn+1)) {
        omega.resize(DOFn, 0);
        c.resize(DOFn, 0);
        req.resize(DOFn, 0);
    }
    std::vector<std::vector<double>*> data = {&omega, &c, &req};
    LoadDataFile(loadfile, data);
}

void HamiltonianFMO::saveModelParameters(const std::string& savefile) {
    // If DOFe sets c/req are generated share one set omega, we should extend
    // omega fistly using same set value.
    // Note that in this case, the total output index is DOFn = DOFe*N.
    if (c.size() == DOFn && omega.size() == N) {
        omega.resize(DOFn, 0);
        for (int i = 1; i < DOFe; i++)
            for (int j = 0; j < N; j++)
                omega[i*N+j] = omega[j];
    }
    std::vector<std::vector<double>*> data = {&omega, &c, &req};
    std::vector<std::string> headers = {"omega/au", "c", "req"};
    SaveDataFile(savefile, headers, data);
}

double HamiltonianFMO::getPotentialEnergy(int index) {
    if (index >= DOFe || index < 0)
        throw std::runtime_error("ERROR: State index out of scope.");
    // Ref: J. Phys. Chem. A, 124, 11006−11016 (2020). eq.22
    // Note that, here, req (D_k in paper) has been extended to N*DOFe with
    // one set of req or DOFe set of req according the list given lambda.
    double PE = 0.0;
    for (int i = 0; i < DOFe; ++i)
        for (int j = 0; j < N; ++j) {
            double X = R[i*N+j] - (i == index ? req[i*N+j] : 0);
            PE += omega[j] * omega[j] * X * X;
        }
    return 0.5 * PE + epsilon[index];
}

void HamiltonianFMO::getForces(int i, int j, std::vector<double>& F) {
    if (i >= DOFe || i < 0 || j >= DOFe || j < 0)
        throw std::runtime_error("ERROR: Called getForces() with wrong state index.");
    F.resize(DOFn, 0);
    if (i == j)  // forces of each state from diagonal Hamiltonian
        for (int s = 0; s < DOFe; ++s)
            for (int k = 0; k < N; ++k)
                F[s*N+k] = -omega[k] * omega[k] * (R[s*N+k] - (s == i ? req[s*N+k] : 0));
    // forces from off-diagonal Hamiltonian (diabatic coupling)
    else if (!Condon_approximation)
        throw std::runtime_error("ERROR: The force from non-Condon diabatic coupling for FMO model is undefined.");
    else // In the case of Condon_approximation, forces are zero
        std::fill(F.begin(), F.end(), 0);
}