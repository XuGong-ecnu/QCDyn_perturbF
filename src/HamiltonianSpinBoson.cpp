/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 2, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "HamiltonianSpinBoson.h"

void HamiltonianSpinBoson::init() {
    // * 0. Initialize data members and Hamiltonian.
    HamiltonianModelBase::init();
    if (model_type != "SB" && model_type != "SB_req" && model_type != "GOA")
        throw std::runtime_error("ERROR: Undefined model_type=" + model_type + " for Spin-Boson model.");
    if (DOFe != 2)
        throw std::runtime_error("ERROR: For Spin-Boson model, DOFe must equal to 2.");
    // For SB model, number of normal modes, N = DOFn.
    N = DOFn;
    // Set minimum energy of each state, epsilon and modify the diagonal element
    // epsilon[0] is epsilon and has been loaded from input in ModelBase::init().
    if (model_type == "SB") {
        H[0][0] = epsilon[0];
        H[1][1] = epsilon[1] = -epsilon[0];
    }
    else { // model_type = SB_req or GOA
        H[0][0] = epsilon[0] *= 2.0/hbar;  // epsilon = 0.5*hbar*omega_DA
        H[1][1] = epsilon[1] = 0.0;
    }

    // * 1. Initialize and get model parameters by building or loading
    omega.resize(N, 0);
    c.resize(N, 0);
    req.resize(N, 0);
    shift.resize(DOFn, 0);
    if (!Condon_approximation) {
        // Currently, for non-Condon SB/SB_req model, the model parameters as
        // including the linear coupling coefficients (gamma) can only be loaded
        // from external file. While for non-Condon GOA model, the gamma can be
        // dertermined from the input value of key "GOA_gamma".
        if (param.getStr("model_load").empty() && model_type != "GOA")
            throw std::runtime_error("ERROR: The non-Condon SB/SB_req model parameters should be loaded from file.");
        gamma.resize(N, 0);
    }
    initializeModelParameters();
    getSamplingShifts();

    // * 2. Save Hamiltonian matrix from file (if required) in input energy unit
    // Do this in subclasses since the subclasses init() may modify it.
    const std::string H_save= param.getStr("H_save");
    if (!H_save.empty())
        saveHamiltonianMatrix(H_save);
}

void HamiltonianSpinBoson::buildModelParameters() {
    // By default, Ohmic bath will be used for SB model.
    std::string spec_density = param.getStr("spec_density");
    if (spec_density.empty())
        spec_density = "Ohmic";
    // the cutoff frequency or characteristic frequency
    const double omega_c = param.getDouble("omega_c") * unit2au;
    // eta is Kondo parameter for Ohmic (dimensionless),
    // It is used to describe the system-bath coupling strength
    const double eta = param.getDouble("eta");
    // lambda is reorganization energy for Debye bath
    const double lambda = param.getDouble("lambda") * unit2au;
    // Firstly, Get omega, c by discretization of spectral density.
    // for SB/SB_req model, this is enough.
    discretizeSpectralDensity(spec_density, N, omega_c, eta, lambda, omega, c);
    // Get req according to req = 2*c/omega^2
    for (int j = 0; j < DOFn; ++j)
        req[j] = 2 * c[j] / (omega[j]*omega[j]);
    // GOA model with bilinear coupled primary and secondary modes
    // Reference: Geva, JCP 148, 102304 (2018) sectioin III
    if (model_type == "GOA") {
        if (spec_density != "Ohmic")
            throw std::runtime_error("ERROR: For GOA model, the Ohmic spectral density should be used.");
        // GOA_Omega is the frequency of the primary mode
        const double GOA_Omega = param.getDouble("GOA_Omega") * unit2au;
        if (GOA_Omega <= 0)
            throw std::runtime_error("ERROR: Illegal value of GOA_Omega (should be > 0).");
        // 2*y0 is the shift in equilibrium geometry of the primary mode
        // between the donor and acceptor states
        // Note: the unit of it should be au or unitless.
        const double GOA_y0    = param.getDouble("GOA_y0");
        if (GOA_y0 == 0)
            throw std::runtime_error("ERROR: Illegal value of GOA_y0 (can't be 0).");
        // The effective freq is obtained by diagonalizing the Hessian matrix
        // using matrix from Eigen lib: MatrixXd, d means double,
        // X means dynamics (unknown size before runnning).
        // The initial values are zero.
        Eigen::MatrixXd D(DOFn, DOFn);
        // Get D matrix according to JCP 148, 102304 (2018) equation 38
        // The above omega, c in SB model aslo can be as freq and coupling
        // of bath modes (from 1 to DOFn-1, index in C++ is 0 ~ DOFn-2) in GOA model
        D(0, 0) = GOA_Omega * GOA_Omega;
        for (int j = 1; j < DOFn; ++j) {
            D(0, 0) += c[j-1]*c[j-1] / (omega[j-1]*omega[j-1]);
            D(0, j) = D(j, 0) = c[j-1];
            D(j, j) = omega[j-1]*omega[j-1];
        }
        // diagonalize D matrix, compute eigenvalues and eigenvectors.
        // Here, D is a real Hermite matrix (or self adjoint matrix)
        // Construct an object of SelfAdjointEigenSolver to compute eigenvalues.
        // For Hermite matrix, SelfAdjointEigenSolver is faster and more accurate
        // than the general purpose eigenvalue algorithms that implemented in
        // EigenSolver and ComplexEigenSolver.
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(D);
        const Eigen::VectorXd& eigenvalues = eig.eigenvalues();
        // T_trans is the unitary diagonalizing transformation matrix
        const Eigen::MatrixXd& T_trans = eig.eigenvectors().transpose();
        // Get omega, eigenvalues are the square of effective freq
        const Eigen::VectorXd& w_tilde = eigenvalues.array().sqrt();
        // Get req, the donor-toacceptor shifts in equilibrium geometry
        // along  the normal mode coordinates (equation 40)
        Eigen::VectorXd C(DOFn);
        C(0) = 1;
        for (int j = 1; j < DOFn; ++j)
            C(j) = -c[j-1] / (omega[j-1]*omega[j-1]);
        const Eigen::VectorXd& R_eq = 2 * GOA_y0 * T_trans * C;
        // Get noeq shifts, used for initial sampling of R
        // positive means shift further away from diabats' crossing
        // JCP 148, 102304 (2018) equation 37, 41-42
        // GOA_shift is noneq shift of the primary mode, 0=equilibrium (default)
        const double GOA_shift = param.getDouble("GOA_shift") * unit2au;
        for (int j = 0; j < DOFn; ++j) {
            shift[j] = T_trans(j, 0) * GOA_shift;
            req[j] = R_eq(j);
            omega[j] = w_tilde(j);
            // Get effective coupling strength normal mode, c_tilde
            c[j] = 0.5/GOA_y0 * req[j] * omega[j] * omega[j];
        }
        // The followings are used for modified non-Condon GOA model.
        // Reference: Xiang Sun, J. Chem. Phys. 144, 244105 (2016). sectioin III
        if (!Condon_approximation) {
            // γ is the the electronic coupling coefficient in non-Condon GOA model,
            // which decides the linear coupling coefficients {γ_j}.
            const double GOA_gamma = param.getDouble("GOA_gamma");
            if (GOA_gamma == 0)
                throw std::runtime_error("ERROR: Illegal value of GOA_gamma (can't be 0) for non-Condon GOA model.");
            // Get linear coupling coefficients, {γ_j} for non-Condon GOA model
            // See Eq. 28 in J. Chem. Phys. 144, 244105 (2016).
            for (int j = 0; j < DOFn; ++j)
                gamma[j] = T_trans(j, 0) * GOA_gamma;
        }
    }
}

void HamiltonianSpinBoson::loadModelParameters(const std::string& loadfile) {
    std::vector<std::vector<double>*> data = {&omega, &c, &req};
    if (model_type == "GOA")
        data.emplace_back(&shift);
    if (!Condon_approximation)
        data.emplace_back(&gamma);
    LoadDataFile(loadfile, data);
}

void HamiltonianSpinBoson::saveModelParameters(const std::string& savefile) {
    std::vector<std::vector<double>*> data = {&omega, &c, &req};
    std::vector<std::string> headers = {"omega", "c", "req"};
    if (model_type == "GOA") {
        data.emplace_back(&shift);
        headers.emplace_back("shift");
    }
    if (!Condon_approximation) {
        data.emplace_back(&gamma);
        headers.emplace_back("gamma");
    }
    SaveDataFile(savefile, headers, data);
}

void HamiltonianSpinBoson::getSamplingShifts() {
    // It is useful to load the equilbrium shifts of R for nuclear sampling.
    // such as, we want to sample on the ground state, but don't want to inlcude
    // it in the dynamics simulation.
    const std::string loadfile = param.getStr("shift_load");
    if (!loadfile.empty()) {
        // Requirments of the shifts file (.csv or .dat file):
        // 1. First line is headers: index,shift
        // 2. The number of lines is DOFn+1.
        // 3. The unit is au.
        std::vector<std::vector<double>*> data = {&shift};
        LoadDataFile(loadfile, data);
    }
    // the shifts are dertermined based on the sample state, which is normal case.
    else {
        // By default, for SB or SB_req model, sampling on Donor state (state=0)
        // However, it can be spcified by user.
        // For GOA model, shift has been setted in build/loadModelParameters, (do nothing here)
        const int sample_state = param.getInt("sample_state");
        // Sample on min of Donor State (sample_state == 0):
        // SB: Donor min at (Rj+c_j/omega_j^2), then shift is c_j/omega_j^2
        // SB_req: Donor min at 0, then shift is 0, (0 is initial value, do nothing here)
        if (sample_state == 0) { // min of donor state
            if (model_type == "SB")
                for (int j = 0; j < DOFn; ++j)
                    shift[j] = c[j] / (omega[j]*omega[j]);
        }
        // Sample on min of acceptor State (sample_state == 1):
        // SB: Accetor min at (Rj-c_j/omega_j^2), then shift is -c_j/omega_j^2
        // SB_req: Accetor min at (Rj-reqj), then shift is -req.
        // ? Must use {} for this branch, otherwise it won't be excuted.
        else if (sample_state == 1) { // min of acceptor state
            if (model_type == "SB")
                for (int j = 0; j < DOFn; ++j)
                    shift[j] = -c[j] / (omega[j]*omega[j]);
            else if (model_type == "SB_req")
                for (int j = 0; j < DOFn; ++j)
                    shift[j] = -req[j];
        }
        // Sample on around 0 (shift is 0) (sample_state == -1)
        else if (sample_state != -1)
            throw std::runtime_error("ERROR: sample_state=" + std::to_string(sample_state) +
                " is illegal for Spin-Boson model.");
    }
}

double HamiltonianSpinBoson::getPotentialEnergy(int index) {
    if (index >= DOFe || index < 0)
        throw std::runtime_error("ERROR: State index out of scope.");
    double PE = 0.0;
    if (model_type == "SB") {
        double cR = 0.0;
        for (int j = 0; j < DOFn; ++j) {
            PE += omega[j] * omega[j] * R[j] * R[j];
            cR += c[j] * R[j];
        }
        PE *= 0.5;
        PE += index == 0 ? cR : -cR;
    }
    else { // model_type = SB_req or GOA
        for (int j = 0; j < DOFn; ++j) {
            double X = R[j] - (index == 0 ? 0 : req[j]);
            PE += omega[j] * omega[j] * X * X;
        }
        PE *= 0.5;
    }
    return PE + epsilon[index];
}

void HamiltonianSpinBoson::getForces(int i, int j, std::vector<double>& F) {
    if (i >= DOFe || i < 0 || j >= DOFe || j < 0)
        throw std::runtime_error("ERROR: Called getForces() with wrong state index.");
    F.resize(DOFn, 0);
    if (i == j) { // froces of each state from diagonal Hamiltonian
        if (model_type == "SB")
            for (int k = 0; k < DOFn; ++k)
                F[k] = -omega[k] * omega[k] * R[k] + (i == 0 ? -c[k] : c[k]);
        else // model_type = SB_req or GOA
            for (int k = 0; k < DOFn; ++k)
                F[k] = -omega[k] * omega[k] * (R[k] - (i == 0 ? 0 : req[k]));
    }
    // forces from off-diagonal Hamiltonian (diabatic coupling)
    // force is F[k] = -dH_ij/dR = - γ[k], this form are same for SB/SB_req/GOA
    // and two-state LVC model, all of them are in the form linear of coordinates.
    else if (!Condon_approximation) {
        for (int k = 0; k < DOFn; ++k)
            F[k] = -gamma[k];
    }
    else // In the case of Condon_approximation, forces are zero
        std::fill(F.begin(), F.end(), 0);
}

double HamiltonianSpinBoson::getNonCondonCoupling(int i, int j) {
    if (i >= DOFe || i < 0 || j >= DOFe || j < 0 || i == j)
        throw std::runtime_error("ERROR: Called getNonCondonCoupling() with wrong state index.");
    // H_ji = H_ij = Σγ * R_k, where γ is the linear coupling coefficient.
    // This form are same for SB/SB_req/GOA and two-state LVC model, all of them
    // are in the form linear of coordinates.
    double coupling = 0;
    for (int k = 0; k < DOFn; ++k)
        coupling += gamma[k] * R[k];
    return coupling;
}

double HamiltonianSpinBoson::getNormalModePE(int state, int mode) {
    if (state >= DOFe || state < 0)
        throw std::runtime_error("ERROR: State index out of scope.");
    if (mode >= N || mode < 0)
        throw std::runtime_error("ERROR: Normal mode index out of scope.");
    double PE = 0.0;
    if (model_type == "SB") {
        double cR = c[mode] * R[mode];
        PE = 0.5 * omega[mode] * omega[mode] * R[mode] * R[mode];
        PE += state == 0 ? cR : -cR;
    }
    else { // model_type = SB_req or GOA
        double X = R[mode] - (state == 0 ? 0 : req[mode]);
        PE = 0.5 * omega[mode] * omega[mode] * X * X;
    }
    return PE;
}

double HamiltonianSpinBoson::getNormalModeKE(int mode) {
    if (mode >= N || mode < 0)
        throw std::runtime_error("ERROR: Normal mode index out of scope.");
    // For model, mass = 1.
    return 0.5 * V[mode] * V[mode];
}

void HamiltonianSpinBoson::updateAdiabaticHamiltonian() {
    // Save current H as H_old before update
    H_old = H;
    Heff_old = Heff;
    // For Condon SB model, the analytical form is used, which doesn't work
    // for SB_req and GOA model.
    if (model_type == "SB" && Condon_approximation)
        updateAdiabaticHamiltonianForSB();
    // While for SB_req and GOA models or non-Condon models, I use the general
    // one in HamiltonianModelBase. Perhaps, someone can write anlytical for them.
    else { // model_type = SB_req or GOA
        updateHAndT();
        updateFallAndNAC();
    }
}

void HamiltonianSpinBoson::updateAdiabaticHamiltonianForSB() {
    // * 1. Compute potential energy of each state
    // Ref1: J. Phys. Chem. Lett. 2018, 9, 5660−5663. eq.S14-S15 in SI
    // Here, cR = Σc_j*R_j + ε, the meaning of it is half of the energy difference
    // between two states. And Vb = Σw_j^2*R_j^2, the meaning of it is the sum
    // of enenrgy of two states
    double cR = 0.0, Vb = 0.0;
    double Gamma = getDiabaticCoupling(0, 1); // diabatic coupling: Γ_01
    for (int j = 0; j < DOFn; ++j) {
        Vb += omega[j] * omega[j] * R[j] * R[j];
        cR += c[j] * R[j];
    }
    cR += epsilon[0];
    double root = sqrt(cR*cR + Gamma*Gamma); // root = sqrt(cR^2 + Γ^2)
    // V± = 0.5*Σw_j^2*R_j^2 ± sqrt(cR^2 + Γ^2), state 0 is ground state (-).
    // Or V± = 0.5*E_sum ± sqrt(E_diff*E_diff/4 + Gamma*Gamma)
    H[0][0] = 0.5*Vb - root;
    H[1][1] = 0.5*Vb + root;
    H[0][1] = H[1][0] = 0;
    // Update effective Hamiltioan matrix (Heff) with removing H_avg
    // this Heff may be used in the symmterical Hamiltonian EOMs.
    Heff = H;
    double H_avg = 0;
    for (int i = 0; i < DOFe; i++)
        H_avg += H[i][i] / DOFe;
    for (int i = 0; i < DOFe; i++)
        Heff[i][i] -= H_avg;
    // * 2. Update the unitary rotation matrix T
    // Ref2: J. Chem. Phys. 116, 2346 (2002). eq.39-40
    // Note that the sign of the varibale, that used here is different from Ref2.
    // The matrix here is same as the matrix got in ModelBase class.
    double G = 1.0/cR * (-Gamma + root); // G = 1/cR * (-Γ + sqrt(cR^2 + Γ^2))
    double factor = 1.0/sqrt(2.0*(1.0+G*G)); // factor = 1/(sqrt(2(1+G^2)))
    T[0][0] = (G-1)*factor; // -sina, T_gD or T_Dg
    T[0][1] = (1+G)*factor; //  cosa, T_eD or T_De
    T[1][0] = (1+G)*factor; //  cosa, T_gA or T_Ag
    T[1][1] = (1-G)*factor; //  sina, T_eA or T_Ae
    // * 3. Update conjugate transpose (Hermitian transpose) matrix of T
    // Here, T is always a real matrix, so T_dagger = T_trans
    Matrix_Transpose(T, T_dagger);
    // * 4. Compute aidabatic forces of each state
    // Ref1: J. Phys. Chem. Lett. 2018, 9, 5660−5663. eq.S16 in SI
    // Note the negative sign: F = -dV/dR
    // Here, F_all is DOFe^2-dimentional-vector, the off-diagnol is 0.
    for (int j = 0; j < DOFn; ++j) {
        F_all[0][j] = -omega[j] * omega[j] * R[j] + c[j]*(cR/root); // F_g
        F_all[3][j] = -omega[j] * omega[j] * R[j] - c[j]*(cR/root); // F_e
    }
    // * 5. Compute nonaidabtic coupling (NAC) vector between ground and excited states
    // Ref1: J. Phys. Chem. Lett. 2018, 9, 5660−5663. eq.S17 in SI
    // Here, NAC[0] is d_01 = -d_10.
    factor = -0.5*Gamma/(root*root); // -0.5*Γ/(cR^2 + Γ^2)
    for (int j = 0; j < DOFn; ++j)
        NAC[0][j] = factor*c[j]; // d_01
}