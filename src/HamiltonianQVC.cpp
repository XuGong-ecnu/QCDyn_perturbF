/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * First created: Nov. 28, 2021                                               *
 * Last updated: Dec. 2, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "HamiltonianQVC.h"

void HamiltonianQVC::init() {
    // * 0. Initialize data members and Hamiltonian.
    HamiltonianModelBase::init();
    if (model_type != "QVC")
        throw std::runtime_error("ERROR: Undefined model_type=" + model_type + " for QVC model.");
    if (Condon_approximation)
        throw std::runtime_error("ERROR: Please set Condon_approximation to false for QVC model.");
    if (DOFe != 2)
        throw std::runtime_error("ERROR: For QVC model, DOFe must equal to 2.");
    if (epsilon[0] != 0)
        throw std::runtime_error("ERROR: For QVC model, the first value of epsilon should be 0.");
    // For QVC model, number of normal modes, N = DOFn, they are donor state.
    N = DOFn;

    // * 1. Initialize and get model parameters by loading from file
    omega.resize(N, 0); // N frequecies for donor normal modes
    omega_A.resize(N, 0); // N frequecies for acceptor normal modes
    req.resize(N, 0); // the shift vector of acceptor state (x_A in the paper)
    shift.resize(N, 0); // used for initial nuclear sampling, they are zero for QVC model
    gamma.resize(N, 0); // linear coupling coefficient
    J.resize(N, std::vector<double>(N, 0)); // Duschinsky matrix
    Theta.resize(N, std::vector<double>(N, 0)); // quadratic coupling coefficient
    initializeModelParameters();

    // * 2. Check the nuclear sampling state
    // If sampling on the donor state (sample_state = 1), there is no shift;
    // If sampling on the acceptor state (sample_state = 0), the nuclear
    // positions/momenta should be rotated and shifted to donor state, which
    // will be done after the initial nuclear sampling of each trajectory.
    const int sample_state = param.getInt("sample_state");
    if (sample_state < 0 || sample_state >= DOFe)
        throw std::runtime_error("ERROR: Illegal sample_state=" + std::to_string(sample_state) +
            " for QVC model. The nuclear sampling should be performed on the donor"
            " (0, default) or acceptor state (1).");
}

void HamiltonianQVC::buildModelParameters() {
    throw std::runtime_error("ERROR: The QVC model parameters should be loaded from file.");
}

void HamiltonianQVC::loadModelParameters(const std::string& loadfile) {
    std::vector<std::vector<double>*> data = {&omega, &omega_A, &req, &gamma};
    for (int j = 0; j < N; ++j)
        data.emplace_back(&Theta[j]);
    for (int j = 0; j < N; ++j)
        data.emplace_back(&J[j]);
    LoadDataFile(loadfile, data);
}

void HamiltonianQVC::saveModelParameters(const std::string& savefile) {
    std::vector<std::vector<double>*> data = {&omega, &omega_A, &req, &gamma};
    std::vector<std::string> headers = {"omega_D", "omega_A", "req_A", "gamma"};
    for (int j = 0; j < N; ++j) {
        data.emplace_back(&Theta[j]);
        headers.emplace_back("Theta_i" + std::to_string(j));
    }
    for (int j = 0; j < N; ++j) {
        data.emplace_back(&J[j]);
        headers.emplace_back("J_i" + std::to_string(j));
    }
    SaveDataFile(savefile, headers, data);
}

double HamiltonianQVC::getPotentialEnergy(int index) {
    if (index >= DOFe || index < 0)
        throw std::runtime_error("ERROR: State index out of scope.");
    // Ref: J. Chem. Phys. 141, 034104 (2014). Eq. 3-7
    // Donor state: H_D = Σ(1/2*P_D^2 + 1/2*omega_D^2*R_D^2)
    // Acceptor state: H_A = Σ(1/2*P_A^2 + 1/2*omega_A^2*R_A^2) + ΔE
    // Since P/R_A can be written as in terms of P/R_D with Duschinsky matrix J:
    // P_A,i = ΣJ_ij (P_D,j); R_A,i = ΣJ_ij (R_D,j + x_A,j)
    // Then we can know: R_D,i = (ΣJ_ji R_A,j) - x_A,i
    // Note that the P_A = P_D, since the J^dag*J is identity matrix. This
    // feature can be verified numerically by the function checkKineticEnergy()
    double PE = 0.0;
    if (index == 0) { // donor state
        for (int i = 0; i < N; ++i)
            PE += 0.5 * omega[i] * omega[i] * R[i] * R[i];
    }
    else { // index == 1, acceptor state
        for (int i = 0; i < N; ++i) {
            // Get normal mode of acceptor using Duschinsky matrix J and shift
            // vector req (the x_A,j in Eq. 7 in literature).
            double R_A = 0;
            for (int j = 0; j < N; ++j)
                R_A += J[i][j] * (R[j] + req[j]);
            PE += 0.5 * omega_A[i] * omega_A[i] * R_A * R_A;
        }
    }
    // Here, epsilon[0] is always 0, and epsilon[1] is the ΔE in Eq. 5
    return PE + epsilon[index];
}

void HamiltonianQVC::getForces(int i, int j, std::vector<double>& F) {
    if (i >= DOFe || i < 0 || j >= DOFe || j < 0)
        throw std::runtime_error("ERROR: Called getForces() with wrong state index.");
    F.resize(DOFn, 0);
    if (i == j)  { // froces of each state from diagonal nuclear Hamiltonian
        if (i == 0) { // froces of donor state
            for (int k = 0; k < N; ++k)
                F[k] = -omega[k] * omega[k] * R[k];
        }
        else { // froces of acceptor state
            // R_A,i = ΣJ_ij (R_D,j + x_A,j), compute -d(R_A,i^2)/d(R_D,i)
            // R_A,i = factor + J_ii (R_D,i + x_A,i)
            // factor = ΣJ_ij (R_D,j + x_A,j), j != i, which does't contain R_D,i
            // R_A,i^2 = factor^2 + J_ii^2*(R_D,i + x_A,i)^2 + 2*factor*J_ii*(R_D,i + x_A,i))
            // -d(R_A,i^2)/d(R_D,i) = J_ii^2(2*R_D,i + 2*x_A,i) + 2*factor*J_ii
            // F_i = -omega_A,i * omega_A,i * (J_ii^2 * (R_D,i + x_A,i) + factor*J_ii)
            for (int k = 0; k < N; ++k) {
                double factor = 0;
                for (int l = 0; l < N; ++l)
                    if (k != l)
                        factor += J[k][l] * (R[l] + req[l]);
                F[k] = -omega_A[k] * omega_A[k] * (J[k][k] * J[k][k] * (R[k] + req[k]) + factor * J[k][k]);
            }
        }
    }
    else { // forces from off-diagonal Hamiltonian (diabatic coupling)
        // V_AD = V_DA = ΣTheta_ij R_D,i *  R_D,j + Σγ_i R_D,i + ΔDA
        // compute force: -dV_DA/d(R_D,i) = - Σ Theta_ij * R_D,j - γ_i
        for (int k = 0; k < N; ++k) {
            F[k] = -gamma[k];
            for (int l = 0; l < N; ++l)
                F[k] -= Theta[k][l] * R[l];
        }
    }
}

double HamiltonianQVC::getNonCondonCoupling(int i, int j) {
    if (i >= DOFe || i < 0 || j >= DOFe || j < 0 || i == j)
        throw std::runtime_error("ERROR: Called getNonCondonCoupling() with wrong state index.");
    // V_AD = V_DA = ΣTheta_ij R_D,i *  R_D,j + Σγ_i R_D,i + ΔDA
    // Here, gamma_DA is the "Δ_DA" in the literature
    double coupling = param.getDouble("gamma_DA");
    for (int k = 0; k < N; ++k) {
        coupling += gamma[k] * R[k];
        for (int l = 0; l < N; ++l)
            coupling += Theta[k][l] * R[k] * R[l];
    }
    return coupling;
}

void HamiltonianQVC::transformNormalModes(const std::vector<double>& X_old, const Real_Matrix& J,
    const std::vector<double>& shift, std::vector<double>& X_new, bool inverse) {
    X_new.resize(N, 0);
    std::fill(X_new.begin(), X_new.end(), 0); // reset all elements to 0
    // * Transfrom donor to acceptor by using X_A,i = ΣJ_ij (X_D,j + shift_j)
    if (!inverse) { // write acceptor in terms of donor (inverse = false) [default]
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                X_new[i] += J[i][j] * (X_old[j] + shift[j]);
    }
    // * Transfrom acceptor to donor by using X_D,i = (ΣJ_ji X_A,j) - shift_i
    else { // write donor in terms of acceptor (inverse = true)
        for (int i = 0; i < N; ++i) {
            X_new[i] = -shift[i];
            for (int j = 0; j < N; ++j)
                X_new[i] += J[j][i] * X_old[j];
        }
    }
}

bool HamiltonianQVC::checkKineticEnergy() {
    double KE_D = 0, KE_A = 0;
    for (int i = 0; i < N; ++i) // donor state
        KE_D += 0.5 * V[i] * V[i];
    for (int i = 0; i < N; ++i) { // acceptor state
        double V_A = 0; // transfrom donor velocities to acceptor
        for (int j = 0; j < N; ++j)
            V_A += J[i][j] * V[j];
        KE_A += 0.5 * V_A * V_A;
    }
    if (abs(KE_D - KE_A) > 1e-10) {
        std::cout << "KE_D = "  << KE_D << "; KE_A = " << KE_A << std::endl;
        return false;
    }
    else
        return true;
}
