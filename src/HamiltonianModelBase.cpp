/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 3, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "HamiltonianModelBase.h"

void HamiltonianModelBase::init() {
    // * 0. Initialize data members and Hamiltonian.
    HamiltonianBase::init();
    if (system_type != "model")
        throw std::runtime_error("ERROR: Unsupported system_type=" + system_type +
            "for model Hamiltonian.");
    if (DOFe < 2 || DOFn < 2)
        throw std::runtime_error("ERROR: For model, DOFe and DOFn must be greater than 1.");

    // * 1. Resize all vectors of in this class and set elements to 0.
    R.resize(DOFn, 0);
    V.resize(DOFn, 0);
    F.resize(DOFn, 0);
    F_avg.resize(DOFn, 0);
    F_all.resize(DOFe * DOFe, std::vector<double>(DOFn, 0));
    NAC.resize((DOFe-1)*DOFe/2, std::vector<double>(DOFn, 0));
    T.resize(DOFe, std::vector<double>(DOFe, 0));
    T_dagger.resize(DOFe, std::vector<double>(DOFe, 0));
    // * 2. Initialize and get model parameters (do this in subclasses).
}

void HamiltonianModelBase::updateDiabaticHamiltonian() {
    // Save current H as H_old before update
    H_old = H;
    Heff_old = Heff;
    // Remember to reset H_avg to 0.
    // the average potential energy of all states, which is removed in
    // effective energy matrix (Heff)
    double H_avg = 0.0;
    for (int i = 0; i < DOFe; ++i) {
        // Get potential energy of each state
        H[i][i] = getPotentialEnergy(i);
        // Compute H_avg
        H_avg += H[i][i]/DOFe;
        // Get forces of each state
        getForces(i, i, F_all[i * DOFe + i]);
        // Get the off-diagonal element of Hamiltonian matrix and forces
        // When using Condon approximation, coupling is constant and force is 0.
        // So, we only need to compute them in the non-Condon case.
        // Here, the Hamiltonian is a real symmetrix matrix, i.e., H_ij = H_ji
        // and F[i][j] = F[j][i] = - dH_ij/dR, i != j.
        if (!Condon_approximation)
            for (int j = i + 1; j < DOFe; ++j) {
                H[j][i] = H[i][j] = getDiabaticCoupling(i, j);
                getForces(i, j, F_all[i * DOFe + j]);
                F_all[j * DOFe + i] = F_all[i * DOFe + j];
            }
    }
    // Get effective Hamiltioan matrix (Heff) with removing H_avg (diagonal only)
    Heff = H;
    for (int i = 0; i < DOFe; ++i)
        Heff[i][i] -= H_avg;
    // Reset all elements to 0, then compute average forces (F_avg) (diagonal only)
    std::fill(F_avg.begin(), F_avg.end(), 0.0);
    for (int j = 0; j < DOFn; ++j) {
        for (int i = 0; i < DOFe; ++i)
            F_avg[j] += F_all[i * DOFe + i][j];
        F_avg[j] /= DOFe;
    }
}

void HamiltonianModelBase::updateAdiabaticHamiltonian() {
    // Save current H as H_old before update
    H_old = H;
    Heff_old = Heff;
    // Transform diabatic to adiabtic Hamiltonian via rotation matrix T
    // This will update data: H, T, T_dagger
    updateHAndT();
    // Compute adiabatic forces and nonadiabatic coupling vectors from
    // diabatic forces via rotation matrix T
    // This will update data: F_all, NAC
    updateFallAndNAC();
}

void HamiltonianModelBase::updateHAndT() {
    // * Get temperoray diabatic Hamiltonian firstly
    Real_Matrix H_diabatic(DOFe, std::vector<double>(DOFe, 0));
    for (int i = 0; i < DOFe; ++i) {
        H_diabatic[i][i] = getPotentialEnergy(i);
        for (int j = i+1; j < DOFe; ++j)
            H_diabatic[j][i] = H_diabatic[i][j] = getDiabaticCoupling(i, j);
    }
    // * Update adiabatic Hamiltonian via the diagonalization of H_diabatic
    // * And update the unitary rotation matrix T
    // Note the order of adiabatic states, here using energy ascending, so, the
    // index 00 in adiabatic Hamiltonian is always ground state.
    // For two-level system, the diagonalization can be done analytically.
    // Ref: Notes: Report 28 Ehrenfest dynamics, Section 2.1 and 2.2
    if (DOFe == 2) {
        // Get Gamma (coupling), average enenrgy, energ differnece from H_diabatic
        const double Gamma = H_diabatic[0][1];
        const double E_avg = 0.5 * (H_diabatic[0][0] + H_diabatic[1][1]);
        const double E_diff = H_diabatic[0][0] - H_diabatic[1][1];
        // Get adiabatic energy (eigenvalue of H_diabatic matrix)  (eq.31-32)
        // Here, E_plus is V_e (excited state), and E_minus is V_g (ground state).
        const double rt = sqrt(E_diff*E_diff + 4*Gamma*Gamma);
        const double E_minus = E_avg - rt * 0.5;
        const double E_plus = E_avg + rt * 0.5;
        // * Update adiabatic Hamiltonian (order of state is energy ascending)
        H[0][0] = E_minus;
        H[1][1] = E_plus;
        H[0][1] = H[1][0] = 0;
        // Get rotation angle alpha in rotation matrix (eq.33)
        // The range of value of atan() in C++ is [-0.5pi, 0.5pi]
        double alpha = 0.5 * atan(2*Gamma / E_diff);
        // if V_D < V_A, alpha ~ [-0.25pi, 0], then add 0.5pi to alpha, then
        // cos(alpha+0.5pi) = -sin(alpha), sin(alpha+0.5pi) = cos(alpha)
        if (E_diff < 0)
            alpha += 0.5*pi;
        // * Update unitary rotation matrix (transformation matrix) (eq.34)
        // Note that, here, T = (c_minu, c_plus), differnent from eq.34
        T[0][0] = -sin(alpha);
        T[1][0] = cos(alpha);
        T[0][1] = cos(alpha);
        T[1][1] = sin(alpha);
    }
    // For multi-level systems, the diagonalization has to be calculated numerically.
    // Here, using Eigen C++ lib to do this.
    else {
        // using matrix from Eigen lib: MatrixXcd, cd means complex double,
        // X means dynamics (unknown size before runnning), d means double,
        // 2d means 2*2 real matrix. The initial values are zero.
        Eigen::MatrixXd H_eigen(DOFe, DOFe);
        // Initialize h by assigning values from diabatic Hamiltonian H
        // (C++ vector<vector>). It is an Hermite matrix or self-adjoint matrix,
        // which means the conjugate transpose matrix of it is itself.
        for (int j = 0; j < DOFe; j++)
            for (int k = 0; k < DOFe; k++)
                H_eigen(j, k) = H_diabatic[j][k];
        // Construct an object of SelfAdjointEigenSolver to compute eigenvalues.
        // For Hermite matrix, SelfAdjointEigenSolver is faster and more accurate
        // than the general purpose eigenvalue algorithms that implemented in
        // EigenSolver and ComplexEigenSolver. And SelfAdjointEigenSolver will
        // sort eigenvalues ascending, since all eigenvalues are real numbers.
        // Moreover, it can make sure the eigenvectors are orthogonal.
        // Note that in EigenSolver or ComplexEigenSolver, the eigenvalues are
        // not sorted in any particular order.
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(H_eigen);
        // Get eigenvalues (ascending) of Hamiltonian h.
        // They are diganol elements of adiabatic Hamiltonian.
        const Eigen::VectorXd& eigenvalues = eig.eigenvalues();
        // * Update adiabatic Hamiltonian (order of state is energy ascending)
        // Note, the off-diagonal elements of it are always zero.
        for (int i = 0; i < DOFe; ++i) {
            H[i][i] = eigenvalues(i);
            for (int j = 0; j < DOFe; ++j)
                if (i != j)
                    H[i][j] = 0;
        }
        // * Update unitary rotation matrix T of Hamiltonian h
        // The columns of T are the coefficients of eigenvectors
        const Eigen::MatrixXd& T_eigen = eig.eigenvectors();
        // update the class data member T by assigning values from Eigen variables
        for (int j = 0; j < DOFe; j++)
            for (int k = 0; k < DOFe; k++)
                T[j][k] = T_eigen(j, k);
    }
    // * Update conjugate transpose (Hermitian transpose) matrix of T
    // Here, T is always a real matrix, so T_dagger = T_trans
    Matrix_Transpose(T, T_dagger);
    // * Update effective Hamiltioan matrix (Heff) with removing H_avg
    // this Heff may be used in the symmterical Hamiltonian EOMs.
    Heff = H;
    double H_avg = 0;
    for (int i = 0; i < DOFe; i++)
        H_avg += H[i][i] / DOFe;
    for (int i = 0; i < DOFe; i++)
        Heff[i][i] -= H_avg;
}

void HamiltonianModelBase::updateFallAndNAC() {
    // This function are modified to consider non-Condon case on Dec. 2, 2021.
    // * Get temperoray diabatic Forces of each state firstly
    Real_Matrix Fall_diabatic(DOFe * DOFe, std::vector<double>(DOFn, 0));
    for (int i = 0; i < DOFe; ++i) {
        getForces(i, i, Fall_diabatic[i * DOFe + i]); // forces of each state
        if (!Condon_approximation) // when it is true, off-diagnoal force is zero
            for (int j = i + 1; j < DOFe; ++j) { // forces from non-Condon diabatic coupling
                getForces(i, j, Fall_diabatic[i * DOFe + j]);
                Fall_diabatic[j * DOFe + i] = Fall_diabatic[i * DOFe + j];
            }
    }
    // * Compute adiabatic forces from diabatic forces via rotation matrix
    // we can get adiabatic force from diabatic forces using the rotation matrix
    // F_adiabatic = T_dagger*F_diabatic*T, then we can derive: the forces for
    // ath adiabtatic state: F_adiabatic[a][a] = (Σ_i (Σ_j T_ja * F_diabatic[j][i])*T_ia),
    // a/ij is index of adiabatic/diabatic state, T is rotation matrix.
    // Here, F_all is DOFe*DOFe matrix (stored as DOFe^2-dimentional vector),
    // but off-diagonal elemnts are 0, since the off-diagonal of adiabatic
    // Hamiltonian is always zero. While, off-diagonal of diabatic force may be
    // not zero in the non-Condon case. The index for adiabatic force of each
    // state is a*DOFe+a, a is the index of adiabatic state.
    for (int a = 0; a < DOFe; ++a) { // a is index of adiabatic state
        std::fill(F_all[a*DOFe+a].begin(), F_all[a*DOFe+a].end(), 0.0); // reset to 0
        for (int k = 0; k < DOFn; ++k) // k is index of nuclear DOF
            for (int i = 0; i < DOFe; ++i) { // i,j is index of diabatic state
                double temp = 0; // Σ_j T_ja * F_diabatic[j][i]
                for (int j = 0; j < DOFe; ++j)
                    temp += T[j][a] * Fall_diabatic[j*DOFe+i][k];
                F_all[a*DOFe+a][k] += temp * T[i][a];
            }
    }
    // Here, we also compute the average adiabatic forces F_avg = (-1/DOFe)F_all
    // this F_avg will be used in the symmterical Hamiltonian EOMs.
    std::fill(F_avg.begin(), F_avg.end(), 0);
    for (int a = 0; a < DOFe; ++a)
        for (int k = 0; k < DOFn; ++k)
            F_avg[k] += F_all[a * DOFe + a][k] / DOFe;
    // * Compute Nonadiabatic Coupling (NAC) vectors from diabatic forces via rotation matrix
    // Nonadiabatic coupling: d_ab = -d_ba = - <a|dH/dR|b>/(E_a-E_b), here, a,b
    // is adiabatic state, a < b. Adiabatic state |a> can be written in terms of
    // diabatic state |i> as |a> = Σ_i T_ia |i>, where, i is diabatic state.
    // diabatic -dH/dR = Σ_ij <i|F_ij|j>, then - <a|dH/dR|b> = (Σ_i T_ia |i>) *
    // (Σ_ij <i|F_ij|j>) * (Σ_i T_ib |i>) = Σ_ij T_ia * T_jb * F_ij
    // here, i, j is the index of diabtaic state, and a,b is adiabatic state.
    // The following loop will update all NAC of pairs.
    for (int a = 0; a < DOFe; ++a) // a, b is index of adiabatic state
        for (int b = a+1; b < DOFe; ++b) {
            int ab = (2*DOFe-a-1)*a/2 + (b-a)-1; // index of d_ab, a < b
            std::fill(NAC[ab].begin(), NAC[ab].end(), 0.0); // reset to 0
            // compute - <a|dH/dR|b> = Σ_ij T_ia * T_jb * F_ij
            for (int i = 0; i < DOFe; ++i)  // i,j is index of diabatic state
                for (int j = 0; j < DOFe; ++j)
                    for (int k = 0; k < DOFn; ++k) // index of nuclear DOF
                        NAC[ab][k] += T[i][a] * T[j][b] * Fall_diabatic[i*DOFe+j][k];
            // compute d_ab = - <a|dH/dR|b>/(E_a-E_b)
            // here, H is adiabtaic H, and has been updated.
            for (int k = 0; k < DOFn; ++k)
                NAC[ab][k] /= H[a][a] - H[b][b];
        }
}

void HamiltonianModelBase::updateQuasiDiabaticHamiltonian() {
    // Copy T of previous step as old one, which will be used later
    Real_Matrix T_old  = T; // transformation matrix at R(t1)
    // Transform diabatic to adiabtic Hamiltonian via rotation matrix T
    // This will update data: H (adiabtic), T, T_dagger at R(t2)
    updateHAndT();
    // Save current adiabatic H as H_old also as quasi-diabatic (QD) H_old
    H_old = H; // V_αβ(R(t1)) (using aidabatic as QD basis)
    // * Compute quasi-diabatic electronic Hamiltonian V_αβ(R(t2))
    // based on Eqs. 19-21 in J. Chem. Phys. 149, 044115 (2018) Section. II. C-D
    // Note that the notation in the following, α, β, μ, v is adiabatic state.
    // And t1 is previous time (old), t2 is current time (new), t2-t1 is nulcear time step.
    // V_αβ(R(t1)) = E_α(R(t1))δ_αβ (at t1, the QD is the adiabatic basis, H_old)
    // V_αβ(R(t2)) = Σ_μv S_αμ * E_μ(R(t2))δ_μv * S_dagger_βv
    // S_αμ = <Φ_α(R0)|Φ_μ(R(t2))>, S_dagger_βv = <Φ_v(R(t2))|Φ_β(R(0))>,
    // (Auctually, S and S_dagger is same for model Hamiltonians).
    // where the QD basis Φ_α(R0) is same as the adiabtic basis Φ_α(R(t1)).
    // Here, S_αμ and {S_βv}_dagger is the wavefunction overlap matrix.
    // For diabatic model, we can compute them by writting the adiabatic states
    // in terms of the diabtic ones: |α> = Σ_i T_iα |i>, where α/i is the index
    // of the adiabatic/diabatic state, and T is rotation (transformation) matrix.
    // Then, S_αμ=Σ_i(T_iα(R(t1))*T_iμ(R(t2))), S_dagger_βv=Σ_i(T_iv(R(t2))*T_iβ(R(t1)))
    // H_QD_new is V_αβ(R(t2)), which is diabatic-like and off-diagonal element are not zero.
    Real_Matrix H_QD_new(DOFe, std::vector<double>(DOFe, 0)); // V_αβ(R(t2))
    Real_Matrix S(DOFe, std::vector<double>(DOFe, 0)); // overlap matrix, S_αμ
    Real_Matrix S_dagger(DOFe, std::vector<double>(DOFe, 0)); // S_dagger_βv
    for (int a = 0; a < DOFe; ++a)
        for (int u = 0; u < DOFe; ++u)
            for (int i = 0; i < DOFe; ++i) // S_αμ = Σ_i (T_iα(R(t1)) * T_iμ(R(t2)))
                S[a][u] += T_old[i][a] * T[i][u];
    for (int b = 0; b < DOFe; ++b)
        for (int v = 0; v < DOFe; ++v)
            for (int i = 0; i < DOFe; ++i) // S_dagger_βv = Σ_i(T_iv(R(t2)) * T_iβ(R(t1)))
                S_dagger[b][v] += T[i][v] * T_old[i][b];
    for (int a = 0; a < DOFe; ++a)
        for (int b = 0; b < DOFe; ++b)
            for (int u = 0; u < DOFe; ++u)
                for (int v = 0; v < DOFe; ++v)
                    if (u == v) // adiabatic H[u][v] is E_μ(R(t2))δ_μv (if u != v, is 0)
                        H_QD_new[a][b] += S[a][u] * H[u][v] * S_dagger[b][v];
    // Save the new quasi-diabatic Hamiltonian V_αβ(R(t2)) to current H.
    H = H_QD_new;
    // * Compute quasi-diabatic forces (negative gradients) -▽V_αβ(R(t2))
    // -▽V_αβ(R(t2)) = Σ_μv S_αμ * -<Φ_μ(R(t2))|▽V(R(r;t2))|Φ_v(R(t2)> * S_dagger_βv
    // Here, V(R(r;t2)) is adiabatic electronic Hamiltonian operator,
    // So -<Φ_μ(R(t2))|▽V(R(r;t2))|Φ_v(R(t2)> is adibatic forces (μ == v)
    // and nonadiabatic coupling d_μv * (E_μ - E_v) (μ != v), which can be computed
    // from diabatic forces and rotation matrix. Please see updateFallAndNAC().
    // -<Φ_μ(R(t2))|▽V(R(r;t2))|Φ_v(R(t2)> = Σ_i T_iμ*T_iv * F_diabatic[i]
    // here, i is the index of diabtaic state. μ, v is adiabatic state.
    // 1. First, get temperoray diabatic Forces of each state
    Real_Matrix Fall_diabatic(DOFe * DOFe, std::vector<double>(DOFn, 0));
    for (int i = 0; i < DOFe; ++i) {
        getForces(i, i, Fall_diabatic[i * DOFe + i]); // forces of each state
        if (!Condon_approximation)
            for (int j = i + 1; j < DOFe; ++j) { // forces from non-Condon diabatic coupling
                getForces(i, j, Fall_diabatic[i * DOFe + j]);
                Fall_diabatic[j * DOFe + i] = Fall_diabatic[i * DOFe + j];
            }
    }
    // 2. Then, compute -<Φ_μ(R(t2))|▽V(R(r;t2))|Φ_v(R(t2)> = Σ_i T_iμ*T_iv * F_diabatic[i]
    // Here, M denotes the -<Φ_μ(R(t2))|▽V(R(r;t2))|Φ_v(R(t2)>, which is stored
    // as DOFe*DOFe dimentional vector, i.e., the index of M_μv matrix is M[μ*DOFe+v].
    Real_Matrix M(DOFe*DOFe, std::vector<double>(DOFn, 0));
    for (int u = 0; u < DOFe; ++u)
        for (int v = 0; v < DOFe; ++v)
            for (int i = 0; i < DOFe; ++i) // index of diabatic state
                for (int k = 0; k < DOFn; ++k) // index of nuclear DOF
                    M[u*DOFe+v][k] += T[i][u]*T[i][v] * Fall_diabatic[i*DOFe+i][k];
    // 3. Last, compute quasi-diabatic forces F_all = -▽V_αβ(R(t2))
    // F_all = -▽V_αβ(R(t2)) = Σ_μv S_αμ * M_μv * S_dagger_βv
    for (int a = 0; a < DOFe; ++a)
        for (int b = 0; b < DOFe; ++b) {
            std::fill(F_all[a*DOFe+b].begin(), F_all[a*DOFe+b].end(), 0); // reset to 0
            for (int u = 0; u < DOFe; ++u)
                for (int v = 0; v < DOFe; ++v) {
                    double factor = S[a][u] * S_dagger[b][v]; // S_αμ * S_dagger_βv
                    for (int k = 0; k < DOFn; ++k) // index of nuclear DOF
                        F_all[a*DOFe+b][k] += factor * M[u*DOFe+v][k];
                }
        }
    // Here, we also compute the average forces F_avg = (-1/DOFe)Σ▽V_αα(R(t2))
    // this F_avg will be used in the symmterical Hamiltonian EOMs.
    std::fill(F_avg.begin(), F_avg.end(), 0);
    for (int a = 0; a < DOFe; ++a)
        for (int k = 0; k < DOFn; ++k)
            F_avg[k] += F_all[a * DOFe + a][k] / DOFe;
}

double HamiltonianModelBase::getKineticEnergy() {
    double KE = 0.0;
    // For model, mass = 1.
    for (int j = 0; j < DOFn; j++)
        KE += V[j] * V[j];
    return 0.5 * KE;
}

void HamiltonianModelBase::initializeModelParameters() {
    // Build model parameters if file isn't provided, else load from file.
    const std::string loadfile = param.getStr("model_load");
    if (loadfile.empty())
        buildModelParameters();
    else
        loadModelParameters(loadfile);
    // Save force filed parameters to file if requested
    const std::string savefile = param.getStr("model_save");
    if (!savefile.empty())
        saveModelParameters(savefile);
}

void HamiltonianModelBase::discretizeSpectralDensity(const std::string& spec_density, int N, double omega_c, double eta,
                                                     double lambda, std::vector<double>& omega, std::vector<double>& c) {
    if (omega_c <= 0)
        throw std::runtime_error("ERROR: Illegal value of cutoff frequency (omega_c) for spectral density (should be > 0).");
    // To obtain the N nuclear mode frequencies and coupling coefficients by
    // discretization the specific spectral denisty. Typical bath is Ohmic and Debye.
    // Discretization can be done in a variety of different ways, all of which lead
    // to the same result when enough bath modes are included.
    // Note that in code, j start from 0, while 1 in paper
    // 3 discretization schemes are implemented for Ohmic bath.
    // and 1 discretization scheme for Debye bath.
    if (spec_density == "Ohmic") {
        if (eta <= 0)
            throw std::runtime_error("ERROR: Illegal value of Kondo parameter (eta) for Ohmic bath (should be > 0).");
        // By default, CMM scheme will be used for Ohmic spectral density.
        std::string Ohmic_discrete = param.getStr("Ohmic_discrete");
        if (Ohmic_discrete.empty())
            Ohmic_discrete = "CMM";
        // 1. This discretization scheme used in Jian's CMM paper (recommand)
        // Ref1: X. He, J. Liu, J. Chem. Phys. 2019, 151, 024105. eq. 39, 44-45
        // Ref2: J. Phys. Chem. Lett. 2021, 12, 2496−2501. in SI, eq.S34
        // Here, if using omega[j] = -omega_c * log((j+1.0)/(1.0+DOFn)), then
        // omega will be from large to small with same value, which is used in
        // J. Chem. Phys. 127, 144503 (2007). eq.3.5
        if (Ohmic_discrete == "CMM")
            for (int j = 0; j < N; j++) {
                omega[j] = -omega_c * log(1.0 - (j+1.0)/(1.0+N));
                c[j] = sqrt(eta * omega_c/(1.0+N)) * omega[j];
            }
        // 2. This discretization scheme is given in Appendix C of
        // Geva, J. Chem. Phys. 150, 034101 (2019), which also is used in Xing's
        // paper and fortran code. The max. freq should be provided.
        // This scheme is also used in most TBSH literatures.
        else if (Ohmic_discrete == "GQME") {
            const double omega_max = param.getDouble("omega_max") * unit2au; // max. freq
            if (omega_max <= 0)
                throw std::runtime_error("ERROR: Illegal value of max frequency (omega_max) for Ohmic bath (should be > 0).");
            double omega_0 = omega_c / N * (1.0 - exp(-omega_max/omega_c));
            for (int j = 0; j < N; j++) {
                omega[j] = -omega_c * log(1.0 - omega_0/omega_c * (j + 1.0));
                c[j] = sqrt(eta * hbar * omega_0) * omega[j];
            }
        }
        // 3. This discretization scheme used in Richardson's spin-maping paper
        // J. E. Runeson and J. O. Richardson, J. Chem. Phys. 151, 044119 (2019)
        // Orignal ref of its discretization procedure:
        // J. Chem. Phys. 122, 084106 (2005) eq. 29, 31-32
        // This scheme is seldom use.
        else if (Ohmic_discrete == "SPM")
            for (int j = 0; j < N; j++) {
                omega[j] = -omega_c * log((j + 0.5) / N);
                c[j] = sqrt(eta * omega_c / N) * omega[j];
            }
        else
           throw std::runtime_error("ERROR: Unsupported Ohmic_discrete=" + Ohmic_discrete + " to discretize Ohmic bath.");
    }
    else if (spec_density == "Debye") {
        if (lambda <= 0)
            throw std::runtime_error("ERROR: Illegal value of reorganization energy (lambda) for Debye bath (should be > 0).");
        // This discretization scheme used in Jian's CMM paper
        // Ref1: X. He, J. Liu, J. Chem. Phys. 2019, 151, 024105. eq. 40, 46, 47
        // Ref2: J. Phys. Chem. Lett. 2021, 12, 2496−2501. in SI, eq.S35-S36
        // And the discretization procedure used in Richardson's generized
        // spin-mapping for FMO model (J. Chem. Phys. 152, 084110 (2020)) is same
        // as this, but omega is from samll to large with same value, formula is
        // omega[j] = omega_c * tan(0.5*pi * ((j+1.0)/(1.0+DOFn)))
        // see ref : J. Chem. Phys. 127, 144503 (2007). eq.3.6
        for (int j = 0; j < N; j++) {
            omega[j] = omega_c * tan(0.5*pi * (1 - (j+1.0)/(1.0+N)));
            c[j] = sqrt(2*lambda/(1.0+N)) * omega[j];
        }
    }
    else if (spec_density == "Brownian") {
        // Ref: Xiang and Geva, J. Phys. Chem. A 2016, 120, 2976−2990
        // J(w) equation 40, and req equation 41, omega_c is Ω.
        // Here, get req according to req = 2*c/omega^2
        // Using this for Dynamics ref: J. Phys. Chem. A 2013, 117, 6196−6204
        // GOA_Omega is the frequency of the primary mode
        const double GOA_Omega = param.getDouble("GOA_Omega") * unit2au;
        if (lambda <= 0 || eta <= 0 || GOA_Omega <= 0)
            throw std::runtime_error("ERROR: Found illegal value in reorganization energy (lambda), "
                "or GOA_Omega, or Kondo parameter (eta) for Brownian spectral density (They should be > 0).");
        const double delta_w = omega_c/N;
        const double factor = delta_w*GOA_Omega*GOA_Omega*eta*lambda/pi;
        for (int j = 0; j < N; j++) {
            omega[j] = (j+1.0) * delta_w;
            double sqr1 = omega[j]*omega[j]-GOA_Omega*GOA_Omega;
            double sqr2 = eta*omega[j];
            c[j] = sqrt(factor/(sqr1*sqr1+sqr2*sqr2)) * omega[j];
        }
    }
    else
        throw std::runtime_error("ERROR: Undefined spec_density=" + spec_density);
}
