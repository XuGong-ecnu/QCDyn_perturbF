/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Oct. 28, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsElec.h"

void DynamicsElec::init() {
    // Parameters used for electronic propagation
    DOFe = param.getInt("DOFe");
    EPN  = param.getInt("EPN");
    if (EPN < 1)
        throw std::runtime_error("ERROR: Illegal value for EPN (requires >= 1)");
    dt = param.getDouble("DT") / EPN;
    // If real units are provided, convert unit of dt from ps to au used for propagation.
    // If the unit is empty or "none", it is a unitless model.
    const std::string unit = param.getStr("energy_unit");
    if (!unit.empty() && unit != "none")
        dt = dt * ps2au;
}

// # The following three functions are used to propagate electronic DOF (q, p)
// # with 4th-order Runge-Kutta method in aidiabatic basis.
void DynamicsElec::MOVE_elec(double t_nucl, const Real_Matrix& H_old, const Real_Matrix& H_new,
                             const std::vector<double>& V_old, const std::vector<double>& V_new,
                             const Real_Matrix& NAC_old, const Real_Matrix& NAC_new,
                             std::vector<double>& q, std::vector<double>& p, std::vector<Complex>& coeff) {
    // cast y[0..2DOFe-1] <- (q[0..DOFe-1], p[0..DOFe-1])
    std::vector<double> y(2*DOFe, 0.0);
    for (int i = 0; i < DOFe; i++) {
        y[i] = q[i];
        y[DOFe + i] = p[i];
    }
    // propagate EPN steps of dt for electronic DOF
    for (int i = 0; i < EPN; i++) {
        // t_nucl is the current time of nuclear propagation (the time of H_old)
        // t_elec is the current time of electronic propagation
        double t_elec = t_nucl + i * dt;
        // explicit classic Runge-Kutta 4th order algorithm for ODE dy/dt = f(t,y)
        rk4(t_nucl, t_elec, H_old, H_new, V_old, V_new, NAC_old, NAC_new, y);
    }
    // cast back y[0..2DOFe-1] -> (q[0..DOFe-1], p[0..DOFe-1])
    for (int i = 0; i < DOFe; i++) {
        q[i] = y[i];
        p[i] = y[DOFe + i];
    }
    // update coeff from new q, p
    // coefficient of electronic wavefunction, coeff[j]=(q[j]+I*p[j])/sqrt(2)
    for (int i = 0; i < DOFe; i++)
        coeff[i] = (q[i] + I*p[i]) / sqrt(2.0);
}

void DynamicsElec::rk4(double t_nucl, double t_elec, const Real_Matrix& H_old, const Real_Matrix& H_new,
                       const std::vector<double>& V_old, const std::vector<double>& V_new,
                       const Real_Matrix& NAC_old, const Real_Matrix& NAC_new, std::vector<double>& y) {
    const int n = y.size();
    std::vector<double> f1(n, 0);
    std::vector<double> f2(n, 0);
    std::vector<double> f3(n, 0);
    std::vector<double> f4(n, 0);
    std::vector<double> yt(n, 0);
    // [1] k1 = dt * f(t,y),  f(t, y) = dydt
    deriv(t_nucl, t_elec, H_old, H_new, V_old, V_new, NAC_old, NAC_new, y, f1);
    // yt = y + k1 / 2
    for (int i = 0; i < n; i++)
        yt[i] = y[i] + 0.5 * dt * f1[i];
    // [2] k2 = dt * f2, where f2 = f(t+dt/2 , y + k1 / 2)
    deriv(t_nucl, t_elec+0.5*dt, H_old, H_new, V_old, V_new, NAC_old, NAC_new, yt, f2);
    // yt = y + k2 / 2
    for (int i = 0; i < n; i++)
        yt[i] = y[i] + 0.5 * dt * f2[i];
    // [3] k3 = dt * f3, where f3 = f(t+dt/2, y + k2 /2)
    deriv(t_nucl, t_elec+0.5*dt, H_old, H_new, V_old, V_new, NAC_old, NAC_new, yt, f3);
    // yt = y + k3
    for (int i = 0; i < n; i++)
        yt[i] = y[i] + dt * f3[i];
    // [4] k4 = dt * f(t+dt, y + k3)
    deriv(t_nucl, t_elec+dt, H_old, H_new, V_old, V_new, NAC_old, NAC_new, yt, f4);
    // At last, accumulate k1-k4 with coefficients
    // y[t+dt] = y[t] + 1/6 * (k1 + 2k2 + 2k3 + k4), k = dt * f
    for (int i = 0; i < n; i++)
        y[i] += dt/6.0 * (f1[i] + 2.0 * f2[i] + 2.0 * f3[i] + f4[i]);
}

void DynamicsElec::deriv(double t_nucl, double t_prime, const Real_Matrix& H_old, const Real_Matrix& H_new,
                         const std::vector<double>& V_old, const std::vector<double>& V_new,
                         const Real_Matrix& NAC_old, const Real_Matrix& NAC_new,
                         const std::vector<double>& y, std::vector<double>& dydt) {
    // In a single nuclear step (DT), electronic DOF can be updated several times,
    // i.e., electronic time step (dt): 10*dt = DT, EPN = DT/dt = 10.
    // Here, use linear interpolation to estimate the enengy matrix, velocites,
    // and nonadiabatic coupling (NAC) in the intermediate time, such as:
    // H' = H_old + (H_new - H_old) * (t'- t)/DT (DT = EPN*dt)
    const double factor = (t_prime - t_nucl) / (EPN * dt);
    // compute enengy matrix at intermediate time
    Real_Matrix H_prime(H_old.size(), std::vector<double>(H_old[0].size(), 0));
    for (int j = 0; j < H_old.size(); j++)
        for (int k = 0; k < H_old[0].size(); k++)
            H_prime[j][k] = H_old[j][k] + (H_new[j][k] - H_old[j][k]) * factor;
    // compute velocities at intermediate time
    std::vector<double> V_prime(V_old.size(), 0);
    for (int j = 0; j < V_old.size(); j++)
        V_prime[j] = V_old[j] + (V_new[j] - V_old[j]) * factor;
    // compute nonadiabatic coupling (NAC) at intermediate time
    Real_Matrix NAC_prime(NAC_old.size(), std::vector<double>(NAC_old[0].size(), 0));
    for (int j = 0; j < NAC_old.size(); j++)
        for (int k = 0; k < NAC_old[0].size(); k++)
            NAC_prime[j][k] = NAC_old[j][k] + (NAC_new[j][k] - NAC_old[j][k]) * factor;
    // * compute force vector (dydt) used to update y (q,p) using H, V, and NAC
    // Here, y[0, ..., 2DOFe-1] <-> (q[0, ..., DOFe-1], p[0, ..., DOFe-1])
    // and, dydt[0, ..., 2DOFe-1] <-> (dqdt[0, ..., DOFe-1], dpdt[0, ..., DOFe-1])
    // Ref: Miller, J. Chem. Phys. 147, 064112 (2017), Eqs. 12a and 12b
    for (int i = 0; i < DOFe; ++i)
        for (int j = 0; j < DOFe; ++j)
            if (i != j) { // if i == j, d_ji is 0, deltaE is 0, and term is 0.
                // index of d_ji in nonadiabatic coupling matrix is
                // (2*DOFe-j-1)*j/2 + (i-j) -1 (j < i).
                // Note, only d_ji (j < i) element is stored in NAC
                // So, to get the d_ji (j > i) by: d_ji = -d_ij
                int sign, index; // sign of d_ji and index of d
                if (j < i) { // j < i: d_ji = d_ji
                    sign = 1;
                    index = (2*DOFe-j-1)*j/2+(i-j)-1; // index of d_ji
                }
                else { // j > i: d_ji = - d_ij
                    sign = -1;
                    index = (2*DOFe-i-1)*i/2+(j-i)-1; // index of d_ij
                }
                double d_dot_V = 0;
                for (int k = 0; k < V_prime.size(); ++k)  // k is index of DOFn
                    d_dot_V += sign * NAC_prime[index][k] * V_prime[k];
                dydt[i] += y[j] * d_dot_V; // y[j] -> q_j
                dydt[DOFe+i] += y[DOFe+j] * d_dot_V; // y[DOFe+j] -> p_j
                double deltaE = (H_prime[i][i] - H_prime[j][j]) / DOFe;
                dydt[i] += y[DOFe+i] * deltaE; // y[DOFe+i] -> p_i
                dydt[DOFe+i] -= y[i] * deltaE; // y[i] -> q_i
            }
}
// =============================================================================

// # The following three functions are used to propagate electronic DOF (q, p)
// # with 4th-order Runge-Kutta method in diabatic basis.
void DynamicsElec::MOVE_elec(double t_nucl, const Real_Matrix& H_old, const Real_Matrix& H_new, std::vector<double>& q, std::vector<double>& p) {
    // cast y[0..2DOFe-1] <- (q[0..DOFe-1], p[0..DOFe-1])
    std::vector<double> y(2*DOFe, 0.0);
    for (int i = 0; i < DOFe; i++) {
        y[i] = q[i];
        y[DOFe + i] = p[i];
    }
    // propagate EPN steps of dt for electronic DOF
    for (int i = 0; i < EPN; i++) {
        // t_nucl is the current time of nuclear propagation (the time of H_old)
        // t_elec is the current time of electronic propagation
        double t_elec = t_nucl + i * dt;
        // explicit classic Runge-Kutta 4th order algorithm for ODE dy/dt = f(t,y)
        rk4(t_nucl, t_elec, H_old, H_new, y);
    }
    // cast back y[0..2DOFe-1] -> (q[0..DOFe-1], p[0..DOFe-1])
    for (int i = 0; i < DOFe; i++) {
        q[i] = y[i];
        p[i] = y[DOFe + i];
    }
}

void DynamicsElec::rk4(double t_nucl, double t_elec, const Real_Matrix& H_old, const Real_Matrix& H_new, std::vector<double>& y) {
    const int n = y.size();
    std::vector<double> f1(n, 0);
    std::vector<double> f2(n, 0);
    std::vector<double> f3(n, 0);
    std::vector<double> f4(n, 0);
    std::vector<double> yt(n, 0);
    // [1] k1 = dt * f(t,y),  f(t, y) = dydt
    deriv(t_nucl, t_elec, H_old, H_new, y, f1);
    // yt = y + k1 / 2
    for (int i = 0; i < n; i++)
        yt[i] = y[i] + 0.5 * dt * f1[i];
    // [2] k2 = dt * f2, where f2 = f(t+dt/2 , y + k1 / 2)
    deriv(t_nucl, t_elec+0.5*dt, H_old, H_new, yt, f2);
    // yt = y + k2 / 2
    for (int i = 0; i < n; i++)
        yt[i] = y[i] + 0.5 * dt * f2[i];
    // [3] k3 = dt * f3, where f3 = f(t+dt/2, y + k2 /2)
    deriv(t_nucl, t_elec+0.5*dt, H_old, H_new, yt, f3);
    // yt = y + k3
    for (int i = 0; i < n; i++)
        yt[i] = y[i] + dt * f3[i];
    // [4] k4 = dt * f(t+dt, y + k3)
    deriv(t_nucl, t_elec+dt, H_old, H_new, yt, f4);
    // At last, accumulate k1-k4 with coefficients
    // y[t+dt] = y[t] + 1/6 * (k1 + 2k2 + 2k3 + k4), k = dt * f
    for (int i = 0; i < n; i++)
        y[i] += dt/6.0 * (f1[i] + 2.0 * f2[i] + 2.0 * f3[i] + f4[i]);
}

void DynamicsElec::deriv(double t_nucl, double t_prime, const Real_Matrix& H_old, const Real_Matrix& H_new,
                         const std::vector<double>& y, std::vector<double>& dydt) {
    // In a single nuclear step, electronic DOF can be updated several times,
    // i.e., electronic time step (dt): 10*dt = DT (nuclear time step). Here,
    // use linear interpolation to estimate this H tilde in the intermediate time.
    // H tilde at t+xDT (H_avg removed)
    Real_Matrix H_eff(DOFe, std::vector<double>(DOFe, 0));
    // Runge Kutta driver: computeDerivatives  f = deriv(t, y[], dydt[])
    for (int j = 0; j < DOFe; j++)
        for (int k = 0; k < DOFe; k++) {
            // Use linear interpolation to estimate H tilde in the intermediate time (tprime).
            // Here, step is the current step of nuclear propagation. (t'- t)/ DT
            H_eff[j][k] = H_old[j][k] + (H_new[j][k] - H_old[j][k]) * (t_prime - t_nucl) / (EPN * dt);
            // dy_j/dt = dq_j/dt = 1/hbar * sum_k{H_eff[j][k] * p_k}, (p_k=y[DOFe+k])
            // dy_(DOFe+j)/dt = dp_j/dt = -1/hbar * sum_k{H_eff[j][k] * q_k}, (q_k=y[k])
            dydt[j] += 1.0/hbar * H_eff[j][k] * y[DOFe + k];
            dydt[DOFe + j] -= 1.0/hbar * H_eff[j][k] * y[k];
        }
}
// =============================================================================


// # The following three functions are used to propagate coefficient of electronic
// # wavefunction (coeff) with 4th-order Runge-Kutta method in diabatic basis.
void DynamicsElec::MOVE_elec(double t_nucl, const Real_Matrix& H_old, const Real_Matrix& H_new,
                             std::vector<double>& q, std::vector<double>& p, std::vector<Complex>& coeff) {
    // propagate EPN steps of dt for electronic DOF
    for (int i = 0; i < EPN; i++) {
        // t_nucl is the current time of nuclear propagation (the time of H_old)
        // t_elec is the current time of electronic propagation
        double t_elec = t_nucl + i * dt;
        // explicit classic Runge-Kutta 4th order algorithm for ODE dy/dt = f(t,y)
        rk4(t_nucl, t_elec, H_old, H_new, coeff);
    }
    // update q, p from coeff
    for (int i = 0; i < DOFe; i++) {
        q[i] = coeff[i].real() * sqrt(2);
        p[i] = coeff[i].imag() * sqrt(2);
    }
}

void DynamicsElec::rk4(double t_nucl, double t_elec, const Real_Matrix& H_old, const Real_Matrix& H_new, std::vector<Complex>& y) {
    const int n = y.size();
    std::vector<Complex> f1(n, 0);
    std::vector<Complex> f2(n, 0);
    std::vector<Complex> f3(n, 0);
    std::vector<Complex> f4(n, 0);
    std::vector<Complex> yt(n, 0);
    // [1] k1 = dt * f(t,y),  f(t, y) = dydt
    deriv(t_nucl, t_elec, H_old, H_new, y, f1);
    // yt = y + k1 / 2
    for (int i = 0; i < n; i++)
        yt[i] = y[i] + 0.5 * dt * f1[i];
    // [2] k2 = dt * f2, where f2 = f(t+dt/2 , y + k1 / 2)
    deriv(t_nucl, t_elec+0.5*dt, H_old, H_new, yt, f2);
    // yt = y + k2 / 2
    for (int i = 0; i < n; i++)
        yt[i] = y[i] + 0.5 * dt * f2[i];
    // [3] k3 = dt * f3, where f3 = f(t+dt/2, y + k2 /2)
    deriv(t_nucl, t_elec+0.5*dt, H_old, H_new, yt, f3);
    // yt = y + k3
    for (int i = 0; i < n; i++)
        yt[i] = y[i] + dt * f3[i];
    // [4] k4 = dt * f(t+dt, y + k3)
    deriv(t_nucl, t_elec+dt, H_old, H_new, yt, f4);
    // At last, accumulate k1-k4 with coefficients
    // y[t+dt] = y[t] + 1/6 * (k1 + 2k2 + 2k3 + k4), k = dt * f
    for (int i = 0; i < n; i++)
        y[i] += dt/6.0 * (f1[i] + 2.0 * f2[i] + 2.0 * f3[i] + f4[i]);
}

void DynamicsElec::deriv(double t_nucl, double t_prime, const Real_Matrix& H_old, const Real_Matrix& H_new,
                         const std::vector<Complex>& coeff, std::vector<Complex>& f) {
    // In a single nuclear step, electronic DOF can be updated several times,
    // i.e., electronic time step (dt): 10*dt = DT (nuclear time step). Here,
    // use linear interpolation to estimate this H tilde in the intermediate time.
    // H tilde at t+xDT (H_avg removed)
    Real_Matrix H_eff(DOFe, std::vector<double>(DOFe, 0));
    // Runge Kutta driver: computeDerivatives  f = deriv(t, y[], dydt[])
    for (int j = 0; j < DOFe; j++) {
        for (int k = 0; k < DOFe; k++) {
            // Use linear interpolation to estimate H tilde in the intermediate time (tprime).
            // Here, step is the current step of nuclear propagation. (t'- t)/ DT
            H_eff[j][k] = H_old[j][k] + (H_new[j][k] - H_old[j][k]) * (t_prime - t_nucl) / (EPN * dt);
            // dc_j/dt = -i/hbar * sum_k{H_eff[j][k] * c_k}
            f[j] += H_eff[j][k] * coeff[k];
        }
        f[j] *= -I/hbar;
    }
}
// =============================================================================


// # The following function is used to propagate coefficient of electronic wavefunction
// # (coeff) via the diagonalization of Hamiltonian in diabatic basis.
void DynamicsElec::MOVE_elec(const Real_Matrix& H, std::vector<double>& q, std::vector<double>& p, std::vector<Complex>& coeff) {
    // Ref: Report 28 Ehrenfest dynamics, Section 2.1 and 2.2
    // Exact propagation of electronic DOF via the diagonalization of Hamiltonian
    // For two-level system, the propagation of the electronic DOF by diagonalization
    // of Hamiltonian can be done analytically.
    if (DOFe == 2) {
        // Get Gamma_DA, average enenrgy, energ differnece from H_diabatic
        const double Gamma = H[0][1];
        const double E_avg = 0.5 * (H[0][0] + H[1][1]);
        const double E_diff = H[0][0] - H[1][1];
        // Get adiabatic energy (eigenvalue of H_diabatic matrix)  (eq.31-32)
        // Here, E_plus is V_e (excited state), and E_minus is V_g (ground state).
        const double rt = sqrt(E_diff*E_diff + 4*Gamma*Gamma);
        const double E_minus = E_avg - rt * 0.5;
        const double E_plus = E_avg + rt * 0.5;
        // Get unitary matrix T from Hamiltonion matrix (equation 23-25)
        Complex_Matrix T(DOFe, std::vector<Complex>(DOFe, 0.0));
        // alpha = 0.5 * theta
        // The value of atan() in C++ is [-0.5pi, 0.5pi]
        double alpha = 0.5 * atan(2*Gamma / E_diff);
        // if V_D < V_A, alpha ~ [-0.25pi, 0], then add 0.5pi to alpha, then
        // cos(alpha+0.5pi) = -sin(alpha), sin(alpha+0.5pi) = cos(alpha)
        if (E_diff < 0)
            alpha += 0.5*pi;
        T[0][0] = cos(alpha);
        T[0][1] = -sin(alpha);
        T[1][0] = sin(alpha);
        T[1][1] = cos(alpha);
        // Get conjugate transpose (Hermitian transpose) matrix of T:
        // T_dagger[j][k] = conj(T[k][j]) (equation 26)
        Complex_Matrix T_dagger(DOFe, std::vector<Complex>(DOFe, 0.0));
        T_dagger[0][0] = T[0][0];
        T_dagger[0][1] = T[1][0];
        T_dagger[1][0] = T[0][1];
        T_dagger[1][1] = T[1][1];
        // Get D matrix (diag{e^{-i*Ep*dt} , e^{-i*Em*dt}}) (equation 27)
        // Here, dt is electronic time step size
        // E_plus and E_minus are eigenvalues.
        Complex_Matrix D(DOFe, std::vector<Complex>(DOFe, 0.0));
        // Here, dt is electronic time step size,
        // EPN = DT/dt, is electronic steps per nuclear step.
        D[0][0] = exp(-I * (double)EPN * dt * E_plus);
        D[1][1] = exp(-I * (double)EPN * dt * E_minus);
        // U(t+dt, t) = T * diag(e^{-i Ep DT} , e^{-i Ep DT}) * T_dagger (equation 27)
        // Compute tmp = T * D
        Complex_Matrix tmp(DOFe, std::vector<Complex>(DOFe, 0.0));
        for (int i = 0; i < DOFe; i++)
            for (int j = 0; j < DOFe; j++)
                for (int k = 0; k < DOFe; k++)
                    tmp[i][j] += T[i][k] * D[k][j];
        // Compute U = tmp * T_dagger
        Complex_Matrix U(DOFe, std::vector<Complex>(DOFe, 0.0));
        for (int i = 0; i < DOFe; i++)
            for (int j = 0; j < DOFe; j++)
                for (int k = 0; k < DOFe; k++)
                    U[i][j] += tmp[i][k] * T_dagger[k][j];
        // Compute coeff_new = U(matrix) * coeff(vector) (equation 40)
        std::vector<Complex> coeff_new(DOFe, 0.0);
        for (int i = 0; i < DOFe; i++)
            for (int j = 0; j < DOFe; j++)
                coeff_new[i] += U[i][j] * coeff[j];
        // Update coeff
        for (int i = 0; i < DOFe; i++)
            coeff[i] = coeff_new[i];
    }
    // For multi-level systems, the diagonalization has to be calculated numerically
    // and similar time evolution operator can be constructed from inverse transformation
    // of diagonal matrix with adiabatic energies.
    else {
        // using matrix from Eigen lib: MatrixXcd, cd means complex double,
        // X means dynamics (unknown size before runnning).
        // others: d means double, 2d means 2*2 real matrix.
        // The initial values are zero.
        Eigen::MatrixXcd h(DOFe, DOFe);
        // Initialize h by assigning values from Hamiltonian H (C++ vector<vector>)
        // Hamiltonian H is Hermite matrix or self-adjoint matrix
        for (int j = 0; j < DOFe; j++)
            for (int k = 0; k < DOFe; k++)
                h(j, k) = H[j][k];
        // Construct an object of SelfAdjointEigenSolver to compute eigenvalues.
        // For Hermite matrix, SelfAdjointEigenSolver is faster and more accurate
        // than the general purpose eigenvalue algorithms that implemented in
        // EigenSolver and ComplexEigenSolver.
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eig(h);
        // Get eigenvalues (real number) of Hamiltonian h
        const Eigen::VectorXd& eigenvalues = eig.eigenvalues();
        // Get unitary matrix T of Hamiltonian h
        // The columns of T are the coefficients of eigenvectors
        const Eigen::MatrixXcd& T = eig.eigenvectors();
        // Compute diagonal values e^{-i*dt/hbar*eigenvalues}
        // Here, dt is electronic time step size,
        // EPN = DT/dt, is electronic steps per nuclear step.
        const Eigen::VectorXcd& D = (-I*(double)EPN*dt/hbar*eigenvalues).array().exp();
        // Compute time evolution operator: U(t+dt, t) = T * D * T_dagger
        // D is diagonal matrix (diag{e^{-I*dt/hbar*eigenvalues[i]},...)
        // Get diagonal matrix D by D.asDiagonal()
        // T_dagger is conjugate transpose matrix of T
        // For unitary matrix, conjugate transpose = inverse = adjoint matrix.
        const Eigen::MatrixXcd U = T * D.asDiagonal() * T.adjoint();
        // Initialize coeff_old by assigning values from coeff (C++ vector<complex>)
        Eigen::VectorXcd coeff_old(DOFe);
        for (int i = 0; i < DOFe; i++)
                coeff_old(i) = coeff[i];
        // Compute coeff_new = U(matrix) * coeff_old(vector<complex>) in diabatic base
        Eigen::VectorXcd coeff_new = U * coeff_old;
        // Cast coeff_new back to coeff (C++)
        for (int i = 0; i < DOFe; i++)
            coeff[i] = coeff_new[i];
    }
    // update q, p from coeff
    for (int i = 0; i < DOFe; i++) {
        q[i] = coeff[i].real() * sqrt(2);
        p[i] = coeff[i].imag() * sqrt(2);
    }
}
// =============================================================================


// # The following function is used to propagate reduced density matrix (sigma)
// # via the diagonalization of Hamiltonian in diabatic basis. (Ehrenfest only)
void DynamicsElec::MOVE_elec(const Real_Matrix& H, Complex_Matrix& sigma) {
    // Ref: Xiang's Report 28 Ehrenfest dynamics, Section 2.1 and 2.2
    // Exact propagation of electronic DOF via the diagonalization of Hamiltonian
    // For two-level system, the propagation of the electronic DOF by diagonalization
    // of Hamiltonian can be done analytically.
    // Note that the order of state is {-, +}, which is different from Report 28,
    // since the order using numerical way to do diagonalization is ascending.
    if (DOFe == 2) {
        // Get Gamma (coupling), average enenrgy, energ differnece from H_diabatic
        const double Gamma = H[0][1];
        const double E_avg = 0.5 * (H[0][0] + H[1][1]);
        const double E_diff = H[0][0] - H[1][1];
        // Get adiabatic energy (eigenvalue of H_diabatic matrix)  (eq.31-32)
        // Here, E_plus is V_e (excited state), and E_minus is V_g (ground state).
        const double rt = sqrt(E_diff*E_diff + 4*Gamma*Gamma);
        const double E_minus = E_avg - rt * 0.5;
        const double E_plus = E_avg + rt * 0.5;
        // Get rotation angle alpha in rotation matrix (eq.33)
        // The range of value of atan() in C++ is [-0.5pi, 0.5pi]
        double alpha = 0.5 * atan(2*Gamma / E_diff);
        // if V_D < V_A, alpha ~ [-0.25pi, 0], then add 0.5pi to alpha, then
        // cos(alpha+0.5pi) = -sin(alpha), sin(alpha+0.5pi) = cos(alpha)
        if (E_diff < 0)
            alpha += 0.5*pi;
        // Get unitary rotation matrix (transformation matrix) (Eq. 34)
        // Note that, here, T = (c_minus, c_plus), differnent order from eq.34
        Complex_Matrix T(DOFe, std::vector<Complex>(DOFe, 0.0));
        T[0][0] = -sin(alpha);
        T[1][0] = cos(alpha);
        T[0][1] = cos(alpha);
        T[1][1] = sin(alpha);
        // Get conjugate transpose (Hermitian transpose) matrix of T:
        // T_dagger[j][k] = conj(T[k][j]). Here, T is real matrix due to
        // coupling is real number, and T[0][1] = T[1][0], so T_dagger = T.
        Complex_Matrix T_dagger = T;
        // Get D matrix (diag{e^{±i*Em*dt}, e^{±i*Ep*dt}}) (Eq. 27,28)
        // Here, dt is electronic time step size,
        // EPN = DT/dt, is electronic steps per nuclear step.
        // Usually, dt is much smaller than nuclear time step DT, e.g., EPN=10
        // However, T is determined by nuclear position R(t), there is no need
        // to propagate multiple times for one cycle of updating nuclear, since
        // the exact propagator is equivalent. Ref: Eq. 29
        Complex_Matrix D(DOFe, std::vector<Complex>(DOFe, 0.0));
        Complex_Matrix D_conj = D;
        D[0][0] = exp(-I * (double)EPN * dt * E_minus);
        D[1][1] = exp(-I * (double)EPN * dt * E_plus);
        D_conj[0][0] = exp(I * (double)EPN * dt * E_minus);
        D_conj[1][1] = exp(I * (double)EPN * dt * E_plus);
        // U(t+dt,t) = T * diag(e^{-i Ep DT} , e^{-i Ep DT}) * T_dagger (Eq. 27)
        // U_dagger(t+dt,t) = T * diag(e^{i Ep DT} , e^{i Ep DT}) * T_dagger (Eq. 28)
        Complex_Matrix tmp, U, U_dagger;
        Matrix_Multiply(T, D, tmp);
        Matrix_Multiply(tmp, T_dagger, U);
        Matrix_Multiply(T, D_conj, tmp);
        Matrix_Multiply(tmp, T_dagger, U_dagger);
        // Update reduced density matrix by: U(t+dt,t)*σ(t)*U_dagger(t+dt,t) (Eq.18)
        Matrix_Multiply(U, sigma, tmp);
        Matrix_Multiply(tmp, U_dagger, sigma);
    }
    // For multi-level systems, the diagonalization has to be calculated numerically
    // and similar time evolution operator can be constructed from inverse transformation
    // of diagonal matrix with adiabatic energies.
    else {
        // using matrix from Eigen lib: MatrixXcd, cd means complex double,
        // X means dynamics (unknown size before runnning), d means double,
        // 2d means 2*2 real matrix. The initial values are zero.
        Eigen::MatrixXd h(DOFe, DOFe);
        // Initialize h by assigning values from diabatic Hamiltonian H
        // (C++ vector<vector>). It is an Hermite matrix or self-adjoint matrix,
        // which means the conjugate transpose matrix of it is itself.
        for (int j = 0; j < DOFe; j++)
            for (int k = 0; k < DOFe; k++)
                h(j, k) = H[j][k];
        // Construct an object of SelfAdjointEigenSolver to compute eigenvalues.
        // For Hermite matrix, SelfAdjointEigenSolver is faster and more accurate
        // than the general purpose eigenvalue algorithms that implemented in
        // EigenSolver and ComplexEigenSolver. And SelfAdjointEigenSolver will
        // sort eigenvalues ascending, since all eigenvalues are real numbers.
        // Moreover, it can make sure the eigenvectors are orthogonal.
        // Note that in EigenSolver or ComplexEigenSolver, the eigenvalues are
        // not sorted in any particular order.
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(h);
        // Get eigenvalues (ascending) of Hamiltonian h.
        // They are diganol elements of adiabatic Hamiltonian.
        // Note the order of eigenvalues is energy ascending.
        const Eigen::VectorXd& eigenvalues = eig.eigenvalues();
        // Get unitary rotation matrix T of Hamiltonian
        // The columns of T are the coefficients of eigenvectors
        const Eigen::MatrixXd& T = eig.eigenvectors();
        // Get diagonal values e^{±i*dt/hbar*eigenvalues} of diagonal matrix
        // Here, dt is electronic time step size,
        // EPN = DT/dt, is electronic steps per nuclear step.  (Eq. 27,28)
        const Eigen::VectorXcd& D = (-I*(double)EPN*dt/hbar*eigenvalues).array().exp();
        const Eigen::VectorXcd& D_conj = (I*(double)EPN*dt/hbar*eigenvalues).array().exp();
        // Compute time evolution operator: U(t+dt, t) = T * D * T_dagger (Eq.27)
        // And U_dagger(t+dt, t) = T * D_conj * T_dagger (Eq.28)
        // D is diagonal matrix (diag{e^{-I*dt/hbar*eigenvalues[i]},...)
        // Get diagonal matrix D by D.asDiagonal()
        // T_dagger is conjugate transpose matrix of T
        // For unitary matrix, conjugate transpose = inverse = adjoint matrix.
        const Eigen::MatrixXcd U = T * D.asDiagonal() * T.adjoint();
        const Eigen::MatrixXcd U_dagger = T * D_conj.asDiagonal() * T.adjoint();
        // Update reduced density matrix by: U(t+dt,t)*σ(t)*U_dagger(t+dt,t) (Eq.18)
        // Initialize sigma_old by assigning values from sigma (C++ vector<Complex>)
        Eigen::MatrixXcd sigma_old(DOFe, DOFe);
        for (int j = 0; j < DOFe; j++)
            for (int k = 0; k < DOFe; k++)
                sigma_old(j, k) = sigma[j][k];
        Eigen::MatrixXcd sigma_new = U * sigma_old * U_dagger;
        // Cast sigma_new back to sigma (C++ vector<Complex>)
        for (int j = 0; j < DOFe; j++)
            for (int k = 0; k < DOFe; k++)
                sigma[j][k] = sigma_new(j, k);
    }
}