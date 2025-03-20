/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Sep. 02, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "Tools.h"
#include "Parameters.h"
#include "HamiltonianElec.h"

/**
 * This class defines the method of electronic propagation.
 *
 * Now the 4th-order Runge-Kutta method and exact propagation of electronic DOF
 * via the diagonalization of Hamiltonian are implemented.
 *
 * It should be noted that all the units used in the electronic propagation are
 * atomic units (a.u.).
 */
class DynamicsElec {
public:
    /**
     * Construct a DynamicsElec object.
     *
     * @param param   the global paramters
     */
    DynamicsElec(Parameters& param) : param(param) {}
    ~DynamicsElec() {}
    /**
     * Initialize data members and check the legality of parameters.
     */
    void init();
    /**
     * Do electronic propagation (electronic mapping variable positions and momenta)
     * with 4th-order Runge-Kutta method.
     *
     * @param t_nucl current nuclear time in au
     * @param H_old  effective energy matrix (H tilde) at time t (Havg removed)
     * @param H_new  effective energy matrix (H tilde) at time t+DT (Havg removed)
     * @param p      electronic mapping variable positions [update]
     * @param q      electronic mapping variable momenta [update]
     */
    void MOVE_elec(double t_nucl, const Real_Matrix& H_old, const Real_Matrix& H_new, std::vector<double>& q, std::vector<double>& p);
    /**
     * Do electronic propagation (coefficient of electronic wavefunction)
     * with 4th-order Runge-Kutta method.
     *
     * It will also update the q, p from new coeff.
     *
     * @param t_nucl current nuclear time in au
     * @param H_old  effective energy matrix (H tilde) at time t (Havg removed)
     * @param H_new  effective energy matrix (H tilde) at time t+DT (Havg removed)
     * @param p      electronic mapping variable positions [update]
     * @param q      electronic mapping variable momenta [update]
     * @param coeff  coefficient of electronic wavefunction, coeff[j]=(q[j] + I * p[j])/sqrt(2) [update]
     */
    void MOVE_elec(double t_nucl, const Real_Matrix& H_old, const Real_Matrix& H_new,
                   std::vector<double>& q, std::vector<double>& p, std::vector<Complex>& coeff);
    /**
     * Do electronic propagation by diagonalize the Hamiltonian.
     * It will also update the q, p from new coeff.
     *
     * For two-level system, the propagation can be done analytically and exactly.
     * For multi-level systems, the diagonalization has to be calculated numerically.
     *
     * Note that the Hamiltonina used here is original Hamiltonian, not effective
     * Hamiltonian (Havg removed) that used in RK4 method. However, using effective
     * Hamiltonian also works (eigenvalues will shift, but eigenvectors are same)
     * and will produce same result for Ehrenfest or FSSH dynamics.
     *
     * @param H      Hamiltonian at time t
     * @param p      electronic mapping variable positions [update]
     * @param q      electronic mapping variable momenta [update]
     * @param coeff  coefficient of electronic wavefunction, coeff[j]=(q[j] + I * p[j])/sqrt(2) [update]
     */
    void MOVE_elec(const Real_Matrix& H, std::vector<double>& q, std::vector<double>& p, std::vector<Complex>& coeff);
    /**
     * Do electronic reduced density matrix propagation by diagonalize the Hamiltonian
     * (Used for Ehrenfest method only).
     *
     * For two-level system, the propagation can be done analytically and exactly.
     * For multi-level systems, the diagonalization has to be calculated numerically.
     *
     * Note that the Hamiltonina used here is original Hamiltonian, not effective
     * Hamiltonian (Havg removed) that used in RK4 method. However, using effective
     * Hamiltonian also works (eigenvalues will shift, but eigenvectors are same)
     * and will produce same result for Ehrenfest or FSSH dynamics.
     *
     * @param H      Hamiltonian at time t
     * @param sigma  electronic reduced density matrix [update]
     */
    void MOVE_elec(const Real_Matrix& H, Complex_Matrix& sigma);
    // # The following function is used for model in adiabatic basis
    /**
     * Do electronic propagation (electronic mapping variable positions and momenta)
     * with 4th-order Runge-Kutta method in aidiabatic basis.
     *
     * It will also update the coeff from q, p.
     *
     * @param t_nucl   current nuclear time in au
     * @param H_old    adibatic energy matrix at time t
     * @param H_new    adibatic energy matrix at time t+DT (DT is nuclear time step)
     * @param V_old    velocites at time t
     * @param V_new    velocites at time t+DT
     * @param NAC_old  nonadiabatic coupling (NAC) at time t
     * @param NAC_new  nonadiabatic coupling (NAC) at time t+DT
     * @param p        electronic mapping variable positions [update]
     * @param q        electronic mapping variable momenta [update]
     * @param coeff    coefficient of electronic wavefunction, coeff[j]=(q[j] + I * p[j])/sqrt(2) [update]
     */
    void MOVE_elec(double t_nucl, const Real_Matrix& H_old, const Real_Matrix& H_new,
                   const std::vector<double>& V_old, const std::vector<double>& V_new,
                   const Real_Matrix& NAC_old, const Real_Matrix& NAC_new,
                   std::vector<double>& q, std::vector<double>& p, std::vector<Complex>& coeff);

private:
    /**
     * Explicit classic Runge-Kutta 4th order algorithm for ODE dy/dt = f(t,y).
     *
     * This is an internal function of MOVE_elec().
     *
     * @param t_nucl current nuclear time in au
     * @param t_elec current electronic time, t_elec = t_nucl + i*dt
     * @param H_old  effective energy matrix (H tilde) at time t (Havg removed)
     * @param H_new  effective energy matrix (H tilde) at time t+DT (Havg removed)
     * @param y      y[2*DOFe] = (q[DOFe], p[DOFe])
     */
    void rk4(double t_nucl, double t_elec, const Real_Matrix& H_old, const Real_Matrix& H_new, std::vector<double>& y);
    /**
     * Runge Kutta driver: computeDerivatives  f = deriv(t', y[], dydt[]).
     * This function will compute H_eff at tprime, and update dydt.
     *
     * This is an internal function of rk4().
     *
     * @param t_nucl  current nuclear time in au
     * @param t_prime intermediate electronic time. In a single nuclear step, electronic
     *                DOF can be updated several times, and use linear interpolation
     *                to estimate the H tilde at t+xDT (Havg removed)
     * @param H_old   effective energy matrix (H tilde) at time t (Havg removed)
     * @param H_new   effective energy matrix (H tilde) at time t+DT (Havg removed)
     * @param y       y[2*DOFe] = (q[DOFe], p[DOFe])
     * @param dydt    the derivative of y, 2*DOFe-dimensional force vector
     */
    void deriv(double t_nucl, double t_prime, const Real_Matrix& H_old, const Real_Matrix& H_new,
               const std::vector<double>& y, std::vector<double>& dydt);
    /**
     * Explicit classic Runge-Kutta 4th order algorithm for ODE dy/dt = f(t,y).
     *
     * The only difference is the data type of y, here is vector<Complex> (coeff).
     *
     * This is an internal function of MOVE_elec() [coeff].
     *
     * @param t_nucl current nuclear time au
     * @param t_elec current electronic time, t_elec = t_nucl + i*dt
     * @param H_old  effective energy matrix (H tilde) at time t (Havg removed)
     * @param H_new  effective energy matrix (H tilde) at time t+DT (Havg removed)
     * @param y      coefficient of electronic wavefunction, coeff[j]=(q[j] + I * p[j])/sqrt(2)
     */
    void rk4(double t_nucl, double t_elec, const Real_Matrix& H_old, const Real_Matrix& H_new, std::vector<Complex>& y);
    /**
     * Runge Kutta driver: computeDerivatives  f = deriv(t', y[], dydt[]).
     * This function will compute H_eff at tprime, and update dydt.
     *
     * This is an internal function of rk4() [coeff].
     *
     * @param t_nucl  current nuclear time au
     * @param t_prime intermediate electronic time. In a single nuclear step, electronic
     *                DOF can be updated several times, and use linear interpolation
     *                to estimate the H tilde at t+xDT (Havg removed)
     * @param H_old   effective energy matrix (H tilde) at time t (Havg removed)
     * @param H_new   effective energy matrix (H tilde) at time t+DT (Havg removed)
     * @param coeff   coefficient of electronic wavefunction, coeff[j]=(q[j] + I * p[j])/sqrt(2)
     * @param f       the derivative of coefficient, DOFe-dimensional force vector
     */
    void deriv(double t_nucl, double t_prime, const Real_Matrix& H_old, const Real_Matrix& H_new,
               const std::vector<Complex>& coeff, std::vector<Complex>& f);
    // # The following function is used for model in adiabatic basis
    /**
     * Explicit classic Runge-Kutta 4th order algorithm for ODE dy/dt = f(t,y).
     *
     * This is an internal function of MOVE_elec(). (aidiabatic basis)
     *
     * @param t_nucl   current nuclear time in au
     * @param t_elec   current electronic time, t_elec = t_nucl + i*dt
     * @param H_old    adibatic energy matrix at time t
     * @param H_new    adibatic energy matrix at time t+DT (DT is nuclear time step)
     * @param V_old    velocites at time t
     * @param V_new    velocites at time t+DT
     * @param NAC_old  nonadiabatic coupling (NAC) at time t
     * @param NAC_new  nonadiabatic coupling (NAC) at time t+DT
     * @param y        y[2*DOFe] = (q[DOFe], p[DOFe])
     */
    void rk4(double t_nucl, double t_elec, const Real_Matrix& H_old, const Real_Matrix& H_new,
             const std::vector<double>& V_old, const std::vector<double>& V_new,
             const Real_Matrix& NAC_old, const Real_Matrix& NAC_new,
             std::vector<double>& y);
    /**
     * Runge Kutta driver: computeDerivatives  f = deriv(t', y[], dydt[]).
     * This function will compute H, V, NAC at tprime, and update dydt.
     *
     * This is an internal function of rk4(). (aidiabatic basis)
     *
     * @param t_nucl   current nuclear time in au
     * @param t_prime  intermediate electronic time.
     * @param H_old    adibatic energy matrix at time t
     * @param H_new    adibatic energy matrix at time t+DT (DT is nuclear time step)
     * @param V_old    velocites at time t
     * @param V_new    velocites at time t+DT
     * @param NAC_old  nonadiabatic coupling (NAC) at time t
     * @param NAC_new  nonadiabatic coupling (NAC) at time t+DT
     * @param y        y[2*DOFe] = (q[DOFe], p[DOFe])
     * @param dydt     the derivative of y, 2*DOFe-dimensional force vector
     */
    void deriv(double t_nucl, double t_prime, const Real_Matrix& H_old, const Real_Matrix& H_new,
               const std::vector<double>& V_old, const std::vector<double>& V_new,
               const Real_Matrix& NAC_old, const Real_Matrix& NAC_new,
               const std::vector<double>& y, std::vector<double>& dydt);

private:
    // Parameters object controls the simulation
    Parameters& param;
    // electronic DOF = number of states/topologies/surfaces.
    // For mapping dynamics, DOFe should be greater than 1.
    int DOFe;
    // electronic steps per nuclear step.
    int EPN;
    // electronic time step size in au: dt = DT / EPN, DT is the nuclear time step
    // In a single nuclear step, electronic DOF can be updated several times,
    // such as 10*dt = DT, which means EPN = 10.
    double dt;
};