/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zengkui Liu @Sun Group @NYU-SH                                     *
 * Last updated: Dec. 2, 2021                                                 *
 * -------------------------------------------------------------------------- */

#pragma once
#include "DynamicsElec.h"
#include "DynamicsMQCBase.h"

/** 
 *  DynamicsECMMCV class defines non-adiabatic dynamics with extented phase space 
 *  classical mapping model with commutator variables (eCMMcv) developed by Jian 
 *  Liu [1], which is a further development of tunning ZPE parameter in the eCMM [2]
 *  family. 
 *  
 *  ZPE parameter, γ, in the eCMM Hamiltonian are treated as a matrix, Γ, in this 
 *  eCMMcv Hamiltonian. The element of Γ, Γ_nm, can be represented as a product of 
 *  a pair of commutator variables 
 *      Γnm = \sum_{k=1}^{DOFe} 0.5*s_k*(xcv_k^n + I*pcv_k^n)*(xcv_k^m - I*pcv_k^m).
 *
 *  Notice: instead of one ZPE parameter, here we employ (2 * DOFe * DOFe) variables to
 *  represent the ZPE parameter participating in the dynamics.
 * 
 *  The eCMMev method employs full-sphere sampling of electronic mapping variables,
 *  and ZPE commutator variables as we employed in the eCMM, CMM, and Spin-Mapping 
 *  method. 
 *    
 *  
 *  Reference: 
 * [1] X. He, B. Wu, Z. Gong, J. Liu, J. Phys. Chem. A. 2021, 125, 6845.
 * [2] X. He, Z. Gong, B. Wu, J. Liu, J. Phys. Chem. Lett. 2021, 12, 2496.
 */
class DynamicsECMMCV : public DynamicsMQCBase {
public: 
    /**
     * Construct a DynamicsECMM object.
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     */
     DynamicsECMMCV(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsMQCBase(param, Ha) {}
     ~DynamicsECMMCV() {}
    /**
     * Initialize data members and check the legality of parameters.
     */
    void init();
    /**
     * Get average electronic reduced density matrix (observables).
     *
     * This is a multi-frame reduced density matrix and each frame of it is
     * averaged by number of trajectories.
     *
     * @param rdm_out vector of reduced density matrix (output)
     */
    // void getDensityMatrix(std::vector<Complex_Matrix>& rdm_out);
    
    /*
    * Here we set up a function to get RDM data member.
    * Notice here RDM here serves as a two-dimensional matrix, RDM
    * Here we take out one of the matrix at given time snapshot, RDM[i][j*DOFe+k],
    * where the first dimensionality gives the position in the initial state(s), 
    * the second dimensionality gives the order of the state which gets correlation with
    * the initial state[i], j*DOFe+k --> |j><k|.
    * 
    */
    const std::vector<std::vector<Complex>>& getRDM() {
        return RDM;
    }


private:
    /**
     * Get initial electronic mapping variables and commutator variables.
     * Same as full-sphere smapling in DynamcisSPM for electric DOF.
     * Same full-sphere sampling for commutator variables.
     * Constraint is taken following the Eq. (43) in the reference
     * J. Phys. Chem. A. 2021, 125, 6845.
     * TODO: foucsed initial sampling
     */
    void samplingElec();
    /**
     * Get initial denisty (population or coherence operator) at time 0.
     * Same as initial denisty in DynamcisSPM.
     *
     * The results of it are stored in data member: init_density
     */
    void getInitialDensity();
    /**
     * Update the electronic reduced density matrix at current time.
     * Note that the values of RDMs will be accumaleated.
     *
     * The results of it are stored in data member: rdm
     */
    void updateDensityMatrix();
        /**
    * Propagate one step for a simulation of model using velocity Verlet integrator.
    *
    * Same as the general one that defined in DynamicsMQCBase() but using
    * different MOVE_elec() only.
    * TODO: apply RESPA algorithm on the electronic mapping variables
    * TODO: apply other multistep algorithm on electronic mapping variables
    */
    void oneStepForDiabaticModel() override;
    /**
     * Update the effective forces which will be used nuclear propagation.
     * This function is used for a simulation of model.
     *
     * For CMMs methods, it is similiar as the general one that defined in
     * DynamicsMQCBase() but using different population/coherence operator.
     */
    void updateDiabaticModelForces() override;
    // $ Model Hamiltonian with proprgation in adiabatic basis
    /**
     * Propagate one step for a simulation of model using velocity Verlet integrator.
     *
     * The EOMs used here are based on Meyer-Miller Hamiltonian with commutator 
     * variables in adiabatic representation proposed by Jian Liu
     *   [J. Phys. Chem. A 2021, 125, 6845.].
     */
    //void oneStepForAdiabaticModel() override;
    /**
     * Update the effective forces which will be used nuclear propagation.
     */
    //void updateAdiabaticModelForces() override;

    /**
     * Do electronic propagation for 2 electronic mapping variables and (2*DOFe*DOFe)
     * commutator variables. 
     * 
     * TODO: Support multi-steps
     */
    void MOVE_elec_eCMMcv();
    // # The following function is used for model in adiabatic basis
    /**
     * Do electronic propagation (electronic mapping variable positions and momenta)
     * with velocity-Verlet method in adiabatic basis.
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
    // void MOVE_elec_eCMMcv1();
    /**
     * Get population (j == k) or coherence (j !=k) operator according to
     * eCMMcv method, current mapping variables, and commutator variables.
     * 
     * @return complex value of population or coherence
     */
    Complex getOperator(int j, int k);
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
    // void rk4(double t_nucl, double t_elec, const Real_Matrix& H_old, const Real_Matrix& H_new,
    //        const std::vector<double>& V_old, const std::vector<double>& V_new,
    //        const Real_Matrix& NAC_old, const Real_Matrix& NAC_new,
    //        std::vector<double>& y);
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
    // void deriv(double t_nucl, double t_prime, const Real_Matrix& H_old, const Real_Matrix& H_new,
     //           const std::vector<double>& V_old, const std::vector<double>& V_new,
     //           const Real_Matrix& NAC_old, const Real_Matrix& NAC_new,
     //           const std::vector<double>& y, std::vector<double>& dydt);

private:
    // std::function is a polymorphic function wrapper. It stores a pointers to
    // member functions: MOVE_elec_CMMX(), X is 1, 3-6. Since the CMMs methods
    // use different electronic propagation. Note that CMM3-6 can use the same
    // algorithm as CMM1, but the performance is not as good as special one.
    // This binding is finished at init().
    std::function<void()> MOVE_elec;

    int traj_debug;
};