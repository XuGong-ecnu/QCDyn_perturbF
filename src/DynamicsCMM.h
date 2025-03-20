/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "DynamicsElec.h"
#include "DynamicsMQCBase.h"

/**
 * DynamicsCMM class defines nonadiabatic dynamics with phase space classical
 * mapping models (CMMs) developed by Jian Liu [1-2].
 *
 * It includes 5 CMM methods (CMM1, and CMM3-6). All of them use 4 electronic
 * mapping variables, so unlike other mapping dynamics with 2 mapping variables
 * they can not use the general functions for that calculation of effective force
 * and electronic propagation.
 *
 * The five methods (CMM1, and CMM3-6) share same initial sampling method,
 * while the definations of population or coherence operator may different.
 * And 5 functions for electronic propagation are spectially wrriten for them.
 * Although CMM3-6 can use the same algorithm for electronic propagation as CMM1,
 * but the performance is not as good as special one. The effective force is
 * defined from its Hamiltonian just like mean-filed method but with different
 * population or coherence operator.
 *
 * Although, these methods use 4 electronic mapping variables, only 2 variables
 * are indenpendent in initial sampling.
 *
 * And the CMM2 method with 2 mapping variables are defined in DynamicsECMM class.
 * Since the CMM2 method a special case (γ=0) of eCMM method [3]. And since they
 * are using general two mapping variables, the general function for force and
 * electronic propagation (RK4) can be applied like other mapping dynamcis.
 *
 * References:
 * [1] J. Liu, J. Chem. Phys. 2016, 145, 204105.
 * [2] X. He, J. Liu, J. Chem. Phys. 2019, 151, 024105.
 * [3] X. He, Z. Gong, B. Wu, J. Liu, J. Phys. Chem. Lett. 2021, 12, 2496.
 */
class DynamicsCMM : public DynamicsMQCBase {
public:
    /**
     * Construct a DynamicsCMM object.
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     */
    DynamicsCMM(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsMQCBase(param, Ha) {}
    ~DynamicsCMM() {}
    /**
     * Initialize data members and check the legality of parameters.
     */
    void init();
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
     * Get initial electronic mapping variables.
     * Similar as full-sphere sampling of spin-maping, but for 4 mapping
     * variables, while only 2 variables are independent.
     */
    void samplingElec();
    /**
     * Get initial denisty (population or coherence operator) at time 0.
     *
     * The results of it are stored in data member: init_density
     */
    void getInitialDensity();
    /**
     * Update the electronic reduced density matrix at current time.
     *
     * Both the current (RDM_current) and average (RDM_average) ones are updated.
     */
    void updateDensityMatrix();
    /**
     * Propagate one step for a simulation of model using velocity Verlet integrator.
     *
     * Same as the general one that defined in DynamicsMQCBase() but using
     * different MOVE_elec() only.
     * TODO: perhaps can use the general one by taking MOVE_elec as parameter
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
    /**
     * Do electronic propagation for 4 electronic mapping variables.
     * Although CMM3-6 can use the same algorithm for electronic propagation as
     * CMM1, but the performance is not as good as special one.
     *
     * For CMM2 method, the original numerical algorithm to do electronic
     * propagation from Jian Liu's paper is also defined here, but we never use.
     * Since the CMM2 method is the special case of eCMM (gamma=0) method, and
     * like other mapping dynamics, the RK4 method will be used.
     * TODO: Support multi-steps
     */
    void MOVE_elec_CMM1();
    void MOVE_elec_CMM2();
    void MOVE_elec_CMM3();
    void MOVE_elec_CMM4();
    void MOVE_elec_CMM5();
    void MOVE_elec_CMM6();
    /**
     * Get population (j == k) or coherence (j !=k) operator according to
     * CMMs method and current mapping variables.
     *
     * @return complex value of population or coherence
     */
    Complex getOperator(int j, int k);

private:
    // std::function is a polymorphic function wrapper. It stores a pointers to
    // member functions: MOVE_elec_CMMX(), X is 1, 3-6. Since the CMMs methods
    // use different electronic propagation. Note that CMM3-6 can use the same
    // algorithm as CMM1, but the performance is not as good as special one.
    // This binding is finished at init().
    std::function<void()> MOVE_elec;
};