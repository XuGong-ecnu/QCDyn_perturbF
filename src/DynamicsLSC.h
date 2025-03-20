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
 * This class defines linearized semiclassical (LSC) mapping dynamics.
 *
 * For LSC mapping dynamics, there are 5 window types[1]: LSC1, LSC2, RI-LSC1-3.
 * Since the window type has no impact on the initial sampling or propagation (only
 * observables are different), all observables (reduced density matrix) corresonding
 * to these types will be computed in one simulation.
 *
 * Note that, the LSC1 method is also called PBME, and LSC2 method is also
 * called LSC-IVR in most papers.
 *
 * For the RI-LSC1/2/3 method are originally proposed by Rochardson, the identity
 * operators are used in the calculation of population observables (the coherence
 * is same as LSC2) to improve the accuracy of population and make the zero-point
 * energy parameter in the equation disappear [2]. RI-LSC1 is the "Single unity"
 * method (the modified verison of LSC1/PBME, sample from Φ) and RI-LSC3 is the
 * "Double unity" (the modified version of LSC2/LSC-IVR, sample from Φ^2) in
 * ref[2]. The methods are used in 7 state FMO model in ref[3].
 *
 * The sample width of initial electronic distribution can be ajusted by paramter
 * LSC_zeta0 [4] (called gamma in Xing's paper).
 *
 * And the reduced density matrices starting from different initial states can be
 * computed in one sumulation, since the initial sampling and propagation has no
 * relation to initial state.
 *
 * Reference:
 * [1] X. Gao, M. A. C. Saller, Y. Liu, A. Kelly, J. O. Richardson, E. Geva,
 *     J. Chem. Theory Comput. 2020, 16, 2883.
 * [2] M. A. C. Saller, A. Kelly, J. O. Richardson,
 *     J. Chem. Phys. 2019, 150, 071101.
 * [3] M. A. C. Saller, A. Kelly, J. O. Richardson,
 *     Faraday Discuss. 2020, 221, 150.
 * [4] X. Gao, E. Geva, J. Phys. Chem. A 2020, 124, 11006.
 */
class DynamicsLSC : public DynamicsMQCBase {
public:
    /**
     * Construct a DynamicsLSC object.
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     */
    DynamicsLSC(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsMQCBase(param, Ha) {}
    ~DynamicsLSC() {}
    /**
     * Initialize data members and check the legality of parameters.
     */
    void init();
    
private:
    /**
     * Get initial Gaussian sampling of electronic mapping variables.
     *
     * All windows types (LSC1/2, RI-LSC1~3) and initial states share same
     * initial Gaussian sampling.
     */
    void samplingElec();
    /**
     * Get initial windows/density of different types (LSC1/2, RI-LSC1~3) for
     * LSC method according to the init_state.
     *
     * The initial window [|j><k|]_W(q0,p0) is the window at time zero in equation
     * of electronic reduced density matrix but without the G(0) which has been
     * adsorbed in sampling.
     *
     * The results of it are stored in data member: init_XXX, (XXX is type)
     */
    void getInitialDensity();
    /**
     * Update the electronic reduced density matrix of different types at current time.
     * Note that the values of RDMs will be accumaleated and averaged.
     *
     * The results of it are stored in data member: RDM_XXX, (XXX is type)
     */
    void updateDensityMatrix();
        /*
    * Here we set up a function to get RDM data member.
    * Notice here RDM here serves as a two-dimensional matrix, RDM
    * Here we take out one of the matrix at given time snapshot, RDM[i][j*DOFe+k],
    * where the first dimensionality gives the position in the initial state(s), 
    * the second dimensionality gives the order of the state which gets correlation with
    * the initial state[i], j*DOFe+k --> |j><k|.
    * 
    */
    std::vector<std::vector<Complex>>& getRDMLSC1() {
        return RDM_LSC1;
    }
    std::vector<std::vector<Complex>>& getRDMLSC2() {
        return RDM_LSC2;
    }
    std::vector<std::vector<Complex>>& getRDMRILSC1() {
        return RDM_RILSC1;
    }
    std::vector<std::vector<Complex>>& getRDMRILSC2() {
        return RDM_RILSC2;
    }
    std::vector<std::vector<Complex>>& getRDMRILSC3() {
        return RDM_RILSC3;
    }

private:
    // zeta is the parameter to adjust width of the window functions used in
    // LSC method. the value of zeta should be greater than 0. And when zeta=1,
    // it is original LSC method. It has influence on the initial electronic
    // sampling and observables. The recommended value is 2.5 according the
    // publication of Xing Gao: J. Phys. Chem. A 2020, 124, 11006−11016, which
    // reproduce the population dynamics rather accurately for the spin-boson model.
    // In the paper of Xing Gao, it is called gamma.
    // Note that for RI-LSC1~3 method, we don't know if it is good or not to
    // introduce this paramter. But we must include it in the code, since we use
    // same initial electronic sampling (with LSC_zeta). For safety, it should be 1.
    double LSC_zeta;
    // Initial windows/density including all types used by LSC method.
    std::vector<Complex> init_LSC1, init_LSC2, init_RILSC1, init_RILSC2, init_RILSC3;
    // Reduced density matrix including all types genenrated by LSC method.
    std::vector<std::vector<Complex>> RDM_current_LSC1, RDM_current_LSC2, RDM_current_RILSC1, RDM_current_RILSC2, RDM_current_RILSC3, RDM_current;
    std::vector<std::vector<Complex>> RDM_average_LSC1, RDM_average_LSC2, RDM_average_RILSC1, RDM_average_RILSC2, RDM_average_RILSC3, RDM_average;
};