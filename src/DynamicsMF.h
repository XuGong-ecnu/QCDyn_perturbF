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
 * This class defines Ehrenfest mean-filed dynamics method (in mapping base).
 *
 * Ehrenfest mean-filed dynamics is not a mapping dynamics, but it is equivalent
 * to propagate classical {R, P, q, p}[1]. Therefore, we implement it here in
 * mapping base just like other mapping methods. Unlike to propagate the reduced
 * density matrix directly in original standard Ehrenfest, we propagtae electronic
 * {q,p} or coefficient of wavefunction with RK4 method like other mapping dynamics.
 * However, in this case, {q0, p0} is uniquely determined by coefficient of initial
 * state (not a distribution).
 *
 * It should be noted that the version of the MF method used here requires that
 * the initial electronic state is described by a pure population state, i.e.,
 * it can not run a simulation starting from coherence state, since the coefficient
 * is unable to be decided.
 *
 * It should be noted that even though the initial electronic state is same for
 * all nuclear trajectories, different nuclear trajectories will give rise to
 * different reduced density matrix.
 *
 * Reference:
 * [1] X. Gao, M. A. C. Saller, Y. Liu, A. Kelly, J. O. Richardson, E. Geva,
 *     J. Chem. Theory Comput. 2020, 16, 2883.
 * [2] X. Gao, Y. Lai, E. Geva, J. Chem. Theory Comput. 2020, 16, 6465.
 */
class DynamicsMF : public DynamicsMQCBase {
public:
    /**
     * Construct a DynamicsMF object.
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     */
    DynamicsMF(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsMQCBase(param, Ha) {}
    ~DynamicsMF() {}
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
     * Get initial sampling of electronic state.
     *
     * For Ehrenfest method, only let coefficient of initial population state is
     * 1, and others are zero (RDM σ_jj(0)=1, others are 0). (starting from
     * coherence state is not supported). And for each trajectory, the initial
     * coefficient of electronic state is same.
     */
    void samplingElec();
    /**
     * No initial window or denisty in Ehrenfest mean-filed method.
     */
    void getInitialDensity() {}
    /**
     * Update the electronic reduced density matrix at current time.
     * For Ehrenfest method, RDM is directly computed by σ_jk=conj(coeff_k)*coeff_j
     *
     * Both the current (RDM_current) and average (RDM_average) ones are updated.
     */
    void updateDensityMatrix();
};