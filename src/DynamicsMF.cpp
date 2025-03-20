/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsMF.h"

void DynamicsMF::init() {
    // 0. Initialize date members of DynamicsBase.
    DynamicsMQCBase::init();
    if (dyn_type != "MF")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for DynamicsMF.");
    // For Ehrenfest method, only one initial state can be computed.
    if (init_state.size() != 1)
        throw std::runtime_error("ERROR: Multi-initial states is not supported for Ehrefest mean-field method.");
    if (init_state[0]/DOFe != init_state[0]%DOFe) // init_state = j*DOFe+k
        throw std::runtime_error("ERROR: For Ehrenfest mean-filed method, starting from coherence is not supported.");
}

void DynamicsMF::samplingElec() {
    // electronic mapping variable: positions, momenta, coefficient of wavefunction
    std::vector<double>& q = Elec->q;
    std::vector<double>& p = Elec->p;
    std::vector<Complex>& coeff = Elec->coeff;
    // Reset all electronic coefficient to zero since we set coeff directely.
    std::fill(coeff.begin(), coeff.end(), 0.0);
    // For Ehrenfest, make sure density matrix elemnet σ_jk = 1 and others is 0.
    // Here, σ_jk = φ*(k)φ(j)= conj(coeff[k]) * coeff[j].
    // For initial population state, coeff = 1.0, otherwise = 0.0
    // Same for all trajectory.
    coeff[init_state[0]/DOFe] = 1.0; // j = init_state[0]/DOFe
    // Get q0, p0 from coeff, coeff[j]=(q[j] + I * p[j])/sqrt(2)
    // They are same for each trajectory in diabatic representation.
    for (int i = 0; i < DOFe; i++) {
        q[i] = coeff[i].real() * sqrt(2);
        p[i] = coeff[i].imag() * sqrt(2);
    }
}

void DynamicsMF::updateDensityMatrix() {
    const std::vector<Complex>& coeff = Elec->coeff;
    // For Ehrenfest method, RDM is directly computed by σ_jk=φ*(k)φ(j)
    int i = 0;
    for (int j = 0; j < DOFe; j++)
        for (int k = 0; k < DOFe; k++) {
            // jk element of RDM matrix at time t: σ_jk(t)
            Complex sigma = std::conj(coeff[k]) * coeff[j];
            // Here, RDM is stored as DOFe^2-dimentional vector (index j*DOFe+k).
            // RDM_current[0][step/RDM_steps][j*DOFe+k] = sigma;
            // For RDM_average, the value of each element is accumaleated (+=)
            // and averaged by number of trajectories.
            // RDM_average[0][step/RDM_steps][j*DOFe+k] += sigma / (double)(ntraj);

            //calculate RDM of the snapshot
            RDM[i][j*DOFe+k] = sigma;
        }
}