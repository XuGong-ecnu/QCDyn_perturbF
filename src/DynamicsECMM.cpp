/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsECMM.h"

void DynamicsECMM::init() {
    // 0. Initialize date members of DynamicsMQCBase.
    DynamicsMQCBase::init();
    if (dyn_type != "eCMM" && dyn_type != "CMM2")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for DynamicsECMM.");
    // 1. Read and check parameters for eCMM mapping dynamics.
    // If run CMM2 method, gamma_MM = 0; for eCMM, default is 1/2, if no input.
    // It can be specified by user, gamma_MM ~ (-1/DOFe, ∞).
    if (dyn_type == "CMM2")
        gamma_MM = 0; // input value ignored for CMM2 method
    else if (param.getStr("gamma_MM").empty())
        gamma_MM = 0.5; // default value for eCMM is 0.5 (standard MMST)
    else // spcified by user
        gamma_MM = param.getDouble("gamma_MM");
    if (gamma_MM <= (-1.0/DOFe))
        throw std::runtime_error("ERROR: Illegal value for gmmma_MM (requires > -1/DOFe) for eCMM method");
}

void DynamicsECMM::samplingElec() {
    // * Same as full-sphere sampling in DynamicsSPM
    // Note the gamma_MM used here is 1/2 of that used in DynamicsSPM
    std::vector<double>& q = Elec->q;
    std::vector<double>& p = Elec->p;
    std::vector<Complex>& coeff = Elec->coeff;
    std::normal_distribution<double> normal_dist(0.0, 1.0);
    double sum = 0.0;
    const double rad2 = 2.0 + DOFe*gamma_MM*2.0; // the squared {q,p}-radius
    for (int i = 0; i < DOFe; i++) {
        q[i] = normal_dist(elec_gen);
        p[i] = normal_dist(elec_gen);
        sum += q[i]*q[i] + p[i]*p[i];
    }
    for (int i = 0; i < DOFe; i++) {
        q[i] *= sqrt(rad2/sum);
        p[i] *= sqrt(rad2/sum);
    }
    for (int i = 0; i < DOFe; i++)
        coeff[i] = (q[i] + I*p[i])/sqrt(2.0);
    // TODO: foucsed initial sampling
}

void DynamicsECMM::getInitialDensity() {
    const std::vector<double>& q = Elec->q;
    const std::vector<double>& p = Elec->p;
    // For eCMM method, the init_density is the A(0) operator in correlation function.
    // Same as getInitialDensity() in DynamicsSPM, but gamma_MM = 1/2*gamma_init
    for (int i = 0; i < init_state.size(); ++i) {
        int j = init_state[i] / DOFe; // initial state i = j*DOFe+k, |j><k|
        int k = init_state[i] % DOFe;
        if (j == k) // population operator
            init_density[i] = 0.5 * (q[j]*q[j] + p[k]*p[k]) - gamma_MM;
        else        // coherence operator (jk)
            init_density[i] = 0.5 * (q[j] - I*p[j]) * (q[k] + I*p[k]);
    }
}

void DynamicsECMM::updateDensityMatrix() {
    const std::vector<double>& q = Elec->q;
    const std::vector<double>& p = Elec->p;
    // The formula can be found in Jian's JPCL paper (eCMM):
    // Ref2: X. He, J. Liu, J Phys. Chem. Lett. 2021, 12, 2496. Eq.20
    // However, they can be used for population-population operator only.
    // The coherence operator can be found in DynamcisSPM.
    // Actually, the formula used in SPM and eCMM are equivalent.
    // γ=0 -> SPM-Q; γ=1 -> SPM-P; γ=1/2 -> SPM-MMST.
    // The population/coherence at t has a scale factor, and for population,
    // there is a shift value. they depend on the DOFe and gamma_MM.
    // And the normalized factor (DOFe) arised from full-sphere sampling will be
    // multiplied in getDensityMatrix() to get final RDMs.
    const double factor = (1.0+DOFe)/(1.0+DOFe*gamma_MM)/(1.0+DOFe*gamma_MM);
    const double shift = (1.0-gamma_MM)/(1.0+DOFe*gamma_MM);
    for (int i = 0; i < init_density.size(); ++i)
        for (int j = 0; j < DOFe; ++j)
            for (int k = 0; k < DOFe; ++k) {
                Complex B = 0 ; // B(q_k, q_j)(t)
                if (j == k)   // population operator
                    B = 0.5 * factor * (q[j]*q[j] + p[k]*p[k]) - shift;
                else // coherence operator (kj) (same as DynamicsSPM)
                    B = 0.5 * factor * (q[k] - I*p[k]) * (q[j] + I*p[j]);
                // jk element of RDM matrix at time t: σ_jk(t)
                // DOFe is normalized factor arised from full-sphere sampling.
                Complex sigma = (double)(DOFe) * init_density[i] * B;
                // Here, RDM is stored as DOFe^2-dimentional vector (index j*DOFe+k).
                // RDM_current[i][step/RDM_steps][j*DOFe+k] = sigma;
                // For RDM_average, the value of each element is accumaleated (+=)
                // and averaged by number of trajectories.
                // RDM_average[i][step/RDM_steps][j*DOFe+k] += sigma / (double)(ntraj);
                
                //calculate RDM of the snapshot
                RDM[i][j*DOFe+k] = sigma;
            }
}