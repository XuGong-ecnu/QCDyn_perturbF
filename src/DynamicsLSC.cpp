/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsLSC.h"

void DynamicsLSC::init() {
    // 0. Initialize date members of DynamicsMQCBase.
    DynamicsMQCBase::init();
    if (dyn_type != "LSC")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for LSC mapping dynamics.");
    // 1. Read and check parameters for LSC mapping dynamics.
    // For original LSC method, there is no publication to discuss the choice of
    // gamma. And the value is 1/2, which means full ZPE effect.
    // For RI-LSC method, the traceless operators do not contain zero-point
    // energy terms arising from the commutators of the harmonic oscillator
    // mapping variables. As a result of this, the difference between
    // RI-LSC1~3 for both operators is simply a factor G(0).
    // Therefore, the potential errors associated with the zero-point energy
    // terms in the traditional LSC methods is avoided. [Ref: Faraday Discuss.,
    // 2020, 221, 150–167].
    // Although the users can specify another value, which is not recommanded.
    if (param.getStr("gamma_MM").empty())
        gamma_MM = 0.5; // defalut value is 0.5
    else // spcified by user
        gamma_MM = param.getDouble("gamma_MM");
    // Get the window width parameter of electronic sampling in LSC method.
    // This is proposed in Xing Gao's paper: J. Phys. Chem. A 2020, 124, 11006−11016
    // By default, it is 1.0, which means the original LSC method.
    LSC_zeta = param.getDouble("LSC_zeta"); // defalut value is 1.0
    if (gamma_MM < 0.0 || gamma_MM > 0.5 || LSC_zeta <= 0.0)
        throw std::runtime_error("ERROR: Illegal value for LSC_zeta (requires > 0), or gamma_MM (requires >=0 and <=1/2)");
    if (gamma_MM != 0.5)
        std::cout << "WARNING: For LSC method, gamma_MM != 1/2 is not well tested." << std::endl;
    // For LSC method, 5 window types will be computed in one simulation.
    // They are LSC1/2, and RI-LSC1~3. So, we don't use the general variables
    // (init_density and RDM) in DynamicsMQCBase, but they have same sizes.
    init_LSC1 = init_LSC2 = init_RILSC1 = init_RILSC2 = init_RILSC3 = init_density;
    RDM_current_LSC1 = RDM_current_LSC2 = RDM_current_RILSC1 = RDM_current_RILSC2 = RDM_current_RILSC3 = RDM_current;
    RDM_average_LSC1 = RDM_average_LSC2 = RDM_average_RILSC1 = RDM_average_RILSC2 = RDM_average_RILSC3 = RDM_average;
}

// const std::vector<Complex_Matrix>& DynamicsLSC::getCurrentDensityMatrix() {
//     // For LSC method, all types of RDM will be inserted to RDM in the order
//     // of LSC1/2, RI-LSC1~3.
//     RDM_current = RDM_current_LSC1;
//     RDM_current.insert(RDM_current.end(), RDM_current_LSC2.begin(), RDM_current_LSC2.end());
//     RDM_current.insert(RDM_current.end(), RDM_current_RILSC1.begin(), RDM_current_RILSC1.end());
//     RDM_current.insert(RDM_current.end(), RDM_current_RILSC2.begin(), RDM_current_RILSC2.end());
//     RDM_current.insert(RDM_current.end(), RDM_current_RILSC3.begin(), RDM_current_RILSC3.end());
//     return RDM_current;
// }

// const std::vector<Complex_Matrix>& DynamicsLSC::getAverageDensityMatrix() {
//     // For LSC method, all types of RDM will be inserted to RDM in the order
//     // of LSC1/2, RI-LSC1~3.
//     RDM_average = RDM_average_LSC1;
//     RDM_average.insert(RDM_average.end(), RDM_average_LSC2.begin(), RDM_average_LSC2.end());
//     RDM_average.insert(RDM_average.end(), RDM_average_RILSC1.begin(), RDM_average_RILSC1.end());
//     RDM_average.insert(RDM_average.end(), RDM_average_RILSC2.begin(), RDM_average_RILSC2.end());
//     RDM_average.insert(RDM_average.end(), RDM_average_RILSC3.begin(), RDM_average_RILSC3.end());
//     return RDM_average;
// }

void DynamicsLSC::samplingElec() {
    // electronic mapping variable: positions, momenta, coefficient of wavefunction
    std::vector<double>& q = Elec->q;
    std::vector<double>& p = Elec->p;
    std::vector<Complex>& coeff = Elec->coeff;
    // * Generate initial electronic sampling according to Gaussian distribution
    // Generates random numbers according to the Normal (Gaussian) random
    // number distribution, mean is 0, standard deviation (sigma) is 1.0
    std::normal_distribution<double> normal_dist(0.0, 1.0);
    // Ref: equation 119 - 120 in Report 35 Mapping
    // Sampling width adjustment parameter LSC_zeta (window width) or 1.
    // LSC_zeta = 1, which means the original LSC method.
    // Sampling using Gaussian distribution in LSC2 Wigner transform
    // density rho = (1/(pi*hbar*LSC_zeta)^F) * exp[-1/hbar * sum_j^DOFe (pj^2+qj^2)/LSC_zeta]
    // standard deviation sigma = sqrt{LSC_zeta*hbar/2} for both q and p
    for (int i = 0; i < DOFe; i++) {
        q[i] = normal_dist(elec_gen) * sqrt(0.5*LSC_zeta*hbar);
        p[i] = normal_dist(elec_gen) * sqrt(0.5*LSC_zeta*hbar);
        // coefficient of electronic wavefunction, coeff[j]=(q[j]+I*p[j])/sqrt(2)
        coeff[i] = (q[i] + I*p[i]) / sqrt(2.0);
    }
}

void DynamicsLSC::getInitialDensity() {
    // Get references to mapping variables (read only)
    const std::vector<double>& q = Elec->q;
    const std::vector<double>& p = Elec->p;
    // Here, initWindow is the window/density at time 0 ([|j><k|]_W(q(0),p(0)))
    // in Eq. 122 in Report 35 Mapping. Note that there is no G(p, q) included,
    // which has been adsorbed in initial sampling. And the sample of them are same
    // Gaussian sampling. The init_density is constant in one trajectory, but is
    // differnent in different trajectory which denpends on the initial sampling.
    // Note that the window operator is expected to be divided by LSC_zeta to make
    // total population normalized to 1, then the factor is 0.5/hbar/LSC_zeta
    // Ref: see equation 127 in Report 35 Mapping
    double sumpop = 0.0; // sumpop will be used in RI-LSC method
    for (int i = 0; i < DOFe; ++i)
        sumpop += (q[i]*q[i] + p[i]*p[i]);
    // G(0) = 2^{DOFe+2} * exp[-1/hbar*sum_j^DOFe(pj^2+qj^2)/LSC_zeta]
    double gpq = pow(2, DOFe+2) * exp(-sumpop / (LSC_zeta*hbar));
    for (int i = 0; i < init_state.size(); ++i) {
        int j = init_state[i] / DOFe; // initial state i = j*DOFe+k, |j><k|
        int k = init_state[i] % DOFe;
        if (j == k) { // population
            // * 1. Original LSC1 (or PBME) method and LSC2 (or LSC-IVR) method
            // Ref: equation 122-129 in Report 35 Mapping
            // Refs: Eitan Geva, et. al., J. Chem. Phys. 151, 074103 (2019).
            // Note that, the LSC_zeta is the parameter to adjust width of the window
            // functions, it should be 1 for original LSC method. And gamma_MM is ZPE
            // paramter, it should be 1/2 for LSC method
            // Both of LSC1 and LSC2 use the LSC2 type Gaussian sampling at time 0,
            // therefore, the type of LSC1/2 initial window, they are same.
            // population: [|k><k|] = 1/(2hbar*LSC_zeta) * (qk^2+pk^2-LSC_zeta*hbar*gamma_MM)
            init_LSC1[i] = 0.5/hbar/LSC_zeta * (q[j]*q[j] + p[j]*p[j] - hbar*gamma_MM*LSC_zeta);
            init_LSC2[i] = init_LSC1[i];
            // * 2. RI-LSC1/2/3 (or mLSC/Φ1Φ1, mLSC/Φ1Φ2, mLSC/Φ2Φ2) method
            // Ref: equation 147 - 149 in Report 35 Mapping
            // Ref: J. Chem. Theory Comput. 2020, 16, 2883.
            // Use Identity and traceless operators for population dynamics, here we
            // call it RI-LSC method.
            // Note that the modifications to the population operator into identity
            // and traceless operators do not change the coherence operator nor the EOM.
            // There are three approaches to choose initial window, but all of
            // them is sampled on G(p(0), q(0)), and have same population operator or
            // coherence operator at time t. The differece is the factor of G(0).
            // The advantage of this method is that the ZPE parameter is vanished.
            // Therefore, in the code of RI-LSC method, no gamma_MM parameter exists.
            // * 2.1 RI-LSC1 (or mLSC/Φ1Φ1) method
            // 1.0/DOFe * Q(0) = 1/(2hbar*LSC_zeta) * [(qm^2+pm^2) - 1.0/DOFe * sum_k^F(qk^2+pk^2)]
            double Q = 0.5/hbar/LSC_zeta * ((q[j]*q[j] + p[j]*p[j]) - sumpop/DOFe);
            // initWindow = 1.0/DOFe * (1 + Q(0))
            // Note that in code, it must be 1.0/DOFe, since 1/DOFe = 0.
            init_RILSC1[i] = 1.0/DOFe + Q;
            // * 2.2 RI-LSC2 (or mLSC/Φ1Φ2) method
            // RI-LSC2 is multiplying a factor G(0) in the part of Q(0) in RI-LSC1.
            // initWindow = 1.0/DOFe * (1 + G(0)*Q(0))
            init_RILSC2[i] = 1.0/DOFe + gpq * Q;
            // * 2.3 RI-LSC3 (or mLSC/Φ2Φ2) method
            // RI-LSC3 is multiplying a factor G(0) in the whole part of RI-LSC1.
            // initWindow = 1.0/DOFe * (1 + Q(0)) * G(0)
            init_RILSC3[i] = gpq * init_RILSC1[i];
        }
        else { // coherence
            // All the initial coherence operators are LSC2 type (without G(0)
            // here, since it has been adsorbed in initial sampling):
            // [|j><k|] = 1/(2hbar*LSC_zeta) * (qj - i*pj) (qk + i*pk)
            init_LSC1[i] = 0.5/hbar/LSC_zeta * (q[j] - I*p[j]) * (q[k] + I*p[k]);
            init_LSC2[i] = init_LSC1[i];
            // But fort RI-LSC method, like population, the G(t)=G(0) also
            // included here, therefore, when compute coherence at t, don't
            // multiply G(t). And the coherence-coherence of RI-LSC1-3 methods
            // are totally same as LSC2 type.
            // ! Note: The RI-LSC methods from coherence are not well defined
            init_RILSC3[i] = init_RILSC2[i] = init_RILSC1[i] = gpq * init_LSC1[i];
        }
    }
}

void DynamicsLSC::updateDensityMatrix() {
    // Get references to needed variables (read only)
    const std::vector<double>& q = Elec->q;
    const std::vector<double>& p = Elec->p;
    // According to the equation 128-129 in Report 35 Mapping. The electronic
    // reduced density matrix elemnet can be divided there parts:
    // 1. prefactor, 2. initWindow/density (time 0), 3. window (time t).
    // Note that In equation 122 the part of nuclear density is adsorbed in
    // initial nuclear sampling, and the electronic initial density rho_e(q,p)
    // is adsorbed in initial electronic sampling. Therefore, we don't include
    // them explicitly here. And the sum_jk{sigma_jk(0)} is the sum of initial
    // density from input, it is 1 for 00 and 0 for others in most case, which
    // means it is 1 (the value used here).
    // 1. The prefactor is left factor of density after the initial density is
    // adsorbed into the initial sampling. For different window_type, they have
    // same prefactor: 2^2 , since they use same Gaussian sample G(p(0), q(0)).
    // 2. The initWindow is the window at time 0 ([|j><k|]_W(q(0),p(0))) without
    // initial density, which can be got by getInitWindow(). The initWindow is
    // constant in one trajectory, but is differnent in different trajectory
    // which denpends on the initial sampling.
    // 3. The window at time t (such as [|n><m|]_W(q(t),p(t))) will be computed
    // at here.
    // Then the density matrix elemnet is pre * initWindow * window, note that
    // the window operator is expected to be divided by LSC_zeta to get total
    // population normalized to 1, then pre-factor is 0.5/hbar/LSC_zeta
    // (Ref: see equation 127 in Report 35 Mapping). In particular,
    // for the population of RI-LSC method, we should add extra 1.0/DOFe.
    // For original (RI)-LSC method, gamma_MM=1/2, LSC_zeta=1.
    Complex window_LSC1, window_LSC2, window_RILSC1, window_RILSC2, window_RILSC3;
    double sumpop = 0.0; // used in RI-LSC method
    for (int i = 0; i < DOFe; ++i)
        sumpop += q[i]*q[i] + p[i]*p[i];
    // G(p, q) = 2^{DOFe+2} * exp[-1/hbar*sum_j^DOFe(pj^2+qj^2)/LSC_zeta]
    double gqp = pow(2, DOFe+2) * exp(-sumpop / (LSC_zeta*hbar));
    for (int j = 0; j < DOFe; ++j)
        for (int k = 0; k < DOFe; ++k) {
            // * 1. Compute [|j><k|]_W(q(t),p(t)) firstly
            if (j == k) { // population
                // * Original LSC1 (or PBME) method
                // Ref: equation 122-129 in Report 35 Mapping
                // population: [|k><k|] = 1/(2hbar) * (qk^2+pk^2-2*LSC_zeta*hbar*gamma_MM)
                // Note that LSC1 initWindow (population) at time 0 is LSC2 type without G(p,q),
                // which is different with LSC1 type here.
                window_LSC1 = 0.5/hbar/LSC_zeta * (q[j]*q[j] + p[j]*p[j] - 2*hbar*gamma_MM*LSC_zeta);
                // * Original LSC2 (or LSC-IVR) method
                // population: [|k><k|] = 1/(2hbar) * (qk^2+pk^2-LSC_zeta*hbar*gamma_MM) * G(p, q)
                // Note that the differece between LSC1 and LSC2, is not only the factor of G(p, q)
                window_LSC2 = 0.5/hbar/LSC_zeta * (q[j]*q[j] + p[j]*p[j] - hbar*gamma_MM*LSC_zeta) * gqp;
                // * RI-LSC1/2/3 (or mLSC/Φ1Φ1, mLSC/Φ1Φ2, mLSC/Φ2Φ2) method
                // Ref: equation 147 - 149 in Report 35 Mapping
                // For RI-LSC method, the initWindow of RI-LSC1/2/3 type is different,
                // but population Qk(t) or coherence [|k><j|]_M(qt,pt) is same, all of them
                // are without G(t), since it has been adsorbed into initWindow (G(t) = G(0)).
                // window at time t: 1.0/DOFe * Qk(t) = 1/(2*hbar) * ((qk^2+pk^2) - 1.0/DOFe * sum_{i}(qi^2+pi^2))
                window_RILSC1 = 0.5/hbar/LSC_zeta * ((q[j]*q[j] + p[j]*p[j]) - sumpop/DOFe);
                window_RILSC3 = window_RILSC2 = window_RILSC1;
            }
            else { // coherence (Note that rho_jk is ... * [|k><j|])
                // * Original LSC1 (or PBME) method
                // coherence: [|k><j|] = 1/(2hbar) * (qk - i*pk) (qj + i*pj)
                window_LSC1 = 0.5/hbar/LSC_zeta * (q[k] - I*p[k]) * (q[j] + I*p[j]);
                // * Original LSC2 (or LSC-IVR) method
                // coherence: [|k><j|] = 1/(2hbar) * (qk - i*pk) (qj + i*pj) * G(p, q)
                window_LSC2 = window_LSC1 * gqp;
                // * RI-LSC1/2/3 (or mLSC/Φ1Φ1, mLSC/Φ1Φ2, mLSC/Φ2Φ2) method
                // Here, the window is same with LSC2 type: [|k><j|]_M(qt,pt)*G(t)
                // But the G(t)=G(0) has been adsorbed into initWindow.
                // So the coherence window of them are using LSC1 type formula.
                window_RILSC3 = window_RILSC2 = window_RILSC1 = window_LSC1;
            }
            // * 2. Compute density matrix element for all possible initial states
            for (int i = 0; i < init_state.size(); ++i) {
                // For LSC1/2, RDM value is prefactor * initWindow * window
                // For all possible initial states, the window at t is same.
                // Therefore they can be computed in one simulation.
                // jk element of RDM matrix at time t: σ_jk(t)
                Complex sigma_LSC1   = 4.0 * init_LSC1[i] * window_LSC1;
                Complex sigma_LSC2   = 4.0 * init_LSC2[i] * window_LSC2;
                Complex sigma_RILSC1 = 4.0 * init_RILSC1[i] * window_RILSC1;
                Complex sigma_RILSC2 = 4.0 * init_RILSC2[i] * window_RILSC2;
                Complex sigma_RILSC3 = 4.0 * init_RILSC3[i] * window_RILSC3;
                // In particular, for the population of RI-LSC method, add extra 1.0/DOFe.
                // while, the coherence of RI-LSC method is same as LSC1/2.
                if (j == k) {
                    sigma_RILSC1 += 1.0 / DOFe;
                    sigma_RILSC2 += 1.0 / DOFe;
                    sigma_RILSC3 += 1.0 / DOFe;
                }
                // RDM store to the results RDM
                RDM_LSC1[i][j*DOFe+k]   = sigma_LSC1;
                RDM_LSC2[i][j*DOFe+k]   = sigma_LSC2;
                RDM_RILSC1[i][j*DOFe+k] = sigma_RILSC1;
                RDM_RILSC2[i][j*DOFe+k] = sigma_RILSC2;
                RDM_RILSC3[i][j*DOFe+k] = sigma_RILSC3;
            }
        }
}