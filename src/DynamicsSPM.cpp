/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsSPM.h"

void DynamicsSPM::init() {
    // 0. Initialize date members of DynamicsBase.
    DynamicsMQCBase::init();
    if (dyn_type != "SPM")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for DynamicsSPM.");
    // 1. Read and check parameters used for spin mapping method.
    // Set representation and gamma_init (gamma_s in paper) and gamma_final
    // Ref: J. Chem. Phys. 152, 084110 (2020), Table I, eq 39
    SPM_type = param.getStr("SPM_type");
    if (SPM_type.empty())
        SPM_type = "W"; // default is W-method, the best one in original paper
    if (SPM_type == "Q")
        gamma_init = 0;
    else if (SPM_type == "P")
        gamma_init = 2;
    else if (SPM_type == "W")
        gamma_init = 2.0/DOFe*(sqrt(DOFe+1.0)-1.0);
    else if (SPM_type == "MMST")
        gamma_init = 1;
    else
        throw std::runtime_error("ERROR: Unsupported SPM_type=" + SPM_type + " for spin-mapping method.");
    // gamma_init (gamma_s) is used Hamiltonian and 2 times of standard gamma_MM
    gamma_MM = 0.5 * gamma_init;
    // Specify full-sphere or foused (default) initial electronic sampling
    elec_sample = param.getStr("elec_sample");
    if (elec_sample.empty())
        elec_sample = "full"; // By default, use full-sphere initial conditions
    if (elec_sample == "focused") {
        gamma_final = gamma_init;
        if (init_state.size() != 1)
            throw std::runtime_error("ERROR: Multi-initial states is not supported for spin-mapping with focused initial condition. "
                "Please use full-sphere sampling (elec_sample=full) in this case.");
        if (init_state[0]/DOFe != init_state[0]%DOFe) // init_state = j*DOFe+k
            throw std::runtime_error("ERROR: For spin-mapping with focused initial condition, starting from coherence is not supported."
                "Please use full-sphere sampling (elec_sample=full) in this case.");
    }
    else if (elec_sample == "full") { // full-sphere
        if (SPM_type == "Q")
            gamma_final = 2;
        else if (SPM_type == "P")
            gamma_final = 0;
        else if (SPM_type == "W")
            gamma_final = gamma_init;
        // For full-sphere MMST case, from Jian's eCMM method.
        // Ref: J. Phys. Chem. Lett. 2021, 12, 2496−2501. eq.20
        else if (SPM_type == "MMST")
            gamma_final = 2.0/(2.0 + DOFe);
    }
    else
        throw std::runtime_error("ERROR: Unsupported elec_sample=" + elec_sample + " for spin-mapping method.");
}

void DynamicsSPM::samplingElec() {
    std::vector<double>& q = Elec->q;
    std::vector<double>& p = Elec->p;
    std::vector<Complex>& coeff = Elec->coeff;
    std::normal_distribution<double> normal_dist(0.0, 1.0);
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    // Ref: J. Chem. Phys. 152, 084110 (2020), Table I, eq. 36
    // For full-sphere initial conditions, requires total population Σ_n{P_n}=1,
    // and 0 << P_n << 1 , here, P_n = 0.5*(q_n^2+p_n^2-gamma_init), which means all
    // {q, p} are uniformly sampled on a surface of hypersphere with a fixed radius=
    // sqrt(2 + DOFe*gamma_init). This is easy to implement by drawing q, p from a standard
    // normal distribution and rescaling them with a common factor radius/sqrt(q^2+p^2).
    // This sampling method is same as used in Jian's CMM2, but gamma_init = 0 for CMM2.
    if (elec_sample == "full") { // like LSC Gaussian sampling, independent of init_state
        double sum = 0.0;
        const double rad2 = 2 + DOFe*gamma_init; // the squared {q,p}-radius
        for (int i = 0; i < DOFe; i++) {
            q[i] = normal_dist(elec_gen);
            p[i] = normal_dist(elec_gen);
            sum += q[i]*q[i] + p[i]*p[i];
        }
        for (int i = 0; i < DOFe; i++) {
            q[i] *= sqrt(rad2/sum);
            p[i] *= sqrt(rad2/sum);
        }
    }
    // Ref: J. Chem. Phys. 152, 084110 (2020), eq. 44
    // For focused initial conditions, if starting from a pure population state i
    // then P_i = 1, and P_n = 0, n!=i. One can interpret this as sampling from
    // the "polar circles" as following:
    else { // elec_sample=focused, depends on the init_state
        const int j = init_state[0] / DOFe; // initial state = j*DOFe+k
        const int k = init_state[0] % DOFe;
        if (j == k) { // initial population state
            // Let the population of of initial state j,j is 1, rad2 is 2+gamma_init
            // And let the population of other states are 0, rad2 is gamma_init
            const double rad_j = sqrt(2+gamma_init);
            const double rad_o = sqrt(gamma_init);
            for (int i = 0; i < DOFe; ++i) {
                double phi = uniform_dist(elec_gen) * 2 * pi;
                double rad = i == j ? rad_j : rad_o;
                q[i] = rad * cos(phi);
                p[i] = rad * sin(phi);
            }
        }
        // ! I don't know how to do focused sampling when starting from coherence.
        // ! This is only a code according to my thought, which may not right.
        // ! Please use full-sphere initial sampling if starting from coherence.
        else { // initial coherence state
            throw std::runtime_error("ERROR: For SPM method with focused initial sampling, "
                "starting from coherence is not supported.");
            // Let the total population of initial states j, k is 1.
            q[j] = normal_dist(elec_gen);
            p[j] = normal_dist(elec_gen);
            q[k] = normal_dist(elec_gen);
            p[k] = normal_dist(elec_gen);
            const double sum = q[j]*q[j] + p[j]*p[j] + q[k]*q[k] + p[k]*p[k];
            const double rad2 = 2 + 2*gamma_init;
            q[j] *= sqrt(rad2/sum);
            p[j] *= sqrt(rad2/sum);
            q[k] *= sqrt(rad2/sum);
            p[k] *= sqrt(rad2/sum);
            // And let the population of other states are 0.
            const double rad_o = sqrt(gamma_init);
            for (int i = 0; i < DOFe; ++i)
                if (i != j && i != k) {
                    double phi = uniform_dist(elec_gen) * 2 * pi;
                    q[i] = rad_o * cos(phi);
                    p[i] = rad_o * sin(phi);
                }
        }
    }
    // Get coefficients of electronic wavefunction according to q, p
    // Since Σ_n{P_n}=1, and P_n can be |coeff_n|^2 - 0.5*gamma_init
    // Here, use standard coeff[j]=(q[j]+I*p[j])/sqrt(2) in MMST model.
    // Then, Σ_n|coeff_n|^2 = 1 + 0.5*DOFe*gamma_init, a scale factor of standard
    // coeff, scale factor = sqrt(1 + 0.5*DOFe*gamma_init). If gamma_init = 0 (Q-method),
    // then it is standard coeff, in this case, it is a equivalence of Ehrenfest.
    // But, here, we only use coeff to store {q, p}, which can be used in the
    // electronic propagation. Therefore, use standard coeff here, too.
    for (int i = 0; i < DOFe; i++)
        coeff[i] = (q[i] + I*p[i])/sqrt(2.0);
}

void DynamicsSPM::getInitialDensity() {
    const std::vector<double>& q = Elec->q;
    const std::vector<double>& p = Elec->p;
    // In spin-maping method, the init_denisty is the A operator (population or
    // coherence) in correlation function in its original paper (A_s(q_0,p_0))
    // Ref: J. Chem. Phys. 152, 084110 (2020) Eq.41-42
    for (int i = 0; i < init_state.size(); ++i) {
        int j = init_state[i] / DOFe; // initial state i = j*DOFe+k, |j><k|
        int k = init_state[i] % DOFe;
        if (j == k) // population operator
            init_density[i] = 0.5 * (q[j]*q[j] + p[k]*p[k] - gamma_init);
        else        // coherence operator (jk)
            init_density[i] = 0.5 * (q[j] - I*p[j]) * (q[k] + I*p[k]);
    }
}

void DynamicsSPM::updateDensityMatrix() {
    const std::vector<double>& q = Elec->q;
    const std::vector<double>& p = Elec->p;
    // Ref: J. Chem. Phys. 152, 084110 (2020) eq.40-42
    // Compute time-dependent density matrix according to correlation function C_AB(t)
    // A(q_n,q_m)(0) has been computed using gamma_init and stored in init_denisty.
    // A|n><m| = init_denisty[n*DOFe+m], here, n,m is indices of initial state.
    // If only one initial state is computed, then init_denisty[0] is the A.
    // B operator is B(q_k, q_j)(t) computed here using gamma_final.
    // reduced density matrix: sigma_jk(t) = prefactor * A_nm(0) * B_kj(t)
    // Note the factor in B: factor = rad2_final/rad2_init,
    // here, rad2 is the squared {q,p}-radius rad2 = 2 + DOFe*gamma.
    // For focused method and W-method, gamma_final=gamma_init, then factor = 1.
    const double factor = (2.0 + DOFe*gamma_final) / (2.0 + DOFe*gamma_init);
    // And prefactor arises from full-sphere sampling over all initial states.
    // RDM should be multiply DOFe, since in this sampling the average population
    // of each state at time 0 is 1/DOFe (like a uniform sampling).
    // However, for focused initial sampling, it is 1 (only one initial state).
    const double prefactor = elec_sample == "full" ? DOFe : 1;
    for (int i = 0; i < init_density.size(); ++i)
        for (int j = 0; j < DOFe; ++j)
            for (int k = 0; k < DOFe; ++k) {
                Complex B = 0; // B(q_k, q_j)(t)
                if (j == k) // population
                    B = 0.5 * (factor * (q[j]*q[j] + p[j]*p[j]) - gamma_final);
                else // coherence
                    B = 0.5 * factor * (q[k] - I*p[k]) * (q[j] + I*p[j]);
                // jk element of RDM matrix at time t: σ_jk(t)
                Complex sigma = prefactor * init_density[i] * B;
                // Here, RDM is stored as DOFe^2-dimentional vector (index j*DOFe+k).
                // RDM_current[i][step/RDM_steps][j*DOFe+k] = sigma;
                // For RDM_average, the value of each element is accumaleated (+=)
                // and averaged by number of trajectories.
                // RDM_average[i][step/RDM_steps][j*DOFe+k] += sigma / (double)(ntraj);
                
                //calculate RDM of the snapshot
                RDM[i][j*DOFe+k] = sigma;
            }
}