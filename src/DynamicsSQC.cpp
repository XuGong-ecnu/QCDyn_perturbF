/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsSQC.h"

void DynamicsSQC::init() {
    // 0. Initialize date members of DynamicsMQCBase.
    DynamicsMQCBase::init();
    if (dyn_type != "SQC")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for SQC mapping dynamics.");
    // 1. Read and check parameters for SQC mapping dynamics.
    window_type = param.getStr("SQC_window");
    if (window_type.empty()) // If not specified by user, use triangle (default)
        window_type = "triangle";
    if (window_type != "square" && window_type != "triangle")
        throw std::runtime_error("ERROR: Unsupported window=" + window_type + " for SQC method.");
    // In Miller's SQC method, it ties the value of gamma to the width of the
    // window which is treated as an adjustable parameter whose value is set
    // based on fitting to several popular benchmark models. [S. J. Cotton and
    // W. H. Miller, J. Chem. Phys., 2013, 139, 234112].
    // For SQC_square method, the optimal value ((sqrt(3)-1)/2=0.366) will be used
    // regardless of input value. And for SQC_triangle method, the optimal value
    // (1.0/3.0) will be used by default if the input value is empty.
    if (param.getStr("gamma_MM").empty()) // use optimal value
        gamma_MM = window_type == "square" ? 0.366 : (1.0/3.0);
    else // spcified by user
        gamma_MM = param.getDouble("gamma_MM");
    if (gamma_MM < 0.0 || gamma_MM > 0.5)
        throw std::runtime_error("ERROR: Illegal value for gamma_MM (requires >=0 and <=1/2)");
    // factor used to renormalize RDM.
    // to-do normalize the RDM in observableRDM class
    //std::vector<Complex_Matrix> RDM1;
    //RDM1.resize(num_init_states, Complex_Matrix(nsteps, std::vector<Complex>(DOFe*DOFe, 0)));
    // norm_pop.resize(nsteps, 0);
    // norm_coh = norm_pop;
    // For SQC method, only one initial state can be computed.
    if (init_state.size() != 1)
        throw std::runtime_error("ERROR: Multi-initial states is not supported for SQC method.");
}
/*
const std::vector<Complex_Matrix>& DynamicsSQC::getAverageDensityMatrix() {
    // 1. For SQC method, all values (inlcuding coherence) of RDMs are renormalized
    // by: value / raw_total_population, this will make sure total population is 1.
    // Ref: Eq. 10 in Cotton,and Miller, J. Chem. Phys. 2016, 145, 144108.
    // This is the popular/original way used in most papers.
    // ! However, in this way, re-normalization from coherence is not well defined.
    // 2. Another way is to average the RDM by the number of trajectoies lie
    // within population/coherence window. (in this case, the normalized factor
    // is different for population and coherence, and can be used starting from
    // coherence. For population-population, this two ways are same, but for
    // population-coherence is minor different according to the test.)
    // Ref: J. Chem. Phys. 148, 181102 (2018) and the code in github:
    // github.com/jprov410/mqds, mqds-master/mqds/src/general_src/windows.f90
    // ! However, also, re-normalization from coherence is not well defined.
    // ! So, when starting from coherence, the re-normalization is unknown.
    const int j = init_state[0] / DOFe; // initial state = j*DOFe+k
    const int k = init_state[0] % DOFe;
    // This is raw averaged RDM, which should be re-normalized.
    Complex_Matrix& RDM = RDM_average[0]; // This is raw averaged RDM.
    // When satrting from population, 1st renormalized way is used
    if (j == k) {
        for (int n = 0; n < RDM.size(); ++n) { // n is step
            double totalpop = 0; // get total population firstly
            for (int i = 0; i < DOFe; ++i)
                totalpop += RDM[n][i*DOFe+i].real();
            for (int a = 0; a < DOFe; a++)
                for (int b = 0; b < DOFe; b++)
                    if (totalpop != 0) // population-population/coherence
                        RDM[n][a*DOFe+b] /= totalpop;
        }
    }
    // When satrting from coherence, 2st renormalized way is used
    // ? The normalization for coherence-population/coherence ?
    // Currently, the raw RDM averaged by ntraj will be outputed.
    else {
        std::cout << "WARNING: For SQC method, the re-normalization of RDM from "
            "coherence state is not well defined, and the raw RDM which is averaged "
            "by number of trajectory will be outputted." << std::endl;
        // ! Note the normalization may not right when starting from coherence
        // ! start from coherence, make sure σ_jk(0) = 1. (only for 2 state)
        // ! for multi-state, here, /2.0 is not right.
        //for (int n = 0; n < RDM.size(); ++n) // n is step
        //    for (int a = 0; a < DOFe; a++)
        //        for (int b = 0; b < DOFe; b++)
        //            if (a == b && norm_pop[n] != 0) // coherence-population
        //                RDM[n][a*DOFe+b] *= (double)(ntraj) / norm_pop[n];
        //            else                            // coherence-coherence
        //                RDM[n][a*DOFe+b] *= (double)(ntraj) / (norm_coh[n] / 2.0);
    }
    // For starting from population state, it is the renormalized RDM, while
    // for starting from coherence state, it is the raw RDM.
    return RDM_average;
}
*/

void DynamicsSQC::samplingElec() {
    // electronic mapping variable: positions, momenta, coefficient of wavefunction
    std::vector<double>& q = Elec->q;
    std::vector<double>& p = Elec->p;
    std::vector<Complex>& coeff = Elec->coeff;
    // Electronic mapping variables n, u are used for SQC method.
    // But we can convert them to p, q to do electronic propagation.
    std::vector<double> e(DOFe, 0); // positive-definite "energy" (action n= e - gamma)
    std::vector<double> u(DOFe, 0); // angle variable
    const int j = init_state[0] / DOFe; // initial state = j*DOFe+k
    const int k = init_state[0] % DOFe;
    // Generate random numbers according to the uniformly distributed on the interval [a, b]
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    // * 1. Generate uniform distribution (square and triangle are different)
    // The angle distribution is same u ~ uniform[0, 2*pi]
    // The energy variable is differnet in square and triangle
    if (window_type == "square") {
        // Ref: J. Chem. Phys. 150, 194110 (2019) Eq. 8-10
        // gamma_MM = 0.366, best for square scheme
        // W_k(e) = w_1(ek) * Prod_{k'!=k} w_0(k') given initial state k
        // w_N(e) = 1 for 0 < e-N < 2*gamma_MM, otherwise =0;
        // N=1 for occupied state, N=0 for unoccupied state
        for (int i = 0; i < DOFe; i++) {
            // angle u ~ uniform[0, 2*pi]
            u[i] = uniform_dist(elec_gen) * 2.0 * pi;
            // energy variable e ~ uniform[0, 2*gamma]
            e[i] = uniform_dist(elec_gen) * 2.0 * gamma_MM;
        }
    }
    else { // window_type == "triangle"
        // Ref: J. Chem. Phys. 150, 194110 (2019) Algorithm 1: Sampling
        // triangle implicitly using gamma = 1/3 here.
        // e[i] = 1/2 * (p[i]^2 + q[i]^2)
        // action n[i] = e[i] - gamma[i]
        // W_k(e) = w_1(ek) * Prod_{k'!=k} w_0(ek,ek') given initial state k
        // w_1(e) = (2-e)^{2-DOFe} for 1<e<2, otherwise =0
        // w_0(e,e') = 1 for e' < 2-e, otherwise =0
        // This initial sampling source code in C++ is provided in paper (Appendix A)
        do { // Weighted triangle sampling for initial occupied state
            e[j] = uniform_dist(elec_gen);
        } while (1 - e[j] < uniform_dist(elec_gen));
        // Unoccupied states
        for (int i = 0; i < DOFe; i++) {
            // angle u ~ uniform[0, 2*pi], same as square
            u[i] = uniform_dist(elec_gen) * 2.0 * pi;
            if (i != j)
                e[i] = uniform_dist(elec_gen) * (1.0 - e[j]);
        }
    }
    // * 2. Shift the occupied state according initial state
    if (j == k)  // Let initial population: W_j(e) = 1
        e[j] += 1;
    else { // Let initial coherence: W_j k(e) = e^{i u_j - i u_k} * 1
        e[j] += 0.5;
        e[k] += 0.5;
    }
    // * 3. Convert action("energy")-angle sampling to MM variable
    // q[i] = + sqrt(2e[i]) * cos(u[i]) = + sqrt(2n[i]+2gamma) * cos(u[i])
    // p[i] = - sqrt(2e[i]) * sin(u[i]) = - sqrt(2n[i]+2gamma) * sin(u[i])
    for (int i = 0; i < DOFe; i++) {
        q[i] =  sqrt(2 * e[i]) * cos(u[i]);
        p[i] = -sqrt(2 * e[i]) * sin(u[i]);
        coeff[i] = (q[i] + I * p[i]) / sqrt(2);
    }
}

void DynamicsSQC::getInitialDensity() {
    // Get references to mapping variables (read only)
    const std::vector<double>& q = Elec->q;
    const std::vector<double>& p = Elec->p;
    const std::vector<Complex>& coeff = Elec->coeff;
    const int j = init_state[0] / DOFe; // initial state = j*DOFe+k
    const int k = init_state[0] % DOFe;
    if (j== k) // start from population
        init_density[0] = 1;
    else { // start from coherence, it is the phase factor e^{i u_j - i u_k}
        // coeff_k = (q_k + i*p_k)/sqrt(2) = sqrt(e_k)*exp(-iu_k)
        // abs(coeff_k) = sqrt((q_k^2 + p_k^2)/2) = sqrt(e_k)
        // e^{i u_k - i u_j} = coeff_j/sqrt(e_j) / (coeff_k/sqrt(e_k))
        const double e_j = 0.5 * (q[j] * q[j] + p[j] * p[j]);
        const double e_k = 0.5 * (q[k] * q[k] + p[k] * p[k]);
        init_density[0] = coeff[k] / sqrt(e_k) / (coeff[j] / sqrt(e_j));
    }
}

void DynamicsSQC::updateDensityMatrix() {
    // Get references to mapping variables (read only)
    const std::vector<double>& q = Elec->q;
    const std::vector<double>& p = Elec->p;
    const std::vector<Complex>& coeff = Elec->coeff;
    Complex window = 1.0; // window at current time
    std::vector<double> e(DOFe, 0);
    for (int i = 0; i < DOFe; i++)
        // "energy" variable e[i] = 1/2 * (p[i]^2 + q[i]^2)
        e[i] = 0.5 * (q[i] * q[i] + p[i] * p[i]);
    if (window_type == "square") {
        // J. Chem. Phys. 150, 194110 (2019) Eq. 8-10
        // gamma_MM = 0.366, best for square scheme
        // W_k(e) = w_1(ek) * Prod_{k'!=k} w_0(k') given initial state k
        // w_N(e) = 1 for 0 < e-N < 2*gamma_MM, otherwise =0;
        // N=1 for occupied state, N=0 for unoccupied state
        for (int j = 0; j < DOFe; j++)
            for (int k = 0; k < DOFe; k++) {
                window = 1.0;
                if (j == k) {
                    // population, binning window for WW_k(e) = w_1(e_k) Prod_i w_0(e_i)
                    for (int i = 0; i < DOFe; i++) {
                        // w_1(e_k) = 1 if 0 < e_k-1 < 2*gamma, = 0 otherwise
                        if (i == k && (e[i] <= 1.0 || e[i] >= (2*gamma_MM+1)))
                            window = 0.0;
                        // w_0(e_i) = 1 if 0 < e_i < 2*gamma, = 0 otherwise
                        // here, e_i always > 0
                        else if (i != k && e[i] >= (2*gamma_MM))
                            window = 0.0;
                    }
                    // Record the number of trajectories lie within window at
                    // current time. Note that the value is accumaleated.
                    // if (window.real() == 1.0)
                        // norm_pop[step/RDM_steps] += 1;
                }
                else {
                    // coherence, binning window for WW_k j(e) =
                    // w_1/2(e_k)w_1/2(e_j) Prod_i w_0(e_i) * e^{i u_j - i u_k}
                    for (int i = 0; i < DOFe; i++) {
                        // w_1/2(e_j) = 1 if 0 < e_j-1/2 < 2 * gamma, = 0 otherwise
                        if ((i == j || i == k) && (e[i] <= 0.5 || e[i] >= (2*gamma_MM+0.5)))
                            window = 0.0;
                        // w_0(e_i) = 1 if 0 < e_i < 2 * gamma, = 0 otherwise
                        else if ((i != j && i != k) && e[i] >= (2*gamma_MM))
                            window = 0.0;
                    }
                    // multiply phase factor e^{i u_k - i u_j}
                    // coeff_k = (q_k + i*p_k)/sqrt(2) = sqrt(e_k)*exp(-iu_k)
                    // abs(coeff_k) = sqrt((q_k^2 + p_k^2)/2) = sqrt(e_k)
                    // e^{i u_k - i u_j} = coeff_j/sqrt(e_j) / (coeff_k/sqrt(e_k))
                    if (window.real() == 1.0) {
                        window = coeff[j] / sqrt(e[j]) / (coeff[k] / sqrt(e[k]));
                        // Record the number of trajectories lie within window at
                        // current time. Note that the value is accumaleated.
                        // Note that for coherence, jk, kj is duplicated so
                        // the value will be accumaleated two times in on step.
                        // norm_coh[step/RDM_steps] += 1;
                    }
                }
                // jk element of RDM matrix at time t: σ_jk(t)
                Complex sigma = init_density[0] * window;
                // Here, RDM is stored as DOFe^2-dimentional vector (index j*DOFe+k).
                // For SQC method, only one RDM will be computed: RDM[0].
                // RDM_current[0][step/RDM_steps][j*DOFe+k] = sigma;
                // For RDM_average, the value of each element is accumaleated (+=)
                // and averaged by number of trajectories.
                // RDM_average[0][step/RDM_steps][j*DOFe+k] += sigma / (double)(ntraj);
                RDM[0][j*DOFe+k] = sigma;
            }
    }
    else { // window_type == "triangle"
        // J. Chem. Phys. 150, 194110 (2019) Algorithm 2: binning
        // WW_k(e) = w_1(e_k) Prod_i w_0(e_k, e_i), where
        // w_1(e) = (2-e)^{2-F} if (1<e<2), =0 otherwise
        // w_0(e_k,e_i) = 1 if e_i < 2 - e_k, =0 otherwise
        // Binning window: assigned to state k, if(e_k>=1, and e_i < 1 for all i!=k)
        // Note that for triangle, e[i] is always between 0 and 2.
        for (int j = 0; j < DOFe; j++)
            for (int k = 0; k < DOFe; k++) {
                window = 1.0;
                if (j == k) { // population
                    // Source code in C++ is provided in paper (Appendix A)
                    for (int i = 0; i < DOFe; i++)
                        if ((i == j && e[i] < 1.0) || (i != j && e[i] >= 1.0))
                            window = 0.0;
                    // Record the number of trajectories lie within window at
                    // current time. Note that the value is accumaleated.
                    // if (window.real() == 1.0) by Cesare, it is muted 
                        // norm_pop[step/RDM_steps] += 1;
                }
                else { // coherence (not shown in paper)
                    // according to the range of window for coherence:
                    // w = 1 if 0.5 < e_j & e_k < 1.5, and 0 < e_l < 1, = 0 otherwise
                    for (int i = 0; i < DOFe; i++) {
                        if ((i == j || i == k) && (e[i] <= 0.5 || e[i] >= 1.5))
                            window = 0.0;
                        else if ((i != j && i != k) && e[i] >= 1.0)
                            window = 0.0;
                    }
                    // multiply phase factor e^{i u_k - i u_j}
                    if (window.real() == 1.0) {
                        window = coeff[j] / sqrt(e[j]) / (coeff[k] / sqrt(e[k]));
                        // Record the number of trajectories lie within window at
                        // current time. Note that the value is accumaleated.
                        // Note that for coherence, jk, kj is duplicated so
                        // the value will be accumaleated two times in on step.
                        // norm_coh[step/RDM_steps] += 1;
                    }
                }
                // jk element of RDM matrix at time t: σ_jk(t)
                Complex sigma = init_density[0] * window;
                // Here, RDM is stored as DOFe^2-dimentional vector (index j*DOFe+k).
                // For SQC method, only one RDM will be computed: RDM[0].
                // RDM_current[0][step/RDM_steps][j*DOFe+k] = sigma;
                // For RDM_average, the value of each element is accumaleated (+=)
                // and averaged by number of trajectories.
                // RDM_average[0][step/RDM_steps][j*DOFe+k] += sigma / (double)(ntraj);
                RDM[0][j*DOFe+k] = sigma;
            }
    }
}
