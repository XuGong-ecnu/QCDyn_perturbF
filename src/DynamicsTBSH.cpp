/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsTBSH.h"

void DynamicsTBSH::init() {
    // * 0. Initialize date members of DynamicsMQCBase and check parameters.
    DynamicsMQCBase::init();
    if (dyn_type != "TBSH")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for DynamicsTBSH.");
    if (DOFe != 2)
        throw std::runtime_error("ERROR: For TBSH-MQCL dynamics, only two-state system is supported now.");
    if (representation != "adiabatic")
        throw std::runtime_error("ERROR: Only adiabatic represention is supported for TBSH-MQCL dynamics.");
    if (system_type != "model")
        throw std::runtime_error("ERROR: Realistic system is not supported for dyn_type=" + dyn_type);
    filter_bound = param.getDouble("TBSH_filter"); // By default it is zero.
    if (filter_bound < 0)
        throw std::runtime_error("ERROR: Illegal value for TBSH_filter. Negative filter bound is not allowed.");
    // Get the initial diabatic density matrix element indices |j><k| specified
    // by user, the input format is "j,k", e.g. 0,0, then the initial diabatic
    // density matrix will be initialized as: rho_00 is 1 and others are 0
    // If init_state=all, then all possible initial states are considered.
    // Therefore, the size of rho_init_diabatic is same as init_state.
    rho_init_diabatic.resize(init_state.size(), Real_Matrix(DOFe, std::vector<double>(DOFe, 0)));
    for (int i = 0; i < init_state.size(); ++i) // index of initial state i = j*DOFe+k
        rho_init_diabatic[i][init_state[0]/DOFe][init_state[0]%DOFe] = 1.0; // j=i/DOFe, k=i%DOFe
    // Initialized the vector of density matrix (stored as DOFe^2-dimentional vector).
    // The size of them are the same as the one in base class, but here we use
    // other names to denotes them explicitly.
    // RDM_current_adiabatic = RDM_current_diabatic = RDM_current;
    // RDM_average_adiabatic = RDM_average_diabatic = RDM_average;
    // Here, we set DT = 0.5DT, since the one full step nuclear propagation in
    // TBSH-MQCL method are divided to 2 half-time step.
    DT *= 0.5;

    // 1. Set random number seed for hopping probability
    const int TBSH_seed = param.getInt("TBSH_seed");
    if (TBSH_seed == 0) { // Use system-time seed for real simulation
        std::random_device rd;
        TBSH_gen.seed(rd());
    }
    else // Use deterministic seed for debug
        TBSH_gen.seed(TBSH_seed);
}

void DynamicsTBSH::samplingElec() {
    // Generate uniform sampling of the initial index for pairs of quantum states
    // it is s_0 in orignal TBSH literature [JPCB,2008,112,424], s = a*DOFe+b,
    // a, b is the index of state, DOFe is number of states (N in TBSH literature).
    // Here, 0 <= a,b <= (DOFe-1), 0 <= s < (DOFe^2-1). a,b,s is interger.
    // Produces random integer values, uniformly distributed on the closed interval [a, b]
    std::uniform_int_distribution<int> uniform_dist(0, DOFe*DOFe-1);
    // Initialize s_init (random), s_prime_init, s_previous, s_current, a, b
    // Here, use the random number seed for electronic sampling.
    s_init = uniform_dist(elec_gen);
    s_previous = s_current = s_init;
    a = s_current / DOFe;
    b = s_current % DOFe;
    s_prime_init = b*DOFe + a;
    // Initial value of phase factor is DOFe^2/ntraj, it is N^2/M in eq.35 in
    // orignal TBSH literature [JPCB,2008,112,424].
    phase_factor = (double)(DOFe*DOFe)/ntraj;
}

void DynamicsTBSH::getInitialDensity() {
    // Link to the orignal data in classes of HamiltonianModel and HamiltonianElec
    std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    const Real_Matrix& T = ha->T; // rotation matrix
    const Real_Matrix& T_dagger = ha->T_dagger; // conjuagte transform of T
    for (int i = 0; i < init_density.size(); ++i) { // index of initial state i = j*DOFe+k
        // Compute T_dagger*rho_diabtic*T to get rho(0) in adiabatic base
        Real_Matrix TEMP, rho_init_adiabatic;
        Matrix_Multiply(T_dagger, rho_init_diabatic[i], TEMP);
        Matrix_Multiply(TEMP, T, rho_init_adiabatic);
        // Get initial denisty matrix elemnt ρ(0) of s_prime_init (s'_0)
        // Here, s_prime_init = b*DOFe+a
        // It is the ρ_W,s0'(R,P) in eq.32 in [JPCB,2008,112,424]
        init_density[i] = rho_init_adiabatic[b][a];
    }
}

void DynamicsTBSH::updateDensityMatrix() {
    // Link to the orignal data in classes of HamiltonianModel and HamiltonianElec
    std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    const Real_Matrix& T = ha->T; // rotation matrix
    const Real_Matrix& T_dagger = ha->T_dagger; // conjuagte transform of T
    // * 1. Transform adiabatic ρ(t) to diabatic firstly
    // Create the adabatic unit matrix of s_current
    Real_Matrix rho_adiabatic;
    rho_adiabatic.resize(DOFe, std::vector<double>(DOFe, 0));
    rho_adiabatic[a][b] = 1.0;
    // And compute T*rho_adiabatic*T_dagger to get rho in diabatic base
    Real_Matrix TEMP, rho_diabatic;
    Matrix_Multiply(T, rho_adiabatic, TEMP);
    Matrix_Multiply(TEMP, T_dagger, rho_diabatic);
    // * 2. Compute the observables B(t) by ρ(t)*ρ(0)
    // See the B(t) eq.32,35 in orignal TBSH literature [JPCB,2008,112,424]
    // Note: ρ(0) may have multi-values if all initial states are considered,
    // in this case the all RDMs corresponding to different ρ(0) will be updated.
    for (int i = 0; i < init_density.size(); ++i) {
        // Accumulate diabatic electronic density matrix of current step.
        Complex factor = phase_factor*init_density[i];
        for (int j = 0; j < DOFe; j++)
            for (int k = 0; k < DOFe; k++) {
                // Note the flip of indice here, RDM_jk (σ_jk) = Tr{ρ(0)ρ_kj(t)}
                // In MQCL paper, obervable B(t) correspond to Tr{ρ(0)ρ_jk(t)}
                // Here, we use reduced denisty matix to store the observable
                // just like other mapping dynamics method.
                // RDM_current_diabatic[i][step/RDM_steps][j*DOFe+k] = factor * rho_diabatic[k][j];
                // It should be noted that the factor (/ntraj) has been adsorbed in factor.
                // RDM_average_diabatic[i][step/RDM_steps][j*DOFe+k] += factor * rho_diabatic[k][j];
                // Here, we let the RDM always equal to the diabatic one, which will be outputed.
                // RDM_current[i][step/RDM_steps][j*DOFe+k] = RDM_current_diabatic[i][step/RDM_steps][j*DOFe+k];
                // RDM_average[i][step/RDM_steps][j*DOFe+k] = RDM_average_diabatic[i][step/RDM_steps][j*DOFe+k];
                                
                //calculate RDM of the snapshot
                RDM[i][j*DOFe+k] = factor * rho_diabatic[k][j];
            }
        // We also update the adiabatic electronic density matrix element of current step.
        // RDM_current_adiabatic[i][step/RDM_steps][b*DOFe+a] = factor;
        // RDM_average_adiabatic[i][step/RDM_steps][b*DOFe+a] += factor;
        AdiabaticRDM[i][b*DOFe+a] = factor;
    }
}

void DynamicsTBSH::updateAdiabaticModelForces() {
    std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    // In TBSH-MQCL method, Hellmann-Feynman force is used for nuclear propagation,
    // ref: eq.8 in in orignal TBSH literature [JPCB,2008,112,424].
    // Here, F_all is adiabatic force of all states in DOFe*DOFe base vector,
    // the index a*DOFe+a is the force of state a in adiabatic base.
    // The off-diagnoal of F_all is zero since diabtic coupling is always 0.
    if (a == b) // if population, use the force of this state (same as FSSH method)
        ha->F = ha->F_all[s_current];
    else // if coherence, use the average force of these two states
        for (int j = 0; j < DOFn; ++j)
            ha->F[j] = 0.5 * (ha->F_all[a*DOFe+a][j] + ha->F_all[b*DOFe+b][j]);
}

void DynamicsTBSH::oneStepForAdiabaticModel() {
    // * Propagate a full step as follows: (See eq.35 in TBSH [JPCB,2008,112,424])
    // * (1) Classical adiabatic nuclear propagation with half-time step
    // The first half-time step: e^{iL_s_{j-1}*δ/2}: from (R,P) to (R',P')
    // and compute W_s_{j-1}(t, t+δ/2) and multiply–accumulate phase_factor
    AdiabaticPropagateForModel();

    // * (2) Nonadiabatic propagation: Decide s_j
    // This will evaluate jump probability (Q1 matrix) and momentum jump
    // operator (J_αβ or J_α->β in M matrix). Then update phase factor.
    // If filter is used, it will be applied after getting jump probability.
    NonadiabaticPropagateForModel();

    // * (3) Classical adiabatic nuclear propagation with another half-time step
    // The second half-time step: e^{iL_s_{j-1}*δ/2}: from (R,P) to (R',P')
    // and compute W_s_{j}(t+δ/2, t+δ) and multiply–accumulate phase_factor
    // This is same as first half-time step, but with new (R, P, s).
    AdiabaticPropagateForModel();
}

void DynamicsTBSH::AdiabaticPropagateForModel() {
    std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    // * Classical adiabatic nuclear propagation with half-time step
    // * Half-time step: e^{iL_s_{j}*δ/2}: from (R,P) to (R',P')
    // Here using Velocity Verlet integrator For nulcear propagation
    // Note that here, variable DT has been 0.5DT, i.e., δ/2
    // Note that the effective forces of step 0 should be updated in beforeOneTraj()
    if (s_previous != s_current) // if jump, update the effective forces
        updateF();    // else, the effective forces are same as last half step
    MOVE_Half_V(ha->F, ha->V);
    MOVE_R(ha->V, ha->R);
    // update quantites: H, F_all, NAC, since R changed
    updateH();
    updateF();
    MOVE_Half_V(ha->F, ha->V);
    // * Compute W_s_{j}(t, t+δ/2) and multiply–accumulate phase_factor
    // here, W_s(t,t') = exp(i ∫t,t' dt w_s) ≈ exp (i*0.5*(w_s(t')+w_s(t))*Δt)
    // (definition of W_s See eq.14 in in TBSH [JPCB,2008,112,424])
    // where w_s = w_ab = (E_a - E_b)/hbar, so, when a=b, W = 0, if a != b
    // W = exp(i*0.5/hbar*(ΔE(t')+ΔE(t))*DT) [Note, here, DT is 0.5DT or δ/2]
    if (a != b) {
        double deltaE_old = ha->H_old[a][a] - ha->H_old[b][b];
        double deltaE_now = ha->H[a][a] - ha->H[b][b];
        Complex W = exp(I*0.5/hbar*(deltaE_old+deltaE_now)*DT);
        phase_factor *= W;
    }
}

void DynamicsTBSH::NonadiabaticPropagateForModel() {
    std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    // * Nonadiabatic propagation: Decide s_j by the jump probability and momentum jump operator
    // ! Note: Q1 and M matrix only works for DOFe = 2 now.
    // Evaluate the probability and momentum jump firstly, if no sufficent momentum,
    // set the probability to 0, then, decide s_j and apply momentum shift.
    // This way is used in Alexei Kananenka's Fortran TBSH-TTM code.[JPCL,2016,7,4809]
    // If doesn't set the probability to 0 for the case of no sufficent momentum
    // The total probability will be overestimated, which will lead to inaccurate result.
    // Note, the nonadiabatic coupling (NAC) vectors has been updated.
    // * 1. Save s_current as s_previous before jump
    s_previous = s_current;
    // * 2. Evaluate the probability of jump to each final state (Q1 matrix)
    // i.e., compute the Q1_{s_j-1,s_j}(δ) matrix, δ is a full step time
    // the angle alpha in Q1 matrix is alpha = (P/M · d_βα) * Δt, M is mass,
    // for model, it is 1, P/M=V, d_αβ = -d_βα is nonadiabatic coupling vector.
    // Here, α,β (α<β) is the index of adiabatic state, if DOFe=2, α=0, β=1.
    // (expliclt definition of Q1 See eq.18 in TBSH [JPCB,2008,112,424])
    const std::vector<double>& d_ab = ha->NAC[0]; // d_αβ
    double V_dot_d_ba = 0.0; // P/M · d_βα
    for (int j = 0; j < DOFn; ++j)
        V_dot_d_ba += ha->V[j] * (-d_ab[j]); // d_βα = -d_αβ
    double alpha = V_dot_d_ba * (DT * 2); // multipy 2 since here DT is half step (0.5DT).
    double cosa = cos(alpha);
    double sina = sin(alpha);
    Real_Matrix Q1(DOFe*DOFe, std::vector<double>(DOFe*DOFe, 0.0));
    Q1[0][0] = Q1[1][1] = Q1[2][2]= Q1[3][3] = cosa*cosa;
    Q1[0][1] = Q1[0][2] = Q1[1][3]= Q1[2][3] = -cosa*sina;
    Q1[1][0] = Q1[2][0] = Q1[3][1]= Q1[3][2] = cosa*sina;
    Q1[0][3] = Q1[3][0] = sina*sina;
    Q1[1][2] = Q1[2][1] = -sina*sina;
    // * 3. Evaluate the momentum jump operator J (M matrix)
    // Here, use M matrix to label (interger) which jump operator should be applied.
    // (For details of this section refer to eq.25-27 in TBSH [JPCB,2008,112,424])
    // If lable = 0, no jump operator, correspond to s_j = s_j-1 (no jump),
    // or both of j-1 and j are coherence state (jump between off-diagnoal);
    // if lable = ±1, use J_αβ (+) or J_βα (-) jump operator, correspond to
    // jump between diagnoal and off-diagnoal adiabatic surface;
    // if lable = ±2, use J_α->β (+)  or J_β->α (-) jump operator,
    // correspond to jump between two diagnoal adiabatic surface.
    std::vector<std::vector<int>> M(DOFe*DOFe, std::vector<int>(DOFe*DOFe, 0));
    M[0][1] = M[0][2] = M[1][3]= M[2][3] = 1;
    M[1][0] = M[2][0] = M[3][1]= M[3][2] = -1;
    M[0][3] = 2;
    M[3][0] = -2;
    // * 4. Check if having sufficent momentum along the drection of nonadiabatic
    // * coupling vector to access new state. [if no, set probability to 0]
    // Here, define factor = sgn(d^_αβ·P) * sqrt((d^_αβ·P)^2 + hbar*w_αβ/M)
    // - (d^_αβ·P) for J_αβ operator (when label = ±1), or factor =
    // sgn(d^_αβ·P) * sqrt((d^_αβ·P)^2 + 2*hbar*w_αβ/M) - (d^_αβ·P) for J_α->β
    // operator (when label = ±2), α,β or β,α denpends on the sign of label.
    // where w_αβ = (E_α - E_β)/hbar, M (mass) is 1 for model, d^_αβ is unit
    // vector along d_αβ, i.e., d^_αβ = d_αβ / |d_αβ|.
    // Then, ΔP_αβ or ΔP_α->β = d^_αβ * factor. (See eq.25-27 in TBSH paper).
    // Note, when α < β (label > 0), w_αβ < 0, here, the index of adibatic
    // state is in enenrgy ascending, the item (d^_αβ·P)^2 + hbar*w_αβ/M
    // may less than 0 leading a imaginary factor, which means no sufficent
    // momentum along the nonadiabatic coupling vector can be provided for
    // the jump, in this case ,the jump is impossible, here we will set its
    // probability to 0 in Q1 matrix before deciding s_j with stochastic
    // probability algorithm.
    // Compute unit vector along d_αβ, i.e., d^_αβ = d_αβ / |d_αβ|
    // But, we don't need d^_αβ explicitly. |d_αβ| is enough.
    double ABS_d = 0.0; // |d_αβ| = |d_βα|
    for (int j = 0; j < DOFn; ++j)
        ABS_d += d_ab[j]*d_ab[j];
    ABS_d = sqrt(ABS_d);
    // Compute d^_αβ · P, here, for model, P = V, M = 1
    // Since we have got P/M · d_βα (V_dot_d_ba) before.
    // Then d^_αβ · P = P/M · -d_βα * M / |d_αβ| = -V_dot_d_ba/ABS_d
    double d_ab_hat_dot_P = -V_dot_d_ba/ABS_d; // d^_αβ · P
    // Compute hbar*w_αβ/M = ΔE_ab, w_αβ = (E_α-E_β)/hbar, for model, M is 1
    double dE_ab = ha->H[0][0] - ha->H[1][1];
    for (int s = 0; s < DOFe*DOFe; ++s) {
        int label = M[s_previous][s];
        if (label > 0) { // α < β (label > 0), w_αβ < 0
            // here, square = (d^_αβ · P)^2 + hbar*w_αβ/M for J_αβ or
            // square = (d^_αβ · P)^2 + 2*hbar*w_αβ/M for J_α->β.
            double square = d_ab_hat_dot_P*d_ab_hat_dot_P + label*dE_ab;
            if (square < 0) { // jump is impossible if no sufficent momentum
                Q1[s_previous][s] = 0; // set probability in Q1 matrix to 0
                M[s_previous][s] = 0;  // set label to 0
            }
        }
    }
    // * 5. Deceide s_j with stochastic probability algorithm
    // Copute each probability |Q1_{s_j-1,s}| and total Σ_s|Q1_{s_j-1,s}|
    // Note that, probabilitiy may be set to 0 due to no sufficent momentum.
    double total_prob = 0.0;
    std::vector<double> probs(DOFe*DOFe, 0.0);
    for (int s = 0; s < DOFe*DOFe; ++s) {
        probs[s] = abs(Q1[s_previous][s]);   // |Q1_{s_j-1,s}|
        total_prob += probs[s];              // Σ_s|Q1_{s_j-1,s}|
    }
    // Generate a random probability (uniform distribution) [0, total_prob]
    // And check its location in which filed, then to decide s_j
    // Note: Here, we use a separated seed for hopping probability.
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    double random = uniform_dist(TBSH_gen);
    double random_prob = random * total_prob;
    for (int s = 0; s < DOFe*DOFe; ++s) {
        double sum = 0.0; // sum = Σ_{0,s}|Q1_{s_j-1,s}|
        for (int i = 0; i <= s; ++i)
            sum += probs[i];
        if (random_prob <= sum) { // if random_prob <= Σ_{0,s}|Q1_{s_j-1,s}|
            s_current = s;        // set s_j = s
            break;
        }
    }
    // the prob_factor is the factor should be multiplied to phase_factor.
    // prob_factor = Σ_s|Q1_{s_j-1,s}|/|Q1_{s_j-1,s_j}|
    double prob_factor = total_prob/probs[s_current];
    // * 6. (optional) Apply filter before apply momentum shift (if any)
    // Filter is another essential component of the algorithm. Large fluctuations
    // will come from the factors in the integrand containing the matrix Q1,
    // These large factors exacerbate the sign problem that comes from the
    // phase factors in the evolution and make it difficult to obtain accurate
    // MC estimates of the observable. These fluctuations can be eliminated
    // in part by using a filter similar to that described in ref [G.; Kapral, R.
    // J. Chem. Phys. 2005, 122, 244505.] that eliminates improperly large biasing
    // fluctuations which should not contribute to the averaged quantity. The
    // simplest form of filter is putting an upper bound on the magnitude of
    // the factor in the square bracket in eq.35. When, at stage j in the
    // calculation of the product in the summand, the running summand exceeds
    // the bound, the factor in the updating of the running product is put
    // to unity and the index s_j is put to s_j-1. The value of the bound
    // depends on the nonadiabaticity of the system and the duration of the
    // simulation. The filter must be used for medium and long time simulations.
    // Also, although the filter is not required for short times, a filter may
    // be employed, and the smallest value of the bound which reproduces the
    // previously obtained converged results may be determined. If needed,
    // for longer-time simulations, one may increase the bound further until
    // convergence is obtained. (See discussion in TBSH [JPCB,2008,112,424])
    if (filter_bound > 0)  { // Apply filter, by default it is 0.
        // Here product is the factor in square bracket in eq.35.
        // Since the factor N^2/M out of square bracket in eq.35 has been
        // included phase_factor, elimate it firstly. Here, phase factor is
        // a Complex, so, use the absolute value for comparision.
        Complex product = (double)(ntraj)/(DOFe*DOFe)*prob_factor*Q1[s_previous][s_current]*phase_factor;
        if (abs(product) >=  filter_bound) { // if exceed bound,
            s_current = s_previous;     // set s_j=s_j-1,
            prob_factor = 1.0;          // and set probability factor = 1.
        }
    }
    // * 7. Apply momentum shift ΔP (if any) after s_j has been deceide.
    int label = M[s_previous][s_current];
    if (label != 0) { // jump is accepted and apply momentum shift
        // Here, using the labe_sgn to label d_αβ (+) or d_βα (-) should be used.
        int label_sgn = label > 0 ? 1 : -1;
        // Get sgn(d^_αβ · P) or sgn(d^_βα · P)
        int d_hat_dot_P_sgn = (label_sgn*d_ab_hat_dot_P) > 0 ? 1 : -1;
        // Here, using the label to label ±1*w_αβ or ±2*w_αβ should be used.
        // here, square = (d^_αβ · P)^2 + hbar*w_αβ/M for J_αβ or
        // square = (d^_αβ · P)^2 + 2*hbar*w_αβ/M for J_α->β.
        // where w_αβ = (E_α - E_β)/hbar, M (mass) is 1 for model.
        double square = d_ab_hat_dot_P*d_ab_hat_dot_P + label*dE_ab;
        // Here, factor = sgn(d^_αβ·P) * sqrt((d^_αβ·P)^2 + hbar*w_αβ/M)
        // - (d^_αβ·P) for J_αβ operator (when label = ±1), or factor =
        // sgn(d^_αβ·P) * sqrt((d^_αβ·P)^2 + 2*hbar*w_αβ/M) - (d^_αβ·P)
        // for J_α->β operator (when label = ±2).
        // Using α,β or β,α denpends on the sign of label (label_sgn)
        // Then, ΔP_αβ or ΔP_α->β = d^_αβ * factor = d_αβ/|d_αβ| * factor.
        // Note, for model, P = V, since mass always = 1.
        // Don't forget to use the lable_sgn to decide d_αβ or d_βα
        double factor = d_hat_dot_P_sgn*sqrt(square)-label_sgn*d_ab_hat_dot_P;
        for (int j = 0; j < DOFn; ++j)
            ha->V[j] += label_sgn * d_ab[j]/ABS_d * factor;
    }
    // * 8. Update phase_factor according to s_j-1,s_j by multiply–accumulate
    // inverse jump probability factor: Σ_s|Q1_{s_j-1,s}|/|Q1_{s_j-1,s_j}|
    // and M_{s_j-1,s_j}. Here, M is same as Q1, but with jump operator.
    // If filter is applied, the prob_factor may be 1.
    phase_factor *= prob_factor * Q1[s_previous][s_current];
    // * 9. Update current indices of adiabatic states, a, b, when jump happens.
    // here, s_current = a*DOFe+b
    if (s_previous != s_current) {
        a = s_current / DOFe;
        b = s_current % DOFe;
    }
}