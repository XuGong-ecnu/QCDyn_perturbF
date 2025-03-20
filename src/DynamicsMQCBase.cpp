/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsMQCBase.h"

void DynamicsMQCBase::init() {
    // * 0. Initialize date members and check parameters.
    DynamicsBase::init();
    if (DOFe < 2)
        throw std::runtime_error("ERROR: For nonadiabatic dynamics, DOFe must be greater than 1.");
    gamma_MM = 0;
    RDM_steps   = param.getInt("RDM_steps");
    if (RDM_steps < 1 || (nsteps % RDM_steps != 0))
        throw std::runtime_error("ERROR: Illegal value for RDM_steps (requires >= 1). "
           "And the value of nsteps must be divided exactly by RDM_steps.");
    // Get the initial reduced density matrix element indices |j><k| specified
    // by user, the input format is "j,k", e.g. 0,0, then the initial density
    // matrix will be initialized as: rho_00 is 1 and others are 0.
    // If init_state=all, then all possible initial states are considered, and
    // all reduced density matrices starting from all initial states wiil be
    // computed and save (works for some methods).
    const std::string init_state_str = param.getStr("init_state");
    if (init_state_str == "all")  { // all possible initial states, DOFe*DOFe
        init_state.resize(DOFe*DOFe, 0);
        for (int i = 0; i < DOFe*DOFe; ++i) // index of initial state i=j*DOFe+k
            init_state[i] = i;
    }
    else { // specfied one initial state with two indices
        std::vector<int> indices;
        SplitString(indices, param.getStr("init_state"));
        if (indices.size() != 2)
            throw std::runtime_error("ERROR: Illegal init_state=" + param.getStr("init_state") +
                ". You must provide 2 indices to represent initial state.");
        for (int i = 0; i < 2; ++i)
            if (indices[i] >= DOFe || indices[i] < 0)
                throw std::runtime_error("ERROR: Illegal init_state index: " + std::to_string(indices[i]) +
                    ", which should be in range [0, DOFe).");
        init_state.resize(1, 0);
        init_state[0] = indices[0]*DOFe + indices[1];
    }

    // * 1. Initialize initial density and reduced density matrix (RDM)
    // If only one initial state is allowed, then the size of init_state,
    // init_denisty, and RDM is one. These three vectors have same size.
    // The size of initial density depends on the initial state (one or DOFe*DOFe)
    init_density.resize(init_state.size(), 0);
    // Initialized the vector of density matrix (stored as DOFe^2-dimentional vector).
    // Plus 1 to the length, since the step 0 and last step will be stored.
    // The size of RDMs depends on the init_denisty (one or DOFe*DOFe).
    const int length = (nsteps/RDM_steps) + 1;
    if (dyn_type != "LSC") 
        RDM.resize(init_density.size(), std::vector<Complex>(DOFe*DOFe, 0));
    else {
        RDM_LSC1.resize(init_density.size(), std::vector<Complex>(DOFe*DOFe, 0));
        RDM_LSC2.resize(init_density.size(), std::vector<Complex>(DOFe*DOFe, 0));
        RDM_RILSC1.resize(init_density.size(), std::vector<Complex>(DOFe*DOFe, 0));
        RDM_RILSC2.resize(init_density.size(), std::vector<Complex>(DOFe*DOFe, 0));
        RDM_RILSC3.resize(init_density.size(), std::vector<Complex>(DOFe*DOFe, 0));
    }
    if (representation == "adiabatic") {
        if (dyn_type == "LSC") {
            AdiabaticRDM_LSC1.resize(init_density.size(), std::vector<Complex>(DOFe*DOFe, 0));
            AdiabaticRDM_LSC2.resize(init_density.size(), std::vector<Complex>(DOFe*DOFe, 0));
            AdiabaticRDM_RILSC1.resize(init_density.size(), std::vector<Complex>(DOFe*DOFe, 0));
            AdiabaticRDM_RILSC2.resize(init_density.size(), std::vector<Complex>(DOFe*DOFe, 0));
            AdiabaticRDM_RILSC3.resize(init_density.size(), std::vector<Complex>(DOFe*DOFe, 0));
        }
        else if (dyn_type == "TBSH") {
            AdiabaticRDM.resize(init_density.size(), std::vector<Complex>(DOFe*DOFe, 0));
        }
    }

    // * 2. Create an object of DynamicsElec and initialize.
    DyElec = std::make_shared<DynamicsElec>(param);
    DyElec->init();

    // * 3. Get smarter pointer to HamiltonianELec.
    Elec = Ha->getElec();

    // * 4. Bind function to update Hamiltonian in diabatic or adiabatic representation.
    // It is bound to the function in Hamiltonian object.
    if (representation == "diabatic")
        updateH = std::bind(&HamiltonianBase::updateDiabaticHamiltonian, Ha);
    else if (representation == "adiabatic")
        updateH = std::bind(&HamiltonianBase::updateAdiabaticHamiltonian, Ha);
    else if (representation == "quasi-diabatic")
        updateH = std::bind(&HamiltonianBase::updateQuasiDiabaticHamiltonian, Ha);
    else
        throw std::runtime_error("ERROR: Unsupported representation=" + representation +
            ". Only diabatic or adiabatic representation is supported now.");

    // * 5. Bind function to update forces and propagation.
    // 1. Model Hamiltonian within different representations.
    if (system_type == "model") {
        if (representation == "diabatic") {
            updateF = std::bind(&DynamicsMQCBase::updateDiabaticModelForces, this);
            oneStep = std::bind(&DynamicsMQCBase::oneStepForDiabaticModel, this);
        }
        else if (representation == "adiabatic") {
            updateF = std::bind(&DynamicsMQCBase::updateAdiabaticModelForces, this);
            oneStep = std::bind(&DynamicsMQCBase::oneStepForAdiabaticModel, this);
        }
        else if (representation == "quasi-diabatic") {
            updateF = std::bind(&DynamicsMQCBase::updateQuasiDiabaticModelForces, this);
            oneStep = std::bind(&DynamicsMQCBase::oneStepForQuasiDiabaticModel, this);
        }
    }
    // 2. onthefly simulation with electronic structure (adiabatic/quasi-diabatic)
    else if (system_type == "onthefly") {
        throw std::runtime_error("ERROR: onthefly dynamics is not implemented.");
    }
    // 3. All-Atom simulation with OpenMM (diabatic only)
    else if (system_type == "allatom") {
        if (representation != "diabatic")
            throw std::runtime_error("ERROR: Dynamics in adiabatic representation is not supported for all-atom simulation.");
        if (allatom_type == "OpenMM") {
           updateF = std::bind(&DynamicsMQCBase::updateAllAtomForces, this);
           oneStep = std::bind(&DynamicsMQCBase::oneStepForAllAtom, this); // Plan B
        }
        else if (allatom_type == "CustomForceField") {
            std::cout<<" Test1 : updateF = std::bind(&DynamicsMQCBase::updateCustomAllAtomForces, this); "<<std::endl;
            updateF = std::bind(&DynamicsMQCBase::updateCustomAllAtomForces, this);
            oneStep = std::bind(&DynamicsMQCBase::oneStepForCustomAllAtom, this); // Plan C
        }
        //oneStep = std::bind(&DynamicsMQCBase::customCPUIntegrator, this); // Plan A
    }

    // * 6. Set random number seed for electronic sampling
    // Note: FSSH, MF/MF-RDM methods don't need this.
    const int elec_seed = param.getInt("elec_seed");
    if (elec_seed == 0) { // Use system-time seed for real simulation
        std::random_device rd;
        elec_gen.seed(rd());
    }
    else // Use deterministic seed for debug
        elec_gen.seed(elec_seed);
}

void DynamicsMQCBase::beforeOneTraj() {
    // Note the order of these functions is important.
    step = 0; // reset current step to 0.
    // 1. Initial nuclear DOF (R, P) sampling
    samplingNucl();
    // 2. Initialize the Hamiltonian and Forces at step 0.
    // For adiabatic represention, the nonadiabtic coupling (NAC) is also updated.
    // And the rotation matrix T (if any) used for transfromation between diabatic
    // and adiabtic representaions is also updated for model Hamiltonian.
    updateH();
    // 3. Initial electronic DOF (q, p) sampling
    samplingElec();
    // 4. Compute effective forces of step 0 used for next propagtaion.
    updateF();
    // 5. Compute initial density in the calculation of reduced density matrix (RDM)
    // It is requried by some methods and is unchanged in one trajectory.
    getInitialDensity();
    // 6. Compute the initial electroinc reduced density matrix (RDM) of step 0.
    updateDensityMatrix();
}

void DynamicsMQCBase::dynamics(int steps) {
    for (int i = 0; i < steps; i++) {
        // Propagate one step
        oneStep();
        // Increment current step (class data member).
        step++;
        // Calculate observable at every RDM_steps * DT
        // Here. step is the data member which records the current step.
        // Note that the first element is step 0.
        if (step % RDM_steps == 0)
            updateDensityMatrix();
    }
}

// # The following two functions are used to propagate nuclear DOF
// # with velocity Verlet integrator for model in diabatic basis.
void DynamicsMQCBase::oneStepForDiabaticModel() {
    // Link to the orignal data in classes of HamiltonianModel and HamiltonianElec
    static std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    static std::vector<Complex>& coeff = Elec->coeff; // coefficient of electronic wavefunction
    static std::vector<double>& q = Elec->q; // electronic mapping coordinates
    static std::vector<double>& p = Elec->p; // electronic mapping momenta
    static std::vector<double>& R = ha->R; // nuclear coordinates
    static std::vector<double>& V = ha->V; // nuclear velecoites/momenta (P=V for model)
    static std::vector<double>& F = ha->F; // effective forces
    static const Real_Matrix&     H = ha->Heff; // diabatic effecive Hamiltonian matrix (Havg removed)
    static const Real_Matrix& H_old = ha->Heff_old; // previous diabatic Hamiltonian
    // * The following is Velocity Verlet integrator
    // Update velocities/momenta with a half step using effective force.
    MOVE_Half_V(F, V);
    // Update positions with a full step
    MOVE_R(V, R);
    // Update H and F_all since R has been changed.
    updateH();
    // TODO: Perhaps take MOVE_elec() as parameter or std::function
    // Update electronic q, p, coeff with RK4 method using effective H
    // The effective H at intermediate time are computed by linear interpolation.
    DyElec->MOVE_elec(step*DT, H_old, H, q, p, coeff);
    // Update effective forces
    updateF();
    // Update velocities/momenta with a half step using new effective forces
    MOVE_Half_V(F, V);
}

void DynamicsMQCBase::updateDiabaticModelForces() {
    // Link to the orignal data in classes of HamiltonianModel and HamiltonianElec
    static std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    static const std::vector<double>& q = Elec->q; // electronic mapping coordinates
    static const std::vector<double>& p = Elec->p; // electronic mapping momenta
    static const std::vector<double>& F_avg = ha->F_avg; // average forces
    static const Real_Matrix&         F_all = ha->F_all; // diabatic forces of each states
    static std::vector<double>&       F     = ha->F; // effective forces (update)
    std::fill(F.begin(), F.end(), 0.0); // Reset all elements to 0 before update
    // This is genenral formula to get the effective forces for mapping based method
    // from definition of diabatic Hamiltonian, which will be used to do nulcear propagation.
    // F = 1.0/hbar * Σ_nm{((-dHnn/dR))*pop_nn + (-dH_nm/dR)*coh_nm (n!=m)}
    // here, off-diagonal Hnm = Hmn* (n!=m) is coupling (in most case, it is real),
    // therfore,(-dH_nm/dR)*coh_nm (n!=m) = 2*(-dH_nm/dR)*coh_nm (n<m)
    // And the imaginary part of coh_nm will be vanished in the above formula.
    // If Condon approximation is used. Hnm (n!=m) is constant, then (-dH_nm/dR) is zero.
    // pop_nn and coh_nm are electronic population and cohenrence operators, repectivelly.
    // They have different definitions in different methods from elctronic mapping variables.
    // Note that the pop and coh is directly from the definition of Hamiltonian,
    // which is different from that last observables in density matrix in some methods.
    // If the total pop is always 1 at any time, i.e., Σ_n{pop_nn}=1), then F can be written as
    // F = (-dH_bar/dR) + 1.0/hbar * Σ_nm{((-dHnn/dR)-(-dH_bar/dR))*pop_nn + 2*(-dH_nm/dR)*coh_nm (n<m)}
    // in which, H_bar is the average of diagonal elemnts of H. And H - H_bar is
    // called effective H (H_tilde). This symmetrized form is almost always
    // employed due to the numerical stability if it can be employed.
    // Ref: J. Chem. Phys. 152, 084110 (2020), eq. 35-36 (SQC paper)
    // Here, we use F_all as -dH/dR, F_avg as -dH_bar/dR, which has been computed
    // based on nulear positions and stored in Hamiltonian object.
    // Moreover, the zero point energy parameter (gamma) (if any) is not included
    // in population operator, since it has no impact on the effective force (constant)
    for (int j = 0; j < DOFn; j++) {
        for (int n = 0; n < DOFe; n++)
            for (int m = 0; m < DOFe; m++)
                if (n == m) // from diagnoal of H (population)
                    F[j] += (F_all[n*DOFe+m][j] - F_avg[j]) * 0.5 * (q[n]*q[n] + p[n]*p[n]);
                else if (n < m) // from off-diagnoal of H (coherence)
                    // factor 0.5 in coherence is ommited since nm is same as mn (multiply 2)
                    // If Condon approximation is used (coupling is contant), it is 0.
                    F[j] += F_all[n*DOFe+m][j] * (q[n]*q[m] + p[n]*p[m]);
        F[j] /= hbar;
        F[j] += F_avg[j];
    }
}
// =============================================================================

// # The following three functions are used to propagate nuclear DOF
// # with for model in quasi-diabatic basis proposed by Pengfei Huo.
void DynamicsMQCBase::oneStepForQuasiDiabaticModel() {
    // Link to the orignal data in classes of HamiltonianModel
    static std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    static std::vector<double>& R = ha->R; // nuclear coordinates
    static std::vector<double>& V = ha->V; // nuclear velecoites/momenta (P=V for model)
    static std::vector<double>& F = ha->F; // effective forces
    static const Real_Matrix&         H = ha->H; // current Hamiltonian matrix V_αβ(R(t2))
    static const Real_Matrix&         T = ha->T; // current transformation matrix at R(t2)
    static const Real_Matrix&         F_all = ha->F_all; // -▽V_αβ(R(t2))
    static const std::vector<double>& F_avg = ha->F_avg; // (-1/DOFe)Σ▽V_αα(R(t2))
    // The Quasi-Diabatic (QD) propagation scheme follows the Pengfei Huo's work
    // Ref. J. Chem. Phys. 149, 044115 (2018) Section. II. D
    // Note, in the following, t1 is previous time, and t2 is current time.
    // The algorithm for the QD-SQC quantum dynamics propagation as following:
    // 0. Create the data will be used in Quasi-Diabatic (QD) propagation scheme
    // Here, the H_old and T_old are always the quantities at R(t1).
    static Real_Matrix H_old; // V_αβ(R(t1)) (using aidabatic as QD basis)
    static Real_Matrix T_old; // transformation matrix at R(t1)
    // The electronic mapping varibales in QD basis used in the propagation
    static std::vector<Complex> coeff(DOFe, 0); // {coeff_α} in the instantaneous QD basis
    static std::vector<double> q(DOFe, 0); // {q_α} in the instantaneous QD basis
    static std::vector<double> p(DOFe, 0); // {p_α} in the instantaneous QD basis
    // 1. Do initial nuclear sampling of nuclear DOF R(t1), P(t1+0.5DT)
    // and initial elctronic mapping variables sampling q, p (diabatic).
    // (initial sampling has been done in beforeOneTraj() ,not included in the loop)
    // At step 0, we need transform the initial diabatic ones into adiabatic basis.
    // Note, if we want to get the adiabatic RDM, we don't need to do this.
    if (step == 0)
        transformCoefficient(Elec->coeff, T, coeff, q, p);
    // 2. Perfrom calculations at t1 to obtain the QD basis |Φ_α(R0)> = |Φ_α(R(t1))>
    // (For step 0, it is done in beforeOneTraj() ,not included in the loop)
    // Note, at time t1, always use the adiabatic basis as QD basis i.e.,
    // V_αβ(R(t1)) = E_α(R(t1))δ_αβ, so we assign the ha->H_old (in this case, it
    // is always the adiabatic Hamiltonian) to H_old as QD basis before update.
    H_old = ha->H_old;
    // And assign the transformation matrix to T_old before update Hamiltonian
    T_old = ha->T;
    // 3. Propagate nuclear positions with DT
    MOVE_R(V, R);
    // 4. Perfrom calculations at t2=t1+DT to get the adiabatic basis |Φ_μ(R(t2))>
    // and this will also compute electronic Hamiltonian V_αβ(R(t2)) and
    // forces -▽V_αβ(R(t2)), the corresponding variables are H and F_all.
    // Note, the ha->H_old is also updated to the current adiabatic Hamiltonina.
    // But the variable H_old in this function is at previous time (V_αβ(R(t1)))
    updateH();
    // 5. Propagate the electronic mapping variables by RK4 method
    // The quasi-diabatic electroinc Hamiltonin element at intermediate time
    //  V_αβ(R(t)) are computed by linear interpolation, i.e.,
    // V_αβ(R(t)) = V_αβ(R(t1)) + (t-t1)/(t2-t1)*[V_αβ(R(t2))-V_αβ(R(t1))]
    // here, H_old is V_αβ(R(t1)) and H is V_αβ(R(t2)).
    // Note, MOVE_elec() is totally same as the diabtaic one, but Here, H is the
    // quasi-diabatic electronic Hamilton V_αβ(R(t) and now q,p is in QD basis.
    DyElec->MOVE_elec(step*DT, H_old, H, q, p, coeff);
    // 6. Compute the effective forces (F(R(t2))) used to update P
    // Note, this is same as the diabtaic one, but not use symmterical one.
    // F = 0.5 * Σ_{αβ} -▽V_αβ(R(t2)) [q_α*q_β + p_α*p_β - 2*γ*δ_αβ]
    // Here, -▽V_αβ(R(t2)) is F_all and γ is ZPE parameter in MMST mapping.
    // Ref. J. Chem. Phys. 149, 044115 (2018) Eq. 7 and the equation in page 6.
    // I found LSC and SQC using the original one, suffer sever numerical problem.
    // This is because 0.5*(q_α^2+p_α^2-γ) may be negative if γ > 0, which means
    // that the system evolves on inverted potentials, which can lead to problems,
    // in particular, for systems with steep potentials. While, the Eherenfest
    // (γ = 0) and eCMM, LSC (full-sphere sampling, q,p radius is constant,totalpop
    // is always 1, q^2+p^2=F*γ+1) doesn't suffer this inverted potential problem.
    // Therefore, we use the symmterical one here, it can verified that the form
    // of it is same as the diabatic one, but here, the ▽V_αβ(R(t2)) is quasi-diabatic.
    std::fill(F.begin(), F.end(), 0.0); // reset to 0
    for (int j = 0; j < DOFn; ++j) { // j is index of nulcer DOF
        for (int a = 0; a < DOFe; ++a)
            for (int b = 0; b < DOFe; ++b) {
                double factor = 0.5 * (q[a]*q[b] + p[a]*p[b]);
                F[j] += (F_all[a*DOFe+b][j] - F_avg[j] * (a == b ? 1 : 0)) * factor;
            }
        F[j] += F_avg[j];
    }
    // 7. Propagate nuclear velecoites/momenta with DT using effective forces
    // P(t1+0.5DT) = P(t1+0.5DT) + F(R(t2))*DT/M (for model, M=1, and P is V)
    MOVE_V(F, V);
    // 8. Transform the mapping variables from the instantaneous QD basis
    // {q_α, p_α} back to the strict diabatic basis {q_i, p_i} with
    // q_i = Σ_α <Φ_α(R(0) | i > q_α and p_i = Σ_α <Φ_α(R(0) | i > p_α
    // Ref: J. Chem. Phys. 149, 044115 (2018) Eq. 23
    // Note, at time t1, the adiabatic basis Φ_α(R(t1)) is used as QD basis Φ_α(R(0).
    // So, transforming QD into diabatic basis is to transform aidabtic to
    // diabatic basis using the transformation matrix at time t1 (T_old)
    // Note Elec->coeff,q,p are always diabatic used for observable calcualtion.
    // Note, if we didn't do this transformation, the computed RDM is adiabatic.
    transformCoefficient(coeff, T_old, Elec->coeff, Elec->q, Elec->p, true);
    // 9. Transform the mapping variables {q_α, p_α} from old QD basis into the
    // new QD basis, |Φ'_μ(R'(0))> = |Φ_μ(R(t2))> for next nuclear propagation
    // step with expression: q_μ = Σ_α <Φ_α(R(t1)|Φ_μ(R(t2)> q_α.
    // Here, <Φ_α(R(t1)|Φ_μ(R(t2)> is the overlap matrix of the adiabatic
    // wavefunction at different time. For diabatic model, it can be computed by
    // writting the adiabatic state in terms of diabatic ones via transformation
    // matrix: |α> = Σ_i T_iα |i>, where α/i is the index of the adiabatic/diabatic
    // state, and T is rotation (transformation) matrix. Then, it can be verified
    // that: S_αμ = <Φ_α(R(t1)|Φ_μ(R(t2)> = Σ_i (T_iα(R(t1)) * T_iμ(R(t2)))
    Real_Matrix S(DOFe, std::vector<double>(DOFe, 0)); // overlap matrix: <Φ_α(R(t1)|Φ_μ(R(t2)>
    for (int a = 0; a < DOFe; ++a)
        for (int u = 0; u < DOFe; ++u)
            for (int i = 0; i < DOFe; ++i)
                S[a][u] += T_old[i][a] * T[i][u]; // S_αμ = Σ_i (T_iα(R(t1)) * T_iμ(R(t2)))
    std::vector<Complex> coeff_new(DOFe, 0); // coeff in new QD basis (temp)
    std::vector<double> q_new(DOFe, 0); // q in new QD basis (temp)
    std::vector<double> p_new(DOFe, 0); // p in new QD basis (temp)
    for (int u = 0; u < DOFe; ++u) {
        for (int a = 0; a < DOFe; ++a) {
            q_new[u] += S[a][u] * q[a]; // q_μ = Σ_α S_αμ q_α
            p_new[u] += S[a][u] * p[a]; // p_μ = Σ_α S_αμ p_α
        }
        coeff_new[u] = (q_new[u]+I*p_new[u])/sqrt(2.0); // coeff=(q+I*p)/sqrt(2)
    }
    // update the mapping variables {q_α, p_α} by assigning new q,p,coeff
    q = q_new;
    p = p_new;
    coeff = coeff_new;
}
// =============================================================================

// # The following three functions are used to propagate nuclear DOF
// # with velocity Verlet integrator for model in aidiabatic basis.
void DynamicsMQCBase::oneStepForAdiabaticModel() {
    // Link to the orignal data in classes of HamiltonianModel
    static std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    static std::vector<double>& R = ha->R; // nuclear coordinates
    static std::vector<double>& V = ha->V; // nuclear velecoites/momenta (P=V for model)
    static std::vector<double>& F = ha->F; // effective forces
    static const Real_Matrix&     H = ha->H; // adiabatic Hamiltonian matrix
    static const Real_Matrix& H_old = ha->H_old; // previous adiabatic Hamiltonian
    static const Real_Matrix&     T = ha->T; // current rotation matrix
    static const Real_Matrix& F_all = ha->F_all; // adiabatic forces of each states
    static const std::vector<double>& F_avg = ha->F_avg; // avearge adiabatic forces
    // The effective Hamiltonian is the one removing the average potential energy
    // which can be used in the equivalent symmetric Hamiltonian (like diabatic one)
    static const Real_Matrix&     Heff = ha->Heff; // adiabatic effective Hamiltonian matrix
    static const Real_Matrix& Heff_old = ha->Heff_old; // previous effective adiabatic Hamiltonian
    // NAC stores nonadiabatic coupling vectors between states i and j (i < j)
    // The number of vectors is (DOFe-1)^DOFe/2 (DOFe > 1). And the size of each
    // vector is DOFn. The index of d_ij is NAC matrix is: index =
    // (2*DOFe-i-1)*i/2 + (j-i)-1, 0 =< i < j < DOFe. And, d_ij = -d_ji
    static const Real_Matrix&           NAC = ha->NAC;
    // The adiabatic propagation scheme follows  the Miller's work.
    // Ref. J. Chem. Phys. 147, 064112 (2017)
    // Create the static data will be used in adiabatic propagation scheme
    // Here, the T_old is always the rotation matrix at previous step.
    static Real_Matrix T_old; // transformation matrix at R(t1)
    static Real_Matrix NAC_old; // NAC at previous step
    static std::vector<double> V_old; // V at previous step
    // The electronic mapping varibales in adiabatic basis used in the propagation
    static std::vector<Complex> coeff(DOFe, 0); // {coeff_α} in the adiabatic basis
    static std::vector<double> q(DOFe, 0); // {q_α} in the adiabatic basis
    static std::vector<double> p(DOFe, 0); // {p_α} in the adiabatic basis
    // Do initial nuclear sampling of nuclear DOF R(t1), P(t1+0.5DT)
    // and initial elctronic mapping variables sampling q, p (diabatic).
    // (initial sampling has been done in beforeOneTraj() ,not included in the loop)
    // At step 0, we need transform the initial diabatic ones into adiabatic basis.
    // And they are always adiabatic in this function (static varibales).
    // But after propation, we transform them back to diabatic representation
    // (the mapping variables in HamiltonianElec) to update the diabatic RDM.
    // we don't need to consider the transformation for the initial sampling and
    // observables, just use the orignal one, which is more clearer and convenient
    // than do the transformation in the initial sampling and observables.
    if (step == 0)
        transformCoefficient(Elec->coeff, T, coeff, q, p);
    // At zero step, we should cast the initial canonical to kinematic momentum
    // Compute kinematic momentum P_kin = P + ΔP, (for model, P == V)
    // ΔP = Σ_{i<j}(q_i p_j - q_j p_i) d_ij(R)
    // Ref: Miller, J. Chem. Phys. 147, 064112 (2017), Eqs. 10a and 3b
    // Note: We only need to do this in step 0, since after that, the momentum
    // in the propagation is always the kinematic momentum.
    if (step == 0)
        for (int i = 0; i < DOFe; ++i)
            for (int j = i+1; j < DOFe; ++j) { // i < j
                double factor = q[i]*p[j] - q[j]*p[i];
                int ij = (2*DOFe-i-1)*i/2 + (j-i)-1; // index of d_ij
                for (int k = 0; k < DOFn; ++k)
                    V[k] += factor * NAC[ij][k];
            }
    // Copy current V, NAC, T as old ones, which will be used in electronic propagation
    V_old   = V;
    NAC_old = NAC;
    T_old   = T;

    // * This is leap-frog Verlet integrator used for nulcear propagation.
    // Update velocities/momenta with a full step using effective forces
    MOVE_V(F, V);
    // Update positions with a full step
    MOVE_R(V, R);
    // Update H, F_all, T, and NAC since R has been changed.
    updateH();
    // Update electronic q, p, coeff with RK4 method using H, V, and NAC
    // The H, V, and NAC at intermediate time are computed by linear interpolation.
    DyElec->MOVE_elec(step*DT, H_old, H, V_old, V, NAC_old, NAC, q, p, coeff);
    // Transfrom updated electronic q, p, coeff back to diabatic.
    // The Elec->coeff, q, p (diabatic) will be used for calculation of observables.
    transformCoefficient(coeff, T, Elec->coeff, Elec->q, Elec->p, true);
    // Update effective forces
    updateF();

    // * The following is Velocity Verlet integrator (similar results)
    // Update velocities/momenta with a half step using effective force.
    //MOVE_Half_V(F, V);
    // Update positions with a full step
    //MOVE_R(V, R);
    // Update H, F_all, T, and NAC since R has been changed.
    //updateH();
    // Update effective forces
    //updateF();
    // Update velocities/momenta with another half step using new effective forces
    //MOVE_Half_V(F, V);
    // Update electronic q, p, coeff with RK4 method using H, V, and NAC
    // The H, V, and NAC at intermediate time are computed by linear interpolation.
    //DyElec->MOVE_elec(step*DT, H_old, H, V_old, V, NAC_old, NAC, q, p, coeff);
    // Transfrom updated electronic q, p, coeff back to diabatic.
    // The Elec->coeff, q, p (dibatic) will be used for calculation of observables.
    //transformCoefficient(coeff, T, Elec->coeff, Elec->q, Elec->p, true);
}

void DynamicsMQCBase::updateAdiabaticModelForces() {
    // Link to the orignal data in classes of HamiltonianModel
    static std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    static const Real_Matrix&         H = ha->H; // adiabatic Hamiltonian matrix
    static const Real_Matrix&       NAC = ha->NAC; // nonadiabatic coupling
    static const Real_Matrix&     F_all = ha->F_all; // adiabatic forces of each states
    static const Real_Matrix&         T = ha->T; // rotation matrix
    static const std::vector<double>& F_avg = ha->F_avg; // avearge adiabatic forces
    static std::vector<double>&       F = ha->F; // effective forces (update)
    static std::vector<Complex> coeff; // coefficient of electronic wavefunction
    static std::vector<double> q;     // electronic mapping coordinates
    static std::vector<double> p;     // electronic mapping momenta
    // * Get adiabatic coeff, q, p from Elec->coeff (diabatic) for propagation
    transformCoefficient(Elec->coeff, T, coeff, q, p);
    // Reset all elements to 0 before update
    std::fill(F.begin(), F.end(), 0.0);
    // Ref: Miller, J. Chem. Phys. 147, 064112 (2017), Eq. 12d
    // On the adiabatic representation of Meyer-Miller electronic-nuclear dynamics
    for (int i = 0; i < DOFe; ++i) // i, j is the index of mapping variables
        for (int j = 0; j < DOFe; ++j)
            if (i == j)
                for (int k = 0; k < DOFn; ++k) // k is the index nuclear DOF
                    F[k] += F_all[i*DOFe+i][k] / DOFe; // F_all[i*DOFe+i] is adiabtic forces of state i
            else if (i < j) {
                double factor1 = 0.5/DOFe * (q[i]*q[i] + p[i]*p[i] - q[j]*q[j] - p[j]*p[j]);
                double factor2 = (q[i]*q[j]+p[i]*p[j]) * (H[i][i]-H[j][j]); // H[i][i] is adiabatic E_i
                int ij = (2*DOFe-i-1)*i/2 + (j-i)-1; // index of d_ij (i<j) in NAC
                for (int k = 0; k < DOFn; ++k)
                    F[k] += factor1*(F_all[i*DOFe+i][k]-F_all[j*DOFe+j][k]) + factor2*NAC[ij][k];
            }
}

void DynamicsMQCBase::transformCoefficient(const std::vector<Complex>& coeff_old, const Real_Matrix& T,
    std::vector<Complex>& coeff, std::vector<double>& q, std::vector<double>& p, bool inverse) {
    coeff.resize(DOFe, 0); // coefficient of electronic wavefunction
    q.resize(DOFe, 0);     // electronic mapping coordinates
    p.resize(DOFe, 0);     // electronic mapping momenta
    std::fill(coeff.begin(), coeff.end(), 0.0); // reset all elements to 0
    // * Transfrom diabatic coefficient to adiabatic by using |i> = Σ_j T_ji |j>,
    // where i/j is the index of the adiabatic/diabatic state, T is rotation matrix.
    if (!inverse) // diabatic -> adiabatic (inverse = false) [default]
        for (int i = 0; i < DOFe; ++i)
            for (int j = 0; j < DOFe; ++j)
                coeff[i] += T[j][i] * coeff_old[j];
    // * Transfrom adiabatic coefficient to diabatic by using |i> = Σ_j T_ij |j>,
    // where i/j is the index of the diabatic/adiabatic state, T is rotation matrix.
    else        // adiabatic -> diabatic (inverse = true)
        for (int i = 0; i < DOFe; ++i)
            for (int j = 0; j < DOFe; ++j)
                coeff[i] += T[i][j] * coeff_old[j];
    // * Get q, p from coeff, coeff[j]=(q[j] + I * p[j])/sqrt(2)
    for (int i = 0; i < DOFe; i++) {
        q[i] = coeff[i].real() * sqrt(2);
        p[i] = coeff[i].imag() * sqrt(2);
    }
}
// =============================================================================

// # The following two functions are used to propagate nuclear DOF using OpenMM
// # internal integrator with providing effective forces in diabatic basis. (Plan B)
void DynamicsMQCBase::oneStepForAllAtom() {
    // This function is modified based on the new strategy on Dec. 3, 2021.
    // New strategy: the data of positions and velocities are stored in the
    // OpenMM Platform and the integration are performed by OpenMM Integrator.
    // But the forces used to do nuclear propagtaion is provided externally.
    // At each step, the force of each state should be downloaded from
    // OpenMM Platform such as CUDA, and weighted by the electronic population
    // and coherence to get the effective forces, then the effective forces will
    // be uploaded to the OpenMM Platform to do the nuclear propagtaion.
    // Link to the orignal data in classes of HamiltonianOpenMM and HamiltonianElec
    static std::shared_ptr<HamiltonianOpenMM> ha = std::static_pointer_cast<HamiltonianOpenMM>(Ha);
    static std::vector<Complex>& coeff = Elec->coeff; // coefficient of electronic wavefunction
    static std::vector<double>& q = Elec->q; // electronic mapping coordinates
    static std::vector<double>& p = Elec->p; // electronic mapping momenta
    static std::vector<OpenMM::Vec3>& F = ha->F; // effective forces used for nuclear propagation
    static const Real_Matrix&     H = ha->Heff; // diabatic effecive Hamiltonian matrix (Havg removed)
    static const Real_Matrix& H_old = ha->Heff_old; // previous effecive diabatic Hamiltonian
    static const std::string integrator = param.getStr("integrator"); // OpenMM integrator
    // The nuclear propagation is perfromed by velocity Verlet integrator.
    // Note that within velocity Verlet integrator, the electronic propagtion
    // is performed after the updating of nuclear position before the second
    // half step of updating velcities.
    if (integrator == "velocityVerlet") {
        // The velocity Verlet integrator is a special integrator, in which the
        // forces used in the updating of velocities at the first half step and
        // at the second half step are different. There is no original OpenMM
        // velocity Verlet Integrator. Here, using a CompoundIntegrator with two
        // CustomIntegrator to do velocity Verlet integrator with external forces.
        static std::shared_ptr<OpenMM::CompoundIntegrator> compound =
            std::static_pointer_cast<OpenMM::CompoundIntegrator>(ha->integrator);
        // At the start of each time step, we need to update OpenMM context state.
        // The thermostat will modify velocities; the barostat will modify positions;
        // the CMMotionRemover will remove motion of centor of mass in this function
        // if they are used in the simulation. It will return true if the positions
        // are modfied, and the potential energy and forces should be re-computed.
        // Otherwise, the previous effective forces can be used directly in the
        // following first half step propagation. At the initial step (step 0),
        // the effective forces have been computed in the beforeOneTraj().
        if (ha->updateContextState()) {
            // Update H and F_all which is defined in the HamiltonianOpenMM.
            updateH();
            // Update effective forces used for nuclear propagation
            updateF();
        }
        // Do nuclear propagation with velocity Verlet integrator
        // Part 1 of velocity Verlet integrator: the first CustomIntegrator
        // update V with first half step using current effective forces and update
        // the R with a full step and then do the distances constraints.
        compound->setCurrentIntegrator(0);
        compound->oneStep(F);
        // Update H and F_all based on the new positions.
        updateH();
        // Update electronic q, p, coeff with RK4 method using effective Hamiltonian
        // of current and previous step. The effective Hamiltonian at intermediate
        // time used for electronic propagaton is computed by linear interpolation.
        // Note, DT is in ps, and it should be converted to au for electronic propagation.
        DyElec->MOVE_elec(step*DT*ps2au, H_old, H, q, p, coeff);
        // Update effective forces with new q,p used for next nuclear propagation
        updateF();
        // Part 2 of velocity Verlet integrator: the second CustomIntegrator
        // update V with second half step with new effective forces: v = v+0.5*dt*f/m
        // and do the velocities constraints so the net velocity of each constrained
        // distance is zero.
        compound->setCurrentIntegrator(1);
        compound->oneStep(F);
        // Finally, we need to modify the time/step in OpenMM platform data since
        // in the above strategy, the simulation was propagated with two steps.
        ha->setTime(ha->getTime() - DT);
        ha->setStep(ha->getStep() - 1);
    }
    // Other OpenMM integrators such as leap-frog/Langevin/LangevinMiddle/NoseHoover
    // Within these integrators, the electronic propagation is performed after the
    // finish of nuclear propagtaion with a full step.
    else {
        // At the start of each time step, we need to update OpenMM context state.
        // The thermostat will modify velocities; the barostat will modify positions;
        // the CMMotionRemover will remove motion of centor of mass in this function
        // if they are used in the simulation. It will return true if the positions
        // are modfied, and the potential energy and forces should be re-computed.
        // Otherwise, the previous effective forces can be used directly in the
        // following first half step propagation. At the initial step (step 0),
        // the effective forces have been computed in the beforeOneTraj().
        if (ha->updateContextState()) {
            // Update H and F_all which is defined in the HamiltonianOpenMM.
            updateH();
            // Update effective forces used for nuclear propagation
            updateF();
        }
        // Do nuclear propagation within OpenMM integrator with the effective forces
        ha->integrator->oneStep(F);
        // Update H and F_all based on the new positions.
        updateH();
        // Update electronic q, p, coeff with RK4 method using effective Hamiltonian
        // of current and previous step. The effective Hamiltonian at intermediate
        // time used for electronic propagaton is computed by linear interpolation.
        // Note, DT is in ps, and it should be converted to au for electronic propagation.
        DyElec->MOVE_elec(step*DT*ps2au, H_old, H, q, p, coeff);
        // Update effective forces with new q,p used for next nuclear propagation
        updateF();
    }
}

void DynamicsMQCBase::oneStepForCustomAllAtom() {
    std::cout<<" oneStepForCustomAllAtom  000"<<std::endl;
    // This function is modified based on the new strategy on Dec. 3, 2021.
    // New strategy: the data of positions and velocities are stored in the
    // OpenMM Platform and the integration are performed by OpenMM Integrator.
    // But the forces used to do nuclear propagtaion is provided externally.
    // At each step, the force of each state should be downloaded from
    // OpenMM Platform such as CUDA, and weighted by the electronic population
    // and coherence to get the effective forces, then the effective forces will
    // be uploaded to the OpenMM Platform to do the nuclear propagtaion.
    // Link to the orignal data in classes of HamiltonianOpenMM and HamiltonianElec
    static std::shared_ptr<HamiltonianForceFieldBase> ha = std::static_pointer_cast<HamiltonianForceFieldBase>(Ha);
    static std::vector<Complex>& coeff = Elec->coeff; // coefficient of electronic wavefunction
    static std::vector<double>& q = Elec->q; // electronic mapping coordinates
    static std::vector<double>& p = Elec->p; // electronic mapping momenta
    static const Real_Matrix&     H = ha->Heff; // diabatic effecive Hamiltonian matrix (Havg removed)
    static const Real_Matrix& H_old = ha->Heff_old; // previous effecive diabatic Hamiltonian
    static const std::string integrator = param.getStr("integrator"); // OpenMM integrator
    // The nuclear propagation is perfromed by velocity Verlet integrator.
    // Note that within velocity Verlet integrator, the electronic propagtion
    // is performed after the updating of nuclear position before the second
    // half step of updating velcities.
    if (integrator == "velocityVerlet") {
        // Update the positions (R) and velocities (V) in Hamiltonian object (not
        // OpenMM Context) with velocity Verlet integrator. Update velocities with a
        // first half step using effective force. Then, update positions with a full step.
        MOVE_Half_V(ha->masses, ha->F, ha->V);
        MOVE_R(ha->masses, ha->V, ha->R);
        // Update H and F_all based on the new positions.
        updateH();
        // Update electronic q, p, coeff with RK4 method using effective Hamiltonian
        // of current and previous step. The effective Hamiltonian at intermediate
        // time used for electronic propagaton is computed by linear interpolation.
        // Note, DT is in ps, and it should be converted to au for electronic propagation.
        DyElec->MOVE_elec(step*DT*ps2au, H_old, H, q, p, coeff);
        // Update effective forces with new q,p used for next nuclear propagation
        updateF();
        // Update velocities with a second half step using new effective force.
        MOVE_Half_V(ha->masses, ha->F, ha->V);
    }
    // Other OpenMM integrators such as leap-frog/Langevin/LangevinMiddle/NoseHoover
    // Within these integrators, the electronic propagation is performed after the
    // finish of nuclear propagtaion with a full step.
    else {
        // Update the positions (R) and velocities (V) in Hamiltonian object (not
        // OpenMM Context) with leap-frog integrator. Update velocities with a full
        // step using effective force. Then, update positions with a full step.
        MOVE_V(ha->masses, ha->F, ha->V);
        MOVE_R(ha->masses, ha->V, ha->R);
        // Update H and F_all based on the new positions.
        updateH();
        // Update electronic q, p, coeff with RK4 method using effective Hamiltonian
        // of current and previous step. The effective Hamiltonian at intermediate
        // time used for electronic propagaton is computed by linear interpolation.
        // Note, DT is in ps, and it should be converted to au for electronic propagation.
        DyElec->MOVE_elec(step*DT*ps2au, H_old, H, q, p, coeff);
        // Update effective forces with new q,p used for next nuclear propagation
        updateF();
    }
    std::cout<<" oneStepForCustomAllAtom  111"<<std::endl;
}

void DynamicsMQCBase::updateAllAtomForces() {
    // This is very simillar as the updateDiabaticModelForces(), but using
    // different data type for forces, here is OpenMM::Vec3. (modified on Dec. 3, 2021)
    // Link to the orignal data in classes of HamiltonianOpenMM and HamiltonianElec
    static std::shared_ptr<HamiltonianOpenMM> ha = std::static_pointer_cast<HamiltonianOpenMM>(Ha);
    static const std::vector<double>& q = Elec->q; // electronic mapping coordinates
    static const std::vector<double>& p = Elec->p; // electronic mapping momenta
    static const std::vector<std::vector<OpenMM::Vec3>>& F_all = ha->F_all; // diabatic forces of each states
    static const std::vector<OpenMM::Vec3>& F_avg = ha->F_avg; // average forces
    static std::vector<OpenMM::Vec3>& F = ha->F; // effective forces to be updated
    // Reset all elements to 0 with keeping same size
    std::fill(F.begin(), F.end(), OpenMM::Vec3());
    // Note: in the following code , the hbar is ommitted since it will be vanished.
    for (int j = 0; j < DOFn; ++j) {
        for (int n = 0; n < DOFe; ++n)
            for (int m = 0; m < DOFe; ++m)
                if (n == m) // from diagnoal of H (weighted by population)
                    F[j] += (F_all[n*DOFe+n][j] - F_avg[j]) * 0.5 * (q[n]*q[n] + p[n]*p[n]);
                else if (n < m) // from off-diagnoal of H (weighted by coherence)
                    // factor 0.5 in coherence is ommited since nm is same as mn (multiply 2)
                    // If Condon approximation is used (coupling is contant), it is 0.
                    F[j] += F_all[n*DOFe+m][j] * (q[n]*q[m] + p[n]*p[m]);
        F[j] += F_avg[j];
    }
}

void DynamicsMQCBase::updateCustomAllAtomForces() {
    // This is very simillar as the updateDiabaticModelForces(), but using
    // different data type for forces, here is Vec3. (modified on Dec. 3, 2021)
    // Link to the orignal data in classes of HamiltonianOpenMM and HamiltonianElec
    static std::shared_ptr<HamiltonianForceFieldBase> ha = std::static_pointer_cast<HamiltonianForceFieldBase>(Ha);
    static const std::vector<double>& q = Elec->q; // electronic mapping coordinates
    static const std::vector<double>& p = Elec->p; // electronic mapping momenta
    std::cout<<"  test Elec q,p: q[0]:  "<<q[0]<<" p[0]:  "<<p[0]<<std::endl;
    static const std::vector<std::vector<Vec3>>& F_all = ha->F_all; // diabatic forces of each states
    static const std::vector<Vec3>& F_avg = ha->F_avg; // average forces
    static std::vector<Vec3>& F = ha->F; // effective forces to be updated
    // Reset all elements to 0 with keeping same size
    std::fill(F.begin(), F.end(), Vec3());
    // Note: in the following code , the hbar is ommitted since it will be vanished.
    for (int j = 0; j < DOFn; ++j) {
        for (int n = 0; n < DOFe; ++n)
            for (int m = 0; m < DOFe; ++m)
                if (n == m) // from diagnoal of H (weighted by population)
                    F[j] += (F_all[n*DOFe+n][j] - F_avg[j]) * 0.5 * (q[n]*q[n] + p[n]*p[n]);
                else if (n < m) // from off-diagnoal of H (weighted by coherence)
                    // factor 0.5 in coherence is ommited since nm is same as mn (multiply 2)
                    // If Condon approximation is used (coupling is contant), it is 0.
                    F[j] += F_all[n*DOFe+m][j] * (q[n]*q[m] + p[n]*p[m]);
        F[j] += F_avg[j];
    }
}

// # The following three functions are used to propagate nuclear DOF with
// # velocity or leapfrog Verlet integrator with OpenMM in diabatic basis.
// # Note, the propagtaion is performed by own CPU code. (Plan A)
// # In this propagation, the virtural sites and constraints are not supported.
void DynamicsMQCBase::customCPUIntegrator() {
    // Link to the orignal data in classes of HamiltonianAllAtom and HamiltonianElec
    static std::shared_ptr<HamiltonianOpenMM> ha = std::static_pointer_cast<HamiltonianOpenMM>(Ha);
    static std::vector<Complex>& coeff = Elec->coeff; // coefficient of electronic wavefunction
    static std::vector<double>& q = Elec->q; // electronic mapping coordinates
    static std::vector<double>& p = Elec->p; // electronic mapping momenta
    static std::vector<OpenMM::Vec3>& R = ha->R; // nuclear coordinates
    static std::vector<OpenMM::Vec3>& V = ha->V; // nuclear velecoites
    static std::vector<OpenMM::Vec3>& F = ha->F; // effective force
    static const std::vector<double>& masses = ha->masses; // masses of atoms
    static const Real_Matrix&     H = ha->Heff; // diabatic effecive Hamiltonian matrix (Havg removed)
    static const Real_Matrix& H_old = ha->Heff_old; // previous diabatic Hamiltonian
    static const std::string integrator = param.getStr("integrator"); // leapfrog or velocityVerlet
    static const bool useBarostat = param.getStr("barostat") == "none" ? false : true;
    static const bool useThermostat = param.getStr("thermostat") == "none" ? false : true;
    static const bool remove_COMMotion = param.getBool("remove_COMMotion");
    static const bool hasVirtualSites = ha->hasVirtualSites();
    static const bool hasConstraints = ha->hasConstraints();
    if (step == 0 && (hasVirtualSites || hasConstraints))
        throw std::runtime_error("ERROR: Virtual sites and constraints are not supported for custom CPU propagation.");
    // At the start of each time step, we need to update OpenMM context state.
    // The thermostat will modify velocities; the barostat will modify positions;
    // the CMMotionRemover will remove motion of centor of mass in this function
    // if they are used in the simulation. It will return true if the positions
    // are modfied, and the potential energy and forces should be re-computed.
    // Otherwise, the previous effective forces can be used directly in the
    // following first half step propagation. At the initial step (step 0),
    // the effective forces have been computed in the beforeOneTraj().
    // Note: only the Andersen thermostat and MonteCarlo barostat can be used.
    if (ha->updateContextState()) {
        // Update H and F_all which is defined in the HamiltonianOpenMM.
        updateH();
        // Update effective forces used for nuclear propagation
        updateF();
    }
    // We need copy positions/velocities from OpenMM Context at the initial step
    // or they changed after ha->updateContextState().
    if (useBarostat || step == 0)
        ha->getPositions();
    if (useThermostat || remove_COMMotion || step == 0)
        ha->getVelocities();
    // The nuclear propagation is perfromed by velocity Verlet integrator.
    // Note that within velocity Verlet integrator, the electronic propagtion
    // is performed after the updating of nuclear position before the second
    // half step of updating velcities.
    if (integrator == "velocityVerlet") {
        // Update the positions (R) and velocities (V) in Hamiltonian object (not
        // OpenMM Context) with velocity Verlet integrator. Update velocities with a
        // first half step using effective force. Then, update positions with a full step.
        MOVE_Half_V(masses, F, V);
        MOVE_R(masses, V, R);
        // ? After the updating of positions, we should apply distance constraints here.
        // ? However, in this case, the OpenMM constraint can not be applied externally.
        // ? Therefore the constraints is not supported in this way.
        // Set the new postions into OpenMM Context.
        ha->setPositions(R);
        // Update H and F_all based on the new positions.
        updateH();
        // Update electronic q, p, coeff with RK4 method using effective Hamiltonian
        // of current and previous step. The effective Hamiltonian at intermediate
        // time used for electronic propagaton is computed by linear interpolation.
        // Note, DT is in ps, and it should be converted to au for electronic propagation.
        DyElec->MOVE_elec(step*DT*ps2au, H_old, H, q, p, coeff);
        // Update effective forces with new q,p used for next nuclear propagation
        updateF();
        // Update velocities with a second half step using new effective force.
        MOVE_Half_V(masses, F, V);
        // ? Make velocity correction (deltaR/DT) here, if constraints are used.
        // ? Apply velocity constraints here, so the net velocity along all constraints is 0,
        // ? if constraints are used. The constraints is not supported in this way.
        // Set the new velocities into OpenMM Context.
        ha->setVelocities(V);
    }
    // Within leapfrog integrator, the electronic propagation is performed after
    // the finish of nuclear propagtaion with a full step.
    else if (integrator == "leapfrog") {
        // Update the positions (R) and velocities (V) in Hamiltonian object (not
        // OpenMM Context) with leap-frog integrator. Update velocities with a full
        // step using effective force. Then, update positions with a full step.
        MOVE_V(masses, F, V);
        MOVE_R(masses, V, R);
        // ? After the updating of positions, we should apply distance constraints here.
        // ? However, in this case, the OpenMM constraint can not be applied externally.
        // ? Therefore the constraints is not supported in this way.
        // Set the new postions and velocities into OpenMM Context.
        ha->setPositions(R);
        ha->setVelocities(V);
        // Update H and F_all based on the new positions.
        updateH();
        // Update electronic q, p, coeff with RK4 method using effective Hamiltonian
        // of current and previous step. The effective Hamiltonian at intermediate
        // time used for electronic propagaton is computed by linear interpolation.
        // Note, DT is in ps, and it should be converted to au for electronic propagation.
        DyElec->MOVE_elec(step*DT*ps2au, H_old, H, q, p, coeff);
        // Update effective forces with new q,p used for next nuclear propagation
        updateF();
    }
    else
        throw std::runtime_error("ERROR: Only leapfrog or velocityVerlet integrator is supported for custom CPU propagation.");
    // Now, the positions and velocities in OpenMM Context and Hamiltonian object
    // are at the time (t+DT). Then update the simluation time (ps) and step in Context.
    // This is useful, since something are dependent on the current step or time,
    // such as the frequnceny of CMMotionRemover.
    ha->setTime(ha->getTime() + DT);
    ha->setStep(ha->getStep() + 1);
}