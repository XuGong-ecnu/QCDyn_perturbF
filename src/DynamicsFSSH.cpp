/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsFSSH.h"

void DynamicsFSSH::init() {
    // 0. Initialize date members of DynamicsMQCBase.
    DynamicsMQCBase::init();
    if (dyn_type != "FSSH")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for DynamicsFSSH.");
    // For surface hopping, only starting from pure population state is allowed.
    if (init_state.size() != 1)
        throw std::runtime_error("ERROR: Multi-initial states is not supported for surface hopping method.");
    if (init_state[0]/DOFe != init_state[0]%DOFe) // init_state = j*DOFe+k
        throw std::runtime_error("ERROR: For surface hopping method, starting from coherence is not supported.");
    if (system_type != "model") // TODO:
        throw std::runtime_error("ERROR: Realistic system is not supported for dyn_type=" + dyn_type);
    if (representation != "diabatic") // TODO: adiabatic
        throw std::runtime_error("ERROR: Do surface hopping in adiabatic is not supported now.");

    // 1. Set random number seed for hopping probability
    const int FSSH_seed = param.getInt("FSSH_seed");
    if (FSSH_seed == 0) { // Use system-time seed for real simulation
        std::random_device rd;
        FSSH_gen.seed(rd());
    }
    else // Use deterministic seed for debug
        FSSH_gen.seed(FSSH_seed);
}

void DynamicsFSSH::samplingElec() {
    // electronic mapping variable: positions, momenta, coefficient of wavefunction
    std::vector<double>& q = Elec->q;
    std::vector<double>& p = Elec->p;
    std::vector<Complex>& coeff = Elec->coeff;
    std::fill(coeff.begin(), coeff.end(), 0.0);
    // For initial population state, coeff = 1.0, otherwise = 0.0.
    // This sampling is same as the Ehrenfest mean-filed method.
    // Note, here the value of init_state is j*DOFe+k (j=k), then the index of
    // the initial state is init_state[0]/DOFe.
    coeff[init_state[0]/DOFe] = 1.0;
    // get q, p from coeff.
    for (int i = 0; i < DOFe; i++) {
        q[i] = coeff[i].real() * sqrt(2);
        p[i] = coeff[i].imag() * sqrt(2);
    }
    // Set the active state at step 0.
    active_state = init_state[0]/DOFe; // Note it is the state index j
}

void DynamicsFSSH::updateDensityMatrix() {
    const std::vector<Complex>& coeff = Elec->coeff;
    int i = 0;
    // For FSSH method, RDM is directly computed by σ_jk=φ*(k)φ(j)
    // It is totally same as the Ehrenfest method.
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

void DynamicsFSSH::oneStepForDiabaticModel() {
    // This is Velocity Verlet integrator
    std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    // * This is special for FSSH method
    // If surface hopping occurs, update active state and scale velocities to
    // keep conserved energy, then compute forces of updated active state
    double scale = switchActiveState();
    if (scale > 0)  { // If factor is 0, then no hop
        for (int j = 0; j < DOFn; j++)
            ha->V[j] *= scale;
        updateF();
    }
    // Update velocities/momenta with a half step using effective force.
    MOVE_Half_V(ha->F, ha->V);
    // Update positions with a full step
    MOVE_R(ha->V, ha->R);
    // Update H and F_all since R has been changed.
    updateH();
    // Update coefficient of electronic wavefunction by diagonalize the Hamiltonian.
    // It will also update q, p. Here, the Hamiltonian is H at time t.
    // Note that for surface hopping method, using original Hamiltonian (H),
    // not the effective Hamiltonian (Heff).
    DyElec->MOVE_elec(ha->H_old, Elec->q, Elec->p, Elec->coeff);
    // Update effective forces
    updateF();
    // Update velocities/momenta with a half step using new effective forces
    MOVE_Half_V(ha->F, ha->V);
}

void DynamicsFSSH::updateDiabaticModelForces() {
    std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    // Here， F_all is DOFe * DOFe base, the i*DOFe+i is the force computed
    // from diagonal element of H. When Condon approximation is used. the
    // force computed from off-diagonal element of H is zero (gamma_DA is constant)
    ha->F = ha->F_all[active_state*DOFe+active_state];
}

double DynamicsFSSH::switchActiveState() {
    const std::vector<Complex>& coeff = Elec->coeff;
    // Refs: JCP 148, 102304 (2018) equation 28 - 33
    // transition probability from active state (a) to target states (b)
    std::vector<double> P_ab(DOFe, 0);
    const int a = active_state;
    // Sum of probability: from active state (a) to other states (index <=b)
    // as condition of hop or not in equation 31
    double P_sum = 0;
    // generate random number (uniform distribution) [0, 1]
    // Note: Here, we use a separated seed for hopping probability.
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    double randnum = uniform_dist(FSSH_gen);
    // prefactor of P_a->b in equation 29 30: = 2/|C_a(t)|^2 * 1/hbar * DT
    // in which, |C_a(t)|^2 is the population of this active state
    double prefactor = 2*DT/hbar / (std::conj(coeff[a])*coeff[a]).real();
    // If surface hopping successful, we need the scale factor of velocities
    double scale = 0.0;
    for (int b = 0; b < DOFe; b++) {
        if (b == a) continue;
        // Ha->H[a][b] is electronic diabatic coupling = imag part of <Φa|h|Φb>
        P_ab[b] = prefactor * Ha->H[a][b] * (std::conj(coeff[a])*coeff[b]).imag();
        if (P_ab[b] < 0)
            P_ab[b] = 0;
        P_sum += P_ab[b];
        // JCP 148, 102304 (2018) equation 31
        if (P_sum > randnum && (P_sum-P_ab[b]) <= randnum) {
            // Ref: JCP 148, 102304 (2018) equation 33
            // a hop to an energetically higher state can only occur if sufficient
            // classical kinetic energy is available. Otherwise, the transition needs
            // to be aborted (frustrated hops)
            double deltaPE = Ha->H[a][a] - Ha->H[b][b];
            double KE = Ha->getKineticEnergy();
            double TotalE_a = KE + Ha->H[a][a];
            double diff = 1 + deltaPE / KE;
            if (diff > 0)  // hop and scale velocities
                scale = sqrt(diff);
            active_state = b; // update active state
            break;
        }
    }
    // If scale is 0, then no hop
    return scale;
}