/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsMFRDM.h"

void DynamicsMFRDM::init() {
    // 0. Initialize date members of DynamicsMQCBase.
    DynamicsMQCBase::init();
    if (dyn_type != "MF-RDM")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for DynamicsMFRDM.");
    if (representation != "diabatic")
            throw std::runtime_error("ERROR: Only diabatic representation is supported for dyn_type=" + dyn_type);
    if (system_type != "model")
        throw std::runtime_error("ERROR: Realistic system is not supported for dyn_type=" + dyn_type);
    // For Ehrenfest method, only one initial state can be computed.
    if (init_state.size() != 1)
        throw std::runtime_error("ERROR: Multi-initial states is not supported for Ehrenfest dynamics.");
    // initialize the initial reduced density matrix (sigma_init, σ(0))
    // The element specified by init_state is 1, and others are 0
    // Here, index init_state = j*DOFe+k, |j><k|
    sigma_init.resize(DOFe, std::vector<Complex>(DOFe, 0));
    sigma_init[init_state[0]/DOFe][init_state[0]%DOFe] = 1.0;
    sigma = sigma_init;
}

void DynamicsMFRDM::samplingElec() {
    // initialize reduced density matrix sigma
    sigma = sigma_init;
}

void DynamicsMFRDM::updateDensityMatrix() {
    int i = 0;
    // For Ehrenfest method, RDM is directly computed by σ_jk=φ*(k)φ(j).
    // Here, the σ_jk is directly propagated.
    for (int j = 0; j < DOFe; j++)
        for (int k = 0; k < DOFe; k++) {
            // Here, RDM is stored as DOFe^2-dimentional vector (index j*DOFe+k).
            // RDM_current[0][step/RDM_steps][j*DOFe+k] = sigma[j][k];
            // For RDM_average, the value of each element is accumaleated (+=)
            // and averaged by number of trajectories.
            // RDM_average[0][step/RDM_steps][j*DOFe+k] += sigma[j][k] / (double)(ntraj);

            //calculate RDM of the snapshot
            RDM[i][j*DOFe+k] = sigma[j][k];
        }
}

// TODO using general one in MappingBase, only MOVE_Elec is differnent
void DynamicsMFRDM::oneStepForDiabaticModel() {
    // This is Velocity Verlet integrator
    std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    // Update velocities/momenta with a half step using effective force.
    MOVE_Half_V(ha->F, ha->V);
    // Update positions with a full step
    MOVE_R(ha->V, ha->R);
    // Update H and F_all since R has been changed.
    updateH();
    // Update current reduced density matrix (sigma) by diagonalize the Hamiltonian.
    // Using the original Hamilton (H), not the effective Hamiltonian (Heff)
    DyElec->MOVE_elec(ha->H, sigma);
    // Update effective forces
    updateF();
    // Update velocities/momenta with a half step using new effective forces
    MOVE_Half_V(ha->F, ha->V);
}

// TODO using general one, only pop and coh is different
void DynamicsMFRDM::updateDiabaticModelForces() {
    std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    // Reset all elements to 0 with keeping same size
    std::fill(ha->F.begin(), ha->F.end(), 0.0);
    // Same as that in DynamcisMappingBase, but population and coherence operators
    // are directly from diagnaol and off-diagnoal elements of current RDM (sigma).
    for (int j = 0; j < DOFn; j++) {
        for (int n = 0; n < DOFe; n++)
            for (int m = 0; m < DOFe; m++)
                if (n == m) // from diagnoal of H (population)
                    ha->F[j] += (ha->F_all[n*DOFe+m][j] - ha->F_avg[j]) * sigma[n][m].real();
                else if (n < m) // from off-diagnoal of H (coherence)
                    // multiply 2 since nm is same as mn
                    // If Condon approximation is used (coupling is contant), it is 0.
                    ha->F[j] += ha->F_all[n*DOFe+m][j] * 2.0 * sigma[n][m].real();
        ha->F[j] /= hbar;
        ha->F[j] += ha->F_avg[j];
    }
}