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
 * This class defines original Ehrenfest mean-filed dynamics method.
 *
 * Unlike the DynamicsMF, the version of the Ehrenfest method defined here is
 * the RDM verison, in which no electronic mapping variables are required and
 * propagate reduced density matrix (RDM, σ) directly, thus can start from the
 * off-diagnoal element (coherence) of initial RDM.
 *
 * It should be noted that even though the initial electronic RDM is same for
 * all trajectories, different initial nuclear phase spase will give rise to
 * different final reduced density matrix.
 */
class DynamicsMFRDM : public DynamicsMQCBase {
public:
    /**
     * Construct a DynamicsMFRDM object.
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     */
    DynamicsMFRDM(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsMQCBase(param, Ha) {}
    ~DynamicsMFRDM() {}
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
     * For Ehrenfest method, just initialize sigma as initial sigma_init.
     */
    void samplingElec();
    /**
     * No initial window or denisty in Ehrenfest method.
     */
    void getInitialDensity() {}
    /**
     * Update the electronic reduced density matrix at current time.
     *
     * Both the current (RDM_current) and average (RDM_average) ones are updated.
     */
    void updateDensityMatrix();
    /**
     * Update the effective forces which will be used nuclear propagation.
     * This function is used for a simulation of model.
     *
     * Same as that in DynamcisMappingBase, but the population or coherence
     * operator is directly from Ehrenfest RDM (sigma).
     */
    void updateDiabaticModelForces() override;
    /**
     * Propagate one step for a simulation of model using velocity Verlet integrator.
     */
    void oneStepForDiabaticModel() override;

private:
    // sigma_init is initial reduced denisty matrix at time 0, σ(0), which is
    // specified by user with the indices of initial state, |j><k|, and sigma is
    // reduced denisty matrix of at time t, σ(t), which will be propagated directly.
    // σ_jk = φ*(k)φ(j), φ is electronic wavefunction, but no explict defined here.
    Complex_Matrix sigma, sigma_init;
};