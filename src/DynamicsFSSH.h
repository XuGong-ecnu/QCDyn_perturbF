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
 * This class defines original fewest switches surface-hopping (FSSH) method.
 *
 * Ref: Alexei A. Kananenka et. al. JCP 148, 102304 (2018).
 */
class DynamicsFSSH : public DynamicsMQCBase {
public:
    /**
     * Construct a DynamicsFSSH object.
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     */
    DynamicsFSSH(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsMQCBase(param, Ha) {}
    ~DynamicsFSSH() {}
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
     * Get initial coefficient of electronic wavefunction aacording to init_state.
     *
     * For surface hopping method, a pure population initial state is
     * accepted only.
     *
     * It is same as the Ehrenfest mean-filed method.
     */
    void samplingElec();
    /**
     * No initial window or denisty in surface hopping method.
     */
    void getInitialDensity() {}
    /**
     * Update the electronic reduced density matrix.
     * rho[j][k] = conj(coeff[k])*coeff[j], same as Ehrefest method.
     *
     * Both the current (RDM_current) and average (RDM_average) ones are updated.
     */
    void updateDensityMatrix();

    // # The forces and propagation for FSSH is different from other mapping dynamics,
    // # So we need explicitly define them to override the one in base class.
    // TODO: Implement adiabatic FSSH.
    /**
     * Propagate one step for a simulation of model using velocity Verlet integrator.
     */
    void oneStepForDiabaticModel() override;
    /**
     * Update the effective forces which will be used nuclear propagation.
     */
    void updateDiabaticModelForces() override;
    /**
     * Internal function to used for oneStepForDiabaticModel().
     * It is to check surface hopping occurs or not.
     * Ref: JCP 148, 102304 (2018) equation 28 - 33
     *
     * If surface hopping occurs, update active state and return the scale
     * factor (> 0) of velocities. If no hop, the return value is 0.
     *
     * @return  the scale factor of velocities (> 0); if no hop, it is 0.
     */
    double switchActiveState();

private:
    // The active state index (starting from 0) for FSSH method.
    // the force of this state is used to do nuclear propagation.
    int active_state;
    // random number generator, which is used hopping probability.
    // Within FSSH method, the surface hopping probability at each step is
    // decided by a stochastic algorithm. So, a random number generator
    // is required.
    std::mt19937 FSSH_gen;
};