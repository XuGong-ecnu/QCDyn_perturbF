/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 19, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "DynamicsBase.h"

/**
 * This subclass is used to run a general MD simulation with OpenMM Integrator.
 * We don't define any dynamics method here, just make it became an interface.
 *
 * Note that all the data (positions, velocities, forces) that being propagated
 * is stored in the OpenMM internal Platform when running a general MD simulation
 * with OpenMM.
 *
 * Note that the internal OpenMM::Integrator can only do a general MD simulation
 * since the forces that are used to propagate is only one state.
 *
 * But It can run a special simulation in which the trajectory is propagated on
 * one arbitrary state, the energies of other any states (as long as the topologies
 * are provided) can be computed and reported on the fly, which can replace the
 * energy recalculation based on propagated trajectory in some cases.
 */
class DynamicsMixPES : public DynamicsBase {
public:
    DynamicsMixPES(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsBase(param, Ha) {}
    ~DynamicsMixPES() {}
    /**
     * Initialize data members and check the legality of parameters.
     */
    void init();
    /**
     * Load initial configurations (positions and/or velocities) from trajectory
     * call samplingNucl(), otherwise do nothing do nothing.
     */
    void beforeOneTraj();
    /**
     * Advance a simulation through time by taking a series of time steps.
     *
     * @param steps the number of time steps to take silently
     */
    void dynamics(int steps);
    /**
     * Do nothing for a general MD simulation with OpenMM.
     */
    void afterOneTraj() {}

private:
    void oneStep();
    void updateF();

private:
    // The smarter pointer to HamiltonianOpenMM object.
    std::shared_ptr<HamiltonianOpenMM> ha;
    // The smarter pointer to OpenMM Integrator object.
    // Integrator is an core object in OpenMM. Integrator implements an algorithm
    // for advancing the simulation through time, such as leapfrog Verlet, velocity
    // Verlet, Langvin integrators.
    std::shared_ptr<OpenMM::Integrator> integrator;
    // updateH() is to update Hamiltonian data that required
    // for dynamics, such as Hamiltonian matrix, forces, and/or nonadiabatic
    // coupling vectors, in diabtaic or adiabatic representation.
    // It stores a pointers to member functions in Hamiltonian object:
    // updateDiabaticHamiltonian() or updateAdiabaticHamiltonian() according to
    // representation. This binding is finished at init().
    std::function<void()> updateH;   
    std::vector<double> PES_weights; // the weight factor of each state to effective force  
};