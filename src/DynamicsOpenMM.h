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
class DynamicsOpenMM : public DynamicsBase {
public:
    DynamicsOpenMM(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsBase(param, Ha) {}
    ~DynamicsOpenMM() {}
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
 //   void perturbDynamics(int steps, int perturbStep, int atomIndex, double fx, double fy, double fz);//add by Xu
    void afterOneTraj() {}
    void perturbDynamics(int steps, int perturbStep, int forceFieldIndex, double scaleFactor);
    
private:
    /**
     * Do nuclear propagation with one step with providing the external forces
     * (instead of the forces which is always computed from force fileds based
     * on current positions) at each step.
     */
    void myOneStep();
    //void applyCustomForce(int forceFieldIndex, double scaleFactor);   
    void applyPerturbForces(int forceFieldIndex, double scaleFactor);
private:
    // The smarter pointer to HamiltonianOpenMM object.
    std::shared_ptr<HamiltonianOpenMM> ha;
    // The smarter pointer to OpenMM Integrator object.
    // Integrator is an core object in OpenMM. Integrator implements an algorithm
    // for advancing the simulation through time, such as leapfrog Verlet, velocity
    // Verlet, Langvin integrators.
    std::shared_ptr<OpenMM::Integrator> integrator;
};
