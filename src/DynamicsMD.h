/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: Jan. 3, 2022                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "DynamicsBase.h"

class DynamicsMD : public DynamicsBase {
public:
    /**
     * Construct a DynamicsElec object.
     *
     * @param param   the global paramters
     */
    DynamicsMD(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsBase(param, Ha) {}
    ~DynamicsMD() {}
    /**
     * Initialize data members and check the legality of parameters.
     */
    void init();
     /**
     * Do nothing for a general MD simulation with ForceFieldBase.
     */
    void beforeOneTraj();
    /**
     * Advance a simulation through time by taking a series of time steps.
     *
     * @param steps the number of time steps to take silently
     */
    void dynamics(int steps);
    /**
     * Do nothing for a general MD simulation.
     */
    void afterOneTraj() {}

private:
    /**
     * Do nuclear propagation with one step with providing the external forces
     * (instead of the forces which is always computed from force fileds based
     * on current positions) at each step.
     */
    void myOneStep();
private:
    // The smarter pointer to HamiltonianForceFieldBase object.
    std::shared_ptr<HamiltonianForceFieldBase> ha;
    // Integrator implements an algorithm for advancing the simulation through time, 
    // such as leapfrog Verlet, velocity Verlet, Langvin integrators.
    std::string integrator;
    // judge the step keep the same.
    bool stepKeep;
    int currentStep;
};