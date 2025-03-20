/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "DynamicsBase.h"
#include "Vec3.h"

class DynamicsReadTraj : public DynamicsBase {
public:
    /**
     * Construct a DynamicsElec object.
     *
     * @param param   the global paramters
     */
    DynamicsReadTraj(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsBase(param, Ha) {}
    ~DynamicsReadTraj() {}
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
};