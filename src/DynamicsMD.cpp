/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: Jan. 16, 2022                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsMD.h"

void DynamicsMD::init() {
    DynamicsBase::init();
    if (system_type != "allatom")
        throw std::runtime_error("ERROR: Unsupported system_type=" + system_type + " for CustomForceField.");
    if (allatom_type != "CustomForceField")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + allatom_type + " for CustomForceField."); 
    if (dyn_type != "MD")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for CustomForceField.");   
    if (ntraj != 1)
        throw std::runtime_error("ERROR: For a classical MD simulation with OpenMM, the value of ntraj must equal to 1.");
    // Create smarter pointer to OpenMM intergrator
    ha = std::static_pointer_cast<HamiltonianForceFieldBase>(Ha);
    integrator = param.getStr("integrator");
    DT = param.getDouble("DT");
}

void DynamicsMD::beforeOneTraj() {
    step = 0; // reset current step to 0.
    ha->updateAllForces(); 
}

void DynamicsMD::dynamics(int steps) {   
    for (int i = 0; i < steps; ++i) {
        myOneStep();  
    }
    if (steps != 1)
        this->step += steps;
}

void DynamicsMD::myOneStep() {
    if (integrator == "velocityVerlet") {
        // Update velocities with a full step using effective force.
        MOVE_Half_V(ha->masses, ha->F, ha->V);
        // Update positions with a full step
        MOVE_R(ha->masses, ha->V, ha->R);
        // Update Forces with new positions
        ha->updateAllForces();
        // Update momenta with a half step using new effective forces
        MOVE_Half_V(ha->masses, ha->F, ha->V);
    }    
    else if (integrator == "leapfrog") {
        // Update velocities with a full step using effective force.
        MOVE_V(ha->masses, ha->F, ha->V);
        // Update positions with a full step
        MOVE_R(ha->masses, ha->V, ha->R);
        // At the start of each time step, we need to update force
        ha->updateAllForces(); 
    }
    else
        throw std::runtime_error("ERROR: Unsupported integrator=" + integrator + " for CustomForceField.");
}