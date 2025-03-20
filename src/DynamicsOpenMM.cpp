/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 19, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsOpenMM.h"

void DynamicsOpenMM::init() {
    DynamicsBase::init();
    if (system_type != "allatom")
        throw std::runtime_error("ERROR: Unsupported system_type=" + system_type + " for OpenMM.");
    if (dyn_type.substr(0,6) != "OpenMM")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for OpenMM.");
    // Create smarter pointer to OpenMM intergrator
    ha = std::static_pointer_cast<HamiltonianOpenMM>(Ha);
    integrator = ha->integrator;
    // Get time step size (in ps) from OpenMM integrator.
    DT = integrator->getStepSize();
}

void DynamicsOpenMM::beforeOneTraj() {
    // Load R,V from trajectory if nucl_load is not empty, otherwise do nothing.
    samplingNucl();
    // Get current time (in ps) and step from OpenMM Context
    // This may not zero, if it is a restart simulation.
    // Here, round() will return the nearest integer value.
    // For multi-traj simulation, they are zero.
    step = round(ha->getTime()/DT);
}

void DynamicsOpenMM::dynamics(int steps) {
    if (dyn_type == "OpenMM") // original OpenMM propagtaion
        integrator->step(steps);
    else // propagation with external forces.
        for (int i = 0; i < steps; ++i)
            myOneStep();
    this->step += steps;
}

//add by Xu
/*
void DynamicsOpenMM::perturbDynamics(int steps, int perturbStep, int atomIndex, double fx, double fy, double fz) {
    
    if (perturbStep < 0 || perturbStep >= steps) {
        std::cerr << "Warning: perturbStep " << perturbStep << " is outside the range [0," << steps-1 << "]" << std::endl;
       
        dynamics(steps);
        return;
    }
    
    
    if (perturbStep > 0) {
        dynamics(perturbStep);
    }

    
    static const int propagate_state = param.getInt("propagate_state");
    
   
    ha->updateContextState();
    ha->getPotentialEnergy(propagate_state, true);
    ha->getForces();
    
    
    if (atomIndex >= 0 && atomIndex < ha->F.size()) {
        std::cout << "Applying custom force at step " << step << " to atom " << atomIndex << std::endl;
        std::cout << "Original force: (" 
                  << ha->F[atomIndex][0] << ", " 
                  << ha->F[atomIndex][1] << ", " 
                  << ha->F[atomIndex][2] << ")" << std::endl;
        
        
        ha->F[atomIndex][0] += fx;
        ha->F[atomIndex][1] += fy;
        ha->F[atomIndex][2] += fz;
        
        std::cout << "Modified force: (" 
                  << ha->F[atomIndex][0] << ", " 
                  << ha->F[atomIndex][1] << ", " 
                  << ha->F[atomIndex][2] << ")" << std::endl;
    } else {
        std::cerr << "Warning: Atom index " << atomIndex 
                  << " is out of bounds (0-" << ha->F.size()-1 << ")" << std::endl;
    }
    
    
    if (dyn_type == "OpenMM") { 
        integrator->oneStep(ha->F);
    } else { 
        myOneStep();
    }
    
    
    this->step += 1;
    
    
    int remainingSteps = steps - perturbStep - 1;
    if (remainingSteps > 0) {
        dynamics(remainingSteps);
    }
}
*/


/*
void DynamicsOpenMM::perturbDynamics(int steps, int perturbStep, int atomIndex, double fx, double fy, double fz) {
    if (perturbStep < 0 || perturbStep >= steps) {
        std::cerr << "Warning: perturbStep " << perturbStep << " is outside the range [0," << steps-1 << "]" << std::endl;
        dynamics(steps);
        return;
    }
    
    
    if (perturbStep > 0) {
        dynamics(perturbStep);
    }

    
    static const int propagate_state = param.getInt("propagate_state");
    
    
    ha->updateContextState();
    
    
    if (atomIndex == -1) {
        std::cout << "Applying ForceFieldPolar calculated perturbF at step " << step << " to all atoms" << std::endl;
        
        
        ForceFieldPolar* ffPolar = dynamic_cast<ForceFieldPolar*>(ha->getForceField());
        if (ffPolar == nullptr) {
            std::cerr << "Error: ForceField is not of type ForceFieldPolar" << std::endl;
            return;
        }
        
        
        ffPolar->calcualtePerturbPolarForce(ha->F);
        
        std::cout << "Applied perturbF to all atoms" << std::endl;
    } 
    
    else if (atomIndex >= 0 && atomIndex < ha->F.size()) {
        std::cout << "Applying custom force at step " << step << " to atom " << atomIndex << std::endl;
        std::cout << "Original force: (" 
                  << ha->F[atomIndex][0] << ", " 
                  << ha->F[atomIndex][1] << ", " 
                  << ha->F[atomIndex][2] << ")" << std::endl;
        
        
        ha->F[atomIndex][0] += fx;
        ha->F[atomIndex][1] += fy;
        ha->F[atomIndex][2] += fz;
        
        std::cout << "Modified force: (" 
                  << ha->F[atomIndex][0] << ", " 
                  << ha->F[atomIndex][1] << ", " 
                  << ha->F[atomIndex][2] << ")" << std::endl;
    } else {
        std::cerr << "Warning: Atom index " << atomIndex 
                  << " is out of bounds (0-" << ha->F.size()-1 << ")" << std::endl;
    }
    
    
    if (dyn_type == "OpenMM") { 
        integrator->oneStep(ha->F);
    } else { 
        myOneStep();
    }
    
    
    this->step += 1;
    
    
    int remainingSteps = steps - perturbStep - 1;
    if (remainingSteps > 0) {
        dynamics(remainingSteps);
    }
}
*/

//new version
void DynamicsOpenMM::perturbDynamics(int steps, int perturbStep, int forceFieldIndex, double scaleFactor) {
    if (perturbStep < 0 || perturbStep >= steps) {
        std::cerr << "Warning: perturbStep " << perturbStep << " is outside the range [0," << steps-1 << "]" << std::endl;
        dynamics(steps);
        return;
    }
    
    std::cout << "Will apply perturbation forces from ForceField " << forceFieldIndex 
              << " at step " << (step + perturbStep) 
              << " with scale factor " << scaleFactor << std::endl;
    
    // Run dynamics up to the perturbation step
    if (perturbStep > 0) {
        dynamics(perturbStep);
    }
    // Apply perturbation forces
    applyPerturbForces(forceFieldIndex, scaleFactor);
                                 
    // Continue with remaining steps
    int remainingSteps = steps - perturbStep - 1;
    if (remainingSteps > 0) {
        dynamics(remainingSteps);
    }
}


void DynamicsOpenMM::myOneStep() {
    static const int propagate_state = param.getInt("propagate_state");
    // Use CompoundIntegrator with two CustomIntegrator to do velocity Verlet
    // integrator with external forces.
    if (param.getStr("integrator") == "velocityVerlet") {
        static std::shared_ptr<OpenMM::CompoundIntegrator> compound =
            std::static_pointer_cast<OpenMM::CompoundIntegrator>(integrator);
        // Part 1 of velocity Verlet integrator:
        // update V with first half step and the R with a full step and the
        // distances constraints.
        // We need compute force at initial step or positions are changed by
        // barostat in ha->updateContextState() [return true if changed]
        // If the positions are the same, then we don't need to compute force here.
        if (ha->updateContextState() || step == 0) {
            ha->getPotentialEnergy(propagate_state, true);
            ha->getForces();
        }
        compound->setCurrentIntegrator(0);
        compound->oneStep(ha->F);
        // Part 2 of velocity Verlet integrator:
        // update V with second half step with new forces: v = v+0.5*dt*f/m
        // and the velocities constraints.
        // Calculate forces based on the updated R
        ha->getPotentialEnergy(propagate_state, true);
        ha->getForces();
        compound->setCurrentIntegrator(1);
        compound->oneStep(ha->F);
        // Note: we need to modify the time/step in platform data since in the
        // above strategy, the simulation was propagated with two steps.
        ha->setTime(ha->getTime() - DT);
        ha->setStep(ha->getStep() - 1);
    }
    else { // The original OpenMM integrator
        ha->updateContextState();
        ha->getPotentialEnergy(propagate_state, true);
        ha->getForces();
        integrator->oneStep(ha->F);
    }
}


// add by xu const f 
/*
void DynamicsOpenMM::applyCustomForce(int atomIndex, double fx, double fy, double fz) {
    
    static const int propagate_state = param.getInt("propagate_state");
    ha->updateContextState();
    ha->getPotentialEnergy(propagate_state, true);
    ha->getForces();
    
    
    if (atomIndex >= 0 && atomIndex < ha->F.size()) {
        std::cout << "Applying custom force at step " << step << " to atom " << atomIndex << std::endl;
        std::cout << "Original force: (" << ha->F[atomIndex][0] << ", " 
                 << ha->F[atomIndex][1] << ", " 
                 << ha->F[atomIndex][2] << ")" << std::endl;
        
        
        ha->F[atomIndex][0] += fx;
        ha->F[atomIndex][1] += fy;
        ha->F[atomIndex][2] += fz;
        
        std::cout << "Modified force: (" << ha->F[atomIndex][0] << ", " 
                 << ha->F[atomIndex][1] << ", " 
                 << ha->F[atomIndex][2] << ")" << std::endl;
    } else {
        std::cerr << "Warning: Atom index " << atomIndex 
                  << " is out of bounds (0-" << ha->F.size()-1 << ")" << std::endl;
    }
    
    if (param.getStr("integrator") == "velocityVerlet") {
        static std::shared_ptr<OpenMM::CompoundIntegrator> compound =
            std::static_pointer_cast<OpenMM::CompoundIntegrator>(integrator);
        
        compound->setCurrentIntegrator(0);
        compound->oneStep(ha->F);
        
        
        ha->getPotentialEnergy(propagate_state, true);
        ha->getForces();
        
        compound->setCurrentIntegrator(1);
        compound->oneStep(ha->F);
        
        ha->setTime(ha->getTime() - DT);
        ha->setStep(ha->getStep() - 1);
    } else {
        integrator->oneStep(ha->F);
    }
}
*/


//add perturb F DKPI



void DynamicsOpenMM::applyPerturbForces(int forceFieldIndex, double scaleFactor) {
    static const int propagate_state = param.getInt("propagate_state");
    
    
    ha->updateContextState();
    ha->getPotentialEnergy(propagate_state, true);
    ha->getForces();
    
   
    std::vector<OpenMM::Vec3> perturbForces;
    
    if (ha->getForceFieldForces(perturbForces)) {
        std::cout << "Applying perturbation forces from ForceField " << forceFieldIndex << std::endl;
        
        
        if (!ha->F.empty() && !perturbForces.empty()) {
            
            std::cout << "First atom original force: (" 
                      << ha->F[0][0] << ", " 
                      << ha->F[0][1] << ", " 
                      << ha->F[0][2] << ") kJ/mol/nm" << std::endl;
            
            
            std::cout << "First atom perturbation force: (" 
                      << perturbForces[0][0] * scaleFactor << ", " 
                      << perturbForces[0][1] * scaleFactor << ", " 
                      << perturbForces[0][2] * scaleFactor << ") kJ/mol/nm" << std::endl;
        }
        
        
        for (int i = 0; i < ha->F.size(); i++) {
           
            ha->F[i][0] += perturbForces[i][0] * scaleFactor;
            ha->F[i][1] += perturbForces[i][1] * scaleFactor;
            ha->F[i][2] += perturbForces[i][2] * scaleFactor;
        }
        
       
        if (!ha->F.empty()) {
            std::cout << "First atom modified force: (" 
                      << ha->F[0][0] << ", " 
                      << ha->F[0][1] << ", " 
                      << ha->F[0][2] << ") kJ/mol/nm" << std::endl;
        }
    } else {
        std::cerr << "Warning: Failed to get perturbation forces from ForceField " << forceFieldIndex << std::endl;
    }
    
    
    ha->uploadForces();
    
    
    if (dyn_type == "OpenMM") { 
        integrator->oneStep(ha->F);
    } else { 
        myOneStep();
    }
    
    this->step += 1;
}

