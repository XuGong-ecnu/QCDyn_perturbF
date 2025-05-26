/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 19, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsOpenMM.h"
#include "ForceFieldPolar.h"

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




//add by Xu //new version
void DynamicsOpenMM::perturbDynamics(int steps, int perturbStep, int forceFieldIndex, double scaleFactor, int pulse_type) {
    const int nsteps = param.getInt("nsteps");
    std::cout << "Will apply perturbation forces from ForceField " << forceFieldIndex 
              << " at step " << (step + perturbStep) 
              << " with scale factor " << scaleFactor << std::endl;
//    std::cout << "Apply force from ForceFieldPolar "<< std::endl;
    applyPerturbForces(forceFieldIndex, scaleFactor,pulse_type);
//    std::cout << "After apply atom force: ("
//                      << ha->F[0][0] << ", "
//                      << ha->F[0][1] << ", "
//                      << ha->F[0][2] << ") kJ/mol/nm" << std::endl;
    
//    std::cout << "Finish apply force from ForceFieldPolar "<< std::endl;                       
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

//add raman perturb F H_int = \Pi_ij Ejk \Pi_ki

void DynamicsOpenMM::applyPerturbForces(int forceFieldIndex, double scaleFactor, int pulse_type) {
    static const int propagate_state = param.getInt("propagate_state");
 
    ha->updateContextState();
    ha->getPotentialEnergy(propagate_state, true);
    ha->getKineticEnergy();
    ha->getForces();
    ha->getPositions(); 
    std::vector<OpenMM::Vec3> perturbForces;
    ha->getVelocities();

    perturbForces.resize(DOFn, OpenMM::Vec3(0.0, 0.0, 0.0));
    double kineticenergy = ha->getKineticEnergy();

    if (ha->getForceFieldForces(perturbForces, pulse_type)) {
        std::cout << "Applying perturbation forces from ForceField " << forceFieldIndex << std::endl;
        
        
        if (!ha->F.empty() && !perturbForces.empty()) {
            std::cout << std::fixed << std::setprecision(10);
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
        

    ha ->uploadForces();
//  ha->updateContextState();
    kineticenergy = ha->getKineticEnergy();
    std::cout << "After added perturb force kinetic energy is: "<<kineticenergy<<std::endl;
//    this->step += 1;
    integrator->oneStep(ha->F);
//    ha->getForces();
//    ha->getVelocity();
//    ha->updateContextState();
//    ha->getVelocities(); 
//    std::cout << "First atom velocity end interation is: ("
//                      << ha->V[0][0] << ", "
//                      << ha->V[0][1] << ", "
//                      << ha->V[0][2] << ") nm/ps" << std::endl;
//    kineticenergy = ha->getKineticEnergy();
//    std::cout << "First atom end kinetic energy is: "<<kineticenergy<<std::endl;
    ha->updateContextState();
    ha->getForces();
    ha->getPositions();
    ha->getVelocities();
/*
    if (ha->getForceFieldForces(perturbForces, pulse_type)) {
        std::cout << "Applying perturbation forces from ForceField " << forceFieldIndex << std::endl;


        if (!ha->F.empty() && !perturbForces.empty()) {
            std::cout << std::fixed << std::setprecision(10);
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
    }
  */  
//    this->step += 1;
}   




/*
void DynamicsOpenMM::applyPerturbForces(int forceFieldIndex, double scaleFactor) {
    static const int propagate_state = param.getInt("propagate_state");
    //ha->updateContextState();
    std::vector<OpenMM::Vec3> perturbForces;
    perturbForces.resize(DOFn, OpenMM::Vec3(0.0, 0.0, 0.0));
        if (ha->getForceFieldForces(perturbForces)) {
        std::cout << "Applying perturbation forces from ForceField " << forceFieldIndex << std::endl;


        if (!ha->F.empty() && !perturbForces.empty()) {
            std::cout << std::fixed << std::setprecision(10);
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

    ha->getForces();
    ha->getPositions();
    ha->getVelocities();
    std::cout << "First atom force to Integrate is : ("
                      << ha->F[0][0] << ", "
                      << ha->F[0][1] << ", "
                      << ha->F[0][2] << ") kJ/mol/nm" << std::endl;

    integrator->oneStep(ha->F);
//    ha->updateContextState();
//    ha->getPotentialEnergy(propagate_state, true);
//    ha->getKineticEnergy();
    ha->getForces();
    ha->getPositions();
//    std::vector<OpenMM::Vec3> perturbForces;
    ha->getVelocities();

//    perturbForces.resize(DOFn, OpenMM::Vec3(0.0, 0.0, 0.0));
    double kineticenergy = ha->getKineticEnergy();

    std::cout<<"ha->getForceFieldForces(perturbForces) begins"<<std::endl;
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "First atom before kinetic energy is: "<<kineticenergy<<std::endl;
    std::cout << "First atom velocity before interation is: ("
                      << ha->V[0][0] << ", "
                      << ha->V[0][1] << ", "
                      << ha->V[0][2] << ") nm/ps" << std::endl;
    if (ha->getForceFieldForces(perturbForces)) {
        std::cout << "Applying perturbation forces from ForceField " << forceFieldIndex << std::endl;


        if (!ha->F.empty() && !perturbForces.empty()) {
            std::cout << std::fixed << std::setprecision(10);
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
        

    ha ->uploadForces();
//    ha->updateContextState();
    std::cout << std::fixed << std::setprecision(10);
            std::cout << "upload atom force: ("
                      << ha->F[0][0] << ", "
                      << ha->F[0][1] << ", "
                      << ha->F[0][2] << ") kJ/mol/nm" << std::endl;

    kineticenergy = ha->getKineticEnergy();
    std::cout << "After added perturb force kinetic energy is: "<<kineticenergy<<std::endl;
    this->step += 1;

    ha->getForces();


//    ha->getVelocity();
//    ha->updateContextState();
    ha->getVelocities(); 
    std::cout << "First atom velocity end interation is: ("
                      << ha->V[0][0] << ", "
                      << ha->V[0][1] << ", "
                      << ha->V[0][2] << ") nm/ps" << std::endl;
    kineticenergy = ha->getKineticEnergy();
    std::cout << "First atom end kinetic energy is: "<<kineticenergy<<std::endl;
    std::cout << "get atom force: ("
                      << ha->F[0][0] << ", "
                      << ha->F[0][1] << ", "
                      << ha->F[0][2] << ") kJ/mol/nm" << std::endl;

} */ 



void DynamicsOpenMM::Pi_update() {
    // Empty implementation that just calculates Pi without storing it
    std::vector<double> temp_pi(9, 0.0);
    this->Pi_update(temp_pi);
}

void DynamicsOpenMM::Pi_update(std::vector<double>& Pi_store) {
    // Make sure Pi_store has the right size
    if (Pi_store.size() != 9) {
        Pi_store.resize(9, 0.0);
    }
    
    // Get position data from the current state
   // ha->getPositions();  // Make sure positions are up to date
    
    // Try to directly access a ForceFieldPolar from the polarForceFields vector
    if (!ha->polarForceFields.empty()) {
        ha->getForceFieldForces(Pi_store);
//	std::cout<<"Pi_store is  " << std::endl
//                 <<Pi_store[0] << "   " << Pi_store[1] << "   " << Pi_store[2]<<std::endl
//                 <<Pi_store[3] << "   " << Pi_store[4] << "   " << Pi_store[5]<<std::endl
//                 <<Pi_store[6] << "   " << Pi_store[7] << "   " << Pi_store[8]<<std::endl;

    } else {
        // If no ForceFieldPolar is available, throw an exception
        throw std::runtime_error("ERROR: Pi_update() requires a ForceFieldPolar, but none is available.");
    }
}

