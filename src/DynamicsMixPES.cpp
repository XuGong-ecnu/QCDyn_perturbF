/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 19, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsMixPES.h"

void DynamicsMixPES::init() {
    DynamicsBase::init();
    if (system_type != "allatom")
        throw std::runtime_error("ERROR: Unsupported system_type=" + system_type + " for OpenMM.");
    if (dyn_type != "MixPES")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for OpenMM.");
    // Create smarter pointer to OpenMM intergrator
    ha = std::static_pointer_cast<HamiltonianOpenMM>(Ha);
    integrator = ha->integrator;
    // Get time step size (in ps) from OpenMM integrator.
    DT = integrator->getStepSize();
    // * Bind function to update Hamiltonian in diabatic or adiabatic representation.
    // It is bound to the function in Hamiltonian object.
    if (representation == "diabatic")
        updateH = std::bind(&HamiltonianBase::updateDiabaticHamiltonian, Ha);
    else
        throw std::runtime_error("ERROR: Unsupported representation=" + representation +
            " for DynamicsMixPES.");
    // load the values from input control file 
    SplitString(PES_weights, param.getStr("PES_weights"));   
    if (PES_weights.size() != DOFe)
        throw std::runtime_error("ERROR: wrong number of PES_weights.");
    // Q: range of weights ??
    // A: we don't concern about the ratio here
}

void DynamicsMixPES::beforeOneTraj() {
    // Note the order of these functions is important.
    step = 0; // reset current step to 0.
    // 1. Initial nuclear DOF (R, P) sampling
    samplingNucl();
    // 2. Initialize the Hamiltonian and Forces at step 0.
    // For adiabatic represention, the nonadiabtic coupling (NAC) is also updated.
    // And the rotation matrix T (if any) used for transfromation between diabatic
    // and adiabtic representaions is also updated for model Hamiltonian.
    updateH();
    // 4. Compute effective forces of step 0 used for next propagtaion.
    updateF();
}

void DynamicsMixPES::dynamics(int steps) {
    for (int i = 0; i < steps; ++i)
        oneStep();
    this->step += steps;
}

void DynamicsMixPES::oneStep() {
    // This function is modified based on the new strategy on Dec. 3, 2021.
    // New strategy: the data of positions and velocities are stored in the
    // OpenMM Platform and the integration are performed by OpenMM Integrator.
    // But the forces used to do nuclear propagtaion is provided externally.
    // At each step, the force of each state should be downloaded from
    // OpenMM Platform such as CUDA, and weighted by the electronic population
    // and coherence to get the effective forces, then the effective forces will
    // be uploaded to the OpenMM Platform to do the nuclear propagtaion.
    // Link to the orignal data in classes of HamiltonianOpenMM and HamiltonianElec
    static std::shared_ptr<HamiltonianOpenMM> ha = std::static_pointer_cast<HamiltonianOpenMM>(Ha);
    static std::vector<OpenMM::Vec3>& F = ha->F; // effective forces used for nuclear propagation
    static const Real_Matrix&     H = ha->H; // diabatic Hamiltonian matrix
    static const std::string integrator = param.getStr("integrator"); // OpenMM integrator
    // The nuclear propagation is perfromed by velocity Verlet integrator.
    // Note that within velocity Verlet integrator, the electronic propagtion
    // is performed after the updating of nuclear position before the second
    // half step of updating velcities.
    if (integrator == "velocityVerlet") {
        // The velocity Verlet integrator is a special integrator, in which the
        // forces used in the updating of velocities at the first half step and
        // at the second half step are different. There is no original OpenMM
        // velocity Verlet Integrator. Here, using a CompoundIntegrator with two
        // CustomIntegrator to do velocity Verlet integrator with external forces.
        static std::shared_ptr<OpenMM::CompoundIntegrator> compound =
            std::static_pointer_cast<OpenMM::CompoundIntegrator>(ha->integrator);
        // At the start of each time step, we need to update OpenMM context state.
        // The thermostat will modify velocities; the barostat will modify positions;
        // the CMMotionRemover will remove motion of centor of mass in this function
        // if they are used in the simulation. It will return true if the positions
        // are modfied, and the potential energy and forces should be re-computed.
        // Otherwise, the previous effective forces can be used directly in the
        // following first half step propagation. At the initial step (step 0),
        // the effective forces have been computed in the beforeOneTraj().
        if (ha->updateContextState()) {
            // Update H and F_all which is defined in the HamiltonianOpenMM.
            updateH();
            // Update effective forces used for nuclear propagation
            updateF();
        }
        // Do nuclear propagation with velocity Verlet integrator
        // Part 1 of velocity Verlet integrator: the first CustomIntegrator
        // update V with first half step using current effective forces and update
        // the R with a full step and then do the distances constraints.
        compound->setCurrentIntegrator(0);
        compound->oneStep(F);
        // Update H and F_all based on the new positions.
        updateH();
        // Update effective forces with new q,p used for next nuclear propagation
        updateF();
        // Part 2 of velocity Verlet integrator: the second CustomIntegrator
        // update V with second half step with new effective forces: v = v+0.5*dt*f/m
        // and do the velocities constraints so the net velocity of each constrained
        // distance is zero.
        compound->setCurrentIntegrator(1);
        compound->oneStep(F);
        // Finally, we need to modify the time/step in OpenMM platform data since
        // in the above strategy, the simulation was propagated with two steps.
        ha->setTime(ha->getTime() - DT);
        ha->setStep(ha->getStep() - 1);
    }
    // Other OpenMM integrators such as leap-frog/Langevin/LangevinMiddle/NoseHoover
    // Within these integrators, the electronic propagation is performed after the
    // finish of nuclear propagtaion with a full step.
    else {
        // At the start of each time step, we need to update OpenMM context state.
        // The thermostat will modify velocities; the barostat will modify positions;
        // the CMMotionRemover will remove motion of centor of mass in this function
        // if they are used in the simulation. It will return true if the positions
        // are modfied, and the potential energy and forces should be re-computed.
        // Otherwise, the previous effective forces can be used directly in the
        // following first half step propagation. At the initial step (step 0),
        // the effective forces have been computed in the beforeOneTraj().
        if (ha->updateContextState()) {
            // Update H and F_all which is defined in the HamiltonianOpenMM.
            updateH();
            // Update effective forces used for nuclear propagation
            updateF();
        }
        // Do nuclear propagation within OpenMM integrator with the effective forces
        ha->integrator->oneStep(F);
        // Update H and F_all based on the new positions.
        updateH();
        // Update effective forces with new q,p used for next nuclear propagation
        updateF();
    }
}

void DynamicsMixPES::updateF() {
    // This is very simillar as the updateDiabaticModelForces(), but using
    // different data type for forces, here is OpenMM::Vec3. (modified on Dec. 3, 2021)
    // Link to the orignal data in classes of HamiltonianOpenMM and HamiltonianElec
    static std::shared_ptr<HamiltonianOpenMM> ha = std::static_pointer_cast<HamiltonianOpenMM>(Ha);
    static const std::vector<std::vector<OpenMM::Vec3>>& F_all = ha->F_all; // diabatic forces of each states
    static std::vector<OpenMM::Vec3>& F = ha->F; // effective forces to be updated
    // Reset all elements to 0 with keeping same size
    std::fill(F.begin(), F.end(), OpenMM::Vec3());
    for (int j = 0; j < DOFn; ++j) 
        for (int i = 0; i < DOFe; ++i)
            F[j] += F_all[i*DOFe+i][j] * PES_weights[i];
}