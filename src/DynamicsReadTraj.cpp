/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsReadTraj.h"
#include "chemfiles.hpp"

void DynamicsReadTraj::init() {
    DynamicsBase::init();
    if (system_type != "allatom")
        throw std::runtime_error("ERROR: Unsupported system_type=" + system_type + " for DynamicsReadTraj.");
    if (dyn_type != "ReadTraj")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for DynamicsReadTraj.");   
    DT = param.getDouble("DT");
}

void DynamicsReadTraj::beforeOneTraj() {
    step = 0; // reset current step to 0.
}

void DynamicsReadTraj::dynamics(int steps) {
    // Get smarter ponter to the HamiltonianreportMDData obejct
    static std::shared_ptr<HamiltonianForceFieldBase> ha = std::static_pointer_cast<HamiltonianForceFieldBase>(Ha);
    std::string loadfile;
    loadfile = param.getStr("nucl_load"); // read initial conditions from file
    // Open a file for reading, automatically guessing the file format from the extension.
    // It will throw an exception if the format is not supported.
    //static auto trajectory = chemfiles::Trajectory(loadfile);
    static chemfiles::Trajectory trajectory(loadfile);
    // Read positions, velocities, and box of current frame.
    static const auto nframes = trajectory.nsteps(); // number of frames
    for (int k = 0; k < steps; k++) {
        this->step += 1;
        const chemfiles::Frame frame = trajectory.read_step((unsigned long)step);
        if (frame.size() != (unsigned long)DOFn)
            throw std::runtime_error("ERROR: The number of atoms in trajectory is not consistent with the system.");
        const std::vector<chemfiles::Vector3D>& positions = frame.positions();
        if (positions.size() != (unsigned long)DOFn)
            throw std::runtime_error("ERROR: Wrong size of positions in trajectory.");
        //const std::vector<chemfiles::Vector3D>& velocities = *frame.velocities();
        //if (velocities.size() != (unsigned long)DOFn)
        //    throw std::runtime_error("ERROR: No velocities or wrong size of velocities in trajectory.");
        // Convert the chemfiles format to OpenMM format
        // Factor 10 because the lengths are in nm in the OpenMM, while
        // in Angstorm in chemfiles.
        std::vector<Vec3> R(DOFn, Vec3()); // positions
        for (int i = 0; i < DOFn; ++i)
            for (int j = 0; j < 3; ++j) {
                R[i][j] = positions[i][j]/10;//nm 单位 
            }
        ha->setPositions(R);
    }
}