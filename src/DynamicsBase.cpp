/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 13, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsBase.h"
#include "chemfiles.hpp"

void DynamicsBase::init() {
    // * Get parameters from Hamilatonian object, the legality of them should
    // be checked by Hamilatonian.
    system_type    = Ha->system_type;
    DOFn           = Ha->DOFn;
    DOFe           = Ha->DOFe;
    representation = Ha->representation;
    allatom_type = Ha->allatom_type;
    // * Read dynamics parameters and check the legality.
    dyn_type = param.getStr("dyn_type");
    PIMD_type = param.getStr("PIMD_type");
    if (PIMD_type == "RPMD") nbeads = param.getInt("nbeads");
    nucl_sample = param.getStr("nucl_sample");
    if (nucl_sample != "Wigner" && nucl_sample != "Classical" && nucl_sample != "none")
        throw std::runtime_error("ERROR: Unsupported nucl_sample=" + nucl_sample);
    temperature = param.getDouble("temperature");
    beta        = param.getDouble("beta");
    DT          = param.getDouble("DT");
    ntraj       = param.getInt("ntraj");
    nsteps      = param.getInt("nsteps");
    step        = 0;
    // If real units are provided, set beta (in au) according to real temperature
    // And convert unit of DT from ps to au used for propagation (model only).
    // If the unit is empty or "none", it is a unitless model.
    // Note: for all-atom simulation with OpenMM, the DT is ps and no convertion.
    const std::string unit = param.getStr("energy_unit");
    if (system_type == "model" && !unit.empty() && unit != "none") {
        // If temperature is zero, then set beta as zero (actually it should be infinity)
        beta = temperature == 0 ? 0 : (1.0 / (temperature * kT2au));
        DT = DT * ps2au;
    }
    if (temperature < 0 || beta < 0 || ntraj < 1 || nsteps < 0 || DT <= 0)
        throw std::runtime_error("ERROR: Illegal value for temperature or beta (requires > 0), "
            "or ntraj (requires >= 1), or nsteps (requires >= 0), or DT (requires > 0).");

    // * Set random number seed for nuclear sampling
    const int nucl_seed = param.getInt("nucl_seed");
    if (nucl_seed == 0) { // Use system-time seed for real simulation
        std::random_device rd;
        nucl_gen.seed(rd());
    }
    else // Use deterministic seed for debug
        nucl_gen.seed(nucl_seed);
}

void DynamicsBase::samplingNucl() {
    static int traj = 0; // it equals the index of current traj
    const std::string loadfile = param.getStr("nucl_load"); // read initial conditions from file
    if (system_type == "model") { // Works for models (SB, GOA, FMO, MSH)
        // The positions, velocities and model parameters can be accessed by HamiltonianModelBase
        std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
        // * Generate initial nuclear sampling according to model parameters
        if (loadfile.empty()) {
            // 1. Initialize sampling distribution width for nuclear DOF
            // sigma is standard deviation, which decides the sampling width
            std::vector<double> sigma_R(DOFn, 0);
            std::vector<double> sigma_V(DOFn, 0);
            if (nucl_sample == "Wigner") {
                // When beta is zero, which means to do initial sampling at
                // zero temeprature (actually beta is infinity, DBL_MAX in C++).
                // Basically, when beta is infinity, tanh(beta hbar omega/2)
                // goes to 1, which should be removed from the Wigner distribution
                // for initial nuclear sampling.
                if (beta == 0) // sample at zero temperature
                    for (int j = 0 ; j < DOFn; j++) {
                        sigma_R[j] = sqrt(hbar / (2 * ha->omega[j]));
                        sigma_V[j] = sigma_R[j] * ha->omega[j];
                    }
                else // sample at finite temperature
                    for (int j = 0 ; j < DOFn; j++) {
                        sigma_R[j] = sqrt(hbar / (2 * ha->omega[j] * tanh(0.5 * beta * hbar * ha->omega[j])));
                        sigma_V[j] = sigma_R[j] * ha->omega[j];
                    }
            }
            else if (nucl_sample == "classical") {
                if (beta == 0)
                    throw std::runtime_error("ERROR: Classical sampling can not be performed at zero temperature.");
                for (int j = 0 ; j < DOFn ; j++) {
                    sigma_R[j] = 1.0 / (ha->omega[j] * sqrt(beta));
                    sigma_V[j] = 1.0 / sqrt(beta);
                }
            }
            else if (nucl_sample == "NCMorse") {
                if (beta == 0) {
                     // see Miller's 2001 papers on 3-state (TBCDI)
                    double omega = 0.005;
                    double M     = 20000;
                    for (int j = 0; j < DOFn; j++) {
                        sigma_R[j] = sqrt(hbar / (2 * omega * M));
                        sigma_V[j] = sqrt(hbar * omega / (2 * M));
                    }
                }
                else 
                    throw std::runtime_error("ERROR: Unsupported nucl_sample=" + nucl_sample 
                        + " for model with finite temperature.");
            }
            else
                throw std::runtime_error("ERROR: Unsupported nucl_sample=" + nucl_sample + " for model.");
            // 2. Generate gaussian random initial conditions (R and V)
            // The noeq or eq shift of initial state has been setted in shift
            // accoring to specfic models. So we can use a same form for R here.
            std::normal_distribution<double> normal_dist(0.0, 1.0);
            for (int j = 0 ; j < DOFn ; j++) {
                ha->R[j] = normal_dist(nucl_gen) * sigma_R[j] - ha->shift[j];
                // Velocities (or momenta) are same for all models
                ha->V[j] = normal_dist(nucl_gen) * sigma_V[j];
            }
            // Note, for the QVC model, if sample on the acceptor state, the initial
            // R/V should be converted to the donor state by Duschinsky matrix
            // and shift vector. [added on Nov. 30, 2021]
            if (ha->model_type == "QVC" && param.getInt("sample_state") == 1) {
                std::shared_ptr<HamiltonianQVC> ha_QVC = std::static_pointer_cast<HamiltonianQVC>(Ha);
                std::vector<double> R_new, V_new; // shift for velocities is 0
                ha_QVC->transformNormalModes(ha_QVC->R, ha_QVC->J, ha_QVC->req, R_new, true);
                ha_QVC->transformNormalModes(ha_QVC->V, ha_QVC->J, std::vector<double>(DOFn, 0), V_new, true);
                ha->R = R_new;
                ha->V = V_new;
            }
        }
        // * Load initial nuclear sampling from file (used for debug)
        else {
            std::vector<std::vector<double>*> data = {&(ha->R), &(ha->V)};
            LoadDataFile(loadfile, data, "", 2, 1, traj);
        }
        // * Save initial nuclear sampling to file if requested (model only)
        const std::string file = param.getStr("nucl_save");
        if (!file.empty()) {
           std::vector<std::string> headers = {"R", "V"};
           std::vector<std::vector<double>*> data = {&(ha->R), &(ha->V)};
           // if traj = 0, create a new file, else, append to this file.
           SaveDataFile(file, headers, data, "", true, "traj=", traj);
        }
    }
    // TODO: On-the-fly electronic structure calculation
    else if (system_type == "onthefly") {

    }
    // To set the Intial velocities for every atom
    else if (system_type == "allatom" && allatom_type == "CustomForceField") {
        static std::shared_ptr<HamiltonianForceFieldBase> ha = std::static_pointer_cast<HamiltonianForceFieldBase>(Ha);
        if (loadfile.empty()) {
             if (nucl_sample == "classical")
                 throw std::runtime_error("ERROR: For all-atom simulation, the classical sampling "
                     "is performed by loading initial positions and velocities from trajectory file.\n"
                     "Please specify nucl_load as a trajectory file to do this.");
             else if (nucl_sample == "Wigner")
                 throw std::runtime_error("ERROR: For all-atom simulation, the Wigner sampling is not supported yet.");
             // nucl_sample=none means to run one traj from specific configuration
             // that loaded from structure file.
             else if (nucl_sample != "none")
                 throw std::runtime_error("ERROR: Unsupported nucl_sample=" + nucl_sample + " for an all-atom simulation.");
        }
         // load initial conditions from trajectory file
         // Here, I use the chemfiles library to read trajectory. It supports
         // both text (xyz, gro, pdb, ...) and binary (nc, xtc, trr, ...) files
         // See: https://github.com/chemfiles/chemfiles
         // Documentation: http://chemfiles.org/chemfiles/latest/cxx-tutorials.html
         // Note, the unit in chemfiles for distance is Angstrom, and for tiem is ps.
        else {
            // Open a file for reading, automatically guessing the file format from the extension.
            // It will throw an exception if the format is not supported.
            //static auto trajectory = chemfiles::Trajectory(loadfile);
            static chemfiles::Trajectory trajectory(loadfile);
            static const auto nframes = trajectory.nsteps(); // number of frames
            // Get parameters to control the reading of trajectory.
            // nucl_start: the first frame to read, 0 <= nucl_start <= total frames
            // nucl_skip: the interval to read, 0 < nucl_skip <= nucl_end,
            // and (nucl_end-nucl_start) must be divided exactly by nucl_skip,
            // 1 means to load each frame unitil the nucl_end.
            // nucl_end: the end frame to read, 0 <= nucl_end <= total frames,
            // -1 means the last frame in the trajectory.
            // The number of frames to be loaded is (nucl_end-nucl_start)/nucl_skip+1
            // which should be not less than the number of trajectories (ntraj) to run.
            // nucl_current: current frame to read, should not greater than nucl_end.
            // Note, they the zero-based index of frame, not the real simulation
            // step recorded in the trajectory.
            static const int nucl_start = param.getInt("nucl_start");
            static const int nucl_skip  = param.getInt("nucl_skip");
            static const int nucl_end   = param.getInt("nucl_end") == -1 ? (nframes-1) : param.getInt("nucl_end");
            static int nucl_current = nucl_start;
            if (traj == 0) { // check the validity of them at the first time only
                if (nframes == 0)
                    throw std::runtime_error("ERROR: No frame in trajectory file: " + loadfile);
                if (nucl_start < 0 || (unsigned long)nucl_start > (nframes-1))
                    throw std::runtime_error("ERROR: nucl_start =" + std::to_string(nucl_start) + " is out of scope.");
                if (nucl_end < 0 || (unsigned long)nucl_end > (nframes-1))
                    throw std::runtime_error("ERROR: nucl_end =" + std::to_string(nucl_end) + " is out of scope.");
                if (nucl_start > nucl_end)
                    throw std::runtime_error("ERROR: nucl_start is greater than nucl_end.");
                if (nucl_skip <= 0 || (nucl_end - nucl_start) % nucl_skip != 0)
                    throw std::runtime_error("ERROR: Illegal nucl_skip = " + std::to_string(nucl_skip));
                if ((nucl_end - nucl_start) / nucl_skip + 1 < ntraj)
                    throw std::runtime_error("ERROR: The number of frames loading from trajectory is less than the ntraj.");
            }
            // * Read the frame from trajectory as initial conditions
            // Read positions, velocities, and box of current frame.
            const chemfiles::Frame frame = trajectory.read_step((unsigned long)nucl_current);
            if (frame.size() != (unsigned long)DOFn)
                throw std::runtime_error("ERROR: The number of atoms in trajectory is not consistent with the system.");
            const std::vector<chemfiles::Vector3D>& positions = frame.positions();
            if (positions.size() != (unsigned long)DOFn)
                throw std::runtime_error("ERROR: Wrong size of positions in trajectory.");
            const std::vector<chemfiles::Vector3D>& velocities = *frame.velocities();
            if (velocities.size() != (unsigned long)DOFn)
                throw std::runtime_error("ERROR: No velocities or wrong size of velocities in trajectory.");
            // Convert the chemfiles format to OpenMM format
            // Factor 10 because the lengths are in nm in the OpenMM, while
            // in Angstorm in chemfiles.
            std::vector<Vec3> R(DOFn, Vec3()); // positions
            std::vector<Vec3> V(DOFn, Vec3()); // velocities
            for (int i = 0; i < DOFn; ++i)
                for (int j = 0; j < 3; ++j) {
                    R[i][j] = positions[i][j] / 10.0;
                    V[i][j] = velocities[i][j] / 10.0;
                }
            // Set the positions and velocities in OpenMM Context
            // Remember to re-initialize the OpenMM Context object before the
            // set positions and velocities, which makes sure the simulation is
            // a brand new simulation for a new trajectory. The reinitialize()
            // will set the simulation time and step to zero, and reset the
            // random number seed. But I found this function takes too much time
            // (sevral seconds).
            // For NVE simulation, the reinitialize() is not necessary, and we only
            // need to set the simulation time and step to zero manully.
            // While for NVT or NPT simulation, it is necessary due to the using
            // of thermostat and barostat with random number generator.
            // ? However, in practice, I found that for one specific initial sampling
            // ? loading from one trajectory, the results in multi-traj simulation
            // ? and in one-single traj simulation are not totally same when the
            // ? thermostat and/or barostat is used, althogh after reinitialize()
            // ? the random number seed is re-initialized and same if it is specfied
            // ? by user, but in my opnion, the difference from random number of
            // ? thermostat/barostat is not problem.
            //static const bool useThermostat = param.getStr("thermostat") != "none" ? true : false;
            //static const bool useBarostat = param.getStr("barostat") != "none" ? true : false;
            //if (traj > 0) // do this for a new trajectory
            //    if (useThermostat || useBarostat) {
            //        ha->reinitialize();
            //    }
            //    else {
            //        ha->setTime(0);
            //        ha->setStep(0);
            //    }
            ha->setPositions(R);
            ha->setVelocities(V);
            // If the PBC is used, we need to set the box too, if the
            // trajectory are sampled from NPT simulation, the box changes.
            //if (ha->usesPeriodicBoundaryConditions()) {
                const chemfiles::UnitCell& box = frame.cell();
                if (box.shape() == chemfiles::UnitCell::CellShape::INFINITE) // no cell
                    throw std::runtime_error("ERROR: PBC is used, while there is no box in trajectory.");
                // Factor 10 because the lengths are in nm in the OpenMM, while
                // in Angstorm in chemfiles. The box in chemfiles is
                // a_x,b_x,c_x;a_y,b_y,c_y;a_z,b_z,c_z, while in OpenMM, it
                // should be a_x,a_y,a_z;b_x,b_y,b_z;c_x,c_y,c_z.
                chemfiles::Matrix3D matrix = box.matrix() / 10.0;
                Vec3 a(matrix[0][0], matrix[1][0], matrix[2][0]);
                Vec3 b(matrix[0][1], matrix[1][1], matrix[2][1]);
                Vec3 c(matrix[0][2], matrix[1][2], matrix[2][2]);
                ha->setPeriodicBoxVectors(a, b, c);
            //}
            // Set the next index of frame to read
            nucl_current += nucl_skip;
        }

        //static std::shared_ptr<HamiltonianForceFieldBase> ha = std::static_pointer_cast<HamiltonianForceFieldBase>(Ha);
        //if (temperature != 0)
        //    throw std::runtime_error("ERROR: Classical sampling can not be performed at zero temperature.");
        //std::vector<double> Sum_V(3, 0);
        //std::vector<double> sigma_V(DOFn, 0);
        //std::normal_distribution<double> normal_dist(0.0, 1.0);
        //double Ke = 0;
        //double Kt = 0;
        //for (int j = 0 ; j < DOFn ; j++) {
        //    sigma_V[j] = sqrt( Rbar * temperature / ha->masses[j]);
        //    ha->V[j][0] = normal_dist(nucl_gen) * sigma_V[j];
        //    ha->V[j][1] = normal_dist(nucl_gen) * sigma_V[j];
        //    ha->V[j][2] = normal_dist(nucl_gen) * sigma_V[j];
        //    Sum_V[0] += ha->V[j][0];
        //    Sum_V[1] += ha->V[j][1];
        //    Sum_V[2] += ha->V[j][2];
        //}
        //Sum_V[0] /= DOFn;
        //Sum_V[1] /= DOFn;
        //Sum_V[2] /= DOFn;
        //for (int j = 0 ; j < DOFn ; j++) {
        //    ha->V[j][0] -= Sum_V[0];
        //    ha->V[j][1] -= Sum_V[1];
        //    ha->V[j][2] -= Sum_V[2];
        //}   
        //for (int j = 0; j < DOFn; j++) {
        //    Ke += 0.5 * ha->masses[j] * ((ha->V[j][0])*(ha->V[j][0]) + (ha->V[j][1])*(ha->V[j][1])  +(ha->V[j][2])*(ha->V[j][2]));
        //}
        //Kt = sqrt(0.5 * Rbar *ha->getSystemDOF() * temperature/Ke);
        //for (int j = 0 ; j < DOFn ; j++) {
        //    ha->V[j][0] *= Kt;
        //    ha->V[j][1] *= Kt;
        //    ha->V[j][2] *= Kt;
        //} 
    }
    else if (system_type == "allatom" && allatom_type == "OpenMM") {
        static std::shared_ptr<HamiltonianOpenMM> ha = std::static_pointer_cast<HamiltonianOpenMM>(Ha);
        if (loadfile.empty()) {
            if (nucl_sample == "classical")
                throw std::runtime_error("ERROR: For all-atom simulation, the classical sampling "
                    "is performed by loading initial positions and velocities from trajectory file.\n"
                    "Please specify nucl_load as a trajectory file to do this.");
            else if (nucl_sample == "Wigner")
                throw std::runtime_error("ERROR: For all-atom simulation, the Wigner sampling is not supported yet.");
            // nucl_sample=none means to run one traj from specific configuration
            // that loaded from structure file.
            else if (nucl_sample != "none")
                throw std::runtime_error("ERROR: Unsupported nucl_sample=" + nucl_sample + " for an all-atom simulation.");
        }
        // load initial conditions from trajectory file
        // Here, I use the chemfiles library to read trajectory. It supports
        // both text (xyz, gro, pdb, ...) and binary (nc, xtc, trr, ...) files
        // See: https://github.com/chemfiles/chemfiles
        // Documentation: http://chemfiles.org/chemfiles/latest/cxx-tutorials.html
        // Note, the unit in chemfiles for distance is Angstrom, and for tiem is ps.
        else {
            // Open a file for reading, automatically guessing the file format from the extension.
            // It will throw an exception if the format is not supported.
            //static auto trajectory = chemfiles::Trajectory(loadfile);
            static chemfiles::Trajectory trajectory(loadfile);
            static const auto nframes = trajectory.nsteps(); // number of frames
            // Get parameters to control the reading of trajectory.
            // nucl_start: the first frame to read, 0 <= nucl_start <= total frames
            // nucl_skip: the interval to read, 0 < nucl_skip <= nucl_end,
            // and (nucl_end-nucl_start) must be divided exactly by nucl_skip,
            // 1 means to load each frame unitil the nucl_end.
            // nucl_end: the end frame to read, 0 <= nucl_end <= total frames,
            // -1 means the last frame in the trajectory.
            // The number of frames to be loaded is (nucl_end-nucl_start)/nucl_skip+1
            // which should be not less than the number of trajectories (ntraj) to run.
            // nucl_current: current frame to read, should not greater than nucl_end.
            // Note, they the zero-based index of frame, not the real simulation
            // step recorded in the trajectory.
            static const int nucl_start = param.getInt("nucl_start");
            static const int nucl_skip  = param.getInt("nucl_skip");
            static const int nucl_end   = param.getInt("nucl_end") == -1 ? (nframes-1) : param.getInt("nucl_end");
            static int nucl_current = nucl_start;
            if (traj == 0) { // check the validity of them at the first time only
                if (nframes == 0)
                    throw std::runtime_error("ERROR: No frame in trajectory file: " + loadfile);
                if (nucl_start < 0 || (unsigned long)nucl_start > (nframes-1))
                    throw std::runtime_error("ERROR: nucl_start =" + std::to_string(nucl_start) + " is out of scope.");
                if (nucl_end < 0 || (unsigned long)nucl_end > (nframes-1))
                    throw std::runtime_error("ERROR: nucl_end =" + std::to_string(nucl_end) + " is out of scope.");
                if (nucl_start > nucl_end)
                    throw std::runtime_error("ERROR: nucl_start is greater than nucl_end.");
                if (nucl_skip <= 0 || (nucl_end - nucl_start) % nucl_skip != 0)
                    throw std::runtime_error("ERROR: Illegal nucl_skip = " + std::to_string(nucl_skip));
                if ((nucl_end - nucl_start) / nucl_skip + 1 < ntraj)
                    throw std::runtime_error("ERROR: The number of frames loading from trajectory is less than the ntraj.");
            }
            // * Read the frame from trajectory as initial conditions
            // Read positions, velocities, and box of current frame.
            const chemfiles::Frame frame = trajectory.read_step((unsigned long)nucl_current);
            if (frame.size() != (unsigned long)DOFn)
                throw std::runtime_error("ERROR: The number of atoms in trajectory is not consistent with the system.");
            const std::vector<chemfiles::Vector3D>& positions = frame.positions();
            if (positions.size() != (unsigned long)DOFn)
                throw std::runtime_error("ERROR: Wrong size of positions in trajectory.");
            const std::vector<chemfiles::Vector3D>& velocities = *frame.velocities();
            if (velocities.size() != (unsigned long)DOFn)
                throw std::runtime_error("ERROR: No velocities or wrong size of velocities in trajectory.");
            // Convert the chemfiles format to OpenMM format
            // Factor 10 because the lengths are in nm in the OpenMM, while
            // in Angstorm in chemfiles.
            std::vector<OpenMM::Vec3> R(DOFn, OpenMM::Vec3()); // positions
            std::vector<OpenMM::Vec3> V(DOFn, OpenMM::Vec3()); // velocities
            for (int i = 0; i < DOFn; ++i)
                for (int j = 0; j < 3; ++j) {
                    R[i][j] = positions[i][j] / 10.0;
                    V[i][j] = velocities[i][j] / 10.0;
                }
            // Set the positions and velocities in OpenMM Context
            // Remember to re-initialize the OpenMM Context object before the
            // set positions and velocities, which makes sure the simulation is
            // a brand new simulation for a new trajectory. The reinitialize()
            // will set the simulation time and step to zero, and reset the
            // random number seed. But I found this function takes too much time
            // (sevral seconds).
            // For NVE simulation, the reinitialize() is not necessary, and we only
            // need to set the simulation time and step to zero manully.
            // While for NVT or NPT simulation, it is necessary due to the using
            // of thermostat and barostat with random number generator.
            // ? However, in practice, I found that for one specific initial sampling
            // ? loading from one trajectory, the results in multi-traj simulation
            // ? and in one-single traj simulation are not totally same when the
            // ? thermostat and/or barostat is used, althogh after reinitialize()
            // ? the random number seed is re-initialized and same if it is specfied
            // ? by user, but in my opnion, the difference from random number of
            // ? thermostat/barostat is not problem.
            static const bool useThermostat = param.getStr("thermostat") != "none" ? true : false;
            static const bool useBarostat = param.getStr("barostat") != "none" ? true : false;
            if (traj > 0) // do this for a new trajectory
                if (useThermostat || useBarostat) {
                    ha->reinitialize();
                }
                else {
                    ha->setTime(0);
                    ha->setStep(0);
                }
            ha->setPositions(R);
            ha->setVelocities(V);
            // If the PBC is used, we need to set the box too, if the
            // trajectory are sampled from NPT simulation, the box changes.
            if (ha->usesPeriodicBoundaryConditions()) {
                const chemfiles::UnitCell& box = frame.cell();
                if (box.shape() == chemfiles::UnitCell::CellShape::INFINITE) // no cell
                    throw std::runtime_error("ERROR: PBC is used, while there is no box in trajectory.");
                // Factor 10 because the lengths are in nm in the OpenMM, while
                // in Angstorm in chemfiles. The box in chemfiles is
                // a_x,b_x,c_x;a_y,b_y,c_y;a_z,b_z,c_z, while in OpenMM, it
                // should be a_x,a_y,a_z;b_x,b_y,b_z;c_x,c_y,c_z.
                chemfiles::Matrix3D matrix = box.matrix() / 10.0;
                OpenMM::Vec3 a(matrix[0][0], matrix[1][0], matrix[2][0]);
                OpenMM::Vec3 b(matrix[0][1], matrix[1][1], matrix[2][1]);
                OpenMM::Vec3 c(matrix[0][2], matrix[1][2], matrix[2][2]);
                ha->setPeriodicBoxVectors(a, b, c);
            }
            // Set the next index of frame to read
            nucl_current += nucl_skip;
        }
    }
    traj++;
}

void DynamicsBase::MOVE_Half_V(const std::vector<double>& masses, const std::vector<OpenMM::Vec3>& F, std::vector<OpenMM::Vec3>& V) {
    for (int j = 0; j < DOFn; j++)
       if (masses[j] != 0.0)
            for (int k = 0; k < 3; k++)
                V[j][k] += 0.5 * DT * F[j][k] / masses[j];
}

void DynamicsBase::MOVE_Half_V(const std::vector<double>& masses, const std::vector<Vec3>& F, std::vector<Vec3>& V) {
    for (int j = 0; j < DOFn; j++)
       if (masses[j] != 0.0)
            for (int k = 0; k < 3; k++)
                V[j][k] += 0.5 * DT * F[j][k] / masses[j];
}

void DynamicsBase::MOVE_V(const std::vector<double>& masses, const std::vector<OpenMM::Vec3>& F, std::vector<OpenMM::Vec3>& V) {
    for (int j = 0; j < DOFn; j++)
       if (masses[j] != 0.0)
            for (int k = 0; k < 3; k++)
                V[j][k] += DT * F[j][k] / masses[j];
}

void DynamicsBase::MOVE_V(const std::vector<double>& masses, const std::vector<Vec3>& F, std::vector<Vec3>& V) {
    std::cout<<" MOVE_V 000"<<std::endl;
    for (int j = 0; j < DOFn; j++)
       if (masses[j] != 0.0)
            for (int k = 0; k < 3; k++)
                V[j][k] += DT * F[j][k] / masses[j];
    std::cout<<" MOVE_V 111"<<std::endl;
}

void DynamicsBase::MOVE_R(const std::vector<double>& masses, const std::vector<OpenMM::Vec3>& V, std::vector<OpenMM::Vec3>& R) {
    for (int j = 0; j < DOFn; j++)
       if (masses[j] != 0.0)
            for (int k = 0; k < 3; k++)
                R[j][k] += DT * V[j][k];
}

void DynamicsBase::MOVE_R(const std::vector<double>& masses, const std::vector<Vec3>& V, std::vector<Vec3>& R) {
    
    for (int j = 0; j < DOFn; j++)
       if (masses[j] != 0.0)
            for (int k = 0; k < 3; k++)
                R[j][k] += DT * V[j][k];  
}

void DynamicsBase::MOVE_Half_V(const std::vector<double>& F, std::vector<double>& V) {
    for (int j = 0; j < DOFn; j++)
        V[j] += 0.5 * DT * F[j];
}

void DynamicsBase::MOVE_V(const std::vector<double>& F, std::vector<double>& V) {
    for (int j = 0; j < DOFn; j++)
        V[j] += DT * F[j];
}

void DynamicsBase::MOVE_R(const std::vector<double>& V, std::vector<double>& R) {
    for (int j = 0; j < DOFn; j++)
        R[j] += DT * V[j];
}

void DynamicsBase::RP_MOVE_Half_V(const std::vector<double>& masses, std::vector<std::vector<Vec3>>& F, std::vector<std::vector<Vec3>>& V) {
    for (int i = 0; i < nbeads; i++)
        for (int j = 0; j < DOFn; j++) 
            if (masses[j] != 0.0)
                for (int k = 0; k < 3; k++) 
                    V[i][j][k] += 0.5 * DT * F[i][j][k] / masses[j];
}

void DynamicsBase::RP_MOVE_V(const std::vector<double>& masses, std::vector<std::vector<Vec3>>& F, std::vector<std::vector<Vec3>>& V) {
    for (int i = 0; i < nbeads; i++)
        for (int j = 0; j < DOFn; j++)
            if (masses[j] != 0.0)
                for (int k = 0; k < 3; k++)
                    V[i][j][k] += DT * F[i][j][k] / masses[j];
}