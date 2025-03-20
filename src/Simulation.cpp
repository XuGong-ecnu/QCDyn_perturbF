/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 19, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "Simulation.h"
#include "chemfiles.hpp"

void Simulation::init(int argc, char *argv[]) {
    // * Load parameters from input control file and command line.
    InitializeParameters(argc, argv, param);
    // * Get job and system type from gloable parameters.
    job_type        = param.getStr("job_type");
    dyn_type        = param.getStr("dyn_type");
    prep_type       = param.getStr("prep_type");
    post_type       = param.getStr("post_type");
    system_type     = param.getStr("system_type");
    model_type      = param.getStr("model_type");
    onthefly_type   = param.getStr("onthefly_type");
    allatom_type    = param.getStr("allatom_type");
    PIMD_type       = param.getStr("PIMD_type");
    RDM_steps       = param.getInt("RDM_steps"); // output frequency
    max_friction_fail_traj = param.getDouble("max_friction_fail_traj"); 
    min_judge_traj_avail   = param.getInt("min_judge_traj_avail");
    polar_obs       = param.getBool("polar_obs");
    // * Create the Hamiltonian, Dynamics objects accordingly
    // The first level is job_type since it will decide whether a Hamiltonian
    // or Dynamics objects is required or not.
    if (job_type == "dynamics") { // dynamics simulation
        std::string name; // name of model or allatom or onthefly interface
        // * 1 Create Hamiltonian objetc and initialize it firtly.
        if (system_type == "model") { // model Hamiltonians
            // 2-level Spin-Boson model
            if (model_type == "SB" || model_type == "SB_req" || model_type == "GOA")
                Ha = std::make_shared<HamiltonianSpinBoson>(param);
            // Multi-State Global Bath Harmonic model
            else if (model_type == "MSH")
                Ha = std::make_shared<HamiltonianMSH>(param);
            // Multi-State FMO model, it also can be standard FMO7, FMO4, FMO5.
            else if (model_type.substr(0, 3) == "FMO")
                Ha = std::make_shared<HamiltonianFMO>(param);
            // linear vibronic coupling (LVC) model (non-Condon)
            else if (model_type == "LVC")
                Ha = std::make_shared<HamiltonianLVC>(param);
            // quadratic vibronic coupling (QVC) model (non-Condon)
            else if (model_type == "QVC")
                Ha = std::make_shared<HamiltonianQVC>(param);
            // non-Condon Morse model (3-state photodissociation model)
            else if (model_type == "NCMorse")
                Ha = std::make_shared<HamiltonianNCMorse>(param);
            else
                throw std::runtime_error("ERROR: Unsupported model_type=" + model_type);
            name = model_type;
        }
        else if (system_type == "onthefly") { // TODO: onthefly noadiabatic dynamics
            throw std::runtime_error("ERROR: Unsupported onthefly_type=" + onthefly_type);
            name = onthefly_type;
        }
        else if (system_type == "allatom") { // allatom simulation
            if (allatom_type == "OpenMM") { // internal interface with OpenMM
                if (PIMD_type == "RPMD") {
                    Ha = std::make_shared<HamiltonianRingPolymer>(param);
                }
                else 
                    Ha = std::make_shared<HamiltonianOpenMM>(param);
            }
            else if (allatom_type == "CustomForceField")
                if (PIMD_type == "RPMD") {
                    Ha = std::make_shared<HamiltonianRingPolymer>(param);
                }
                else 
                    Ha = std::make_shared<HamiltonianForceFieldBase>(param);
            else
                throw std::runtime_error("ERROR: Unsupported allatom_type=" + allatom_type);
            name = allatom_type;
        }
        Ha->init(); // initialize the Hamiltonian object
        // * 2 Then create Dynamics object and initialize it.
        // (1) classical MD with OpenMM
        if (dyn_type == "OpenMM") 
            if (PIMD_type == "RPMD") {
                // This is to Do
                //Dy = std::make_shared<DynamicsRPOpenMM>(param, Ha);
            }
            else 
                Dy = std::make_shared<DynamicsOpenMM>(param, Ha);
        // (2) class force field molecule dynamics 
        else if (dyn_type == "MD")
            if (PIMD_type == "RPMD") {
                Dy = std::make_shared<DynamicsRPMD>(param, Ha);
            }
            else
                Dy = std::make_shared<DynamicsMD>(param, Ha);
        // (3) Read Trr file 
        else if (dyn_type == "ReadTraj")
            Dy = std::make_shared<DynamicsReadTraj>(param, Ha); 
        // (4) Ehrenfest mean-filed dynamics (RDM version)
        else if (dyn_type == "MF-RDM")
            Dy = std::make_shared<DynamicsMFRDM>(param, Ha);
        // (5) Ehrenfest mean-filed dynamics (in mapping base)
        else if (dyn_type == "MF")
            Dy = std::make_shared<DynamicsMF>(param, Ha);
        // (6) fewest switches surface-hopping (FSSH) dynamics
        else if (dyn_type == "FSSH")
            Dy = std::make_shared<DynamicsFSSH>(param, Ha);
        // (7) linearized semiclassical (LSC) mapping dynamics
        else if (dyn_type == "LSC")
            Dy = std::make_shared<DynamicsLSC>(param, Ha);
        // (8) symmetrical quasi-classical (SQC) mapping dynamics
        else if (dyn_type == "SQC")
            Dy = std::make_shared<DynamicsSQC>(param, Ha);
        // (9) spin-mapping (SPM) dynamics
        else if (dyn_type == "SPM")
            Dy = std::make_shared<DynamicsSPM>(param, Ha);
        // (10) classical mapping model dynamics (CMM1 and CMM3-6)
        else if (dyn_type.substr(0,3) == "CMM" && dyn_type[3] != '2')
            Dy = std::make_shared<DynamicsCMM>(param, Ha);
        // (11) extend classical mapping model dynamics (eCMM and CMM2)
        else if (dyn_type == "CMM2" || dyn_type == "eCMM")
            Dy = std::make_shared<DynamicsECMM>(param, Ha);
        // (12) Trotter-Based Surface-Hopping (TBSH-MQCL) dynamics
        else if (dyn_type == "TBSH")
            Dy = std::make_shared<DynamicsTBSH>(param, Ha);
        // (13)  extend classical mapping model dynamics with commutator var (eCMMcv)
        else if (dyn_type == "eCMMcv") 
            Dy = std::make_shared<DynamicsECMMCV>(param,Ha);    
        else
            throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type);
        
        Dy->init(); // initialize the Dynamics object
        std::cout << "The " << dyn_type << " dynamics for " << system_type <<
            " (" << name << ") will be performed.\n" << std::endl;
        // * 3 Finally create Observables object using Ha and Dy and call init()
        Obs = std::make_shared<Observables>(param, Ha, Dy);
        Obs->init();
    }
    else if (job_type == "preprocess") {
        if (prep_type == "model_gen") { // generate model parameters
            std::string model_load = param.getStr("model_load");
            if (!model_load.empty())
                throw std::runtime_error("ERROR: When prep_type=" + prep_type +
                ", it is not resonable to load model parameters from file spcified by model_load.");
            std::string model_save = param.getStr("model_save");
            if (model_save.empty())
                throw std::runtime_error("ERROR: When prep_type=" + prep_type +
                ", the filename used to save model parameters should be provided by model_save.");
            if (system_type == "model")  // model Hamiltonians
                // 2-level Spin-Boson model
                if (model_type == "SB" || model_type == "SB_req" || model_type == "GOA")
                    Ha = std::make_shared<HamiltonianSpinBoson>(param);
                // Multi-State Global Bath Harmonic model
                else if (model_type == "MSH")
                    Ha = std::make_shared<HamiltonianMSH>(param);
                // Multi-State FMO model, typical is 7-state FMO.
                else if (model_type.substr(0, 3) == "FMO")
                    Ha = std::make_shared<HamiltonianFMO>(param);
                else
                    throw std::runtime_error("ERROR: Unsupported model_type=" + model_type);
            else
                throw std::runtime_error("ERROR: Unsupported system_type=" + system_type
                    + "when prep_type=" + prep_type);
            Ha->init(); // Generate model parameters at the initilization
            std::cout << "Generated " << model_type << " model parameters (in a.u.) file: " <<
                model_save << ".\n" << std::endl;
        }
        else
            throw std::runtime_error("ERROR: Unsupported prep_type=" + prep_type);
    }
    else if (job_type == "postprocess") {
        // Transfer tensor method (TTM), which is used to predict the long time
        // reduced denisty matrix (RDM) from short time RDM results.
        // In this case, Hamiltonian is not necessary. So It is more like a
        // post-processing method.
        if (post_type == "TTM") { // TODO: perhap move this part to run()
            std::cout << "Start to run transfer tensor method (TTM) since " << CurrentTime() << ".\n" << std::endl;
            DynamicsTTM ttm(param);
            ttm.init();
            ttm.buildDynamicalMaps();
            ttm.buildTransferTensors();
            ttm.propagate();
            std::cout << "The post-processing job with TTM is finished at " << CurrentTime() << ".\n" << std::endl;
            const std::vector<Complex_Matrix>& RDM = ttm.getDensityMatrix();
            // TODO: use the general one
            // Output RMD results.
            const int DOFe = param.getInt("DOFe");
            std::string filename = param.getStr("RDM_save");
            // Get initial state with two indices j,k
            std::vector<std::string> init_state;
            SplitString(init_state, param.getStr("init_state"));
            if (filename.empty()) { // get the default name
                filename = "RDM_TTM";
                if (init_state[0] == "all") // the index of initial state will be added later
                    filename = filename + ".csv";
                else
                    filename = filename + "_" + init_state[0] + init_state[1] + ".csv";
            }
            // Make sure the format of file is csv file.
            else if (GetFileSuffix(filename) != "csv") {
                std::cout << "WARNING: The extension name of RDM_save is not csv, but csv file will be generated." << std::endl;
                filename = GetFilePrefix(filename) + ".csv";
            }
            const std::string prefix = GetFilePrefix(filename); // final prefix of filename
            for (int n = 0; n < RDM.size(); ++n) {
                // if init_state = all, the number of RDMs must be larger than 1
                if (RDM.size() > 1)
                    filename = prefix + "_" + std::to_string(n/DOFe) + std::to_string(n%DOFe) + ".csv";
                writeTTMFile(RDM[n], filename);
                // Obs->output();
            }
            std::cout << std::endl;
        }
        else
            throw std::runtime_error("ERROR: Unsupported post_type=" + post_type);
    }
    else
        throw std::runtime_error("ERROR: Unsupported job_type=" + job_type);
}

void Simulation::run() {
    if (job_type != "dynamics") return; // only for dynamics jobs

    // * 1. Start simulation and record start time.
    std::cout << "Run dynamics simulation since " << CurrentTime() << "." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    // * 2. Get global common parameters that are used to control the simulation
    // The parameters has been verified in dynamics object.
    const int ntraj  = param.getInt("ntraj");  // number of trajectories
    const int nsteps = param.getInt("nsteps"); // number of steps

    // * 3. Run simulation according to different dynamics.
    // -- Clssical MD simulation with OpenMM
    // -- Nonadiabatic dynamics simulation for realistic system with OpenMM
    if (system_type == "allatom" && allatom_type == "OpenMM") {
        // Both classical MD and noadiabatic MD support multi-traj simulation.
        // beforeOneTraj() will load the positions and velocities from a given
        // trajectory if nucl_load is not empty (and set simulation time/step to
        // zero for new trajectory).
        for (int traj = 0; traj < ntraj; ++traj) {
            Dy->beforeOneTraj();
            static const bool perturb = param.getBool("perturb");
            // const int skip_steps = 1; Plan A
            if (dyn_type == "OpenMM"){
                if (perturb) {
                    std::cout << "Running OpenMM dynamics with perturbation" << std::endl;
                    runPerturbOpenMM(traj+1, nsteps);
                } else {
                runOpenMM(traj+1, nsteps);
                }
            }
            else { 
                // Plan A for NAMD
                // Obs->updateOneSnapshot();
                // for (int step = 0; step < nsteps; ) {
                //     runOpenMM(traj+1, skip_steps);
                //     int traj_skip_index = Obs->updateOneSnapshot();
                //     // std::cout << step <<" dynamics steps is done" << std::endl;
                //     // if there are a not-a-number, breaks; else the traj will normally propagates
                //     if (traj_skip_index == 0)  
                //         step += skip_steps;
                //     else 
                //         break;
                // Plan B for NAMD
                Obs->updateOneSnapshot();
                runOpenMM(traj+1, nsteps); // run a complete trajectory
            }
            Dy->afterOneTraj();
            // Save reduced density matrix (RDM) of current single trajectory
            // Do this only when ntraj > 1, since when ntraj = 1, the averaged
            // RDM saved after all trajectory is same. (nonadiabatic dynamics)
            if (dyn_type != "OpenMM" && dyn_type != "ReadTraj" && ntraj > 1) {
                // There are more than one RDMs in some simulation, so the filenames is
                // a vector, whose size will equal to the size of RDM got from dynamcis.
                // int traj_skip_index = Obs->updateOneTraj();
                int traj_avail_index = Obs->updateOneTraj();
                std::cout << "traj_avail_index " << traj_avail_index << std::endl;
                traj += traj_avail_index;
                
                // const std::vector<Complex_Matrix>& RDM = Obs->getCurrentDensityMatrix();
                // static std::vector<std::string> filenames;
                // if (traj == 0) // do this at first traj is enough
                //     getRDMFileName(filenames);
                //for (int i = 0; i < RDM.size(); ++i) {
                    // For multi-traj simulation, add a suffix "_traj?" to filename, where "?" is the index of traj (starting from 1).
                    // std::string filename = GetFilePrefix(filenames[i]) + "_traj" + std::to_string(traj+1) + "." + GetFileSuffix(filenames[i]);
                    // writeRDMFile(RDM[i], filename);
                 // to-do all-atom RDM class
                //}
            }
        }
        if (dyn_type != "OpenMM") {
            Obs->output();
        }
    }
    // -- Clssical MD simulation with custom force field
    else if (system_type == "allatom" && allatom_type == "CustomForceField") {
        // Both classical MD and noadiabatic MD support multi-traj simulation.
        // beforeOneTraj() will load the positions and velocities from a given
        // trajectory if nucl_load is not empty (and set simulation time/step to
        // zero for new trajectory).
        static const bool perturb = param.getBool("perturb");
        for (int traj = 0; traj < ntraj; ++traj) {
            Dy->beforeOneTraj();
            if (dyn_type == "MD") {
                if (perturb) {
                    std::cout <<"begin to run the perturb MD" <<std::endl;
                    runPerturbMD(traj+1, nsteps);
                }
                else {
                    std::cout <<"begin to run MD" <<std::endl;
                    runMD(traj+1, nsteps); // run a complete trajectory
                }
            }
            else if (dyn_type == "ReadTraj") {
                runTrajMD(traj+1, nsteps);
            }
            else if (dyn_type != "MD" ||dyn_type != "ReadTraj" ){ 
                // for NAMD
                Obs->updateOneSnapshot();
                std::cout<<" run NAMD"<<std::endl;
                runMD(traj+1, nsteps); // run a complete trajectory
            }
            Dy->afterOneTraj();
            // Save reduced density matrix (RDM) of current single trajectory
            // Do this only when ntraj > 1, since when ntraj = 1, the averaged
            // RDM saved after all trajectory is same. (nonadiabatic dynamics)
            if (dyn_type != "CustomForceField" && dyn_type != "ReadTraj" && ntraj > 1) {
                // There are more than one RDMs in some simulation, so the filenames is
                // a vector, whose size will equal to the size of RDM got from dynamcis.
                // int traj_skip_index = Obs->updateOneTraj();
                int traj_avail_index = Obs->updateOneTraj();
                std::cout << "traj_avail_index " << traj_avail_index << std::endl;
                traj += traj_avail_index;
            }
        }
        if (dyn_type != "CustomForceField" && dyn_type != "ReadTraj") {
            Obs->output();
        }
    }
    // -- Nonadiabatic dynamics simulation for model system
    else if (system_type == "model") {
        // This is special simulation used to analysis the potential energy and
        // and kinetic energy of each nromal mode along the trajectory. Therefore,
        // the ensemble-averaged energies of each mode will be outputed to a file.
        // The frequency of output is controlled by the key "energy_steps".
        // And the filename is controlled by the key "data_file".
        if (param.getBool("energy_decompose")) {
            // TODO: Currently, this function only supports SB model.
            if (model_type != "SB" && model_type != "SB_req" && model_type != "GOA")
                throw std::runtime_error("ERROR: The energy decompose for model only "
                    "works for SB/SB_req/GOA models now.");
            // The "energy_steps" controls the report frequency (number of steps)
            const int energy_steps = param.getInt("energy_steps");
            if (energy_steps <= 0 || energy_steps > nsteps || nsteps % energy_steps != 0)
                throw std::runtime_error("ERROR: Illeagal value for energy_steps when "
                "performing energy decompose for model.");
            // Here, use skip_trajs to control the frequency to report simulation
            // progress in stdout. By default, the progress interval is 10%,
            // then skip_trajs is 10%*ntraj (at least 1).
            const int skip_trajs = (ntraj*10/100) == 0 ? 1 : (ntraj*10/100);
            for (int traj = 0; traj < ntraj; ++traj) {
                Dy->beforeOneTraj();
                reportModelData(traj+1, 0); // report the initial step always
                for (int step = 0; step < nsteps; ) {
                    Dy->dynamics(energy_steps);
                    step += energy_steps;
                    // get and report kinetic and potential energy of each mode
                    reportModelData(traj+1, step);
                }
                Dy->afterOneTraj();
                if ((traj+1) % skip_trajs == 0 || (traj+1) == ntraj)
                    report(traj+1, nsteps);
            }
        } // TODO: Perhaps merge it to the above one
        else { // Normal simulation without energy components input
            // For model system, the dynamics is taken with one traj by one traj.
            // Here, use skip_trajs to control the frequency to report simulation
            // progress in stdout. By default, the progress interval is 10%,
            // then skip_trajs is 10%*ntraj (at least 1).
            const int skip_trajs = (ntraj*10/100) == 0 ? 1 : (ntraj*10/100);
            // Added on Dec. 31, 2021
            // TODO: 
            // control the frequency to call functions of observables objects
            // const int skip_steps = 1;
            // debug codes by Cesare
            // std::cout << "simulation is on" << std::endl;
            for (int traj = 0; traj < ntraj; ++traj) {
                // debug codes by Cesare
                // std::cout << "simulation is on" << std::endl;
                Dy->beforeOneTraj();
                // debug codes by Cesare
                // std::cout << "initialization is done" << std::endl;
                // update observables always at the step 0
                Obs->updateOneSnapshot();
                for (int step = 0; step < nsteps; ) {
                    Dy->dynamics(RDM_steps);
                    int traj_skip_index = Obs->updateOneSnapshot();
                    // std::cout << step <<" dynamics steps is done" << std::endl;
                    // if there are a not-a-number, breaks; else the traj will normally propagates
                    if (traj_skip_index == 0)  
                        step += RDM_steps;
                    else 
                        break;
                }
                Dy->afterOneTraj();
                // debug codes by Cesare
                // std::cout << "Dynamics is done for 1 traj "<< std::endl;
                // if there is not-a-number or infinity appears, we drop this trajectories
                int traj_avail_index = Obs->updateOneTraj();
                // std::cout << "traj_avail_index " << traj_avail_index << std::endl;
                traj += traj_avail_index;
                if ((traj+1) % skip_trajs == 0 || (traj+1) == ntraj)
                    report(traj+1, nsteps);
                
                double fail_traj_ratio = Obs->getCurrentFailRatio();
                // const int begin_judge_traj_avail = 5000; ////Notice same index in ObservableRDM.cpp
                if (fail_traj_ratio > max_friction_fail_traj && traj > min_judge_traj_avail) {
                   break;
                }
            }
            // debug codes by Cesare
            Obs->output();
            // std::vector<std::string> filenames;
            // getRDMFileName(filenames);
            // for (int file_count = 0; file_count < filenames.size(); ++file_count) 
        }
    }
    // * 4. Finish the simulation and record complete time
    std::cout << "The dynamics simulation is finished at " << CurrentTime() << ".\n";
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedTime = end - start;
    std::cout << "The dynamics elapsed time is " << std::fixed << std::setprecision(3) <<
        elapsedTime.count()  << " seconds.\n" << std::endl;
}

void Simulation::report(int traj, int step) {
    // For simulation with OpenMM including classical MD and nonadiabatic
    // dynamics simultion, report progress, energies and trajectory and so on.
    if (system_type == "allatom" && allatom_type == "OpenMM")
        reportOpenMMData(traj, step);
    // For CustomForceField, report report progress, energies and trajectory and so on.
    else if (system_type == "allatom" && allatom_type == "CustomForceField")

        if (dyn_type == "ReadTraj") {
            reportReadTrajData(traj, step);
        }
        else {
            if (PIMD_type == "RPMD") 
                reportRPMDData(traj, step);
            else 
                reportMDData(traj, step);
        }
    // For model, report the simulation progress only
    else if (system_type == "model")
        reportProgress(traj, step);
}

void Simulation::reportReadTrajData(int traj, int step) {
    static int count = 0; // record the number of times to eneter here
    static int traj_prev = traj; // the traj of previous call
    if (traj_prev != traj) { // if it is not a same traj,
        count = 0;           // we reset the count to 0,
        traj_prev = traj;    // and set current traj as previous one.
    }
    static const double                DT = param.getDouble("DT");
    static const std::string default_name = param.getStr("default_name");
    static const bool useDefault = (default_name.empty() || default_name == "none") ? false : true;
    static const int                ntraj = param.getInt("ntraj");
    static const int               nsteps = param.getInt("nsteps");
    static const int progress_interval = param.getInt("progress_interval");
    static const int progress_steps    = ntraj * nsteps * progress_interval / 100;
    double                          time  = step * DT;
    if (count == 0 || progress_steps == 0 || ((traj-1) * nsteps + step) % progress_steps == 0 || step == nsteps)
    // when traj > 1, the step 0 of current traj and last step of las traj is repetitive
    // which should be aviod to report the step 0 in this case.
    if (traj == 1 || (traj > 1 && step != 0))
        reportProgress(traj, step);
    static const std::string data_file = useDefault ? (default_name + ".csv") : param.getStr("data_file");
    std::string filename = ntraj == 1 ? data_file : GetFilePrefix(data_file) + "_traj" + std::to_string(traj) + "." + GetFileSuffix(data_file);
    static std::vector<double> data;
    getPiState(traj, data);
    writePiState(filename, data, step, time);
    count++;
}

void Simulation::getPiState(int traj, std::vector<double>& data) {
    // Get smarter ponter to the HamiltonianreportMDData obejct
    static std::shared_ptr<HamiltonianForceFieldBase> ha = std::static_pointer_cast<HamiltonianForceFieldBase>(Ha);
    static int count = 0;
    static int traj_prev = traj; // the traj of previous call
    if (traj_prev != traj) { // if it is not a same traj,
        count = 0;           // we reset the count to 0,
        traj_prev = traj;    // and set current traj as previous one.
    }
    int DOFe = param.getInt("DOFe");
    data.resize(9 * DOFe, 0.0);
    std::fill(data.begin(), data.end(), 0);
    std::vector<double> PiCurrentTensor;
    PiCurrentTensor.resize(9, 0.0);
    for (int i = 0; i < DOFe; i++) { // loop from state 1
        ha->getPiTensor(i, PiCurrentTensor); 
        data[0 + i * 9] = PiCurrentTensor[0]; 
        data[1 + i * 9] = PiCurrentTensor[1]; 
        data[2 + i * 9] = PiCurrentTensor[2]; 
        data[3 + i * 9] = PiCurrentTensor[3];
        data[4 + i * 9] = PiCurrentTensor[4]; 
        data[5 + i * 9] = PiCurrentTensor[5]; 
        data[6 + i * 9] = PiCurrentTensor[6];
        data[7 + i * 9] = PiCurrentTensor[7];
        data[8 + i * 9] = PiCurrentTensor[8];
    }
    count++;
}

void Simulation::writePiState(const std::string& file, std::vector<double>& data, int step, double time) {
    // * Get required parameters.
    static int count = 0;
    static std::string file_prev = file; // the filename of previous call
    if (file_prev != file) { // if it is not a same file,
        count = 0;           // we reset the count to 0,
        file_prev = file;    // and set current file name as previous one.
    }
    const std::string type = GetFileSuffix(file); // file format type
    static const double                DT = param.getDouble("DT");
    int DOFe = param.getInt("DOFe");
    FILE* dataFile = CheckFile(file, count);
    if (type == "dat") { // * Using fixed space separated format (.dat)
        for (int state = 0; state < DOFe; ++state) { // Write data of each state.
            fprintf(dataFile, " %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f", data[0+state*9], data[1+state*9], 
                data[2+state*9], data[3+state*9], data[4+state*9], data[5+state*9], 
                data[6+state*9], data[7+state*9],data[8+state*9]);
            fprintf(dataFile, "\n");
        }
    }
    else if (type == "csv") { // * comma (,) separated value (CSV) file
        for (int state = 0; state < DOFe; ++state) { // Write data of each state
            fprintf(dataFile, " %.24g  %.24g  %.24g  %.14g  %.24g  %.24g   %.24g   %.24g   %.24g", data[0+state*9], data[1+state*9], 
                data[2+state*9], data[3+state*9], data[4+state*9], data[5+state*9], 
                data[6+state*9], data[7+state*9], data[8+state*9]);
            fprintf(dataFile, "\n");
        }
    }
    else
        throw std::runtime_error("ERROR: Unknown file type for output: " + type);
    fclose(dataFile);
    count++;
}

void Simulation::runTrajMD(int traj, int nsteps) {
    // Get current simulation step, which may not zero if it is a restart
    // simulation for classical MD simulation, while it is always zero for
    // nonadiabatic dynamics simulations.
    int step = Dy->getStep();
    if (nsteps == 0) // if total step is 0, compute energy only
        return;      // and skip the following dynbamics
    // * Read the frame from trajectory as initial conditions
    std::string loadfile;
    loadfile = param.getStr("nucl_load"); // read initial conditions from file
    // Open a file for reading, automatically guessing the file format from the extension.
    // It will throw an exception if the format is not supported.
    //static auto trajectory = chemfiles::Trajectory(loadfile);
    static chemfiles::Trajectory trajectory(loadfile);
    const int ntraj  = param.getInt("ntraj");  // number of trajectories
    static const auto nframes = trajectory.nsteps(); // number of frames
    // Read positions, velocities, and box of current frame.
    static const int nucl_start = param.getInt("nucl_start");
    static const int nucl_skip  = param.getInt("nucl_skip");
    static const int nucl_end   = param.getInt("nucl_end") == -1 ? (nframes-1) : param.getInt("nucl_end");
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

    for (; step < nsteps ;) {
        Dy->dynamics(nucl_skip);
        int traj_skip_index = Obs->updateOneSnapshot();
        if (traj_skip_index == 0)  
                step += nucl_skip;
            else 
                break;   

        report(traj, step);
    }
    Obs->output(); 
}

int Simulation::getReportFrequency() {
    // Get the report frequency of each item to control simulaton.
    // If XXX_steps = 0, that is to say no output of this item.
    const int nsteps = param.getInt("nsteps");
    std::vector<int> report_frequency;
    report_frequency.push_back(param.getInt("energy_steps"));
    report_frequency.push_back(param.getInt("position_steps"));
    report_frequency.push_back(param.getInt("chk_steps"));
    // Check the validity of them:
    // 1) cannot < 0
    // 2) If they are larger than nsteps, then set them to nsteps.
    // 3) nsteps is a multiply of them
    for (auto& i : report_frequency) // use reference to write
        if (i < 0) // < 0 is not allowed
            throw std::runtime_error("ERROR: Found negative frequency in energy_steps or "
            " position_steps or chk_steps, please check your input.");
        else if (i > nsteps) // if > nsteps, reset to nsteps
            i = nsteps;
        else if (i != 0 && nsteps % i != 0) // nsteps is a multiply of them
            throw std::runtime_error("ERROR: The value of nsteps must be divided "
                "exactly by reprot frequency (energy_steps, position_steps, chk_steps).");
    // Get the skip_steps, the minimum value of them (non-zero), which means the
    // frequency to stop to call report(), also denotes the slient steps to
    // run dynamics in interal OpenMM integrator.
    std::sort(report_frequency.begin(), report_frequency.end()); // ascending
    int skip_steps = 0;
    for (auto i : report_frequency) // get the non-zero minimum value
        if (i == 0)
            continue;
        else {  // the minimum value is the first non-zero value
            skip_steps = i;
            break;
        }
    if (skip_steps == 0) // if all of them are zero, then set it to nsteps.
        skip_steps = nsteps;
    // Check if the skip_steps is the common divisor of them.
    for (auto i : report_frequency)
        if (i > skip_steps && (i % skip_steps != 0))
            throw std::runtime_error("ERROR: Please check the report frequecny "
                "(energy_steps, position_steps, chk_steps), the non-zero minimum "
                "value of them should be a common divisor of non-zero others.");
    // Get and check the progress_interval (in percentage) which decides the
    // output frequecny of simulation progress in the STDOUT.
    const int progress_interval = param.getInt("progress_interval");
    if (progress_interval < 0 || progress_interval > 100)
        throw std::runtime_error("ERROR: Illegal value of progress_interval: " + param.getStr("progress_interval"));
    return skip_steps;
}

void Simulation::reportMDData(int traj, int step) {
    static int count = 0; // record the number of times to eneter here
    static int traj_prev = traj; // the traj of previous call
    if (traj_prev != traj) { // if it is not a same traj,
        count = 0;           // we reset the count to 0,
        traj_prev = traj;    // and set current traj as previous one.
    }
    // Get smarter ponter to the HamiltonianreportMDData obejct
    static std::shared_ptr<HamiltonianForceFieldBase> ha = std::static_pointer_cast<HamiltonianForceFieldBase>(Ha);
    // Get the total number of trajectory (for classical MD, ntraj is always 1)
    // and steps and the frequency (number of steps) to report.
    // Note, the data of step 0 (initial configuration) and last step will be
    // reported all the time unless 0, which means the file won't be genenrated.
    static const int                ntraj = param.getInt("ntraj");
    static const int               nsteps = param.getInt("nsteps");
    static const int         energy_steps = param.getInt("energy_steps");
    static const int       position_steps = param.getInt("position_steps");
    // Get current MD simulation time.
    static const double                DT = param.getDouble("DT");
    double                          time  = step * DT;
    // If default_name is specified, it will be set to the file name of output
    static const std::string default_name = param.getStr("default_name");
    static const bool useDefault = (default_name.empty() || default_name == "none") ? false : true;
    // Get the progress_steps to decide report simulation progress at which step
    // Here, progress_interval (in percentage) decides the output frequecny of
    // simulation progress in the STDOUT.
    // Note, for multi-traj nonadiabatic simulation, the ntraj can be > 1.
    static const int progress_interval = param.getInt("progress_interval");
    static const int progress_steps    = ntraj * nsteps * progress_interval / 100;
    // 1. Report simulation progress
    // If progress_steps is 0, the simulation progress will be reported always.
    // And at the initial step and final steps, it will be reported always.
    // Here, ((traj-1) * nsteps + step) is the current step in the all steps of
    // multi-traj simulation. For classical MD simulation, it is just step.
    if (count == 0 || progress_steps == 0 || ((traj-1) * nsteps + step) % progress_steps == 0 || step == nsteps)
        // when traj > 1, the step 0 of current traj and last step of las traj is repetitive
        // which should be aviod to report the step 0 in this case.
        if (traj == 1 || (traj > 1 && step != 0)) {
            reportProgress(traj, step);
        }
    // 2. Report energies, and other data, such as temeprature
    // For multi-traj simultaion (ntraj > 1), the data filename of each trajectory
    // will add a suffix "_traj?", where "?" is the index of traj (starting from 1).
    if (energy_steps > 0 && (step % energy_steps == 0 || step == nsteps)) {
        static const std::string data_file = useDefault ? (default_name + ".csv") : param.getStr("data_file");
        std::string filename = ntraj == 1 ? data_file : GetFilePrefix(data_file) + "_traj" + std::to_string(traj) + "." + GetFileSuffix(data_file);
        static const bool energy_decompose = param.getBool("energy_decompose");
        static std::vector<double> data;
        getMDData(traj, data, false, energy_decompose);
        writeMDData(filename, data, step, time);
    }
    // 3. Report trajectory: coordinates, (optional velocities and forces)
    // TODO: using structure object in Hamiltonina directly
    // TODO with function updateStructureData(includeForces)
    if (position_steps > 0 && (step % position_steps == 0 || step == nsteps)) {
        static int frame = 0; // record the frame in trajectory file
        static int traj_prev2 = traj; // the traj of previous call
        if (traj_prev2 != traj) { // if it is not a same traj,
            frame = 0;            // we reset the count to 0,
            traj_prev2 = traj;    // and set current traj as previous one.
        }
        static std::string traj_file, traj_type;
        static const std::string traj_group = param.getStr("traj_group");
        static const int      xtc_precision = param.getInt("xtc_precision");
        static const bool     velocity_save = param.getBool("velocity_save");
        static const bool        force_save = param.getBool("force_save");
        static const int    propagate_state = param.getInt("propagate_state");
        static const bool useBarostat = param.getStr("barostat") != "none" ? true : false;
        static const bool usePBC = param.getStr("PBC") != "none" ? true : false;
        // decide the trajectory file name and check the type of trajectory at
        // first time if default_name is not used.
        if (frame == 0 && traj == 1) { // do this only at first time for first trajectory
            if (useDefault) {
                if (velocity_save || force_save)
                    traj_file = default_name + ".trr";
                else
                    traj_file = default_name + ".xtc";
                traj_type = GetFileSuffix(traj_file);
            }
            else {
                traj_file = param.getStr("traj_file");
                traj_type = GetFileSuffix(traj_file);
                if (traj_type == "xtc" && (velocity_save || force_save))
                    throw std::runtime_error("ERROR: xtc trajectory file doesn't support "
                        "to store velocities or forces, please use trr file.");
                if (traj_type == "gro" && force_save)
                    throw std::runtime_error("ERROR: gro trajectory file doesn't support "
                        "to store forces, please use trr file.");
            }
        }
        // For multi-traj simulation, add a suffix "_traj?" to filename, where "?" is the index of traj (starting from 1).
        std::string filename = ntraj == 1 ? traj_file : GetFilePrefix(traj_file) + "_traj" + std::to_string(traj) + "." + GetFileSuffix(traj_file);
        // Create a Structure object, and get the atominfo, positions,
        // velocities, boxVectors of current step from HamiltonianOpenMM
        // object, then use method in Structure to save trajectory file.
        static StructureVec3 structure; // empty Structure object
        // Get writable reference to structure data member
        static std::vector<std::string>&  atominfo   = structure.getAtomInfo();
        static std::vector<Vec3>& positions  = structure.getPositions();
        static std::vector<Vec3>& velocities = structure.getVelocities();
        static std::vector<Vec3>& forces     = structure.getForces();
        static Vec3         (&boxVectors)[3] = structure.getBoxVectors();
        // Set the atominfo and natoms, and never change
        if (frame == 0 && traj == 1) { // do this at first time for first trajectory
            atominfo = ha->getAtomInfo();
            structure.setNumAtoms(atominfo.size());
        }
        // boxVectors changes when barostat is used or for a new trajectory
        // For a non-PBC simulation, the boxVectors are always zero.
        if (usePBC && (frame == 0 || useBarostat))
            ha->getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        // Get positions and always needed
        ha->getPositions(positions);
        // If a subsystem is required to save to tajectory file, and the trajectory
        // format is not gro file. The initial stucture of this subsystem will be
        // saved to a Gromacs structure file with a filename with suffix "_subsystem".
        // This subsystem structure file is necessary when you want to visualize the
        // subsystem binary trajectory file with suah as VMD, since in binary
        // trajectory file, the atom infomation is not included.
        if (frame == 0 && traj == 1 && !traj_group.empty() && traj_group != "system" && traj_type != "gro")
            structure.saveTrajectory(GetFilePrefix(traj_file) + "_subsystem.gro", step, time, traj_group);
        if (velocity_save)
            ha->getVelocities(velocities);
        if (force_save) {
            // For classical M, the forces of current positions should be
            // recalculated. And only the state used to propagate is saved.
            // TODO: How to aviod this calculation for classical MD
            if (dyn_type == "MD") { // classical MD, recompute the forces
                ha->getPotentialEnergy(propagate_state);
                ha->getForces(forces);
            }
            else { // nonadiabatic dynamics, copy effective forces from Hamiltonian
                ha->getForces(forces);
            }
        }
        // Save current sturcture data to trajectory file
        structure.saveTrajectory(filename, step, time, traj_group, xtc_precision);
        frame++;
    }
    // 4. Create checkpoint file (the initial step is skipped)
    //if (chk_steps > 0 && (step % chk_steps == 0 || step == nsteps)) {
    //    static const std::string chk_file = useDefault ? (default_name + ".chk") : param.getStr("chk_file");
    //    // For multi-traj simulation, add a suffix "_traj?" to filename, where "?" is the index of traj (starting from 1).
    //    std::string filename = ntraj == 1 ? chk_file : GetFilePrefix(chk_file) + "_traj" + std::to_string(traj) + "." + GetFileSuffix(chk_file);
    //    // At the initial step, the check point file will be saved with a suufix
    //    // "_init", which can be used to restart a simulation in the future.
    //    if (count == 0)
    //        ha->createCheckpoint(GetFilePrefix(filename) + "_init.chk");
    //    else
    //        ha->createCheckpoint(filename);
    //}
    //5. Save structure file at last step, containing positions and velocities
    //The structure file will always save whole system (traj_group doesn't
    //infulence it) at last step.
    if (step == nsteps && param.getStr("conf_file") != "none") {
        // TODO: using structure object in Hamiltonian directly
        // Create a structure object and get data from Hamiltonian obejct
        StructureVec3 structure;
        std::vector<std::string>&  atominfo   = structure.getAtomInfo();
        std::vector<Vec3>& positions  = structure.getPositions();
        std::vector<Vec3>& velocities = structure.getVelocities();
        Vec3         (&boxVectors)[3] = structure.getBoxVectors();
        atominfo = ha->getAtomInfo();
        structure.setNumAtoms(atominfo.size());
        ha->getPositions(positions);
        ha->getVelocities(velocities);
        if (param.getStr("PBC") != "none")
            ha->getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        // If same name as this trajectory file is used. Then the filename
        // with a suffix "_last" will be used.
        std::string conf_file = useDefault ? (default_name + ".gro") : param.getStr("conf_file");
        if (conf_file == param.getStr("traj_file"))
            conf_file = GetFilePrefix(conf_file) + "_last." + GetFileSuffix(conf_file);
        // For multi-traj simulation, add a suffix "_traj?" to filename, where "?" is the index of traj (starting from 1).
        if (ntraj > 1)
            conf_file = GetFilePrefix(conf_file) + "_traj" + std::to_string(traj) + "." + GetFileSuffix(conf_file);
        structure.saveStructure(conf_file, step, time);
    }
    count++;
}

void Simulation::reportRPMDData(int traj, int step) {
    static int count = 0; // record the number of times to eneter here
    static int traj_prev = traj; // the traj of previous call
    if (traj_prev != traj) { // if it is not a same traj,
        count = 0;           // we reset the count to 0,
        traj_prev = traj;    // and set current traj as previous one.
    }
    // Get smarter ponter to the HamiltonianreportMDData obejct
    static std::shared_ptr<HamiltonianRingPolymer> ha = std::static_pointer_cast<HamiltonianRingPolymer>(Ha);
    // Get the total number of trajectory (for classical MD, ntraj is always 1)
    // and steps and the frequency (number of steps) to report.
    // Note, the data of step 0 (initial configuration) and last step will be
    // reported all the time unless 0, which means the file won't be genenrated.
    static const int                ntraj = param.getInt("ntraj");
    static const int               nsteps = param.getInt("nsteps");
    static const int         energy_steps = param.getInt("energy_steps");
    static const int       position_steps = param.getInt("position_steps");
    // Get current MD simulation time.
    static const double                DT = param.getDouble("DT");
    double                          time  = step * DT;
    // If default_name is specified, it will be set to the file name of output
    static const std::string default_name = param.getStr("default_name");
    static const bool useDefault = (default_name.empty() || default_name == "none") ? false : true;
    // Get the progress_steps to decide report simulation progress at which step
    // Here, progress_interval (in percentage) decides the output frequecny of
    // simulation progress in the STDOUT.
    // Note, for multi-traj nonadiabatic simulation, the ntraj can be > 1.
    static const int progress_interval = param.getInt("progress_interval");
    static const int progress_steps    = ntraj * nsteps * progress_interval / 100;
    // 1. Report simulation progress
    // If progress_steps is 0, the simulation progress will be reported always.
    // And at the initial step and final steps, it will be reported always.
    // Here, ((traj-1) * nsteps + step) is the current step in the all steps of
    // multi-traj simulation. For classical MD simulation, it is just step.
    if (count == 0 || progress_steps == 0 || ((traj-1) * nsteps + step) % progress_steps == 0 || step == nsteps)
        // when traj > 1, the step 0 of current traj and last step of las traj is repetitive
        // which should be aviod to report the step 0 in this case.
        if (traj == 1 || (traj > 1 && step != 0)) {
            reportProgress(traj, step);
        }
    // 2. Report energies, and other data, such as temeprature
    // For multi-traj simultaion (ntraj > 1), the data filename of each trajectory
    // will add a suffix "_traj?", where "?" is the index of traj (starting from 1).
    if (energy_steps > 0 && (step % energy_steps == 0 || step == nsteps)) {
        static const std::string data_file = useDefault ? (default_name + ".csv") : param.getStr("data_file");
        std::string filename = ntraj == 1 ? data_file : GetFilePrefix(data_file) + "_traj" + std::to_string(traj) + "." + GetFileSuffix(data_file);
        static const bool energy_decompose = param.getBool("energy_decompose");
        static std::vector<double> data;
        getRPMDData(traj, data, false, energy_decompose);
        writeMDData(filename, data, step, time);
    }
    // 3. Report trajectory: coordinates, (optional velocities and forces)
    // TODO: using structure object in Hamiltonina directly
    // TODO with function updateStructureData(includeForces)
    if (position_steps > 0 && (step % position_steps == 0 || step == nsteps)) {
      static int frame = 0; // record the frame in trajectory file
      static int traj_prev2 = traj; // the traj of previous call
        if (traj_prev2 != traj) { // if it is not a same traj,
          frame = 0;            // we reset the count to 0,
          traj_prev2 = traj;    // and set current traj as previous one.
        }
        static std::string traj_file, traj_type;
        static const std::string traj_group = param.getStr("traj_group");
        static const int      xtc_precision = param.getInt("xtc_precision");
        static const bool     velocity_save = param.getBool("velocity_save");
        static const bool        force_save = param.getBool("force_save");
        static const int    propagate_state = param.getInt("propagate_state");
        static const bool useBarostat = param.getStr("barostat") != "none" ? true : false;
        static const bool usePBC = param.getStr("PBC") != "none" ? true : false;
        
        // decide the trajectory file name and check the type of trajectory at
        // first time if default_name is not used.
        if (frame == 0 && traj == 1) { // do this only at first time for first trajectory
            if (useDefault) {
                if (velocity_save || force_save)
                    traj_file = default_name + ".trr";
                else
                    traj_file = default_name + ".xtc";
                traj_type = GetFileSuffix(traj_file);
            }
            else {
                  traj_file = param.getStr("traj_file");
                  traj_type = GetFileSuffix(traj_file);
                  if (traj_type == "xtc" && (velocity_save || force_save))
                      throw std::runtime_error("ERROR: xtc trajectory file doesn't support "
                          "to store velocities or forces, please use trr file.");
                  if (traj_type == "gro" && force_save)
                      throw std::runtime_error("ERROR: gro trajectory file doesn't support "
                          "to store forces, please use trr file.");
              }
        }
        // For multi-traj simulation, add a suffix "_traj?" to filename, where "?" is the index of traj (starting from 1).
        std::string filename = ntraj == 1 ? traj_file : GetFilePrefix(traj_file) + "_traj" + std::to_string(traj) + "." + GetFileSuffix(traj_file);
        // Create a Structure object, and get the atominfo, positions,
        // velocities, boxVectors of current step from HamiltonianOpenMM
        // object, then use method in Structure to save trajectory file.
        static StructureVec3 structure; // empty Structure object
        // Get writable reference to structure data member
        std::vector<std::string> atominfo;
        std::vector<std::vector<Vec3>> positions_RP;
        std::vector<std::vector<Vec3>> velocities_RP;
        std::vector<std::vector<Vec3>> forces_RP;
        Vec3         (&boxVectors)[3] = structure.getBoxVectors();
        // Set the atominfo and natoms, and never change
        if (frame == 0 && traj == 1) { // do this at first time for first trajectory
            atominfo = ha->getAtomInfo();
            structure.setNumAtoms(atominfo.size());
        }
        // boxVectors changes when barostat is used or for a new trajectory
        // For a non-PBC simulation, the boxVectors are always zero.
        if (usePBC && (frame == 0 || useBarostat))
            ha->getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        // Get positions and always needed
        ha->getPositions(positions_RP);
        // If a subsystem is required to save to tajectory file, and the trajectory
        // format is not gro file. The initial stucture of this subsystem will be
        // saved to a Gromacs structure file with a filename with suffix "_subsystem".
        // This subsystem structure file is necessary when you want to visualize the
        // subsystem binary trajectory file with suah as VMD, since in binary
        // trajectory file, the atom infomation is not included.
        if (frame == 0 && traj == 1 && !traj_group.empty() && traj_group != "system" && traj_type != "gro")
            structure.saveTrajectory(GetFilePrefix(traj_file) + "_subsystem.gro", step, time, traj_group);
        if (velocity_save)
            ha->getVelocities(velocities_RP);
        if (force_save) {
            // For classical M, the forces of current positions should be
            // recalculated. And only the state used to propagate is saved.
            // TODO: How to aviod this calculation for classical MD
            if (dyn_type == "MD") { // classical MD, recompute the forces
                ha->getPotentialEnergy(propagate_state);
                ha->getForces(forces_RP);
            }
            else { // nonadiabatic dynamics, copy effective forces from Hamiltonian
                ha->getForces(forces_RP);
            }
        }
        // Save current sturcture data to trajectory file
        structure.saveTrajectory(filename, step, time, traj_group, xtc_precision);
        frame++;
    }
    //5. Save structure file at last step, containing positions and velocities
    //The structure file will always save whole system (traj_group doesn't
    //infulence it) at last step.
    if (step == nsteps && param.getStr("conf_file") != "none") {
        // TODO: using structure object in Hamiltonian directly
        // Create a structure object and get data from Hamiltonian obejct
        StructureVec3 structure;
        std::vector<std::string>&  atominfo   = structure.getAtomInfo();
        std::vector<std::vector<Vec3>> positions_RP;
        std::vector<std::vector<Vec3>> velocities_RP;
        Vec3         (&boxVectors)[3] = structure.getBoxVectors();
        atominfo = ha->getAtomInfo();
        structure.setNumAtoms(atominfo.size());
        ha->getPositions(positions_RP);
        ha->getVelocities(velocities_RP);
        if (param.getStr("PBC") != "none")
            ha->getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        // If same name as this trajectory file is used. Then the filename
        // with a suffix "_last" will be used.
        std::string conf_file = useDefault ? (default_name + ".gro") : param.getStr("conf_file");
        if (conf_file == param.getStr("traj_file"))
            conf_file = GetFilePrefix(conf_file) + "_last." + GetFileSuffix(conf_file);
        // For multi-traj simulation, add a suffix "_traj?" to filename, where "?" is the index of traj (starting from 1).
        if (ntraj > 1)
            conf_file = GetFilePrefix(conf_file) + "_traj" + std::to_string(traj) + "." + GetFileSuffix(conf_file);
        structure.saveStructure(conf_file, step, time);
    }
    count++;
}

void Simulation::reportModelData(int traj, int step) {
    // This function only works for HamiltonianSpinBoson class now.
    // Get smarter ponter to the HamiltonianSpinBoson obejct
    static std::shared_ptr<HamiltonianSpinBoson> ha = std::static_pointer_cast<HamiltonianSpinBoson>(Ha);
    static const int N = ha->getNumNormalModes(); // number of normal mdoes
    static const int DOFe = param.getInt("DOFe"); // number of states
    static const int ntraj  = param.getInt("ntraj");  // number of trajectories
    static const int nsteps = param.getInt("nsteps"); // number of steps
    static const int energy_steps = param.getInt("energy_steps"); // report frequency
    static const double DT  = param.getDouble("DT"); // nuclear time step size
    // PE contains the potential energy of each noraml mode of each state at each
    // report step. The size of it is the step, and the element is the PE of nomral
    // modes of first state, sencond state, ..., so the length of it is N*DOFe.
    static std::vector<std::vector<double>> PE(nsteps/energy_steps+1, std::vector<double>(N*DOFe, 0));
    // Similar as PE, KE contains the kinetic energy of each noraml mode, but the
    // length of element is N, since the velocities are same for each state.
    static std::vector<std::vector<double>> KE(nsteps/energy_steps+1, std::vector<double>(N, 0));
    // Get PE and KE of each normal mode from Hamiltonian object.
    // Note: here, the PE and KE has been accumulated and averaged.
    for (int j = 0; j < N; ++j) { // index of normal modes
        for (int i = 0; i < DOFe; ++i) // index of state
            PE[step/energy_steps][i*N+j] += ha->getNormalModePE(i, j)/ntraj;
        KE[step/energy_steps][j] += ha->getNormalModeKE(j)/ntraj;
    }
    // Save PE and KE to csv file when this is the final step of last trajectory.
    // Note: traj is 1, ..., ntraj and step is 0, 1, ..., nsteps.
    if (traj == ntraj && step == nsteps) {
        const std::string file = param.getStr("data_file");
        std::ofstream out(file);
        if (!out)
            throw std::runtime_error("ERROR: Can't create the file: " + file);
        // 1st line is headers, here, index of mode/state starts from 1 in output
        // They are time, KE of each mode, PE of each mode of 1st state, 2nd state...
        std::string headers = "time";
        for (int j = 0; j < N; ++j)
            headers += ",KE_" + std::to_string(j+1);
        for (int i = 0; i < DOFe; ++i)
            for (int j = 0; j < N; ++j)
                headers += ",PE" + std::to_string(i+1) + "_" + std::to_string(j+1);
        out << headers << "\n";
        // Then write the energy data.
        out.precision(12); // the precision of energy, default is 6
        for (int n = 0; n < nsteps/energy_steps+1; ++n) {
            out << DT*n*energy_steps; // this is time
            for (int j = 0; j < N; ++j)
                out << "," << KE[n][j]; // this is KE
            for (int i = 0; i < DOFe; ++i)
                for (int j = 0; j < N; ++j)
                    out << "," << PE[n][i*N+j]; // this is PE
            out << "\n";
        }
        out.close();
    }
}

void Simulation::runPerturbMD(int traj, int nsteps) {
    // Get current simulation step, which may not zero if it is a restart
    // simulation for classical MD simulation, while it is always zero for
    // class dynamics simulations.
    int step = Dy->getStep();
    report(traj, step); // report the initial step always,here the traj starts 1.
    if (nsteps == 0) // if total step is 0, compute energy only
        return;      // and skip the following dynamics
    // For simulation with QCDyn, the dynamics is taken by several silent
    // internal dynamics steps following the reporting of simulation data,
    // such as temeperature, energies, positions, velocities, and so on.
    // Note, at the step 0 (initial configuration) and last step, it will
    // call report() always. Get the report frequency (skip_steps) to stop the
    // propagtion to call report(). Do this once is enough in a multi-traj
    // nonadiabatic dynamics simultion.
    static const int steps_for_config = param.getDouble("steps_for_config");
    static const int LEN_TRAJ = param.getDouble("LEN_TRAJ");
    static const int LEN_SKIP = param.getDouble("LEN_SKIP");
    std::vector<int> perturb_parameters;
    SplitString(perturb_parameters, param.getStr("perturb_parameters"));
    static std::shared_ptr<HamiltonianForceFieldBase> ha = std::static_pointer_cast<HamiltonianForceFieldBase>(Ha);
    // This part is the noneq Pi_traj_plus result
    // Normal integration with sevral skip_steps unitil reach the nsteps.
    for (; step < nsteps ;) {
        if (step % steps_for_config == 0) {
            //remember the eq CONFIG ,R , V, F.set the same step in the Dynamics.
            ha->saveTheEqConfig();
            // This is the positive perturb MD simulation
            ha->setPositiveDirection();
            for (int i = 0; i < LEN_TRAJ; i++) {
                // This is to update the Pi tensor on the Observable.
                if (i % LEN_SKIP == 0) {
                    Obs->updateOneSnapshot();
                    // add one step to commit the Obs and ha
                    ha->addTheOneStep();
                }
                // Judge the MD perturb true or false
                ha->setFalseMdPerturb();
                if(i == perturb_parameters[0]) ha->setTrueMdPerturb();
                Dy->dynamics(1);
            }
            //set the eq CONFIG ,R , V, F.set the same step in the Dynamics.
            ha->setTheEqConfig();
            // This is the negtive perturb MD simulation
            ha->setNegativeDirection();
            for (int i = 0; i < LEN_TRAJ; i++) {
                if (i % LEN_SKIP == 0) {
                    // add one step to commit the Obs and ha
                    Obs->updateOneSnapshot();
                    ha->addTheOneStep();
                }
                ha->setFalseMdPerturb();
                if(i == perturb_parameters[0]) ha->setTrueMdPerturb();
                Dy->dynamics(1);
            }
            ha->setTheEqConfig();
        }
        ha->setFalseMdPerturb();
        Dy->dynamics(steps_for_config);
        step += steps_for_config;
        report(traj, step);
    }
    Obs->output(); 
}

void Simulation::runMD(int traj, int nsteps) {
    // Get current simulation step, which may not zero if it is a restart
    // simulation for classical MD simulation, while it is always zero for
    // class dynamics simulations.
    int step = Dy->getStep();
    report(traj, step); // report the initial step always,here the traj starts 1.
    if (nsteps == 0) // if total step is 0, compute energy only
        return;      // and skip the following dynamics
    // For simulation with QCDyn, the dynamics is taken by several silent
    // internal dynamics steps following the reporting of simulation data,
    // such as temeperature, energies, positions, velocities, and so on.
    // Note, at the step 0 (initial configuration) and last step, it will
    // call report() always. Get the report frequency (skip_steps) to stop the
    // propagtion to call report(). Do this once is enough in a multi-traj
    // nonadiabatic dynamics simultion.
    static const int skip_steps = getReportFrequency();
    // Perform a simulated annealing if required. This is valid for classical MD
    // simulation only, not work for nonadiabatic dynamics simulation.
    // Make sure the current step is an end of report interval before
    // enter the next normal intergation which do intergration with
    // several skip_steps. This is possible in a restart of classical
    // MD simulation, while not possible in the case of nonadiabatic simulation.
    // const int obs_skip_steps = RDM_steps;
    if (step % skip_steps != 0) {
        int run_steps = skip_steps - (step % skip_steps);
        // Plan B for NAMD
        if (dyn_type == "ReadTraj" ) {
            for ( ;step < run_steps; ) {
                Dy ->dynamics(RDM_steps);
                int traj_skip_index = Obs->updateOneSnapshot();
                // if there are a not-a-number, breaks; else the traj will normally propagates
                if (traj_skip_index == 0)  
                    step += RDM_steps;
                else 
                    break; 
            }
        }
        else {
            Dy->dynamics(run_steps);
            step += run_steps;
        }
        report(traj, step);
    }
    // Normal integration with sevral skip_steps unitil reach the nsteps.
    for (; step < nsteps ;) {
        if (dyn_type == "MD") {
            Dy->dynamics(skip_steps);
            if (polar_obs) {
                int traj_skip_index = Obs->updateOneSnapshot();
                if (traj_skip_index == 0)  
                    step += skip_steps;
                else 
                    break; 
            }
            else {
                step += skip_steps;
            }
            report(traj, step);
        }
        else { // Non-adiabatic MD
            // Dy->dynamics(skip_steps);
            Dy ->dynamics(RDM_steps);
            int traj_skip_index = Obs->updateOneSnapshot();
            // if there are a not-a-number, breaks; else the traj will normally propagates
            if (traj_skip_index == 0)  
                step += RDM_steps;
            else 
                break; 
            // step += skip_steps;
            if (step % skip_steps == 0)
                report(traj, step);
        }
    }
}

void Simulation::runOpenMM(int traj, int nsteps) {
    // Get current simulation step, which may not zero if it is a restart
    // simulation for classical MD simulation, while it is always zero for
    // nonadiabatic dynamics simulations.
    int step = Dy->getStep();
    report(traj, step); // report the initial step always,here the traj starts 1.
    if (nsteps == 0) // if total step is 0, compute energy only
        return;      // and skip the following dynamics
    // For simulation with OpenMM, the dynamics is taken by several silent
    // internal dynamics steps following the reporting of simulation data,
    // such as temeperature, energies, positions, velocities, and so on.
    // Note, at the step 0 (initial configuration) and last step, it will
    // call report() always. Get the report frequency (skip_steps) to stop the
    // propagtion to call report(). Do this once is enough in a multi-traj
    // nonadiabatic dynamics simultion.
    static const int skip_steps = getOpenMMReportFrequency();
    // Perform a simulated annealing if required. This is valid for classical MD
    // simulation only, not work for nonadiabatic dynamics simulation.
    if (param.getStr("annealing") != "none" && dyn_type == "OpenMM")
        runOpenMMAnnealing(step, skip_steps);
    // Make sure the current step is an end of report interval before
    // enter the next normal intergation which do intergration with
    // several skip_steps. This is possible in a restart of classical
    // MD simulation, while not possible in the case of nonadiabatic simulation.
    if (step % skip_steps != 0) {
        int run_steps = skip_steps - (step % skip_steps);
        // Plan B for NAMD        
        if (dyn_type != "OpenMM") {
            for ( ;step < run_steps; ) {
                Dy ->dynamics(RDM_steps);
                int traj_skip_index = Obs->updateOneSnapshot();
                // if there are a not-a-number, breaks; else the traj will normally propagates
                if (traj_skip_index == 0)  
                    step += RDM_steps;
                else 
                    break; 
            }
        }
        else { //Plan A for MD 
            Dy->dynamics(run_steps);
            step += run_steps;
        }
        report(traj, step);
    }
    // Normal integration with sevral skip_steps unitil reach the nsteps.
    for (; step < nsteps ;) {
        if (dyn_type == "OpenMM") {
            Dy->dynamics(skip_steps);
            step += skip_steps;
            report(traj, step);
        }
        else { // Non-adiabatic MD
            // Dy->dynamics(skip_steps);
            Dy ->dynamics(RDM_steps);
            int traj_skip_index = Obs->updateOneSnapshot();
            // if there are a not-a-number, breaks; else the traj will normally propagates
            if (traj_skip_index == 0)  
                step += RDM_steps;
            else 
                break; 
            // step += skip_steps;
            if (step % skip_steps == 0)
                report(traj, step);
        }
    }
}

//add by xu f is a const
/*
void Simulation::runPerturbOpenMM(int traj, int nsteps) {
    
    int step = Dy->getStep();
    report(traj, step); 
    
    if (nsteps == 0) 
        return;      
  
    static const int skip_steps = getOpenMMReportFrequency();
    
  
    std::vector<int> perturb_parameters;
    SplitString(perturb_parameters, param.getStr("perturb_parameters"));
    
   
    if (perturb_parameters.size() < 4) {
        throw std::runtime_error("ERROR: perturb_parameters should contain at least 4 values: step,atom,fx,fy,fz");
    }
    
    int perturbStep = perturb_parameters[0];
    int atomIndex = perturb_parameters[1];
    double fx = perturb_parameters[2];
    double fy = perturb_parameters[3];
    double fz = perturb_parameters.size() > 4 ? perturb_parameters[4] : 0.0;
    
    
    std::shared_ptr<DynamicsOpenMM> dynOpenMM = std::static_pointer_cast<DynamicsOpenMM>(Dy);
    
   
    if (step % skip_steps != 0) {
        int run_steps = skip_steps - (step % skip_steps);
        
      
        if (perturbStep >= step && perturbStep < step + run_steps) {
         
            dynOpenMM->perturbDynamics(run_steps, perturbStep - step, atomIndex, fx, fy, fz);
        } else {
        
            dynOpenMM->dynamics(run_steps);
        }
        
        step += run_steps;
        report(traj, step);
    }
  
    for (; step < nsteps;) {
        int stepsToRun = std::min(skip_steps, nsteps - step);
        
       
        if (perturbStep >= step && perturbStep < step + stepsToRun) {
       
            dynOpenMM->perturbDynamics(stepsToRun, perturbStep - step, atomIndex, fx, fy, fz);
        } else {
        
            dynOpenMM->dynamics(stepsToRun);
        }
        
        step += stepsToRun;
        report(traj, step);
    }
}
*/


void Simulation::runPerturbOpenMM(int traj, int nsteps) {
    int step = Dy->getStep();
    report(traj, step); 
    
    if (nsteps == 0) 
        return;      
  
    static const int skip_steps = getOpenMMReportFrequency();
    
    std::vector<int> perturb_parameters;
    SplitString(perturb_parameters, param.getStr("perturb_parameters"));
    
    if (perturb_parameters.size() < 2) {
        throw std::runtime_error("ERROR: perturb_parameters should contain at least 2 values: step,forceFieldIndex,scaleFactor");
    }
    
    int perturbStep = perturb_parameters[0];
    int forceFieldIndex = perturb_parameters[1];
    double scaleFactor = perturb_parameters.size() > 2 ? perturb_parameters[2] : 1.0;
    
    std::shared_ptr<DynamicsOpenMM> dynOpenMM = std::static_pointer_cast<DynamicsOpenMM>(Dy);
    
    if (step % skip_steps != 0) {
        int run_steps = skip_steps - (step % skip_steps);
        
        if (perturbStep >= step && perturbStep < step + run_steps) {
            
            dynOpenMM->perturbDynamics(run_steps, perturbStep - step, forceFieldIndex, scaleFactor);
        } else {
            dynOpenMM->dynamics(run_steps);
        }
        
        step += run_steps;
        report(traj, step);
    }
  
    for (; step < nsteps;) {
        int stepsToRun = std::min(skip_steps, nsteps - step);
        
        if (perturbStep >= step && perturbStep < step + stepsToRun) {
            
            dynOpenMM->perturbDynamics(stepsToRun, perturbStep - step, forceFieldIndex, scaleFactor);
        } else {
            dynOpenMM->dynamics(stepsToRun);
        }
        
        step += stepsToRun;
        report(traj, step);
    }
}










int Simulation::getOpenMMReportFrequency() {
    // Get the report frequency of each item to control simulaton.
    // If XXX_steps = 0, that is to say no output of this item.
    const int nsteps = param.getInt("nsteps");
    std::vector<int> report_frequency;
    report_frequency.push_back(param.getInt("energy_steps"));
    report_frequency.push_back(param.getInt("position_steps"));
    report_frequency.push_back(param.getInt("chk_steps"));
    // Check the validity of them:
    // 1) cannot < 0
    // 2) If they are larger than nsteps, then set them to nsteps.
    // 3) nsteps is a multiply of them
    for (auto& i : report_frequency) // use reference to write
        if (i < 0) // < 0 is not allowed
            throw std::runtime_error("ERROR: Found negative frequency in energy_steps or "
            " position_steps or chk_steps, please check your input.");
        else if (i > nsteps) // if > nsteps, reset to nsteps
            i = nsteps;
        else if (i != 0 && nsteps % i != 0) // nsteps is a multiply of them
            throw std::runtime_error("ERROR: The value of nsteps must be divided "
                "exactly by reprot frequency (energy_steps, position_steps, chk_steps).");
    // Get the skip_steps, the minimum value of them (non-zero), which means the
    // frequency to stop to call report(), also denotes the slient steps to
    // run dynamics in interal OpenMM integrator.
    std::sort(report_frequency.begin(), report_frequency.end()); // ascending
    int skip_steps = 0;
    for (auto i : report_frequency) // get the non-zero minimum value
        if (i == 0)
            continue;
        else {  // the minimum value is the first non-zero value
            skip_steps = i;
            break;
        }
    if (skip_steps == 0) // if all of them are zero, then set it to nsteps.
        skip_steps = nsteps;
    // Check if the skip_steps is the common divisor of them.
    for (auto i : report_frequency)
        if (i > skip_steps && (i % skip_steps != 0))
            throw std::runtime_error("ERROR: Please check the report frequecny "
                "(energy_steps, position_steps, chk_steps), the non-zero minimum "
                "value of them should be a common divisor of non-zero others.");
    // Get and check the progress_interval (in percentage) which decides the
    // output frequecny of simulation progress in the STDOUT.
    const int progress_interval = param.getInt("progress_interval");
    if (progress_interval < 0 || progress_interval > 100)
        throw std::runtime_error("ERROR: Illegal value of progress_interval: " + param.getStr("progress_interval"));
    return skip_steps;
}

void Simulation::runOpenMMAnnealing(int& step, int skip_steps) {
    const int nsteps = param.getInt("nsteps");
    const double DT  = param.getDouble("DT"); // time step size in ps
    const std::string annealing = param.getStr("annealing");
    if (annealing != "single") // TODO: periodic ?
        throw std::runtime_error("ERROR: Currently only single sequence annealing is supported.");
    if (dyn_type != "OpenMM")
        throw std::runtime_error("ERROR: The annealing is used for classical MD simulation with OpenMM only.");
    // Get smarter ponter to the HamiltonianOpenMM obejct
    std::shared_ptr<HamiltonianOpenMM> ha = std::static_pointer_cast<HamiltonianOpenMM>(Ha);
    // Get parameters to control the annealing process.
    // annealing_npoints denotes the number of annealing reference/control points used
    // annealing_time/step denotes the time/step at the annealing reference/control points
    // annealing_temperature denotes the temperatures at the annealing reference/control points
    // The number of annealing_step and annealing_temperature should be equal
    // to the value of annealing_npoints and each entry of them is separated
    // by a comma (,).
    const int annealing_npoints = param.getInt("annealing_npoints");
    if (annealing_npoints != 2) // TODO: more than 2
        throw std::runtime_error("ERROR: Currently only annealing_points = 2 is supported.");
    std::vector<double> annealing_time, annealing_temperature;
    SplitString(annealing_time, param.getStr("annealing_time"));
    SplitString(annealing_temperature, param.getStr("annealing_temperature"));
    // Make sure the values of them are valid.
    if (annealing_time.size() != annealing_npoints || annealing_temperature.size() != annealing_npoints)
        throw std::runtime_error("ERROR: The number of entries of annealing_time or "
            "annealing_temperature should equal to the number given in annealing_npoints");
    for (auto temp : annealing_temperature)
        if (temp < 0)
            throw std::runtime_error("ERROR: Negative temperature for annealing_temperature.");
    // Get step form time, round() will return the nearest integer value.
    std::vector<int> annealing_step(annealing_npoints, 0);
    for (int i = 0; i < annealing_npoints; ++i)
        annealing_step[i] = round(annealing_time[i]/DT);
    for (auto a : annealing_step)
        if (a < 0 || a > nsteps)
            throw std::runtime_error("ERROR: Negative time or larger than total time for annealing_time.");
    // In simulated annealing, the reference temperature is a piecewise
    // linear function. The actual annealing is performed by dynamically
    // changing the reference temperature used in the thermostat algorithm
    // selected, therefore we need to decided the change interval of
    // the reference temperature and the integration steps for one
    // specific reference temperature.
    // The interval of reference temperature to keep contant in integration
    // This is arbitrary, but should not be large. In Gromacs, I
    // found it is 0.1, but I think perhaps too small.
    double temperature_increment = param.getDouble("temperature_increment");
    if (temperature_increment <= 0 || temperature_increment >= annealing_temperature[1])
        throw std::runtime_error("ERROR: Illegal temperature_increment: " + param.getStr("temperature_increment"));
    double diff_temp = annealing_temperature[1] - annealing_temperature[0];
    if (diff_temp < 0) // < 0 means to decrease temperature
        temperature_increment *= -1;
    if (abs(diff_temp) < 1)
        throw std::runtime_error("ERROR: The difference of two annealing_temperature is too small.");
    // the number of segments of reference temperature to keep
    // + 1 since the initial temperature also is included
    // for example, 0 -> 10 K, the number of segments is 21.
    int nsegments = diff_temp / temperature_increment;
    int diff_step = annealing_step[1] - annealing_step[0];
    if (diff_step < 1)
        throw std::runtime_error("ERROR: The difference of two annealing_step is too small.");
    // the number of steps at one segment reference temperature to
    // propagate by integrator (sometimes not divisible)
    int segment_nsteps = diff_step / nsegments;
    // the remaing number of steps if not divisible
    int left_nsteps = annealing_step[1] - segment_nsteps * nsegments;
    // Start to do integration step by step.
    // The following procedure also works for restart (current step is not zero)
    // Perhap it is not to do anealing from step 0, in this case
    // run several steps unitil reach the first annealing_step.
    for (; step < annealing_step[0]; ) {
        Dy->dynamics(1);
        step++;
        if (step % skip_steps == 0)
            report(0, step);
    } // now current step equal to annealing_step[0]
    // In the following, change the reference temperature of heat
    // bath, if reaching the segment_nsteps of the previous temperature.
    for (; step < annealing_step[1]; ) {
        int shift_step = step - annealing_step[0]; // the start step of one segment
        if (shift_step % segment_nsteps == 0) {
            // index of segment, 1, 2, ..., nsegments.
            int index = shift_step/segment_nsteps + 1;
            // reference temperature of heat bath for current segment
            double temp = annealing_temperature[0]+temperature_increment*index;
            // When this is the last segment, we use the second reference
            // temperature point directly, since sometimes diff_temp/temperature_increment
            // is not divisible. And it will keep this temperature until
            // reach annealing_step[1].
            if (index == nsegments)
                temp =  annealing_temperature[1];
            ha->setTemperature(temp);
        }
        Dy->dynamics(1);
        step++;
        if (step % skip_steps == 0)
            report(0, step);
    } // now current step equal to annealing_step[1]
}

void Simulation::reportProgress(int traj, int step) {
    static int count = 0; // record the number of times to eneter here
    static const int ntraj  = param.getInt("ntraj");
    static const int nsteps = param.getInt("nsteps");
    // It is used to record the elpased time during the progress interval.
    static std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
    end = std::chrono::high_resolution_clock::now();
    // The elapsedTime in seconds.
    std::chrono::duration<double> elapsedTime = count > 0 ? (end-start) : (end-end);
    start = end;
    double progress = 0.0;
    // TODO: add ellapsed time, remaing time, performance
    // TODO: and they can use the same one.
    // -- Nonadiabatic dynamics simulation for model system
    if (system_type == "model") {
        if (nsteps == 0)
            progress = 100 * (double)(traj) / (double)(ntraj);
        else
            progress = 100 * (double)((traj-1) * nsteps + step) / (double)(ntraj * nsteps);
        //printf("    step: %10d/%-10d traj: %10d/%-10d progress: %5.1f%\n", step, nsteps, traj, ntraj, progress);
        // Here, the setw is one time manipulation, while others are valid always.
        std::cout << "    step: " << std::right << std::setw(10) << step << "/" <<
            std::left << std::setw(10) << nsteps << " traj: " << std::right <<
            std::setw(10) << traj << "/" <<  std::left << std::setw(10)  << ntraj <<
            " progress: " << std::fixed << std::setprecision(1) << std::right <<
            std::setw(5) << progress << "%" << std::endl;
    }
    // -- Clssical MD simulation with OpenMM and QCDyn (ntraj is always 1)
    // -- Nonadiabatic dynamics simulation for realistic system with OpenMM
    // -- custom force field with QCDyn
    // Here, step starts from 0, step 0 is the initial configuration.
    // Here, when count = 0, the elapsed time is zero and not print.
    else if (system_type == "allatom" && count > 0) 
        if (allatom_type == "OpenMM" || allatom_type == "CustomForceField") {
        if (nsteps == 0)
            progress = 100 * (double)(traj) / (double)(ntraj);
        else
            progress = 100 * (double)((traj-1) * nsteps + step) / (double)(ntraj * nsteps);
        // Here, the setw is one time manipulation, while others are valid always.
        std::cout << "    step: " << std::right << std::setw(8) << step << "/" <<
            std::left << std::setw(8) << nsteps << " traj: " << std::right <<
            std::setw(6) << traj << "/" <<  std::left << std::setw(6)  << ntraj <<
            " progress: " << std::fixed << std::setprecision(1) << std::right <<
            std::setw(5) << progress << "%" <<  "   elapsed time: " << std::setprecision(2) <<
            elapsedTime.count() << " seconds" << std::endl;
    }
    count++;
}

void Simulation::reportOpenMMData(int traj, int step) {
    static int count = 0; // record the number of times to eneter here
    static int traj_prev = traj; // the traj of previous call
    if (traj_prev != traj) { // if it is not a same traj,
        count = 0;           // we reset the count to 0,
        traj_prev = traj;    // and set current traj as previous one.
    }
    // Get smarter ponter to the HamiltonianOpenMM obejct
    static std::shared_ptr<HamiltonianOpenMM> ha = std::static_pointer_cast<HamiltonianOpenMM>(Ha);
    // Get the total number of trajectory (for classical MD, ntraj is always 1)
    // and steps and the frequency (number of steps) to report.
    // Note, the data of step 0 (initial configuration) and last step will be
    // reported all the time unless 0, which means the file won't be genenrated.
    static const int                ntraj = param.getInt("ntraj");
    static const int               nsteps = param.getInt("nsteps");
    static const int         energy_steps = param.getInt("energy_steps");
    static const int       position_steps = param.getInt("position_steps");
    static const int            chk_steps = param.getInt("chk_steps");
    // Get current MD simulation time.
    static const double                DT = param.getDouble("DT");
    double                          time  = step * DT;
    // If default_name is specified, it will be set to the file name of output
    static const std::string default_name = param.getStr("default_name");
    static const bool useDefault = (default_name.empty() || default_name == "none") ? false : true;
    // Get the progress_steps to decide report simulation progress at which step
    // Here, progress_interval (in percentage) decides the output frequecny of
    // simulation progress in the STDOUT.
    // Note, for multi-traj nonadiabatic simulation, the ntraj can be > 1.
    static const int progress_interval = param.getInt("progress_interval");
    static const int progress_steps    = ntraj * nsteps * progress_interval / 100;
    // 1. Report simulation progress
    // If progress_steps is 0, the simulation progress will be reported always.
    // And at the initial step and final steps, it will be reported always.
    // Here, ((traj-1) * nsteps + step) is the current step in the all steps of
    // multi-traj simulation. For classical MD simulation, it is just step.
    if (count == 0 || progress_steps == 0 || ((traj-1) * nsteps + step) % progress_steps == 0 || step == nsteps)
        // when traj > 1, the step 0 of current traj and last step of las traj is repetitive
        // which should be aviod to report the step 0 in this case.
        if (traj == 1 || (traj > 1 && step != 0))
            reportProgress(traj, step);
    // 2. Report energies, and other data, such as temeprature
    // For multi-traj simultaion (ntraj > 1), the data filename of each trajectory
    // will add a suffix "_traj?", where "?" is the index of traj (starting from 1).
    if (energy_steps > 0 && (step % energy_steps == 0 || step == nsteps)) {
        static const std::string data_file = useDefault ? (default_name + ".csv") : param.getStr("data_file");
        std::string filename = ntraj == 1 ? data_file : GetFilePrefix(data_file) + "_traj" + std::to_string(traj) + "." + GetFileSuffix(data_file);
        static const bool energy_decompose = param.getBool("energy_decompose");
        static std::vector<double> data;
        getStateData(traj, data, false, energy_decompose);
        writeStateData(filename, data, step, time);
    }
    // 3. Report trajectory: coordinates, (optional velocities and forces)
    // TODO: using structure object in Hamiltonina directly
    // TODO with function updateStructureData(includeForces)
    if (position_steps > 0 && (step % position_steps == 0 || step == nsteps)) {
        static int frame = 0; // record the frame in trajectory file
        static int traj_prev2 = traj; // the traj of previous call
        if (traj_prev2 != traj) { // if it is not a same traj,
            frame = 0;            // we reset the count to 0,
            traj_prev2 = traj;    // and set current traj as previous one.
        }
        static std::string traj_file, traj_type;
        static const std::string traj_group = param.getStr("traj_group");
        static const int      xtc_precision = param.getInt("xtc_precision");
        static const bool     velocity_save = param.getBool("velocity_save");
        static const bool        force_save = param.getBool("force_save");
        static const int    propagate_state = param.getInt("propagate_state");
        static const bool useBarostat = param.getStr("barostat") != "none" ? true : false;
        static const bool usePBC = param.getStr("PBC") != "none" ? true : false;
        // decide the trajectory file name and check the type of trajectory at
        // first time if default_name is not used.

        if (frame == 0 && traj == 1) { // do this only at first time for first trajectory
            if (useDefault) {
                if (velocity_save || force_save)
                    traj_file = default_name + ".trr";
                else
                    traj_file = default_name + ".xtc";
                traj_type = GetFileSuffix(traj_file);
            }
            else {
                traj_file = param.getStr("traj_file");
                traj_type = GetFileSuffix(traj_file);
                if (traj_type == "xtc" && (velocity_save || force_save))
                    throw std::runtime_error("ERROR: xtc trajectory file doesn't support "
                        "to store velocities or forces, please use trr file.");
                if (traj_type == "gro" && force_save)
                    throw std::runtime_error("ERROR: gro trajectory file doesn't support "
                        "to store forces, please use trr file.");
            }
        }
        // For multi-traj simulation, add a suffix "_traj?" to filename, where "?" is the index of traj (starting from 1).
        std::string filename = ntraj == 1 ? traj_file : GetFilePrefix(traj_file) + "_traj" + std::to_string(traj) + "." + GetFileSuffix(traj_file);
        // Create a Structure object, and get the atominfo, positions,
        // velocities, boxVectors of current step from HamiltonianOpenMM
        // object, then use method in Structure to save trajectory file.
        static Structure structure; // empty Structure object
        // Get writable reference to structure data member
        static std::vector<std::string>&  atominfo   = structure.getAtomInfo();
        static std::vector<OpenMM::Vec3>& positions  = structure.getPositions();
        static std::vector<OpenMM::Vec3>& velocities = structure.getVelocities();
        static std::vector<OpenMM::Vec3>& forces     = structure.getForces();
        static OpenMM::Vec3         (&boxVectors)[3] = structure.getBoxVectors();
        // test for print force data for first atom. March 14 2025
      //  std::cout<<"This is step: "<<step<<std::endl;

        //std::cout << "First force: ("
       //           << forces[0] << ", "  // x
       //           << forces[1] << ", "  // y
       //           << forces[2] << ")"    // z
       //           << std::endl;


        // Set the atominfo and natoms, and never change
        if (frame == 0 && traj == 1) { // do this at first time for first trajectory
            atominfo = ha->getAtomInfo();
            structure.setNumAtoms(atominfo.size());
        }
        // boxVectors changes when barostat is used or for a new trajectory
        // For a non-PBC simulation, the boxVectors are always zero.
        if (usePBC && (frame == 0 || useBarostat))
            ha->getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        // Get positions and always needed
        ha->getPositions(positions);
        // If a subsystem is required to save to tajectory file, and the trajectory
        // format is not gro file. The initial stucture of this subsystem will be
        // saved to a Gromacs structure file with a filename with suffix "_subsystem".
        // This subsystem structure file is necessary when you want to visualize the
        // subsystem binary trajectory file with suah as VMD, since in binary
        // trajectory file, the atom infomation is not included.
        if (frame == 0 && traj == 1 && !traj_group.empty() && traj_group != "system" && traj_type != "gro")
            structure.saveTrajectory(GetFilePrefix(traj_file) + "_subsystem.gro", step, time, traj_group);
        if (velocity_save)
            ha->getVelocities(velocities);
        if (force_save) {
            // For classical M, the forces of current positions should be
            // recalculated. And only the state used to propagate is saved.
            // TODO: How to aviod this calculation for classical MD
            if (dyn_type == "OpenMM") { // classical MD, recompute the forces
                ha->getPotentialEnergy(propagate_state, true);
                ha->getForces(forces);
            }
            else { // nonadiabatic dynamics, copy effective forces from Hamiltonian
                const HamiltonianOpenMM& omm = *ha; // get a constant reference
                forces = omm.getForces(); // return the data member F
            }
            // Print the forces of the first atom
            if (!forces.empty()) { // Ensure there are forces to print
                const OpenMM::Vec3& firstAtomForce = forces[0];
                std::cout << "Force on the first atom: ("
                          << firstAtomForce[0] << ", "
                          << firstAtomForce[1] << ", "
                          << firstAtomForce[2] << ")"
                          << std::endl;
            }


        }
        // Save current sturcture data to trajectory file
        structure.saveTrajectory(filename, step, time, traj_group, xtc_precision);
        frame++;
    }
    // 4. Create checkpoint file (the initial step is skipped)
    if (chk_steps > 0 && (step % chk_steps == 0 || step == nsteps)) {
        static const std::string chk_file = useDefault ? (default_name + ".chk") : param.getStr("chk_file");
        // For multi-traj simulation, add a suffix "_traj?" to filename, where "?" is the index of traj (starting from 1).
        std::string filename = ntraj == 1 ? chk_file : GetFilePrefix(chk_file) + "_traj" + std::to_string(traj) + "." + GetFileSuffix(chk_file);
        // At the initial step, the check point file will be saved with a suufix
        // "_init", which can be used to restart a simulation in the future.
        if (count == 0)
            ha->createCheckpoint(GetFilePrefix(filename) + "_init.chk");
        else
            ha->createCheckpoint(filename);
    }
    // 5. Save structure file at last step, containing positions and velocities
    // The structure file will always save whole system (traj_group doesn't
    // infulence it) at last step.
    if (step == nsteps && param.getStr("conf_file") != "none") {
        // TODO: using structure object in Hamiltonian directly
        // Create a structure object and get data from Hamiltonian obejct
        Structure structure;
        std::vector<std::string>&  atominfo   = structure.getAtomInfo();
        std::vector<OpenMM::Vec3>& positions  = structure.getPositions();
        std::vector<OpenMM::Vec3>& velocities = structure.getVelocities();
        OpenMM::Vec3         (&boxVectors)[3] = structure.getBoxVectors();
        atominfo = ha->getAtomInfo();
        structure.setNumAtoms(atominfo.size());
        ha->getPositions(positions);
        ha->getVelocities(velocities);
        if (param.getStr("PBC") != "none")
            ha->getPeriodicBoxVectors(boxVectors[0], boxVectors[1], boxVectors[2]);
        // If same name as this trajectory file is used. Then the filename
        // with a suffix "_last" will be used.
        std::string conf_file = useDefault ? (default_name + ".gro") : param.getStr("conf_file");
        if (conf_file == param.getStr("traj_file"))
            conf_file = GetFilePrefix(conf_file) + "_last." + GetFileSuffix(conf_file);
        // For multi-traj simulation, add a suffix "_traj?" to filename, where "?" is the index of traj (starting from 1).
        if (ntraj > 1)
            conf_file = GetFilePrefix(conf_file) + "_traj" + std::to_string(traj) + "." + GetFileSuffix(conf_file);
        structure.saveStructure(conf_file, step, time);
    }
    count++;
}

void Simulation::getMDData(int traj, std::vector<double>& data, bool includeForces, bool energy_decompose) {
    // Get smarter ponter to the HamiltonianForceFieldBase obejct
    static std::shared_ptr<HamiltonianForceFieldBase> ha = std::static_pointer_cast<HamiltonianForceFieldBase>(Ha);
    static int count = 0;
    static int traj_prev = traj; // the traj of previous call
    if (traj_prev != traj) { // if it is not a same traj,
        count = 0;           // we reset the count to 0,
        traj_prev = traj;    // and set current traj as previous one.
    }
    static const int nvalues = 13; // currently, 13 quantities/values for each state
    static const int DOFe = param.getInt("DOFe");
    static const double systemMass = ha->getSystemMass()*amu2kg; // mass in kg
    static const int systemDOF = ha->getSystemDOF();
    static const int propagate_state = param.getInt("propagate_state");
    static const bool useBarostat = param.getStr("barostat") != "none" ? true : false;
    static const bool usePBC = param.getStr("PBC") != "none" ? true : false;
    static const bool needForces = true;
    static double volume = 0, density = 0, kineticEnergy, temperature;
    static double potentialEnergy = 0; // the energy of propagate_state
    // If this is a non-PBC simulation, the volume and density is 0.
    // Only when using barostat, the volume/density changes during simulation.
    if (usePBC && (count == 0 || useBarostat)) {
        volume = ha->getPeriodicBoxVolume(0); // Get volume (in nm^3)
        density = systemMass * 1e27 / volume; // Get density (in kg/m^3)
    }
    // Get kinetic energy (in kj/mol)
    // For some integrators, such as leap-frog, Langevin, the current velocities
    // is half step "delayed", and we need a pre-computed forces based on current
    // R to get shift velocities at time t to compute kinetic energy. The re-calculation
    // of force is necessary, since now the froces in Context is the previous step.
    // This is the case of classical MD simulation with OpenMM. While in the case
    // of nonadiabatic dynamics simulation, we should upload the effective forces
    // of current step (the data F in Hamiltonian object) to OpenMM Context before
    // to compute kinetic energy.
    // TODO: How to avoid the recomputation of foreces when saving forces to
    // trajectory for classical MD simulation with OpenMM.
    kineticEnergy   = ha->getKineticEnergy();
    // Get instantaneous temperature (in K)
    // temperature is computed by the formula: T = 2 * Ek / DOF / Rbar,
    // Rbar is gas constant, Rbar = 8.31446261815324e-3 kj/mol/K.
    //temperature = 2.0 * kineticEnergy / systemDOF / Rbar;
    temperature = 2.0 * kineticEnergy / systemDOF / Rbar; 
    // The data in vector should be (energies in kj/mol):
    // data[0]: Temperature (in K), data[1]: Volume (in nm^3),
    // data[2]: Density (in kg/m^3), data[3]: TotalEnergy, data[4]: KineticEnergy,
    // data[5-12]: total potential energy of state 0, and the energy components
    // of potential energy: Nonbond, Bond, Angle, Dihedral, Others, Coulomb, vdW.
    // Here, the Others means from Custom Force such as potsition restraint.
    // data[13-25]: same as above, but for second state (state 1).
    // ...
    // if energy_decompose = true, then energy components of potential energy
    // will be computed, otherwise, they are zero.
    data.resize(nvalues * DOFe, 0.0);
    std::fill(data.begin(), data.end(), 0); // reset to 0 for existing values (if any)
    // Note: the data of Temperature, Volume, Density, KineticEnergy,
    // are the same for all states.
    for (int i = 0; i < DOFe; i++) { // loop from state 1
        data[0 + i * nvalues] = temperature; // data[0]: Temperature (in K)
        data[1 + i * nvalues] = volume; // data[1]: Volume (in nm^3)
        data[2 + i * nvalues] = density; // data[2 ]: Density (in kg/m^3)
        data[4 + i * nvalues] = kineticEnergy; // data[3]: KineticEnergy (in kj/mol)
    }
    // Get potential energy (in kj/mol) of each state
    if (!energy_decompose) { 
        // for classical MD simulation, we should compute the energy.
        if (dyn_type == "MD") {
            for (int i = 0; i < DOFe; i++) {
                // The potential energy of propagate_state has been computed
                // when computing forces to compute shift kinetic energy.
                if (i == propagate_state && needForces)
                    data[5 + i * nvalues] = potentialEnergy;
                else
                    data[5 + i * nvalues] = ha->getPotentialEnergy(i);
        
                } 
        }
        // for nonadiabatic dyanmics, the potential energy of each state have
        // been computed when update Hamiltonian matrix and stored in the data
        // member PE (in kj/mol) in Hamiltonian object.
        else {
            const std::vector<double>& PE = ha->getPE();
            for (int i = 0; i < DOFe; i++)
                data[5 + i * nvalues] = PE[i];
        }
    }
    // Get enegry components of each state. In this case, both for classical MD
    // and nonadibatic dynamics, the energy should be computed.
    else {
        for (int i = 0; i < DOFe; i++) {
            // The potential energy of propagate_state has been computed
            // when computing forces to compute shift kinetic energy.
            std::vector<double> energygroups;
            energygroups.resize(7, 0);
            energygroups = ha->getEnergyGroups(i);
            data[6 + i * nvalues] = energygroups[4] + energygroups[5]; // data[6]: Non-Bond (in kj/mol)
            data[7 + i * nvalues] = energygroups[1]; // data[7]: BondEnergy (in kj/mol)
            data[8 + i * nvalues] = energygroups[2]; // data[8]: AngleEnergy (in kj/mol)
            data[9 + i * nvalues] = energygroups[3]; // data[9]: DihedralEnergy (in kj/mol)
            data[10 + i * nvalues] = energygroups[5]; // data[10 ]: VdwEnergy (in kj/mol)
            data[11 + i * nvalues] = energygroups[4]; // data[11]: ElecEnergy (in kj/mol)  
            data[5 + i * nvalues] =  energygroups[1] + energygroups[2] + 
                energygroups[3] + energygroups[4] + energygroups[5];
        } 
    }
    // Get total energy (data[3]) of each state: sum of kinetic energy data[4]
    // and potential energy data[5].
    for (int i = 0; i < DOFe; i++)
        data[3 + i * nvalues] = data[4 + i * nvalues] + data[5 + i * nvalues];
    count++;
}

void Simulation::getRPMDData(int traj, std::vector<double>& data, bool includeForces, bool energy_decompose) {
    // Get smarter ponter to the HamiltonianForceFieldBase obejct
    static std::shared_ptr<HamiltonianRingPolymer> ha = std::static_pointer_cast<HamiltonianRingPolymer>(Ha);
    static int count = 0;
    static int traj_prev = traj; // the traj of previous call
    if (traj_prev != traj) { // if it is not a same traj,
        count = 0;           // we reset the count to 0,
        traj_prev = traj;    // and set current traj as previous one.
    }
    static const int nvalues = 13; // currently, 13 quantities/values for each state
    static const int DOFe = param.getInt("DOFe");
    static const double systemMass = ha->getSystemMass()*amu2kg; // mass in kg
    static const int systemDOF = ha->getSystemDOF();
    static const int propagate_state = param.getInt("propagate_state");
    static const bool useBarostat = param.getStr("barostat") != "none" ? true : false;
    static const bool usePBC = param.getStr("PBC") != "none" ? true : false;
    static const bool needForces = true;
    static double volume = 0, density = 0, kineticEnergy, temperature;
    static double potentialEnergy = 0; // the energy of propagate_state
    // If this is a non-PBC simulation, the volume and density is 0.
    // Only when using barostat, the volume/density changes during simulation.
    if (usePBC && (count == 0 || useBarostat)) {
        volume = ha->getPeriodicBoxVolume(0); // Get volume (in nm^3)
        density = systemMass * 1e27 / volume; // Get density (in kg/m^3)
    }
    // Get kinetic energy (in kj/mol)
    // For some integrators, such as leap-frog, Langevin, the current velocities
    // is half step "delayed", and we need a pre-computed forces based on current
    // R to get shift velocities at time t to compute kinetic energy. The re-calculation
    // of force is necessary, since now the froces in Context is the previous step.
    // This is the case of classical MD simulation with OpenMM. While in the case
    // of nonadiabatic dynamics simulation, we should upload the effective forces
    // of current step (the data F in Hamiltonian object) to OpenMM Context before
    // to compute kinetic energy.
    // TODO: How to avoid the recomputation of foreces when saving forces to
    // trajectory for classical MD simulation with OpenMM.
    kineticEnergy   = ha->getKineticEnergy();
    // Get instantaneous temperature (in K)
    // temperature is computed by the formula: T = 2 * Ek / DOF / Rbar,
    // Rbar is gas constant, Rbar = 8.31446261815324e-3 kj/mol/K.
    //temperature = 2.0 * kineticEnergy / systemDOF / Rbar;
    temperature = 2.0 * kineticEnergy / systemDOF / Rbar; 
    // The data in vector should be (energies in kj/mol):
    // data[0]: Temperature (in K), data[1]: Volume (in nm^3),
    // data[2]: Density (in kg/m^3), data[3]: TotalEnergy, data[4]: KineticEnergy,
    // data[5-12]: total potential energy of state 0, and the energy components
    // of potential energy: Nonbond, Bond, Angle, Dihedral, Others, Coulomb, vdW.
    // Here, the Others means from Custom Force such as potsition restraint.
    // data[13-25]: same as above, but for second state (state 1).
    // ...
    // if energy_decompose = true, then energy components of potential energy
    // will be computed, otherwise, they are zero.
    data.resize(nvalues * DOFe, 0.0);
    std::fill(data.begin(), data.end(), 0); // reset to 0 for existing values (if any)
    // Note: the data of Temperature, Volume, Density, KineticEnergy,
    // are the same for all states.
    for (int i = 0; i < DOFe; i++) { // loop from state 1
        data[0 + i * nvalues] = temperature; // data[0]: Temperature (in K)
        data[1 + i * nvalues] = volume; // data[1]: Volume (in nm^3)
        data[2 + i * nvalues] = density; // data[2 ]: Density (in kg/m^3)
        data[4 + i * nvalues] = kineticEnergy; // data[3]: KineticEnergy (in kj/mol)
        data[5 + i * nvalues] = ha->getPotentialEnergy(i); // data[4]: PotentialEnergy (in kj/mol)
    }
    if (energy_decompose) { 
        for (int i = 0; i < DOFe; i++) {
            std::vector<double> energygroups;
            energygroups.resize(7, 0);
            energygroups = ha->getEnergyGroups(i);
            data[6 + i * nvalues] = energygroups[4] + energygroups[5]; // data[6]: Non-Bond (in kj/mol)
            data[7 + i * nvalues] = energygroups[1]; // data[7]: BondEnergy (in kj/mol)
            data[8 + i * nvalues] = energygroups[2]; // data[8]: AngleEnergy (in kj/mol)
            data[9 + i * nvalues] = energygroups[3]; // data[9]: DihedralEnergy (in kj/mol)
            data[10 + i * nvalues] = energygroups[5]; // data[10 ]: VdwEnergy (in kj/mol)
            data[11 + i * nvalues] = energygroups[4]; // data[11]: ElecEnergy (in kj/mol)   
        } 
    }
    // Get total energy (data[3]) of each state: sum of kinetic energy data[4]
    // and potential energy data[5].
    for (int i = 0; i < DOFe; i++)
        data[3 + i * nvalues] = data[4 + i * nvalues] + data[5 + i * nvalues];
    count++;
}

void Simulation::writeMDData(const std::string& file, std::vector<double>& data, int step, double time) {
    // * Get required parameters.
    static int count = 0;
    static std::string file_prev = file; // the filename of previous call
    if (file_prev != file) { // if it is not a same file,
        count = 0;           // we reset the count to 0,
        file_prev = file;    // and set current file name as previous one.
    }
    static const bool energy_decompose = param.getBool("energy_decompose");
    static const int nvalues = 13; // 13 values for each state in data
    const int nstates = data.size()/nvalues; // number of states
    const std::string type = GetFileSuffix(file); // file format type
    // * Save current state data to file
    // The data in vector are (energies in kj/mol):
    // data[0]: Temperature (in K), data[1]: Volume (in nm^3),
    // data[2]: Density (in kg/m^3), data[3]: TotalEnergy, data[4]: KineticEnergy,
    // data[5-12]: total potential energy of state 0, and the energy components
    // of potential energy: Nonbond, Bond, Angle, Dihedral, Others, Coulomb, vdW.
    // Here, the Others means from Custom Force such as potsition restraint.
    // data[13-25]: same as above, but for second state (state 1).
    // data[26-38]: same as above, but for third state (state 2).
    FILE* dataFile = CheckFile(file, count);
    if (type == "dat") { // * Using fixed space separated format (.dat)
        if (count == 0) { // Print headers and units at the first time.
            fprintf(dataFile, "%10s %12s %5s %12s %12s %12s %16s %16s %16s", "Step", "Time", "State",
                "Temperature", "Volume", "Density","TotalEnergy", "KineticEnergy",  "PotentialEnergy");
           if (energy_decompose)
                fprintf(dataFile, "%16s %16s %16s %16s %16s %16s", "Nonbonded_E" ,"Bond_E",
                    "Angle_E", "Dihedral_E",  "Elec_E","LJ_E");
            fprintf(dataFile, "\n");
            fprintf(dataFile, "%10s %12s %5s %12s %12s %12s %16s %16s %16s", " ", "(ps)",
                " ", "(K)", "(nm^3)", "(kg/m^3)", "(kJ/mol)", "(kJ/mol)", "(kJ/mol)");
            if (energy_decompose)
                fprintf(dataFile, "%16s %16s %16s %16s %16s %16s", "(kJ/mol)",
                    "(kJ/mol)", "(kJ/mol)", "(kJ/mol)", "(kJ/mol)", "(kJ/mol)");
            fprintf(dataFile, "\n");
        }
        for (int state = 0; state < nstates; ++state) { // Write data of each state.
            fprintf(dataFile, "%10d %12.4f %5d %12.4f %12.4f %12.4f %16.6f %16.6f",
                step, time, state, data[0+state*nvalues], data[1+state*nvalues], data[2+state*nvalues],
                data[3+state*nvalues], data[4+state*nvalues]);
            if (energy_decompose) // energy components: data[6-12]
                for (int i = 6; i < 13; ++i)
                    fprintf(dataFile, " %16.6f", data[i + state*nvalues]);
            fprintf(dataFile, "\n");
        }
    }
    else if (type == "csv") { // * comma (,) separated value (CSV) file
        if (count == 0) { // Print headers at the first time.
            fprintf(dataFile, "Step,Time,State,Temperature,Volume,Density,TotalEnergy, KineticEnergy,PotentialEnergy");
            if (energy_decompose)
                fprintf(dataFile, ",Nonbonded_E, Bond_E, Angle_E, Dihedral_E, Elec_E, LJ_E ");
            fprintf(dataFile, "\n");
        }
        for (int state = 0; state < nstates; ++state) { // Write data of each state
            fprintf(dataFile, "%d,%g,%d,%g,%g,%g,%.12g,%.12g,%.12g", step, time, state,
                data[0+state*nvalues], data[1+state*nvalues], data[2+state*nvalues],
                data[3+state*nvalues], data[4+state*nvalues], data[5+state*nvalues]);
            if (energy_decompose) // energy components: data[6-12]
                for (int i = 6; i < 12; ++i)
                    fprintf(dataFile, ",%.12g", data[i + state*nvalues]);
            fprintf(dataFile, "\n");
        }
    }
    else
        throw std::runtime_error("ERROR: Unknown file type for output: " + type);
    fclose(dataFile);
    count++;
}

void Simulation::getStateData(int traj, std::vector<double>& data, bool includeForces, bool energy_decompose) {
    // Get smarter ponter to the HamiltonianOpenMM obejct
    static std::shared_ptr<HamiltonianOpenMM> ha = std::static_pointer_cast<HamiltonianOpenMM>(Ha);
    static int count = 0;
    static int traj_prev = traj; // the traj of previous call
    if (traj_prev != traj) { // if it is not a same traj,
        count = 0;           // we reset the count to 0,
        traj_prev = traj;    // and set current traj as previous one.
    }
    static const int nvalues = 13; // currently, 13 quantities/values for each state
    static const int DOFe = ha->getNumStates();
    static const int systemDOF = ha->getSystemDOF();
    static const double systemMass = ha->getSystemMass()*amu2kg; // mass in kg
    static const int propagate_state = param.getInt("propagate_state");
    static const bool useBarostat = param.getStr("barostat") != "none" ? true : false;
    static const bool usePBC = param.getStr("PBC") != "none" ? true : false;
    static const bool needForces = ha->kineticEnergyRequiresForce();
    static double volume = 0, density = 0, kineticEnergy, temperature;
    static double potentialEnergy = 0; // the energy of propagate_state
    // If this is a non-PBC simulation, the volume and density is 0.
    // Only when using barostat, the volume/density changes during simulation.
    if (usePBC && (count == 0 || useBarostat)) {
        volume = ha->getPeriodicBoxVolume(); // Get volume (in nm^3)
        density = systemMass * 1e27 / volume; // Get density (in kg/m^3)
    }
    // Get kinetic energy (in kj/mol)
    // For some integrators, such as leap-frog, Langevin, the current velocities
    // is half step "delayed", and we need a pre-computed forces based on current
    // R to get shift velocities at time t to compute kinetic energy. The re-calculation
    // of force is necessary, since now the froces in Context is the previous step.
    // This is the case of classical MD simulation with OpenMM. While in the case
    // of nonadiabatic dynamics simulation, we should upload the effective forces
    // of current step (the data F in Hamiltonian object) to OpenMM Context before
    // to compute kinetic energy.
    // TODO: How to avoid the recomputation of foreces when saving forces to
    // trajectory for classical MD simulation with OpenMM.
    if (needForces)
        if (dyn_type == "OpenMM") // this will update the forces in Context
            potentialEnergy = ha->getPotentialEnergy(propagate_state, needForces);
        else // For nonadiabatic dynamics, upload the effective forces (F) to Context
            ha->uploadForces();
    kineticEnergy = ha->getKineticEnergy();
    // Get instantaneous temperature (in K)
    // temperature is computed by the formula: T = 2 * Ek / DOF / Rbar,
    // Rbar is gas constant, Rbar = 8.31446261815324e-3 kj/mol/K.
    temperature = 2.0 * kineticEnergy / systemDOF / Rbar;
    // The data in vector should be (energies in kj/mol):
    // data[0]: Temperature (in K), data[1]: Volume (in nm^3),
    // data[2]: Density (in kg/m^3), data[3]: TotalEnergy, data[4]: KineticEnergy,
    // data[5-12]: total potential energy of state 0, and the energy components
    // of potential energy: Nonbond, Bond, Angle, Dihedral, Others, Coulomb, vdW.
    // Here, the Others means from Custom Force such as potsition restraint.
    // data[13-25]: same as above, but for second state (state 1).
    // ...
    // if energy_decompose = true, then energy components of potential energy
    // will be computed, otherwise, they are zero.
    data.resize(nvalues * DOFe, 0.0);
    std::fill(data.begin(), data.end(), 0); // reset to 0 for existing values (if any)
    // Note: the data of Temperature, Volume, Density, KineticEnergy,
    // are the same for all states.
    for (int i = 0; i < DOFe; i++) { // loop from state 1
        data[0 + i * nvalues] = temperature; // data[0]: Temperature (in K)
        data[1 + i * nvalues] = volume; // data[1]: Volume (in nm^3)
        data[2 + i * nvalues] = density; // data[2]: Density (in kg/m^3)
        data[4 + i * nvalues] = kineticEnergy; // data[4]: KineticEnergy (in kj/mol)
    }
    // Get potential energy (in kj/mol) of each state
    if (!energy_decompose) { // get the total potential energy directly
        // for classical MD simulation, we should compute the energy.
        if (dyn_type == "OpenMM") {
            for (int i = 0; i < DOFe; i++)
                // The potential energy of propagate_state has been computed
                // when computing forces to compute shift kinetic energy.
                if (i == propagate_state && needForces)
                    data[5 + i * nvalues] = potentialEnergy;
                else
                    data[5 + i * nvalues] = ha->getPotentialEnergy(i);
        }
        // for nonadiabatic dyanmics, the potential energy of each state have
        // been computed when update Hamiltonian matrix and stored in the data
        // member PE (in kj/mol) in Hamiltonian object.
        else {
            const std::vector<double>& PE = ha->getPE();
            for (int i = 0; i < DOFe; i++)
                data[5 + i * nvalues] = PE[i];
        }
    }
    // Get enegry components of each state. In this case, both for classical MD
    // and nonadibatic dynamics, the energy should be computed.
    else {
        // Note that it only works for AMBER-like force filed now.
        // For each state, the components of potential energy are Nonbond, Bond,
        // Angle, Dihedral, Others, Coulomb, vdW. They belong to different force groups.
        // Others is from common group 0 for all states. group 1 is nonbond (total),
        //  2 is bond, 3 is angel, 4 is dihedral, and 5 is a fake LennardJonesForce.
        // (vdW), Coulomb is got by nonbond - vdW.
        std::vector<int> common_group(1, 1 << 0); // common group 0 for Others
        double others = ha->getPotentialEnergy(common_group, false);
        for (int i = 0; i < DOFe; i++) { // state index
            data[10 + i * nvalues] = others; // Others: data[10]
            data[5 + i * nvalues] = others; // Others is same for all states
            // The forces from group i will be included if (groups&(1<<i)) != 0.
            std::vector<int> groups((i + 1), 0);
            for (int j = 0; j < 4; j++) { // Nonbond, Bond, Angle, Dihedral
                groups[i] = 1 << (j + 1); // data[6-9]
                data[6 + j + i * nvalues] =  ha->getPotentialEnergy(groups, false);
                // data[5] Potential energy = sum of Nonbond, Bond, Angle, Dihedral
                data[5 + i * nvalues] += data[6 + j + i * nvalues];
            }
            groups[i] = 1 << 5; //  LennardJonesForce: vdW, data[12]
            data[12 + i * nvalues] =  ha->getPotentialEnergy(groups, false);
            // data[11] Coulomb = nonbond - vdW
            data[11 + i * nvalues] = data[6 + i * nvalues] - data[12 + i * nvalues];
        }
    }
    // Get total energy (data[3]) of each state: sum of kinetic energy data[4]
    // and potential energy data[5].
    for (int i = 0; i < DOFe; i++)
        data[3 + i * nvalues] = data[4 + i * nvalues] + data[5 + i * nvalues];
    count++;
}

void Simulation::writeStateData(const std::string& file, std::vector<double>& data, int step, double time) {
    // * Get required parameters.
    static int count = 0;
    static std::string file_prev = file; // the filename of previous call
    if (file_prev != file) { // if it is not a same file,
        count = 0;           // we reset the count to 0,
        file_prev = file;    // and set current file name as previous one.
    }
    static const bool energy_decompose = param.getBool("energy_decompose");
    static const int nvalues = 13; // 13 values for each state in data
    const int nstates = data.size()/nvalues; // number of states
    const std::string type = GetFileSuffix(file); // file format type
    // * Save current state data to file
    // The data in vector are (energies in kj/mol):
    // data[0]: Temperature (in K), data[1]: Volume (in nm^3),
    // data[2]: Density (in kg/m^3), data[3]: TotalEnergy, data[4]: KineticEnergy,
    // data[5-12]: total potential energy of state 0, and the energy components
    // of potential energy: Nonbond, Bond, Angle, Dihedral, Others, Coulomb, vdW.
    // Here, the Others means from Custom Force such as potsition restraint.
    // data[13-25]: same as above, but for second state (state 1).
    // data[26-38]: same as above, but for third state (state 2).
    FILE* dataFile = CheckFile(file, count);
    if (type == "dat") { // * Using fixed space separated format (.dat)
        if (count == 0) { // Print headers and units at the first time.
            fprintf(dataFile, "%10s %12s %5s %12s %12s %12s %16s %16s %16s", "Step", "Time", "State",
                "Temperature", "Volume", "Density", "TotalEnergy", "KineticEnergy", "PotentialEnergy");
            if (energy_decompose)
                fprintf(dataFile, " %16s %16s %16s %16s %16s %16s %16s", "Nonbonded_E", "Bond_E",
                    "Angle_E", "Dihedral_E", "Others_E", "Coulomb_E", "Lennard-Jones_E");
            fprintf(dataFile, "\n");
            fprintf(dataFile, "%10s %12s %5s %12s %12s %12s %16s %16s %16s", " ", "(ps)",
                " ", "(K)", "(nm^3)", "(kg/m^3)", "(kJ/mol)", "(kJ/mol)", "(kJ/mol)");
            if (energy_decompose)
                fprintf(dataFile, " %16s %16s %16s %16s %16s %16s %16s", "(kJ/mol)",
                    "(kJ/mol)", "(kJ/mol)", "(kJ/mol)", "(kJ/mol)", "(kJ/mol)", "(kJ/mol)");
            fprintf(dataFile, "\n");
        }
        for (int state = 0; state < nstates; ++state) { // Write data of each state.
            fprintf(dataFile, "%10d %12.4f %5d %12.4f %12.4f %12.4f %16.6f %16.6f %16.6f",
                step, time, state, data[0+state*nvalues], data[1+state*nvalues], data[2+state*nvalues],
                data[3+state*nvalues], data[4+state*nvalues], data[5+state*nvalues]);
            if (energy_decompose) // energy components: data[6-12]
                for (int i = 6; i < 13; ++i)
                    fprintf(dataFile, " %16.6f", data[i + state*nvalues]);
            fprintf(dataFile, "\n");
        }
    }
    else if (type == "csv") { // * comma (,) separated value (CSV) file
        if (count == 0) { // Print headers at the first time.
            fprintf(dataFile, "Step,Time,State,Temperature,Volume,Density,TotalEnergy,KineticEnergy,PotentialEnergy");
            if (energy_decompose)
                fprintf(dataFile, ",Nonbonded_E,Bond_E,Angle_E,Dihedral_E,Others_E,Coulomb_E,Lennard-Jones_E");
            fprintf(dataFile, "\n");
        }
        for (int state = 0; state < nstates; ++state) { // Write data of each state
            fprintf(dataFile, "%d,%g,%d,%g,%g,%g,%.12g,%.12g,%.12g", step, time, state,
                data[0+state*nvalues], data[1+state*nvalues], data[2+state*nvalues],
                data[3+state*nvalues], data[4+state*nvalues], data[5+state*nvalues]);
            if (energy_decompose) // energy components: data[6-12]
                for (int i = 6; i < 13; ++i)
                    fprintf(dataFile, ",%.12g", data[i + state*nvalues]);
            fprintf(dataFile, "\n");
        }
    }
    else
        throw std::runtime_error("ERROR: Unknown file type for output: " + type);
    fclose(dataFile);
    count++;
}

void Simulation::writeTTMFile(const Complex_Matrix& RDM, const std::string& file) {
    // Get the needed parameters for output.e
    const double DT         = param.getDouble("DT"); // nuclear time step in ps
    const int    DOFe       = param.getInt("DOFe"); // number of states
    const int    RDM_steps  = param.getInt("RDM_steps"); // output frequency
    // Save reduced density matrix to file.
    FILE* outfile = CheckFile(file);
    // CSV file, fields are separated by one comma (,).
    fprintf(outfile, "%s", "time");
    for (int j = 0; j < DOFe; j++)
        for (int k = 0; k < DOFe; k++) {
                fprintf(outfile, ",%s%d%d%s", "rho", j, k, "re");
                fprintf(outfile, ",%s%d%d%s", "rho", j, k, "im");
            }
    fprintf(outfile, ",%s", "totalpop");
    if (DOFe == 2) // sigmaz = pop_0-pop_1 (DOFe = 2 only)
        fprintf(outfile, ",%s", "sigmaz");
    fprintf(outfile, "\n");
    // Here the size of RDM is nsteps/RDM_steps + 1, since both the values
    // of step 0 and last step are stored. And, denisty matrix of each step is
    // stroed as DOFe^2-dimentional vector. The tranformation between index s in
    // vector and index (a,b) in matrix is: s = a*DOFe+b, 0 <= a,b < DOFe, and
    // a = s / DOFe, b = s % DOFe.
    for (int i = 0; i < RDM.size(); i++ ) {
        double sum = 0.0;
        for (int n = 0; n < DOFe; n++)
            sum += RDM[i][n*DOFe+n].real(); // get total population firstly
        fprintf(outfile, "%.12g", i*DT*RDM_steps); // time
        for (int j = 0; j < DOFe; j++)
            for (int k = 0; k < DOFe; k++)
                fprintf(outfile, ",%.12g,%.12g", RDM[i][j*DOFe+k].real(), RDM[i][j*DOFe+k].imag());
        fprintf(outfile, ",%.12g", sum);
        // Output sigmaz = pop_0-pop_1 (DOFe = 2 only)
        if (DOFe == 2)
            fprintf(outfile, ",%.12g", (RDM[i][0].real() - RDM[i][3].real()));
        fprintf(outfile, "\n");
    }
    fclose(outfile);
}
