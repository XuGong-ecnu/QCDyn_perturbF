/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 23, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "Parameters.h"

void Parameters::loadParamters(int argc, char *argv[], bool isNew) {
    if (argc > 1) {
        // Load parameters form input file firstly.
        int n = 0;
        for (int i = 1; i < argc; i++)
            if (argv[i][0] != '-') {
                std::cout << "Loaded the following parameters from input control file: " << argv[i] << std::endl;
                loadParamFromFile(argv[i], isNew);
                n = i; // records the postion of control file.
                break;
            }
        // Overwrite parameters from command line like "--jobid=123".
        if (!(n == 1 && argc == 2)) {
            std::cout << "Overwrote the following parameters from command line: " << std::endl;
            for (int i = 1; i < argc; i++)
                if (argv[i][0] == '-')
                    loadParamFromString(argv[i]);
        }
    }
    // Print the usage of program and all parameters with default values
    else {
        std::cout << "Usage: " << argv[0] << " [control.inp] [--key1=value1 --Key2=value2...]\n";
        std::cout << "Example 1: " << argv[0] << " control.inp\n";
        std::cout << "Example 2: " << argv[0] << " --structure==bent_GR.gro --topology=bent_GR.top\n";
        std::cout << "Example 3: " << argv[0] << " control.inp --temperature=500\n";
        std::cout << "Example 4: " << argv[0] << " --temperature=500 control.inp\n";
        if (!parameters.empty()) {
            std::cout << "The supported keys and default values of parameters are as following: " << std::endl;
            printAllParam();
        }
        exit(0);
    }
    std::cout << std::endl;
}

void Parameters::loadParamFromFile(const std::string& file, bool isNew) {
    std::ifstream inpFile(file.c_str());
    if (!inpFile)
        throw std::runtime_error("ERROR: Unable to open input file: " + file);
    for (std::string line; getline(inpFile, line); ) {
        // A line starting with (#), or empty (including space, \t), will be ignored.
        if (line[0] == '#' || IsBlankString(line)) continue;
        // Read content from input file, and separate each column to list.
        // Each column should be divided with at least a space.
        std::vector<std::string> list;
        SplitLine(list, line);
        if (isNew) {
            if (list.size() < 3)
                throw std::runtime_error("ERROR: Invalid line for key: " + list[0]);
            addParam(list[0], list[1], list[3]);
        }
        else {
            if (list.size() < 2)
                throw std::runtime_error("ERROR: Invalid line for key: " + list[0]);
            setValue(list[0], list[1]);
        }
        printParam(list[0]);
    }
    inpFile.close();
}

void Parameters::loadParamFromString(const std::string& command) {
    if (command.substr(0, 2) != "--")
        throw std::runtime_error("ERROR: Command line overwrite parameter should be starting with --, such as --id=123.");
    // Split "id=123" style into key and value to list.
    std::vector<std::string> list;
    SplitString(list, command.substr(2), '=');
    if (list.size() != 2)
        throw std::runtime_error("ERROR: Command line overwrite parameter should be like --id=123, indicating both key and value!");
    setValue(list[0], list[1]);
    printParam(list[0]);
}

void Parameters::addParam(const std::string& key, const std::string& value, const std::string& type) {
    if (type != "str" && type != "int" && type != "double" && type != "bool")
        throw std::runtime_error("ERROR: Unsupported type: " + type);
    if (parameters.find(key) != parameters.end())
        throw std::runtime_error("ERROR: The parameter: " + key + " already exists.");
    parameters[key] = ParamEntry(key, value, type);
    keys.emplace_back(key);
}

void Parameters::setValue(const std::string& key, const std::string& value) {
    if (parameters.find(key) == parameters.end())
        throw std::runtime_error("ERROR: Called setValue() for a nonexistent key: " + key);
    parameters[key] = ParamEntry(key, value, parameters.at(key).type);
}

std::string Parameters::getType(const std::string& key) const {
    if (parameters.find(key) == parameters.end())
        throw std::runtime_error("ERROR: Called getType() for a nonexistent key: " + key);
    else
        return parameters.at(key).type;
}

std::string Parameters::getStr(const std::string& key) const {
    if (parameters.find(key) == parameters.end())
        throw std::runtime_error("ERROR: Called getStr() for a nonexistent key: " + key);
    else
        return parameters.at(key).str_value;
}

int Parameters::getInt(const std::string& key) const {
    if (parameters.find(key) == parameters.end())
        throw std::runtime_error("ERROR: Called getInt() for a nonexistent key: " + key);
    else
        return parameters.at(key).int_value;
}

double Parameters::getDouble(const std::string& key) const {
    if (parameters.find(key) == parameters.end())
        throw std::runtime_error("ERROR: Called getDouble() for a nonexistent key: " + key);
    else
        return parameters.at(key).double_value;
}

bool Parameters::getBool(const std::string& key) const {
    if (parameters.find(key) == parameters.end())
        throw std::runtime_error("ERROR: Called getBool() for a nonexistent key: " + key);
    else
        return parameters.at(key).bool_value;
}

void Parameters::printParam(const std::string& key) const {
    if (parameters.find(key) == parameters.end())
        throw std::runtime_error("ERROR: Called printParam() for a nonexistent key: " + key);
    std::cout << std::setw(24) << std::left << key << " " << parameters.at(key).str_value << std::endl;
}

void Parameters::printAllParam() const {
    for (const auto& key : keys)
        printParam(key);
}

Parameters::ParamEntry::ParamEntry(const std::string& key, const std::string& value, const std::string& type) : name(key), type(type), str_value(value) {
    // If value is not a number, the int/double_value will be zero.
    std::stringstream ss;
    ss << value;
    ss >> double_value;
    int_value = double_value;
    if (type == "str")
        bool_value = (str_value == "true" ? true : false);
    else if (type == "int" || type == "double")
        bool_value = (int_value == 0 ? false : true);
    else if (type == "bool") {
        if (str_value == "true")  
            bool_value = true;
        else if (str_value == "false")
            bool_value = false;
        else 
            throw std::runtime_error("ERROR: parameter <"+ key + "> must be true or false! ");
        }
}

void InitializeParameters(int argc, char *argv[], Parameters& param) {
    // Add parameters with keys and default values.
    AddDefaultParameters(param);
    // Load parameters from input control file and command line.
    param.loadParamters(argc, argv);
    // Print all the parameters will be used.
    // std::cout << "The following parameters will be used in this simulation: " << std::endl;
    // param.printAllParam();
}

void AddDefaultParameters(Parameters& param) {
    // Create a tuple to store default parameter: key, value, type
    std::vector<std::tuple<std::string, std::string, std::string>> defaultParam;
    // * Global Settings (System and Dynamics type)
    // Settings for simulation job type: dynamics, preprocess, postprocess
    // dyn_type: OpenMM, LSC, SQC, MF, TBSH, FSSH, ...
    // prep_type: model_gen
    // post_type: TTM
    defaultParam.emplace_back("job_type",                "dynamics",   "str"   );
    defaultParam.emplace_back("dyn_type",                "",           "str"   );
    defaultParam.emplace_back("prep_type",               "model_gen",  "str"   );
    defaultParam.emplace_back("post_type",               "TTM",        "str"   );
    defaultParam.emplace_back("md_type",                  "",          "str"   );
    defaultParam.emplace_back("N_beads",                  "",          "int"   );
    // This is for PIMD simulation parameters
    defaultParam.emplace_back("PIMD_type",                "",          "str"   ); 
    defaultParam.emplace_back("nbeads",                   "",          "int"   );
    defaultParam.emplace_back("PIMD_thermostat",          "",          "str"   );// none, Langevin
    defaultParam.emplace_back("ilngv_flag",               "",          "str"   );// 0, none, 1: not on centroid; TO DO-> 2: on centroid
    defaultParam.emplace_back("tau_lngv_ps",              "",       "double"   );
    // This is to compute the TCF functions
    defaultParam.emplace_back("VibSpec_save",            "",           "str"   );
    defaultParam.emplace_back("LEN_TRAJ",                "",           "int"   );
    defaultParam.emplace_back("LEN_TCF",                 "",           "int"   );
    defaultParam.emplace_back("LEN_SKIP",                "",           "int"   );
    defaultParam.emplace_back("LEN_CORR",                "",           "int"   );
    defaultParam.emplace_back("OKE_type",                "",           "str"   );
    defaultParam.emplace_back("polar_obs",          "false",           "bool"  ); 
    defaultParam.emplace_back("forcefield_type",         "",           "str"   ); 
    defaultParam.emplace_back("DampingParameter1",       "",           "double");
    defaultParam.emplace_back("DampingParameter2",       "",           "double");
    defaultParam.emplace_back("DampingParameter3",       "",           "double");
    defaultParam.emplace_back("iter_max",                "",           "int"   );
    defaultParam.emplace_back("TOLERANCE",               "",           "double");
    defaultParam.emplace_back("damping_type",            "",           "str"   );
    defaultParam.emplace_back("dampingParam1",            "",       "double"   );
    defaultParam.emplace_back("dampingParam2",            "",       "double"   );
    defaultParam.emplace_back("dampingParam3",            "",       "double"   );
    defaultParam.emplace_back("xxxxxx",            "",       "double"   ); //gongxu
    // This is for C60 model vdw model
    defaultParam.emplace_back("LJ_type",                 "",           "str"   );
    // Settings for MixPES dyn_type
    defaultParam.emplace_back("PES_weights",             "",           "str"   );
    // Settings for perturb MD production
    // This is the perturb elec part paramenters
    defaultParam.emplace_back("perturb_parameters",      "",           "str"   );
    // This is to judge the perturb
    defaultParam.emplace_back("perturb",                 "false",      "bool"  );
    defaultParam.emplace_back("steps_for_config",          "",          "int"  );
    defaultParam.emplace_back("LEN_steps",                 "",          "int"  );
    // This is for  DES system COM calculation.
    defaultParam.emplace_back("DES_type",                 "",           "str" );
    defaultParam.emplace_back("MolI",                     "",          "int"  );
    defaultParam.emplace_back("MolN",                     "",          "int"  );
    defaultParam.emplace_back("MolS",                   "",          "int"  );
    defaultParam.emplace_back("DirecA",                 "",          "int"  );
    defaultParam.emplace_back("DirecB",                 "",          "int"  );
    defaultParam.emplace_back("OUTDES",             "",      "str"   );
    defaultParam.emplace_back("DES_save",             "",      "str"   );
    // Settings for simulation system type: model, onthefly, allatom
    // model_type: SB (default), GOA, FMO, MSH, ...
    // onthefly_type: not implemented now
    // allatom: OpenMM and CustomForceField
    defaultParam.emplace_back("system_type",             "model",      "str"   );
    defaultParam.emplace_back("model_type",              "SB",         "str"   );
    defaultParam.emplace_back("onthefly_type",           "",           "str"   );
    defaultParam.emplace_back("turn_off_force_contribution",  "",      "str"   );
    defaultParam.emplace_back("allatom_type",            "",           "str"   );
    //TO DO: defaultParam.emplace_back("GreenKuboTransport", "Diffusion,ShearViscosity,BulkViscosity,Friction,ThermalConductivity", "str" );
    // * OpenMM Settings
    // Settings for Platform and Platform-specific properties
    defaultParam.emplace_back("platform",                "",           "str"   );
    defaultParam.emplace_back("precision",               "",           "str"   );
    defaultParam.emplace_back("deterministic_forces",    "",           "str"   );
    defaultParam.emplace_back("CPU_PME",                 "",           "str"   );
    defaultParam.emplace_back("disable_PMEStream",       "",           "str"   );
    defaultParam.emplace_back("OpenCL_index",            "",           "str"   );
    defaultParam.emplace_back("CUDA_compiler",           "",           "str"   );
    defaultParam.emplace_back("temp_directory",          "",           "str"   );
    defaultParam.emplace_back("device_index",            "",           "str"   );
    defaultParam.emplace_back("blocking_sync",           "",           "str"   );
    defaultParam.emplace_back("CPU_threads",             "",           "str"   );
    defaultParam.emplace_back("plugins_directory",       "",           "str"   );
    // Settings for input files
    defaultParam.emplace_back("structure",               "",           "str"   );
    defaultParam.emplace_back("topology",                "",           "str"   );
    defaultParam.emplace_back("restart",                 "none",       "str"   );
    // Settings for output control
    defaultParam.emplace_back("default_name",            "none",       "str"   );
    defaultParam.emplace_back("data_file",               "data.csv",   "str"   );
    defaultParam.emplace_back("energy_steps",            "100",        "int"   );
    defaultParam.emplace_back("energy_decompose",        "false",      "bool"  );
    defaultParam.emplace_back("conf_file",               "confout.gro","str"   );
    defaultParam.emplace_back("traj_file",               "traj.xtc",   "str"   );
    defaultParam.emplace_back("traj_group",              "system",     "str"   );
    defaultParam.emplace_back("xtc_precision",           "4",          "int"   );
    defaultParam.emplace_back("enforce_box",             "false",      "bool"  );
    defaultParam.emplace_back("position_steps",          "0",          "int"   );
    defaultParam.emplace_back("velocity_save",           "false",      "bool"  );
    defaultParam.emplace_back("force_save",              "false",      "bool"  );
    defaultParam.emplace_back("chk_file",                "state.chk",  "str"   );
    defaultParam.emplace_back("chk_steps",               "10000",      "int"   );
    defaultParam.emplace_back("progress_interval",       "10",         "int"   );
    // Settings for integrator used for propagation
    defaultParam.emplace_back("integrator",          "velocityVerlet", "str"   );
    defaultParam.emplace_back("init_time",               "0",          "double");
    defaultParam.emplace_back("nsteps",                  "1000",       "int"   );
    defaultParam.emplace_back("DT",                      "0.001",      "double");
    defaultParam.emplace_back("propagate_state",         "0",          "int"   );
    // Settings for initial velocities
    defaultParam.emplace_back("initial_velocity",        "auto",       "str"   );
    defaultParam.emplace_back("Boltzmann_temperature",   "300",        "double");
    defaultParam.emplace_back("Boltzmann_seed",          "0",          "int"   );
    defaultParam.emplace_back("velocity_correction",     "false",      "bool"  );
    // Settings for energy minimization
    defaultParam.emplace_back("energy_minimization",     "false",      "bool"  );
    defaultParam.emplace_back("force_tolerance",         "10",         "double");
    defaultParam.emplace_back("max_iterations",          "10000",      "int"   );
    defaultParam.emplace_back("minimized_structure",     "em.gro",     "str"   );
    // Settings for PBC and nonbonded method
    defaultParam.emplace_back("PBC",                     "xyz",        "str"   );
    defaultParam.emplace_back("nonbonded_method",        "PME",        "str"   );
    defaultParam.emplace_back("cutoff",                  "1.0",        "double");
    defaultParam.emplace_back("dispersion_correction",   "true",       "bool"  );
    defaultParam.emplace_back("Ewald_tolerance",         "0.0005",     "double");
    defaultParam.emplace_back("PME_parameters",          "0",          "str"   );
    defaultParam.emplace_back("LJPME_parameters",        "0",          "str"   );
    defaultParam.emplace_back("switching_function",      "false",      "bool"  );
    defaultParam.emplace_back("switching_distance",      "0.9",        "double");
    defaultParam.emplace_back("reaction_field",          "false",      "bool"  );
    defaultParam.emplace_back("RF_dielectric",           "78.3",       "double");
    // Settings for constrains
    defaultParam.emplace_back("constraints",             "none",       "str"   );
    defaultParam.emplace_back("constraint_tolerance",    "0.00001",    "double");
    defaultParam.emplace_back("constraint_start",        "false",      "bool"  );
    defaultParam.emplace_back("rigid_water",             "true",       "bool"  );
    defaultParam.emplace_back("repart_HMass",            "false",      "bool"  );
    // Settings for center of mass (COM) motion removal
    defaultParam.emplace_back("remove_COMMotion",        "true",       "bool"  );
    defaultParam.emplace_back("remove_steps",            "1",          "int"   );
    // Settings for thermostat
    defaultParam.emplace_back("thermostat",              "none",       "str"   );
    defaultParam.emplace_back("temperature",             "300.0",      "double");
    defaultParam.emplace_back("friction_coefficient",    "1.0",        "double");
    defaultParam.emplace_back("collision_frequency",     "1.0",        "double");
    defaultParam.emplace_back("NHC_parameters",          "3,3,7",      "str"   );
    defaultParam.emplace_back("thermostat_seed",         "0",          "int"   );
    // Settings for simulated annealing
    defaultParam.emplace_back("annealing",               "none",       "str"   );
    defaultParam.emplace_back("annealing_npoints",       "2",          "int"   );
    defaultParam.emplace_back("annealing_time",          "0,50",       "str"   );
    defaultParam.emplace_back("annealing_temperature",   "0,300",      "str"   );
    defaultParam.emplace_back("temperature_increment",   "1",          "double");
    // Settings for barostat (pressure coupling)
    defaultParam.emplace_back("barostat",                "none",       "str"   );
    defaultParam.emplace_back("pressure",                "1.0",        "double");
    defaultParam.emplace_back("pressure_steps",          "25",         "int"   );
    defaultParam.emplace_back("barostat_seed",           "0",          "int"   );
    // Settings for position restraint with harmonic force
    defaultParam.emplace_back("position_restraint",      "false",      "str"   );
    defaultParam.emplace_back("restraint_atoms",         "",           "str"   );
    defaultParam.emplace_back("restraint_force",         "1000",       "double");
    defaultParam.emplace_back("reference_structure",     "",           "str"   );
    // Settings for implicit solvation model
    defaultParam.emplace_back("implicit_solvation",      "false",      "bool"  );
    defaultParam.emplace_back("solvent_dielectric",      "78.5",       "double");
    defaultParam.emplace_back("solute_dielectric",       "1.0",        "double");
    // * Settings for Model System
    // Settings of model system
    defaultParam.emplace_back("DOFn",                    "100",        "int"   );
    defaultParam.emplace_back("DOFe",                    "2",          "int"   );
    defaultParam.emplace_back("energy_unit",             "",           "str"   );
    defaultParam.emplace_back("gamma_DA",                "1.0",        "double");
    defaultParam.emplace_back("epsilon",                 "1.0",        "double");
    defaultParam.emplace_back("Condon_approximation",    "true",       "bool"  );
    defaultParam.emplace_back("H_load",                  "",           "str"   );
    defaultParam.emplace_back("H_save",                  "",           "str"   );
    // Settings of model paramters
    defaultParam.emplace_back("spec_density",            "Ohmic",      "str"   );
    defaultParam.emplace_back("Ohmic_discrete",          "CMM",        "str"   );
    defaultParam.emplace_back("eta",                     "0",          "double");
    defaultParam.emplace_back("lambda",                  "0",          "double");
    defaultParam.emplace_back("omega_c",                 "0",          "double");
    defaultParam.emplace_back("omega_max",               "0",          "double");
    // The following 3 parameters are used by GOA model only
    defaultParam.emplace_back("GOA_Omega",               "0",          "double");
    defaultParam.emplace_back("GOA_y0",                  "0",          "double");
    defaultParam.emplace_back("GOA_shift",               "0",          "double");
    defaultParam.emplace_back("GOA_gamma",               "0",          "double");
    defaultParam.emplace_back("model_load",              "",           "str"   );
    defaultParam.emplace_back("model_save",              "",           "str"   );
    // Settings for contructing multi-state harmonic (MSH) model
    defaultParam.emplace_back("MSH_type",                "I",          "str"   );
    defaultParam.emplace_back("energy_load",             "",           "str"   );
    defaultParam.emplace_back("energy_correction",       "",           "str"   );
    defaultParam.emplace_back("TCF_prefix",              "MSH",        "str"   );
    defaultParam.emplace_back("TCF_XY",                  "1,2",        "str"   );
    defaultParam.emplace_back("TCF_ntraj",               "1",          "int"   );
    defaultParam.emplace_back("TCF_corsteps",            "800",        "int"   );
    defaultParam.emplace_back("TCF_DT",                  "0.005",      "double");
    defaultParam.emplace_back("TCF_maxiter",             "1000",       "int"   );
    defaultParam.emplace_back("TCF_tol",                 "1e-8",       "double");
    // Settings for constructing non-Condon Morse model by Miller
    defaultParam.emplace_back("Morse_preset",            "",           "str"   );
    // * Settings for Nonadiabtic Dynamics
    defaultParam.emplace_back("ntraj",                   "1",          "int"   );
    defaultParam.emplace_back("gamma_MM",                "",           "str"   );
    defaultParam.emplace_back("representation",          "diabatic",   "str"   );
    defaultParam.emplace_back("gamma_in_EOM",            "false",      "bool"  );
    defaultParam.emplace_back("adjust_gamma_per_traj",   "false",      "bool"  );
    // Settings for specific dynamics method
    defaultParam.emplace_back("LSC_zeta",                "1.0",        "double");
    defaultParam.emplace_back("SQC_window",              "triangle",   "str"   );
    defaultParam.emplace_back("SPM_type",                "W",          "str"   );
    defaultParam.emplace_back("TBSH_filter",             "0",          "double");
    defaultParam.emplace_back("FSSH_seed",               "0",          "int"   );
    defaultParam.emplace_back("TBSH_seed",               "0",          "int"   );
    // Settings for the initial nuclear sampling
    defaultParam.emplace_back("nucl_sample",             "none",       "str"   );
    defaultParam.emplace_back("beta",                    "1",          "double");
    defaultParam.emplace_back("sample_state",            "0",          "int"   );
    defaultParam.emplace_back("shift_load",              "",           "str"   );
    defaultParam.emplace_back("nucl_seed",               "0",          "int"   );
    defaultParam.emplace_back("nucl_load",               "",           "str"   );
    defaultParam.emplace_back("nucl_start",              "0",          "int"   );
    defaultParam.emplace_back("nucl_end",                "-1",         "str"   );
    defaultParam.emplace_back("nucl_skip",               "1",          "str"   );
    defaultParam.emplace_back("nucl_save",               "",           "str"   );
    // Settings for the electronic sampling and propagation
    defaultParam.emplace_back("elec_sample",             "",           "str"   );
    defaultParam.emplace_back("elec_seed",               "0",          "int"   );
    defaultParam.emplace_back("elec_load",               "",           "str"   );
    defaultParam.emplace_back("elec_save",               "",           "str"   );
    defaultParam.emplace_back("init_state",              "0,0",        "str"   );
    defaultParam.emplace_back("EPN",                     "1",          "int"   );
    // Settings for the output
    defaultParam.emplace_back("RDM_steps",               "1",          "int"   );
    defaultParam.emplace_back("RDM_save",                "",           "str"   );
    // Settings for transfer tensor method (TTM) only
    defaultParam.emplace_back("RDM_load",                "",           "str"   );
    defaultParam.emplace_back("TTM_kmax",                "0",          "int"   );
    // Settings for Observable class
    defaultParam.emplace_back("observables",             "",        "str"   );
    defaultParam.emplace_back("max_friction_fail_traj",  "0.3",        "double"); 
    defaultParam.emplace_back("min_judge_traj_avail",    "1000",       "int"   ); 

    // Add all default parameters to Paramters.
    for (const auto& i : defaultParam)
        param.addParam(std::get<0>(i), std::get<1>(i), std::get<2>(i));
}
