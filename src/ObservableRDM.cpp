/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu and Zengkui Liu @Sun Group @NYU-SH                       *
 * Last updated: Jan. 6, 2022                                                 *
 * -------------------------------------------------------------------------- */

#include "ObservableRDM.h"

void ObservableRDM::init() {
    dyn_type        = param.getStr("dyn_type");
    post_type       = param.getStr("post_type");
    system_type     = param.getStr("system_type");
    model_type      = param.getStr("model_type");
    onthefly_type   = param.getStr("onthefly_type");
    allatom_type    = param.getStr("allatom_type");
    representation  = param.getStr("representation");
    max_friction_fail_traj = param.getDouble("max_friction_fail_traj");
    min_judge_traj_avail   = param.getInt("min_judge_traj_avail");
    // initalization of the counting parameters 
    traj_rejected = 0;
    traj_accepted = 0;
    ObservableBase::init();
    RDM_steps = param.getInt("RDM_steps"); // default is 1
    // Get the initial reduced density matrix element indices |j><k| specified
    // by user, the input format is "j,k", e.g. 0,0, then the initial density
    // matrix will be initialized as: rho_00 is 1 and others are 0.
    // If init_state=all, then all possible initial states are considered, and
    // all reduced density matrices starting from all initial states wiil be
    // computed and save (works for some methods).
    const std::string init_state_str = param.getStr("init_state");
    if (init_state_str == "all")  { // all possible initial states, DOFe*DOFe
        init_state.resize(DOFe*DOFe, 0);
        for (int i = 0; i < DOFe*DOFe; ++i) // index of initial state i=j*DOFe+k
            init_state[i] = i;
    }
    else { // specfied one initial state with two indices
        std::vector<int> indices;
        SplitString(indices, param.getStr("init_state"));
        if (indices.size() != 2)
            throw std::runtime_error("ERROR: Illegal init_state=" + param.getStr("init_state") +
                ". You must provide 2 indices to represent initial state.");
        for (int i = 0; i < 2; ++i)
            if (indices[i] >= DOFe || indices[i] < 0)
                throw std::runtime_error("ERROR: Illegal init_state index: " + std::to_string(indices[i]) +
                    ", which should be in range [0, DOFe).");
        init_state.resize(1, 0);
        init_state[0] = indices[0]*DOFe + indices[1];
    }
    if (RDM_steps < 1 || (nsteps % RDM_steps != 0))
        throw std::runtime_error("ERROR: Illegal value for RDM_steps (requires >= 1). "
           "And the value of nsteps must be divided exactly by RDM_steps.");   
    // Initialized the vector of density matrix (stored as DOFe^2-dimentional vector).
    // Plus 1 to the length, since the step 0 and last step will be stored.
    // The size of RDMs depends on the init_denisty (one or DOFe*DOFe).
    LEN_TRAJ = nsteps/RDM_steps + 1; // + 1 since the step 0
    // Get number of initial states
    num_init_states = 1;
    if (param.getStr("init_state") == "all") // inlcude all initial states
        num_init_states = DOFe * DOFe;
    RDM_current.resize(num_init_states, Complex_Matrix(LEN_TRAJ, std::vector<Complex>(DOFe*DOFe, 0)));
    RDM_average.resize(num_init_states, Complex_Matrix(LEN_TRAJ, std::vector<Complex>(DOFe*DOFe, 0)));  
    // For LSC method, 5 window types will be computed in one simulation. They are LSC1/2, and RI-LSC1~3. 
    if (dyn_type == "LSC") {
        RDM_current_LSC1 = RDM_current_LSC2 = RDM_current_RILSC1 = RDM_current_RILSC2 = RDM_current_RILSC3 = RDM_current;
        RDM_average_LSC1 = RDM_average_LSC2 = RDM_average_RILSC1 = RDM_average_RILSC2 = RDM_average_RILSC3 = RDM_average;      
    }
    else if (dyn_type == "SQC") {
        // Normalization factors for population and coherences
        norm_pop.resize(nsteps, 0);
        norm_coh.resize(nsteps, 0);
    }
    if (representation == "adiabatic") {
        if (dyn_type == "LSC") {
            AdiabaticRDM_current_LSC1 = AdiabaticRDM_current_LSC2 = AdiabaticRDM_current_RILSC1 = AdiabaticRDM_current_RILSC2 = AdiabaticRDM_current_RILSC3 = RDM_current;
            AdiabaticRDM_average_LSC1 = AdiabaticRDM_average_LSC2 = AdiabaticRDM_average_RILSC1 = AdiabaticRDM_average_RILSC2 = AdiabaticRDM_average_RILSC3 = RDM_average;
        }
        else if (dyn_type == "TBSH") {
            AdiabaticRDM_current = AdiabaticRDM_average = RDM_current;
        }

    }
}

int ObservableRDM::updateOneSnapshot() {
    int step = Dy -> getStep();
    // Fur LSC method there are 5 kinds of observables fur dynamics, which should be identified by 
    if (dyn_type == "LSC") {
        for (int i = 0; i < num_init_states; ++i) {
            for (int j = 0; j < DOFe; ++j)
                for (int k = 0; k < DOFe; ++k) { 
                    // get values in the class
                    std::vector<std::vector<Complex> > RDM_LSC1 = Dy ->getRDMLSC1();
                    flag_traj_avail = std::isfinite(RDM_LSC1[i][j*DOFe+k].real());
                    if (flag_traj_avail) { // return true for valid value
                        std::vector<std::vector<Complex> > RDM_LSC2   = Dy->getRDMLSC2();
                        std::vector<std::vector<Complex> > RDM_RILSC1 = Dy->getRDMRILSC1();
                        std::vector<std::vector<Complex> > RDM_RILSC2 = Dy->getRDMRILSC2();
                        std::vector<std::vector<Complex> > RDM_RILSC3 = Dy->getRDMRILSC3();

                        RDM_current_LSC1[i][step/RDM_steps][j*DOFe+k]   = RDM_LSC1[i][j*DOFe+k];
                        RDM_current_LSC2[i][step/RDM_steps][j*DOFe+k]   = RDM_LSC2[i][j*DOFe+k];
                        RDM_current_RILSC1[i][step/RDM_steps][j*DOFe+k] = RDM_RILSC1[i][j*DOFe+k];
                        RDM_current_RILSC2[i][step/RDM_steps][j*DOFe+k] = RDM_RILSC2[i][j*DOFe+k];
                        RDM_current_RILSC3[i][step/RDM_steps][j*DOFe+k] = RDM_RILSC3[i][j*DOFe+k];
                    }
                    else // NAN or INF
                        return -1;
                }
        }
    }
    else { // Other methods
        // std::cout << "updated RDM at step " << step << std::endl;
        // RDM_current[i][Dy->step/RDM_steps][j*DOFe+k] = Dy->RDM[i][j*DOFe+k];
        std::vector<std::vector<Complex> > RDM = Dy->getRDM();
        if (dyn_type == "SQC") {
           int i = 0;
            for (int j = 0; j < DOFe; ++j){
                for (int k = 0; k < DOFe; ++k) {                        
                    flag_traj_avail = std::isfinite(RDM[i][j*DOFe+k].real());
                    //std::cout << "updated RDM at step " << step << ", ijk " << i << j << k << ", flag " << flag_traj_avail <<", RDM elements" << RDM[i][j*DOFe+k] << std::endl;
                    if (flag_traj_avail)  // return true for valid value
                         RDM_current[i][step/RDM_steps][j*DOFe+k] = RDM[i][j*DOFe+k]; // if there is no not-a-number or infinite
                    else // NAN or INF
                        return -1;
                }
            }
            
        }
        else {
            for (int i = 0; i < num_init_states; ++i) {
                for (int j = 0; j < DOFe; ++j){
                    for (int k = 0; k < DOFe; ++k) {                   
                        flag_traj_avail = std::isfinite(RDM[i][j*DOFe+k].real());
                        // std::cout << "updated RDM at step " << step << ", ijk " << i << j << k << ", flag " << flag_traj_avail <<", RDM elements" << RDM[i][j*DOFe+k] << std::endl;
                        if (flag_traj_avail)   // return true for valid value
                             RDM_current[i][step/RDM_steps][j*DOFe+k] = RDM[i][j*DOFe+k]; // if there is no not-a-number or infinite
                        else // NAN or INF
                            return -1;
                    }
                }
            }
        }
    }
    if (representation == "adiabatic" && dyn_type == "TBSH") {
        std::vector<std::vector<Complex> > AdiabaticRDM = Dy->getAdiabaticRDM();
        for (int i = 0; i < num_init_states; ++i) {
            for (int j = 0; j < DOFe; ++j) {
                for (int k = 0; k < DOFe; ++k) {
                    flag_traj_avail = std::isfinite(AdiabaticRDM[i][j*DOFe+k].real());
                    if (flag_traj_avail) // return true for valid value
                         AdiabaticRDM_current[i][step/RDM_steps][j*DOFe+k] = AdiabaticRDM[i][j*DOFe+k];
                    else // NAN or INF
                        return -1;
                }
            }
        }
    }
    return 0; 
}

int ObservableRDM::updateOneTraj() {
    if (flag_traj_avail) {
        traj_accepted += 1;
	// std::cout << "traj accepted " << traj_accepted << std::endl; 
        if (dyn_type == "LSC") {
            for (int i = 0; i < num_init_states; ++i) {
                for (int step_count = 0; step_count < (nsteps/RDM_steps + 1); ++step_count) {
                    for (int j = 0; j < DOFe; ++j) {
                        for (int k = 0; k < DOFe; ++k) { 
                            RDM_average_LSC1[i][step_count][j*DOFe+k]   += RDM_current_LSC1[i][step_count][j*DOFe+k]; //(double) (ntraj);
                            RDM_average_LSC2[i][step_count][j*DOFe+k]   += RDM_current_LSC2[i][step_count][j*DOFe+k]; //(double) (ntraj);
                            RDM_average_RILSC1[i][step_count][j*DOFe+k] += RDM_current_RILSC1[i][step_count][j*DOFe+k]; //(double) (ntraj);
                            RDM_average_RILSC2[i][step_count][j*DOFe+k] += RDM_current_RILSC2[i][step_count][j*DOFe+k]; //(double) (ntraj);
                            RDM_average_RILSC3[i][step_count][j*DOFe+k] += RDM_current_RILSC3[i][step_count][j*DOFe+k]; //(double) (ntraj);
                        }
                    }
                }   
            }
        }
        else if (dyn_type == "TBSH") { // TBSH 
            for (int i = 0; i < num_init_states; ++i) 
                for (int step_count = 0; step_count < (nsteps/RDM_steps + 1); ++step_count) 
                    for (int j = 0; j < DOFe; ++j) 
                        for (int k = 0; k < DOFe; ++k)  
                            RDM_average[i][step_count][j*DOFe+k] += RDM_current[i][step_count][j*DOFe+k];
        }
        else if (dyn_type == "SQC") { // Other methods
            for (int i = 0; i < num_init_states; ++i) 
                for (int step_count = 0; step_count < (nsteps/RDM_steps + 1); ++step_count) 
                    for (int j = 0; j < DOFe; ++j) 
                        for (int k = 0; k < DOFe; ++k)  
                            RDM_average[i][step_count][j*DOFe+k] += RDM_current[i][step_count][j*DOFe+k]; //(double) (ntraj);
            
            if (traj_accepted == ntraj) {
                // 1. For SQC method, all values (inlcuding coherence) of RDMs are renormalized
                // by: value / raw_total_population, this will make sure total population is 1.
                // Ref: Eq. 10 in Cotton,and Miller, J. Chem. Phys. 2016, 145, 144108.
                // This is the popular/original way used in most papers.
                // ! However, in this way, re-normalization from coherence is not well defined.
                // 2. Another way is to average the RDM by the number of trajectoies lie
                // within population/coherence window. (in this case, the normalized factor
                // is different for population and coherence, and can be used starting from
                // coherence. For population-population, this two ways are same, but for
                // population-coherence is minor different according to the test.)
                // Ref: J. Chem. Phys. 148, 181102 (2018) and the code in github:
                // github.com/jprov410/mqds, mqds-master/mqds/src/general_src/windows.f90
                // ! However, also, re-normalization from coherence is not well defined.
                // ! So, when starting from coherence, the re-normalization is unknown.
                const int j = init_state[0] / DOFe; // initial state = j*DOFe+k
                const int k = init_state[0] % DOFe;
                // This is raw averaged RDM, which should be re-normalized.
                // Complex_Matrix& RDM = RDM_average[0]; // This is raw averaged RDM. muted by cesare
                // When satrting from population, 1st renormalized way is used
                if (j == k) {
                    for (int n = 0; n < LEN_TRAJ; ++n) {// for (int n = 0; n < RDM_average[0].size(); ++n) { // n is step
                        double totalpop = 0; // get total population firstly
                        for (int i = 0; i < DOFe; ++i)
                            totalpop += RDM_average[0][n][i*DOFe+i].real();// totalpop += RDM[n][i*DOFe+i].real();
                        for (int a = 0; a < DOFe; a++)
                            for (int b = 0; b < DOFe; b++)
                                if (totalpop != 0) // population-population/coherence
                                    RDM_average[0][n][a*DOFe+b] /= totalpop; // RDM[n][a*DOFe+b] /= totalpop;
                    }
                }
                else {
                    std::cout << "WARNING: For SQC method, the re-normalization of RDM from "
                        "coherence state is not well defined, and the raw RDM which is averaged "
                        "by number of trajectory will be outputted." << std::endl;
                    // ! Note the normalization may not right when starting from coherence
                    // ! start from coherence, make sure σ_jk(0) = 1. (only for 2 state)
                    // ! for multi-state, here, /2.0 is not right.
                    //for (int n = 0; n < RDM.size(); ++n) // n is step
                    //    for (int a = 0; a < DOFe; a++)
                    //        for (int b = 0; b < DOFe; b++)
                    //            if (a == b && norm_pop[n] != 0) // coherence-population
                    //                RDM[n][a*DOFe+b] *= (double)(ntraj) / norm_pop[n];
                    //            else                            // coherence-coherence
                    //                RDM[n][a*DOFe+b] *= (double)(ntraj) / (norm_coh[n] / 2.0);
                    // This is muted by Zhubin
                }
            }             
        }
        else if (dyn_type == "eCMMcv") { // eCMMcv 
            double pop_check = 0; //check final total population
            for (int j = 0; j < DOFe; ++j) 
                pop_check += RDM_current[0][nsteps/RDM_steps][j*DOFe+j].real();
            double pop_tol_max = 10;
            if (pop_check > pop_tol_max) {
                traj_accepted -= 1;
                traj_rejected += 1;
                return -1;
            }
            for (int i = 0; i < num_init_states; ++i) 
                for (int step_count = 0; step_count < (nsteps/RDM_steps + 1); ++step_count) {
                    for (int j = 0; j < DOFe; ++j) {
                        for (int k = 0; k < DOFe; ++k) {
                            RDM_average[i][step_count][j*DOFe+k] += RDM_current[i][step_count][j*DOFe+k];
                        }
                    }
                }
        }
        else { // Other methods
            for (int i = 0; i < num_init_states; ++i) 
                for (int step_count = 0; step_count < (nsteps/RDM_steps + 1); ++step_count) 
                    for (int j = 0; j < DOFe; ++j) 
                        for (int k = 0; k < DOFe; ++k)  
                            RDM_average[i][step_count][j*DOFe+k] += RDM_current[i][step_count][j*DOFe+k]; //(double) (ntraj);        
        }
        if (representation == "adiabatic" && dyn_type == "TBSH") {
           for (int i = 0; i < num_init_states; ++i) {
               for (int step_count = 0; step_count < (nsteps/RDM_steps + 1); ++step_count) {
                   for (int j = 0; j < DOFe; ++j) {
                       for (int k = 0; k < DOFe; ++k) { 
                           AdiabaticRDM_average[i][step_count][j*DOFe+k] += AdiabaticRDM_current[i][step_count][j*DOFe+k]; //(double) (ntraj);
                       }
                   }
               }
           }
        }
        // check if this is the last update. If it is the last update, average the accumulated value. 
        int traj_total_now = traj_accepted + traj_rejected;
        double fail_traj_ratio = (double) (traj_rejected) / (double) traj_total_now;
        // const int begin_judge_traj_avail = 5000;  //Notice same index in simulation.cpp
        if (system_type == "model" && (traj_accepted + traj_rejected) > min_judge_traj_avail && dyn_type != "TBSH" && (traj_accepted == ntraj || fail_traj_ratio > max_friction_fail_traj )) {
            if (dyn_type == "LSC") {
                std::cout << "RDM Accepted final " << traj_accepted << "  RDM Rejected final " << traj_rejected << "  fail_traj_ratio " << fail_traj_ratio << std::endl;
                for (int i = 0; i < num_init_states; ++i) {
                    for (int step_count = 0; step_count < (nsteps/RDM_steps + 1); ++step_count) {
                        for (int j = 0; j < DOFe; ++j) {
                            for (int k = 0; k < DOFe; ++k) { 
                                RDM_average_LSC1[i][step_count][j*DOFe+k]   /= (double) (traj_accepted);
                                RDM_average_LSC2[i][step_count][j*DOFe+k]   /= (double) (traj_accepted);
                                RDM_average_RILSC1[i][step_count][j*DOFe+k] /= (double) (traj_accepted);
                                RDM_average_RILSC2[i][step_count][j*DOFe+k] /= (double) (traj_accepted);
                                RDM_average_RILSC3[i][step_count][j*DOFe+k] /= (double) (traj_accepted);
                            }
                        }
                    }   
                }
            }
            else { // Other methods
		        std::cout << "RDM Accepted final " << traj_accepted << "  RDM Rejected final " << traj_rejected << "  fail_traj_ratio " << fail_traj_ratio << std::endl;
                for (int i = 0; i < num_init_states; ++i) 
                    for (int step_count = 0; step_count < (nsteps/RDM_steps + 1); ++step_count) 
                        for (int j = 0; j < DOFe; ++j) 
                            for (int k = 0; k < DOFe; ++k)  
                                RDM_average[i][step_count][j*DOFe+k] /= static_cast<double> (traj_accepted); //(double) (ntraj);  
            }
        }
        if (system_type == "allatom" && dyn_type != "TBSH" && (traj_accepted == ntraj || fail_traj_ratio > max_friction_fail_traj )) {
            if (dyn_type == "LSC") {
                for (int i = 0; i < num_init_states; ++i) {
                    for (int step_count = 0; step_count < (nsteps/RDM_steps + 1); ++step_count) {
                        for (int j = 0; j < DOFe; ++j) {
                            for (int k = 0; k < DOFe; ++k) {
                                RDM_average_LSC1[i][step_count][j*DOFe+k]   /= (double) (traj_accepted);
                                RDM_average_LSC2[i][step_count][j*DOFe+k]   /= (double) (traj_accepted);
                                RDM_average_RILSC1[i][step_count][j*DOFe+k] /= (double) (traj_accepted);
                                RDM_average_RILSC2[i][step_count][j*DOFe+k] /= (double) (traj_accepted);
                                RDM_average_RILSC3[i][step_count][j*DOFe+k] /= (double) (traj_accepted);
                            }
                        }
                    }
                }
            }
            else { // Other methods
		        // std::cout<<traj_accepted<<" , traj_accepted "<< std::endl;
                for (int i = 0; i < num_init_states; ++i)
                    for (int step_count = 0; step_count < (nsteps/RDM_steps + 1); ++step_count)
                        for (int j = 0; j < DOFe; ++j)
                            for (int k = 0; k < DOFe; ++k)
                                RDM_average[i][step_count][j*DOFe+k] /= (double) (traj_accepted); //(double) (ntraj);
            }
        }
        if (system_type == "allatom")
            writeSingleTrajRDMFile();
        return 0;
    }
    else {
        traj_rejected += 1;
        // check if this is the last update. If it is the last update, average the accumulated value. 
        int traj_total_now = traj_accepted + traj_rejected;
        double fail_traj_ratio = (double) (traj_rejected) / (double) traj_total_now;
        // const int begin_judge_traj_avail = 5000;
        if (system_type == "model" && (traj_accepted + traj_rejected) >   min_judge_traj_avail && dyn_type != "TBSH" && (traj_accepted == ntraj || fail_traj_ratio > max_friction_fail_traj )) {
            if (dyn_type == "LSC") {
                std::cout << "RDM Accepted final " << traj_accepted << "  RDM Rejected final " << traj_rejected << "  fail_traj_ratio " << fail_traj_ratio << std::endl;
                for (int i = 0; i < num_init_states; ++i) {
                    for (int step_count = 0; step_count < (nsteps/RDM_steps + 1); ++step_count) {
                        for (int j = 0; j < DOFe; ++j) {
                            for (int k = 0; k < DOFe; ++k) { 
                                RDM_average_LSC1[i][step_count][j*DOFe+k]   /= (double) (traj_accepted);
                                RDM_average_LSC2[i][step_count][j*DOFe+k]   /= (double) (traj_accepted);
                                RDM_average_RILSC1[i][step_count][j*DOFe+k] /= (double) (traj_accepted);
                                RDM_average_RILSC2[i][step_count][j*DOFe+k] /= (double) (traj_accepted);
                                RDM_average_RILSC3[i][step_count][j*DOFe+k] /= (double) (traj_accepted);
                            }
                        }
                    }   
                }
            }
            else { // Other methods
		        std::cout << "RDM Accepted final " << traj_accepted << "  RDM Rejected final " << traj_rejected << "  fail_traj_ratio " << fail_traj_ratio << std::endl;
                for (int i = 0; i < num_init_states; ++i) 
                    for (int step_count = 0; step_count < (nsteps/RDM_steps + 1); ++step_count) 
                        for (int j = 0; j < DOFe; ++j) 
                            for (int k = 0; k < DOFe; ++k)  
                                RDM_average[i][step_count][j*DOFe+k] /= static_cast<double> (traj_accepted); //(double) (ntraj);            
            }
        }
        if (system_type == "allatom" && dyn_type != "TBSH" && (traj_accepted == ntraj || fail_traj_ratio > max_friction_fail_traj )) {
            if (dyn_type == "LSC") {
                for (int i = 0; i < num_init_states; ++i) {
                    for (int step_count = 0; step_count < (nsteps/RDM_steps + 1); ++step_count) {
                        for (int j = 0; j < DOFe; ++j) {
                            for (int k = 0; k < DOFe; ++k) {
                                RDM_average_LSC1[i][step_count][j*DOFe+k]   /= (double) (traj_accepted);
                                RDM_average_LSC2[i][step_count][j*DOFe+k]   /= (double) (traj_accepted);
                                RDM_average_RILSC1[i][step_count][j*DOFe+k] /= (double) (traj_accepted);
                                RDM_average_RILSC2[i][step_count][j*DOFe+k] /= (double) (traj_accepted);
                                RDM_average_RILSC3[i][step_count][j*DOFe+k] /= (double) (traj_accepted);
                            }
                        }
                    }
                }
            }
            else { // Other methods
		        // std::cout<<traj_accepted<<" , traj_accepted "<< std::endl;
                for (int i = 0; i < num_init_states; ++i)
                    for (int step_count = 0; step_count < (nsteps/RDM_steps + 1); ++step_count)
                        for (int j = 0; j < DOFe; ++j)
                            for (int k = 0; k < DOFe; ++k)
                                RDM_average[i][step_count][j*DOFe+k] /= (double) (traj_accepted); //(double) (ntraj);
            }
        }
        return -1;
    }
}

const std::vector<Complex_Matrix>& ObservableRDM::getAverageDensityMatrices() {
    if (dyn_type != "LSC") {
        return RDM_average;
    }
    else {
        RDM_average = RDM_average_LSC1;
        RDM_average.insert(RDM_average.end(), RDM_average_LSC2.begin(), RDM_average_LSC2.end());
        RDM_average.insert(RDM_average.end(), RDM_average_RILSC1.begin(), RDM_average_RILSC1.end());
        RDM_average.insert(RDM_average.end(), RDM_average_RILSC2.begin(), RDM_average_RILSC2.end());
        RDM_average.insert(RDM_average.end(), RDM_average_RILSC3.begin(), RDM_average_RILSC3.end());
        return RDM_average;
    }
}

const std::vector<Complex_Matrix>& ObservableRDM::getCurrentDensityMatrices() {
    if (dyn_type != "LSC") {
        return RDM_current;
    }
    else {
        RDM_current = RDM_current_LSC1;
        RDM_current.insert(RDM_current.end(), RDM_current_LSC2.begin(), RDM_current_LSC2.end());
        RDM_current.insert(RDM_current.end(), RDM_current_RILSC1.begin(), RDM_current_RILSC1.end());
        RDM_current.insert(RDM_current.end(), RDM_current_RILSC2.begin(), RDM_current_RILSC2.end());
        RDM_current.insert(RDM_current.end(), RDM_current_RILSC3.begin(), RDM_current_RILSC3.end());
        return RDM_current;
    }
}

const double& ObservableRDM::getCurrentFailRatio() {
    fail_traj_ratio = (double) (traj_rejected) / static_cast<double> (traj_accepted + traj_rejected);
    return fail_traj_ratio; 
}

int ObservableRDM::output() {
    if (dyn_type != "OpenMM") {
        // There are more than one RDMs in some simulation, so the filenames is
        // a vector, whose size will equal to the size of RDM got from dynamcis.
        const std::vector<Complex_Matrix>& RDM = getAverageDensityMatrices();
        std::vector<std::string> filenames;
        getRDMFileName(filenames);
        for (int i = 0; i < RDM.size(); ++i) {
            std::cout << "Generated electronic reduced density matrix file: " << filenames[i] << std::endl;
            writeRDMFile(RDM[i], filenames[i]);
        }
        total_traj = traj_accepted + traj_rejected;
        double reject_ratio = (double) (traj_rejected) / (double) (total_traj);
        std::cout << "trajectory acceptance-rejection ratio information  " << std::endl;
        std::cout << "trajectory accepted " << traj_accepted << ", trajectory rejected " << traj_rejected 
            << ", reject_ratio " << reject_ratio << std::endl;
        return 0;
    }
}

void ObservableRDM::getRDMFileName(std::vector<std::string>& filenames) {
    // * Get the file names of RDM file.
    // If the name is not specified by user, it is deduced from the model and dynamics
    // type and initial state: RDM_[model_type]_[dynamics_type]_[init_state].csv
    // Note: Here, the dynamics_type in file name is the name of specific version.
    // For example, for LSC dynamics, the dynamics_type will be LSC1/2, RI-LSC1/2/3.
    std::string basename = param.getStr("RDM_save"); // the basename of filenames
    std::string dynamics_type = dyn_type;
    const int DOFe = param.getInt("DOFe");
    // If the simulation successfully ended, we get the normal RDM file names
    if (traj_accepted == ntraj) {
        // For SQC/SPM methods, add the window or type name to dynamics type.
        if (dyn_type == "SQC")
            dynamics_type = dyn_type + "-" + param.getStr("SQC_window");
        else if (dyn_type == "SPM")
            dynamics_type = dyn_type + "-" + param.getStr("SPM_type") + "-" + param.getStr("elec_sample");
        // For LSC methods, there are 5 versions RDM files to generate.
        // The string "LSCREPLACE" will be replaced by specific version name later.
        else if (dyn_type == "LSC")
            dynamics_type = "LSCREPLACE";
        // Get initial state with two indices j,k
        std::vector<std::string> init_state;
        SplitString(init_state, param.getStr("init_state"));
        if (basename.empty()) { // get the default name
            if (system_type == "model")
                basename = "RDM_" +  model_type + "_" + dynamics_type;
            else if (system_type == "allatom")
                basename = "RDM_" +  allatom_type + "_" + dynamics_type;
            else
                basename = "RDM_" +  onthefly_type + "_" + dynamics_type;
            // For the case of init_state = all, the index of initial state will be added later
            if (init_state[0] != "all")
                basename = basename + "_" + init_state[0] + init_state[1];
        }
        else { // specified by user
            if (GetFileSuffix(basename) != "csv") // Make sure the format of file is csv file.
                std::cout << "WARNING: The extension name of RDM_save is not csv, but csv file will be generated." << std::endl;
            basename = GetFilePrefix(basename);
            // For LSC methods, we need to add a suffix to distinguish them even
            // if the filename is specified by user. The string "LSCREPLACE" will be
            // replaced by specific version name later.
            if (dyn_type == "LSC")
                basename = basename + "_" + "LSCREPLACE";
        }
        // * Generate the final file names of each RDM
        // the number of RDMs depends on the number of initial states and dynamcis methods.
        if (dyn_type == "LSC") { // generate RDMs with 5 LSC types
            std::string names[5] = {"LSC1", "LSC2", "RI-LSC1", "RI-LSC2", "RI-LSC3"};
            // If all initial state is computed, and we need add the suffix
            // of the initial state jk to the filename. Here, i = j*DOFe+k, i is the
            // index in the RDM vectors. The number of RDMs of each LSC type is DOFe*DOFe.
            // And total is DOFe*DOFe * 5.
            if (init_state[0] == "all") { // all initial state
                filenames.resize(DOFe*DOFe * 5);
                for (int n = 0; n < 5; ++n)  // 5 types LSC methods
                    for (int i = 0; i < DOFe*DOFe; ++i) // initial state
                        // We replace the "LSCREPLACE" with LSC type and add the
                        // suffix of initial state for them
                        filenames[n*DOFe*DOFe + i] = regex_replace(basename, std::regex("LSCREPLACE"), names[n]) +
                            "_" + std::to_string(i/DOFe) + std::to_string(i%DOFe) + ".csv";
            }
            else { // one initial state, 5 RDMs
                filenames.resize(5);
                for (int n = 0; n < 5; ++n) // 5 types LSC methods
                    filenames[n] = regex_replace(basename, std::regex("LSCREPLACE"), names[n]) + ".csv";
            }
        }
        // These methods also support to produce all RDM files in one simulation.
        // And they only produce one RDM file for one initial state.
        else if (dyn_type == "SPM" || dyn_type.substr(0, 3) == "CMM" || dyn_type == "eCMM" || dyn_type == "TBSH") {
            if (init_state[0] == "all") { // all initial state, have DOFe*DOFe RDMs
                filenames.resize(DOFe*DOFe);
                for (int i = 0; i < DOFe*DOFe; ++i)
                    filenames[i] = basename + "_" + std::to_string(i/DOFe) + std::to_string(i%DOFe) + ".csv";
            }
            else { // one initial state, one RDM.
                filenames.resize(1);
                filenames[0] = basename + ".csv";
            }
        }
        // One RDM file with one initial state: SQC, FSSH, MF, MF-RDM
        // These methods doesn't support init_state = all.
        else {
            filenames.resize(1);
            filenames[0] = basename + ".csv";
        }
    }
    // If the simulation is not successfully ended, we get the RDM file names 
    // averaged trajectory number.
    else { 
        // For SQC/SPM methods, add the window or type name to dynamics type.
        if (dyn_type == "SQC")
            dynamics_type = dyn_type + "-" + param.getStr("SQC_window");
        else if (dyn_type == "SPM")
            dynamics_type = dyn_type + "-" + param.getStr("SPM_type") + "-" + param.getStr("elec_sample");
        // For LSC methods, there are 5 versions RDM files to generate.
        // The string "LSCREPLACE" will be replaced by specific version name later.
        else if (dyn_type == "LSC")
            dynamics_type = "LSCREPLACE";
        // Get initial state with two indices j,k
        std::vector<std::string> init_state;
        SplitString(init_state, param.getStr("init_state"));
        if (basename.empty()) { // get the default name
            if (system_type == "model")
                basename = "Avg_" + std::to_string(traj_accepted) + "_RDM_" +  model_type + "_" + dynamics_type;
            else if (system_type == "allatom")
                basename = "Avg_" + std::to_string(traj_accepted) + "_RDM_" +  allatom_type + "_" + dynamics_type;
            else
                basename = "Avg_" + std::to_string(traj_accepted) + "_RDM_" +  onthefly_type + "_" + dynamics_type;
            // For the case of init_state = all, the index of initial state will be added later
            if (init_state[0] != "all")
                basename = basename + "_" + init_state[0] + init_state[1];
        }
        else { // specified by user
            if (GetFileSuffix(basename) != "csv") // Make sure the format of file is csv file.
                std::cout << "WARNING: The extension name of RDM_save is not csv, but csv file will be generated." << std::endl;
            basename = GetFilePrefix(basename);
            // For LSC methods, we need to add a suffix to distinguish them even
            // if the filename is specified by user. The string "LSCREPLACE" will be
            // replaced by specific version name later.
            if (dyn_type == "LSC")
                basename = basename + "_" + "LSCREPLACE";
        }
        // * Generate the final file names of each RDM
        // the number of RDMs depends on the number of initial states and dynamcis methods.
        if (dyn_type == "LSC") { // generate RDMs with 5 LSC types
            std::string names[5] = {"LSC1", "LSC2", "RI-LSC1", "RI-LSC2", "RI-LSC3"};
            // If all initial state is computed, and we need add the suffix
            // of the initial state jk to the filename. Here, i = j*DOFe+k, i is the
            // index in the RDM vectors. The number of RDMs of each LSC type is DOFe*DOFe.
            // And total is DOFe*DOFe * 5.
            if (init_state[0] == "all") { // all initial state
                filenames.resize(DOFe*DOFe * 5);
                for (int n = 0; n < 5; ++n)  // 5 types LSC methods
                    for (int i = 0; i < DOFe*DOFe; ++i) // initial state
                        // We replace the "LSCREPLACE" with LSC type and add the
                        // suffix of initial state for them
                        filenames[n*DOFe*DOFe + i] = regex_replace(basename, std::regex("LSCREPLACE"), names[n]) +
                            "_" + std::to_string(i/DOFe) + std::to_string(i%DOFe) + ".csv";
            }
            else { // one initial state, 5 RDMs
                filenames.resize(5);
                for (int n = 0; n < 5; ++n) // 5 types LSC methods
                    filenames[n] = regex_replace(basename, std::regex("LSCREPLACE"), names[n]) + ".csv";
            }
        }
        // These methods also support to produce all RDM files in one simulation.
        // And they only produce one RDM file for one initial state.
        else if (dyn_type == "SPM" || dyn_type.substr(0, 3) == "CMM" || dyn_type == "eCMM" || dyn_type == "TBSH") {
            if (init_state[0] == "all") { // all initial state, have DOFe*DOFe RDMs
                filenames.resize(DOFe*DOFe);
                for (int i = 0; i < DOFe*DOFe; ++i)
                    filenames[i] = basename + "_" + std::to_string(i/DOFe) + std::to_string(i%DOFe) + ".csv";
            }
            else { // one initial state, one RDM.
                filenames.resize(1);
                filenames[0] = basename + ".csv";
            }
        }
        // One RDM file with one initial state: SQC, FSSH, MF, MF-RDM
        // These methods doesn't support init_state = all.
        else {
            filenames.resize(1);
            filenames[0] = basename + ".csv";
        }
    }
}

void ObservableRDM::writeRDMFile(const Complex_Matrix& RDM, const std::string& file) {
    // Get the needed parameters for output.
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

void ObservableRDM::writeSingleTrajRDMFile() {
    static std::vector<std::string> filenames;
    getRDMFileName(filenames);
    for (int i = 0; i < RDM_current.size(); ++i) {
        // For multi-traj simulation, add a suffix "_traj?" to filename, where "?" is the index of traj (starting from 1).
        std::string filename = GetFilePrefix(filenames[i]) + "_traj" + std::to_string(traj_accepted) + "." + GetFileSuffix(filenames[i]);
        writeRDMFile(RDM_current[i], filename);
    }
    // Output RDM files for one single trajectories
}