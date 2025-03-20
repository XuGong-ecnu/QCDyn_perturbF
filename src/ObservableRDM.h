/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "ObservableBase.h"

class ObservableRDM : public ObservableBase {
public:
    /**
     * Construct an object of ObservableRDM.
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     * @param Dy      Dynamics object
     */
    ObservableRDM(Parameters& param, std::shared_ptr<HamiltonianBase> Ha, 
        std::shared_ptr<DynamicsBase> Dy) : ObservableBase(param, Ha, Dy) {}
    ~ObservableRDM() {}
    /**
     * Initialize data members and check the legality of parameters.
     */
    void init();
    /**
     * Here data members in RDM_current are updated for one RDM snapshot,
     * where RDM information is taken out of the temper storage from
     * the on-the-fly calculated RDM data member in the dynamics base.
     * If we get not-a-number of infinity in the RDM, it will return
     * -1 and this trajectories will not be recorded; else it will return 
     * zero. 
     */
    int updateOneSnapshot();
    /**
     *  If RDM_current is not abandoned, it will be recorded into 
     *  RDM_average. If system is all_atom and dynamic type is 
     *  non-adiabatic, RDM of each trajectory will be written. 
     */
    int updateOneTraj();
    /**
     * Write RDM output file in the end of QCDyn.
     */
    int output();
    /**
     * @brief Get the Average Density Matrix object
     * 
     */
    const std::vector<Complex_Matrix>& getAverageDensityMatrices() override;
    /**
     * @brief Get the Current Density Matrix object
     * 
     */
    const std::vector<Complex_Matrix>& getCurrentDensityMatrices() override;
    /*
    *  Get maximum tolerance percentage of current failed trajectories
    */
    const double& getCurrentFailRatio() override;
    
private:
    int RDM_steps;  // the frequecny (in time step) to obtain RDM
    int num_init_states;    // number of initial states
    int LEN_TRAJ;   // the number of obsevable (RDM) snapshots 
    int total_traj; // total trajectory numbers during the simulation (including abandoned and restored)
    int traj_rejected;    // internal count of the trajectories which are rejected 
    int traj_accepted;         // total number of accepted trajectories 
    // minimun trajectories begin to decide if the trajectories is counted 
    int min_judge_traj_avail; 
    // RDM stores the reduced density matrix along the time (final results).
    // the density matrix of one step is stored as DOFe^2-dimentional vector.
    // RDM_current represents the RDM of current trajectory, and RDM_average
    // represents the ensemble avearge RDM, in which the values have been
    // accumaleated and averaged by number of trajectories.
    // Since the density matrice that starting from differnent initial states
    // can be computed at one simulation for some methods, therefore it is a
    // vector of RDM. e.g., RDM[i][t][f], i is init_state, t is the
    // step (time = DT*step), and f is the index of element in rho(t).
    // Here, i and f is, i/f = a*DOFe+b, where a, b is the subscripts of rho.
    std::vector<Complex_Matrix> RDM_current, RDM_average;
    std::vector<Complex_Matrix> AdiabaticRDM_current, AdiabaticRDM_average;
    /*Normalization factor for the SQC method, which can be executed for each traj or once for all*/ 
    std::vector<double>  norm_pop, norm_coh;
    // Reduced density matrix including all types genenrated by LSC method.
    // Note that they will be used only for LSC dynamics methods.
    std::vector<Complex_Matrix> RDM_current_LSC1, RDM_current_LSC2, RDM_current_RILSC1, RDM_current_RILSC2, RDM_current_RILSC3;
    std::vector<Complex_Matrix> RDM_average_LSC1, RDM_average_LSC2, RDM_average_RILSC1, RDM_average_RILSC2, RDM_average_RILSC3; 
    std::vector<Complex_Matrix> AdiabaticRDM_current_LSC1, AdiabaticRDM_current_LSC2, AdiabaticRDM_current_RILSC1, AdiabaticRDM_current_RILSC2, AdiabaticRDM_current_RILSC3;
    std::vector<Complex_Matrix> AdiabaticRDM_average_LSC1, AdiabaticRDM_average_LSC2, AdiabaticRDM_average_RILSC1, AdiabaticRDM_average_RILSC2, AdiabaticRDM_average_RILSC3; 
    std::string dyn_type, post_type, system_type, onthefly_type, model_type, allatom_type, representation;   
    std::vector<int> init_state;
    /*
    * name RDM Files 
    */
    void getRDMFileName(std::vector<std::string>& filenames);
    /* 
    * write RDM Files
    */
    void writeRDMFile(const Complex_Matrix& RDM, const std::string& file);
    /*
    * write RDM Files for one single trajectories for all-atom simulation
    */
    void writeSingleTrajRDMFile();
    /*
    * Flag of trajectory numerical is finite: if current member in this trajectory is finite, return true; else return false
    */
    bool flag_traj_avail;
    /**
     * @brief Flag of doing average now: if current fail ratio reaches maximum or the traj_accepts reachs ntraj 
     * 
     */
    bool taking_average_now;
    // maximum friction of failed trajectories
    double max_friction_fail_traj;
    // friction of failed trajectories
    double fail_traj_ratio;
};
