/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Jan. 6, 2022                                                 *
 * -------------------------------------------------------------------------- */

#pragma once
#include "DynamicsMF.h"
#include "DynamicsLSC.h"
#include "DynamicsSQC.h"
#include "DynamicsSPM.h"
#include "DynamicsCMM.h"
#include "DynamicsTTM.h"
#include "DynamicsECMM.h"
#include "DynamicsFSSH.h"
#include "DynamicsTBSH.h"
#include "DynamicsMFRDM.h"
#include "DynamicsECMMCV.h"
#include "DynamicsOpenMM.h"
#include "HamiltonianFMO.h"
#include "HamiltonianMSH.h"
#include "HamiltonianLVC.h"
#include "HamiltonianQVC.h"
#include "HamiltonianSpinBoson.h"
#include "HamiltonianNCMorse.h"
// #include "Observables.h"

/**

 */
class ObservableBase {
public:

    virtual ~ObservableBase() {}
    /**
     * Initialize data members and check the legality of parameters.
     */
    virtual void init();
    /**
     * update observables in one snapshot
     * 
     */
    virtual int updateOneSnapshot() = 0;
    /**
     *  
     * update observables in one trajectory
     */
    virtual int updateOneTraj() = 0;
    /**
     * output observables
     * 
     */
    virtual int output() = 0;
    /**
     * @brief Get the Average Density Matrix object
     * 
     */
    virtual const std::vector<Complex_Matrix>& getAverageDensityMatrices() {
        throw std::runtime_error("ERROR: This dynamics method doesn't have electronic reduced density matrix.");
    }
    /**
     * @brief Get the Current Density Matrix object
     * 
     */
    virtual const std::vector<Complex_Matrix>& getCurrentDensityMatrices() {
        throw std::runtime_error("ERROR: This dynamics method doesn't have electronic reduced density matrix.");
    }
    /*
    * Virtual function defined in the RDM class
    */
    virtual const double& getCurrentFailRatio() {
        throw std::runtime_error("ERROR: This dynamics method doesn't have electronic reduced density matrix.");
    }
    /*
    * Virtual function defined in the VibSpec class
    */
    virtual const std::vector<double>& getCurrentPiTensors() {
        throw std::runtime_error("ERROR: This dynamics method doesn't have the vibrational spectroscopy.");
    }   
    /*
    * Virtual function defined in the DES class
    */
    virtual const std::vector<double>& getCurrentDipoleTensors() {
        throw std::runtime_error("ERROR: This dynamics method doesn't have the DES obs.");
    }   
    /*
    * Virtual function defined in the VibSpec class
    */
    virtual const Vec3& getCurrentMuTensors() {
        throw std::runtime_error("ERROR: This dynamics method doesn't have the vibrational spectroscopy.");
    }  
    
protected:
    /**
     * This function is only used in the constructor of subclass, since you can
     * not construct an object for an abstract class.
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     * @param Dy      Dynamics object
     */
    ObservableBase(Parameters& param, std::shared_ptr<HamiltonianBase> Ha, 
        std::shared_ptr<DynamicsBase> Dy) : param(param), Ha(Ha), Dy(Dy) {}

protected:
    // Parameters object controls the simulation
    Parameters& param;
    // Hamiltonian object defines the system to be simulated
    std::shared_ptr<HamiltonianBase> Ha;
    // Dynamics object defines a method for simulating the system
    std::shared_ptr<DynamicsBase> Dy; 
    // dyn_type: dynamics method, OpenMM, LSC, SQC, ...
    std::string dyn_type;      
    // electronic DOF = number of states/topologies/surfaces.
    // For nonadiabatic dynamics, DOFe should be greater than 1.
    int DOFe;
    // the number of trajectories should be propagated.
    // For a general MD simultion, ntraj always equals 1 regardless of input value.
    // For an all-atom mapping dynamics, ntraj should be equal to the number of
    // nuclear sampling (if nucl_sample=none, then ntraj=1).
    int ntraj;
    // total steps for one trajectory
    int nsteps;
    // the time interval of observables. same unit as DT.
    double DeltaT;
    // the steps for the simualtion steps
    int nucl_start, nucl_skip, nucl_end;
    // number of traj has been accumulated
    int traj_count;
    // step count of MD simulation
    int step;
    // std::string dyn_type, post_type, system_type, onthefly_type, model_type;
};
