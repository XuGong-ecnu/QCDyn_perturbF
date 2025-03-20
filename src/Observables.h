/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Nov. 26, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "Tools.h"
#include "ObservableBase.h"
#include "ObservableRDM.h"
#include "ObservableVibSpec.h"
#include "ObservableDES.h"

class Observables {
public:
    /**
     * Construct a Observables object.
     *
     * @param param the global simulation parameters
     */
    Observables(Parameters& param, std::shared_ptr<HamiltonianBase> Ha, 
        std::shared_ptr<DynamicsBase> Dy) : param(param), Ha(Ha), Dy(Dy) {}
    ~Observables() {}  
    void init();
    int updateOneSnapshot();
    int updateOneTraj();
    int output();
    // For electronic reduced density matrix.
    const std::vector<Complex_Matrix>& getAverageDensityMatrix();
    const std::vector<Complex_Matrix>& getCurrentDensityMatrix();
    const double& getCurrentFailRatio();
    // For vibrational spectroscopy.
    const std::vector<double>& getCurrentPiTensor();
    const std::vector<double>& getCurrentPiDipoleTensor();
    const Vec3& getCurrentMuTensor();
 
private:
    // Parameters object controls the simulation
    Parameters& param;
    // Hamiltonian object defines the system to be simulated
    std::shared_ptr<HamiltonianBase> Ha;
    // Dynamics object defines a method for simulating the system
    std::shared_ptr<DynamicsBase> Dy; 
    // dyn_type: dynamics method, OpenMM, LSC, SQC, ...
    // post_type: postprocess type, TTM
    std::string dyn_type, post_type, system_type, onthefly_type, model_type, allatom_type;
    // The map is used to store the keys and values of Observables.
    // The vector is used to stores the keys only with the adding order of Observables.  
    std::map<std::string, std::shared_ptr<ObservableBase>> observables;
    std::vector<std::string> keys;
    int DOFe;
};
