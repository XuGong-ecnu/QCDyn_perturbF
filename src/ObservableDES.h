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

class ObservableDES : public ObservableBase {
public:
    /**
     * Construct an object of ObservableDES.
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     * @param Dy      Dynamics object
     */
    ObservableDES(Parameters& param, std::shared_ptr<HamiltonianBase> Ha, 
        std::shared_ptr<DynamicsBase> Dy) : ObservableBase(param, Ha, Dy) {}
    ~ObservableDES() {}
    /**
     * Initialize data members and check the legality of parameters.
     */
    void init();
    /**
     * Here data members in Spectra_current are updated for one Spectra snapshot,
     * where Spectra information is taken out of the temper storage from
     * the on-the-fly calculated Spectra data member in the dynamics base.
     * If we get not-a-number of infinity in the Spectra, it will return
     * -1 and this trajectories will not be recorded; else it will return 
     * zero. 
     */
    int updateOneSnapshot();
    /**
     *  If Spectra_current is not abandoned, it will be recorded into 
     *  Spectra_average. If system is all_atom and dynamic type is 
     *  non-adiabatic, Spectra of each trajectory will be written. 
     */
    int updateOneTraj();
    /**
     * Write Spectra output file in the end of QCDyn.
     */
    int output();
        /**
     * @brief Get the Average Density Matrix object
     * 
     */
    /**
     * Uing the eq method to get the writeDESdipoleResult result file.
     */
    void writeDESdipoleResult(std::string& file);
    const std::vector<double>& getCurrentDipoleTensors() override;
protected:
    
private:
    std::vector<double> DipoleCurrent;
    std::vector<double> DipoleCurrentTensor;
    std::shared_ptr<HamiltonianForceFieldBase> ha;
    std::string DES_type;
    int NPiSteps;
    double cosleft;
    double cosright;
    double Rdipoleleft;
    double Rdipoleright;
    double DT;
};