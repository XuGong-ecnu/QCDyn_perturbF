/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu and Zengkui Liu @Sun Group @NYU-SH                       *
 * Last updated: Jan. 6, 2022                                                 *
 * -------------------------------------------------------------------------- */

#include "ObservableVibSpec.h"

void ObservableVibSpec::init() {
    ObservableBase::init();
    PiCurrentTensor.resize(9);
    PiCurrent.resize(9);
    NPiSteps = nsteps;
    PiAll.resize(9 * NPiSteps);
    MuAll.resize(3 * NPiSteps);
    OKE_type  = param.getStr("OKE_type");
    DT = param.getDouble("DT");
    ha = std::static_pointer_cast<HamiltonianForceFieldBase>(Ha);
    // organize DOFe ObservableTCF class, put the param to the ObservableTCF
//    OKETCF.resize(DOFe, ObservableTCF(param));
    // Initialize data members.
//    for (int i = 0; i < DOFe; i++) {
//        OKETCF[i].init();
//    }
    OKETCF.resize(DOFe);
for (int i = 0; i < DOFe; i++) {
    OKETCF[i] = std::make_shared<ObservableTCF>(param);
    OKETCF[i]->init();
}

}

int ObservableVibSpec::updateOneSnapshot() {
    // get the every steps polarizability and calculate the Rxzxz
    if (OKE_type == "eqRxzxz") {
        int step = Dy->getStep();
        std::vector<double> Pi_tensor;
        Pi_tensor.resize(9, 0.0);
        for (int index = 0; index < DOFe; index++) {
            ha->updatePiTensor(index, Pi_tensor);
       //     OKETCF[index].setEachStepPiTensor(step, Pi_tensor);
            OKETCF[index]->setEachStepPiTensor(step, Pi_tensor); 
        }
    }
    else if (OKE_type == "noneqRxzxz") {
        int step = Dy->getStep();
        std::vector<double> Pi_tensor;
        Pi_tensor.resize(9, 0.0);
        bool direction;
        for (int index = 0; index < DOFe; index++) {
            ha->updatePiTensor(index, Pi_tensor);
            bool direction = ha->getDirection();
//            if (direction) 
//                OKETCF[index].setEachPositiveStepPiTensor(step, Pi_tensor);
//            else {
//                OKETCF[index].setEachNegtiveStepPiTensor(step, Pi_tensor);
//            }
        }
    }
    else
        throw std::runtime_error("ERROR: Unsupported OKE_type=" + OKE_type + " for ObservableVibSpec.");
       
    return 0; 
}

int ObservableVibSpec::updateOneTraj() {
    return 0; 
}

int ObservableVibSpec::output() {
    std::vector<std::string> filenames;
    std::string basename = param.getStr("VibSpec_save"); 
    const int DOFe = param.getInt("DOFe");
    filenames.resize(DOFe);
    for (int i = 0; i < DOFe; i++) {
        filenames[i] = basename + "_" + OKE_type + "_" + std::to_string(i) + ".csv";
        std::cout << "The VibSpec part output file: " << filenames[i] << std::endl;
    }
    // Save OKE result to file.
    for (int i = 0; i < DOFe; i++) {
        if (OKE_type == "eqRxzxz") {
//           OKETCF[i].writeEqRxzxzResult(filenames[i]);
           OKETCF[i]->writeEqRxzxzResult(filenames[i]);
        }
        if (OKE_type == "noneqRxzxz") {
//           OKETCF[i].writeNoneqRxzxzResult(filenames[i]);
        }
    }
    std::cout << "End of the output file result" << std::endl;
    return 0;
}

const Vec3& ObservableVibSpec::getCurrentMuTensors() {
    if (dyn_type != "ReadTraj") {
        return MuCurrent;
    }
    else {
        MuCurrent = MuCurrentTensor;
        return MuCurrent;
    }
}

const std::vector<double>& ObservableVibSpec::getCurrentPiTensors() {
    if (dyn_type != "ReadTraj") {
        return PiCurrent;
    }
    else {
        PiCurrent = PiCurrentTensor;
        return PiCurrent;
    }
}
