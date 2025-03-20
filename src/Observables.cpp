/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 23, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "Observables.h"

void Observables::init() {
    // Get parameters to control the computation of observables
    dyn_type        = param.getStr("dyn_type");
    post_type       = param.getStr("post_type");
    system_type     = param.getStr("system_type");
    model_type      = param.getStr("model_type");
    onthefly_type   = param.getStr("onthefly_type");
    allatom_type    = param.getStr("allatom_type");
    DOFe            = param.getInt("DOFe");
    // get "," separated keys of observables
    std::string observables_names = param.getStr("observables");
    SplitString(keys, observables_names);
    // Add calculation of RDM for nonadiabatic dynamics if not exist 
    bool hasRDM = std::find(keys.begin(), keys.end(), "RDM") != keys.end();
    if ((dyn_type !=  "OpenMM" && dyn_type !=  "MD" && DOFe > 1) && !hasRDM) 
        keys.emplace_back("RDM");
    // construct object and call init() of each observable subclass
    for (int i = 0; i < keys.size(); i++) {
        // construct the object of specific observable
        if (keys[i] == "RDM") {
            observables[keys[i]] = std::make_shared<ObservableRDM>(param, Ha, Dy);
        }
        else if (keys[i] == "OKE") {
            //Enable this when writing this actual case by Xiaofang
            observables[keys[i]] = std::make_shared<ObservableVibSpec>(param, Ha, Dy); 
        }
        else if (keys[i] == "DES") {
            //Enable this when writing this actual case by Xiaofang
            observables[keys[i]] = std::make_shared<ObservableDES>(param, Ha, Dy); 
        }
        // TO DO
        //else if (keys[i] == "GreenKubo")
            //Enable this line by Xianfang
            //observables[keys[i]] = std::make_shared<ObservableGreenKubo>(param, Ha, Dy);            
        else 
            throw std::runtime_error("ERROR: Unsupported observable=" + keys[i]); 
        // call init() of current observable object
        observables[keys[i]]->init();   
    }
}

int Observables::updateOneSnapshot() {
    std::vector<int> tag_valid(keys.size(), 0);
    for (int i = 0; i < keys.size(); i++) {
        tag_valid[i] = observables[keys[i]]->updateOneSnapshot();
        if (tag_valid[i] != 0)
            return -1; // failed
    }
    return 0; // success
}

int Observables::updateOneTraj() {
    std::vector<int> tag_valid(keys.size(), 0);
    for (int i = 0; i < keys.size(); i++) {
        tag_valid[i] = observables[keys[i]]->updateOneTraj();
        if (tag_valid[i] != 0)
            return -1; // failed
    }
    return 0; // success
}

int Observables::output() {
    std::vector<int> tag_valid(keys.size(), 0);
    for (int i = 0; i < keys.size(); i++) 
        tag_valid[i] = observables[keys[i]]->output();
    for (int i = 0; i < tag_valid.size(); i++) 
        if (tag_valid[i] != 0)
            return -1; // failed
    return 0; // success
}

const std::vector<Complex_Matrix>& Observables::getAverageDensityMatrix() {
    // std::vector<Complex_Matrix> tag_valid;
    for (int i = 0; i < keys.size(); i++) 
        if (keys[i] == "RDM")
            return observables[keys[i]]->getAverageDensityMatrices();
}

const std::vector<Complex_Matrix>& Observables::getCurrentDensityMatrix() {
    // std::vector<Complex_Matrix> tag_valid;
    for (int i = 0; i < keys.size(); i++) 
        if (keys[i] == "RDM")
            return observables[keys[i]]->getCurrentDensityMatrices();
}

const double& Observables::getCurrentFailRatio() {
    // std::vector<Complex_Matrix> tag_valid;
    for (int i = 0; i < keys.size(); i++) 
        if (keys[i] == "RDM")
            return observables[keys[i]]->getCurrentFailRatio();
}

const std::vector<double>& Observables::getCurrentPiTensor() {
    // std::vector<Complex_Matrix> tag_valid;
    for (int i = 0; i < keys.size(); i++) 
        if (keys[i] == "OKE")
            return observables[keys[i]]->getCurrentPiTensors();
}

const std::vector<double>& Observables::getCurrentPiDipoleTensor() {
    // std::vector<Complex_Matrix> tag_valid;
    for (int i = 0; i < keys.size(); i++) 
        if (keys[i] == "DES")
            return observables[keys[i]]->getCurrentDipoleTensors();
}

const Vec3& Observables::getCurrentMuTensor() {
    // std::vector<Complex_Matrix> tag_valid;
    for (int i = 0; i < keys.size(); i++) 
        if (keys[i] == "SPSP")
            return observables[keys[i]]->getCurrentMuTensors();
}