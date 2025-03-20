/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu and Zengkui Liu @Sun Group @NYU-SH                       *
 * Last updated: Jan. 6, 2022                                                 *
 * -------------------------------------------------------------------------- */

#include "ObservableDES.h"

void ObservableDES::init() {
    ObservableBase::init();
    DES_type  = param.getStr("DES_type");
    DT = param.getDouble("DT");
    ha = std::static_pointer_cast<HamiltonianForceFieldBase>(Ha);
    DipoleCurrent.resize(3, 0.0);
    DipoleCurrentTensor.resize(3, 0.0);
    cosleft = 0;
    cosright = 0;
    Rdipoleleft = 0;
    Rdipoleright = 0;
}


int ObservableDES::updateOneSnapshot() {
    // get the every steps dipole size and direction
    std::cout<<"Test: update 2"<<std::endl;
    if (DES_type == "dipole") {
        int step = Dy->getStep();
        for (int index = 0; index < DOFe; index++) {
            //ha->updateDipoleRleftAright(index, cosleft, cosright, Rdipoleleft, Rdipoleright);
        }
    }
    else
        throw std::runtime_error("ERROR: Unsupported DES_type=" + DES_type + " for ObservableVibSpec.");

    return 0; 
}

int ObservableDES::updateOneTraj() {
    return 0; 
}

int ObservableDES::output() {
    std::vector<std::string> filenames;
    std::string basename = param.getStr("DES_save"); 
    const int DOFe = param.getInt("DOFe");
    filenames.resize(DOFe);
    for (int i = 0; i < DOFe; i++) {
        filenames[i] = basename + "_" + DES_type + "_" + std::to_string(i) + ".csv";
        std::cout << "The DES part output file: " << filenames[i] << std::endl;
    }
    // Save OKE result to file.
    for (int i = 0; i < DOFe; i++) {
        if (DES_type == "dipole") {
           writeDESdipoleResult(filenames[i]);
        }
    }
    std::cout << "End of the output file result" << std::endl;
    return 0;
}


const std::vector<double>& ObservableDES::getCurrentDipoleTensors() {
    if (dyn_type != "ReadTraj") {
        return DipoleCurrent;
    }
    else {
        DipoleCurrent = DipoleCurrentTensor;
        return DipoleCurrent;
    }
}

void ObservableDES::writeDESdipoleResult(std::string& file) {
    // Save Rxzxz to file.
    //FILE* outfile = CheckFile(file);
    //fprintf(outfile, "%.12g  %.12g  %.12g   %.12g", cosleft, cosright, Rdipoleleft, Rdipoleright);
    //fprintf(outfile, "\n");
    //fclose(outfile);
}