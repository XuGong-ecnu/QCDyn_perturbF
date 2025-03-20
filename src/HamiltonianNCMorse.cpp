/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zengkui Liu @Sun Group @NYU-SH                                     *
 * First created: Dec. 16, 2021                                               *
 * Last updated: Dec. 24, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "HamiltonianNCMorse.h"

void HamiltonianNCMorse::init() {
    // * 0. Initialization of data members and preset Hamiltonians
    HamiltonianModelBase::init();
    if (model_type.substr(0,7) != "NCMorse")
        throw std::runtime_error("ERROR: Undefined model_type=" + model_type + " for non-Condon Morse model.");
    Morse_preset = param.getStr("Morse_preset");
    if (Morse_preset != "preset1" & Morse_preset != "preset2" &  Morse_preset != "preset3" & Morse_preset != "read" ) 
        throw std::runtime_error("ERROR: Undefined Morse preset model parameter sets: " + Morse_preset );

    if (param.getBool("Condon_approximation"))
        throw std::runtime_error("ERROR: Please set Condon_approximation to false for NCMorse model.");
    if (DOFe != 3) 
        throw std::runtime_error("ERROR: Please set DOFe to 3 for this non-Condon Morse Oscillator Model");
    D    .resize(DOFe, 0);
    b    .resize(DOFe, 0);
    Req  .resize(DOFe, 0);
    C    .resize(DOFe, 0);
    A    .resize(DOFe, 0);
    a    .resize(DOFe, 0);
    Rcp  .resize(DOFe, 0);

    std::cout << "Das ist Ver. 2021 1224 " << std::endl;

    // * 1. give model parameters values with respect to the input parameters
    initializeModelParameters(); 
    // DOFn = 1;
    // Debug codes by Cesare
    // std::cout << " Notice here for non-Condon Morse potential DOFn has been redirected to one "  << std::endl;
    
    // reshape the Nuclear DOF variables
    // R.resize(DOFn, 0);
    // V.resize(DOFn, 0);
    // F.resize(DOFn, 0);
    // // masses.resize(DOFn, 0);
    // F_avg.resize(DOFn, 0);
    // F_all.resize(DOFe * DOFe, std::vector<double>(DOFn, 0));
    // NAC.resize((DOFe-1)*DOFe/2, std::vector<double>(DOFn, 0));
    // T.resize(DOFe, std::vector<double>(DOFe, 0));
    // T_dagger.resize(DOFe, std::vector<double>(DOFe, 0));

    shift.resize(DOFn, 0);
    omega.resize(DOFn, 0);
      
    // Initial sampling parameters for nuclear degree of freedom
    if (param.getStr("Morse_preset") == "preset1") {
        omega[0] = 0.005;   // in the unit of Hartree
        shift[0] = -2.9;    // in the unit of Hartree
    }
    else if (param.getStr("Morse_preset") == "preset2") {
        omega[0] = 0.005;   // in the unit of Hartree
        shift[0] = -3.3;    // in the unit of Hartree
    }
    else if (param.getStr("Morse_preset") == "preset3") {
        omega[0] = 0.005;   // in the unit of Hartree
        shift[0] = -2.1;    // in the unit of Hartree
    }
    // std::cout.precision(16);
    // debug codes by Cesare
    // std::cout << "unit 2 au " << unit2au << std::endl;
}

void HamiltonianNCMorse::buildModelParameters() {
    if (param.getStr("Morse_preset") == "preset1") {
        // Model parameters is taken from Table 1 of W.H. Miller's original papers
        // Chem. Phys. Lett. 349, 521 (2001)
        // Morse potential of 3 electronic states, |0> ,|1>, |2>.
        D[0] = 0.003;   b[0] = 0.65;    Req[0] = 5.0;   C[0] = 0.0;
        D[1] = 0.004;   b[1] = 0.60;    Req[1] = 4.0;   C[1] = 0.01;
        D[2] = 0.003;   b[2] = 0.65;    Req[2] = 6.0;   C[2] = 0.006;
        // Coupling parameters of 6 diagonal elements
        // Here we have index of 0 points to the parameters corresponding to
        // the coherence |0> <1| and |1><0|, index 1 points to |0><2| 
        // and |2><0|; index 3 points to |1><2| and |2><1|.
        A[0] = 0.002;   a[0] = 16.0;    Rcp[0] = 3.40;
        A[1] = 0.0;     a[1] = 0.0;     Rcp[1] = 0.0;
        A[2] = 0.002;   a[2] = 16.0;    Rcp[2] = 4.80;
    } 
    else if (param.getStr("Morse_preset") == "preset2") {
        // Model parameters is taken from Table 2 of W.H. Miller's original papers:
        // Chem. Phys. Lett. 349, 521 (2001)
        // The physical properties are all in atomic unit (Hartree unit system)
        // Morse potential of 3 electronic states, |0> ,|1>, |2>.
        D[0] = 0.02;    b[0] = 0.65;    Req[0] = 4.5;   C[0] = 0.0;
        D[1] = 0.01;    b[1] = 0.40;    Req[1] = 4.0;   C[1] = 0.01;
        D[2] = 0.003;   b[2] = 0.65;    Req[2] = 4.4;   C[2] = 0.02;
        // Coupling parameters of 6 diagonal elements
        // Here we have index of 0 points to the parameters corresponding to
        // the coherence |0> <1| and |1><0|, index 1 points to |0><2| 
        // and |2><0|; index 3 points to |1><2| and |2><1|.
        A[0] = 0.005;   a[0] = 32.0;    Rcp[0] = 3.66;
        A[1] = 0.005;   a[1] = 32.0;    Rcp[1] = 3.34;
        A[2] = 0.0;     a[2] = 0.0;     Rcp[2] = 0.0;
    }
    else if (param.getStr("Morse_preset") == "preset3") {
        // Model parameters is taken from Table 3 of W.H. Miller's original papers
        // Chem. Phys. Lett. 349, 521 (2001)
        // Morse potential of 3 electronic states, |0> ,|1>, |2>.
        D[0] = 0.02;    b[0] = 0.40;    Req[0] = 4.0;   C[0] = 0.02;
        D[1] = 0.02;    b[1] = 0.65;    Req[1] = 4.5;   C[1] = 0.0;
        D[2] = 0.003;   b[2] = 0.65;    Req[2] = 6.0;   C[2] = 0.02;
        // Coupling parameters of 6 diagonal elements
        // Here we have index of 0 points to the parameters corresponding to
        // the coherence |0> <1| and |1><0|, index 1 points to |0><2| 
        // and |2><0|; index 3 points to |1><2| and |2><1|.
        A[0] = 0.005;   a[0] = 32.0;    Rcp[0] = 3.4;
        A[1] = 0.005;   a[1] = 32.0;    Rcp[1] = 4.97;
        A[2] = 0.0;     a[2] = 0.0;     Rcp[2] = 0.0;
    }
}

void HamiltonianNCMorse::loadModelParameters(const std::string& loadfile) {
    // std::vector
    throw std::runtime_error("ERROR: Constructing non-Condon Morse Oscillator Model is not supported yet");
}

void HamiltonianNCMorse::saveModelParameters(const std::string& savefile) {
    std::vector <std::vector<double>*> data = {&D, &b, &Req, &C, &A, &a, &Rcp};
    std::vector<std::string> headers = {"D/au", "b", "req", "C", "A", "a", "rcp"};
    SaveDataFile(savefile, headers, data);
}


double HamiltonianNCMorse::getPotentialEnergy(int index){
    if (index >= DOFe || index < 0)
        throw std::runtime_error("ERROR: State index out of scope.");
    double part = 0.0;
    double X = 0.0;
    double PE = 0.0;
    X  = R[0] - Req[index];

    // debug code by Cesare
    // std::cout << "nuclear position in potential update " << R[0] << std::endl;
    // std::cout << "index " << index << ", X " << X << ", V " << V[0] ;
    
    part = 1.0 - exp(- X * b[index]);
    // std::cout << ", part " << part ;
    PE = D[index] * part * part + C[index];
    // std::cout << ", C " << C[index] ;
    // debug code by Cesare
    // std::cout << ", PE " << PE << std::endl;
    // std::cout << std::endl;
    return PE;
}

void HamiltonianNCMorse::getForces(int i, int j, std::vector<double>& F){
    if (i >= DOFe || i < 0 || j >= DOFe || j < 0)
        throw std::runtime_error("ERROR: Called getForces() with wrong state index.");
    F.resize(DOFn, 0);
    double X = 0.0;
    double e1 = 0.0;
    if (i == j) { // Diagonal force
        for (int k = 0; k < DOFn; k++) {
            X    = R[k] - Req[i];
            e1   = exp (- b[i] * X);
            F[k] = 2.0 * D[i] * b[i] * (e1 * e1 - e1) / M;
        }
    }
    else {
        for (int k = 0; k < DOFn; k++) { // Force provided by the coupling term 
            int index = i < j ? ((2*DOFe-i-1)*i/2 + j-i - 1) : ((2*DOFe-j-1)*j/2 + i-j - 1);
            X    = R[k] - Rcp[index];
            e1   = exp (- a[index] * X * X);
            F[k] = 2.0 * A[index] * a[index] * X * e1 / M;
        }
    }
    
    // debug code by Cesare
    // std::cout << "ij " << i << ", " << j << ", " << F[0] << std::endl;
    // std::cout << std::endl;
    //So here we return the acceleration of the nuclear position
}

double HamiltonianNCMorse::getNonCondonCoupling(int i, int j) {
    if (i >= DOFe || i < 0 || j >= DOFe || j < 0 || i == j)
        throw std::runtime_error("ERROR: Called getNonCondonCoupling() with wrong state index.");
    // H_ji = H_ij = Aij exp(-a_ij (x-Rcp[i][j])) 
    // where A_ij is the coupling coefficient.
    int index = i < j ? ((2*DOFe-i-1)*i/2 + j-i - 1) : ((2*DOFe-j-1)*j/2 + i-j - 1);
    double coupling = 0.0;
    double X = 0.0;
    double e1 = 0.0;
    X =  R[0] - Rcp[index];
    e1 =  - X * X * a[index];
    coupling = A[index] * exp(e1);

    // Debug codes by Cesare
    // std::cout << "index_ij " << i << j << " coupling " << coupling << ", A " << A[index] << ", a " << a[index] << ", Rcp " << Rcp[index]  << ", X = " << X << std::endl;
    // std::cout << "position nucl at present " << R[0] << ", Δ(x-R_ij) = " << X << std::endl;
    // std::cout << std::endl;
    return coupling;
}

