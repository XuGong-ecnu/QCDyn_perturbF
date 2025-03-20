/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsTTM.h"

void DynamicsTTM::init() {
    // * 0. Initialize date members from input paramters.
    DOFe   = param.getInt("DOFe");
    DT     = param.getDouble("DT");
    nsteps = param.getInt("nsteps");
    k_max  = param.getInt("TTM_kmax");
    if (DOFe < 2 || DT <= 0 || k_max < 1 || nsteps < k_max)
        throw std::runtime_error("ERROR: For transfer tensor method (TTM), found illegal "
            "value for DOFe (requires > 1), or nsteps (requires > 0), or DT (requires > 0), or k_max (requires > 1)");
    // Get the initial diabtic density matrix element indices |j><k| specified
    // by user, the input format is "j,k", e.g. 0,0, then the initial diabatic
    // density matrix will be initialized as: rho_00 is 1 and others are 0
    // If init_state=all, then all possible initial states are considered.
    // Note, this the init_state rho(0) used for propagation by TTM.
    const std::string init_state_str = param.getStr("init_state");
    if (init_state_str == "all") { // all possible initial states, DOFe*DOFe
        init_state.resize(DOFe*DOFe, 0);
        for (int i = 0; i < DOFe*DOFe; ++i) // index of initial state i = j*DOFe+k
            init_state[i] = i;
    }
    else { // specfied one initial state
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
    // * 1. Initialize the vectors or matrice.
    Complex_Matrix tmp(DOFe*DOFe, std::vector<Complex>(DOFe*DOFe, 0));
    DM.resize(k_max+1, tmp); // k = 0,1,2,...,k_max, so the size is k_max+1
    TT = DM;
    // The size of RDM_output depends on the initial states (one or all)
    // and nsteps+1 is from time 0.
    RDM_output.resize(init_state.size(), Complex_Matrix(nsteps+1, std::vector<Complex>(DOFe*DOFe, 0.0)));
    // * 2. Load input reduced density matrice (if not assigned by constructor)
    if (RDM_input.empty())
        loadDensityMatrix();
    else { // simple check the valid of input density matrice (assigned by constructor)
        if (RDM_input.size() != DOFe*DOFe)
            throw std::runtime_error("ERROR: For transfer tensor method (TTM), all "
                "density matrice strating from all possible inital states (rho(0)) should be provided.");
        if (RDM_input[0].size() < k_max)
            throw std::runtime_error("ERROR: For transfer tensor method (TTM), the steps of"
                "input density matrice used for learning should be greater than k_max.");
    }
}

void DynamicsTTM::loadDensityMatrix() {
    // * Initialize RDM_input firstly
    RDM_input.resize(DOFe*DOFe, Complex_Matrix(k_max+1, std::vector<Complex>(DOFe*DOFe, 0)));
    // * Load input density matrice from files "basename_ab.csv"
    // the names of input reduced density matrix file should be: basename_ab.csv,
    // where ab is the indice of init_state (rho(0)). The files should follow rule:
    // 1) 1st column is headers; 2) 1st row is time and other rows is the values
    // of real and imag part of density matrix elements. 3) each column should be
    // separated by a comma ",". 4) So the number of columns of each row should be
    // 1 + DOFe*DOFe*2, (*2 is from the real and imag part).
    std::string basename = param.getStr("RDM_load");
    if (basename.empty())
        throw std::runtime_error("ERROR: The basename of input density matrice for TTM method is not provided.");
    const int cols = 1 + DOFe*DOFe*2; // the number of columns of each row
    for (int i = 0; i < DOFe*DOFe; ++i) { // init_state=a*DOFe+b, a=i/DOFe, b= i%/DOFe
        std::string file = basename + "_" + std::to_string(i/DOFe) + std::to_string(i%DOFe) + ".csv";
        std::ifstream inpfile(file);
        if (!inpfile)
                throw std::runtime_error("ERROR: Can't open input reduced density matrix file: " + file);
        std::string line;
        getline(inpfile, line); // skip the 1st row (headers) firstly
        int index = 0; // index = 0,1,2,...,k_max
        for ( ; getline(inpfile, line); ++index) {
            std::vector<double> fields;
            SplitString(fields, line);
            if (fields.size() < cols)
                throw std::runtime_error("ERROR: Too few columns at line " + std::to_string(index+2));
            for (int f = 0; f < DOFe*DOFe; ++f) // final_state (rho(t)) = a*DOFe+b
                // the index of columns are 1:rho00re, 2:rho00im, 3:rho01re, 4:rho01im, ....
                // construct a nameless complex, then assign it to RDM_input
                RDM_input[i][index][f] = Complex(fields[f*DOFe+1], fields[f*DOFe+2]);
            if (index == k_max) break;
        }
        if (index != k_max) // check if suficent steps is provied
            throw std::runtime_error("ERROR: Too few rows of data in " + file);
    }
}

void DynamicsTTM::buildDynamicalMaps() {
    // # Construct dynamical maps (DM) from input density matrix (RDM_input)
    // Satisfy: ρ(k) = ε_k*ρ(0) Ref: J. Phys. Chem. Let. 2016, 7, 4809. Eqs. (1)
    // Here, DOFe^2-dimentional vector ρ(0), ρ(k) is the ρ from input denisty
    // matrice at time 0 and step k (t=DT*k). ε_k is a DOFe^2*DOFe^2 matrix.
    // And, the value (complex) of DM[k][f][i] is the element of input density
    // matrix at step k, with final state f (rho(t)) and inital state i (rho(0)),
    // in which i/f = a*DOFe+b, where a, b is the subscripts of rho.
    for (int k = 0; k <= k_max; ++k)  // index of step, k = 0,1,..,k_max
        for (int i = 0; i < DOFe*DOFe; ++i) // index of inital state
            for (int f = 0; f < DOFe*DOFe; ++f) // index of final state
                DM[k][f][i] = RDM_input[i][k][f];
}

void DynamicsTTM::buildTransferTensors() {
    // # Construct transfer tensors (TT) from dynmacal maps (DM):
    // # T_n = ε_n - Σ_{m=1,n-1}T_{n-m}*ε_m
    // Ref: J. Phys. Chem. Let. 2016, 7, 4809. Eq. (2)
    // Note: n start from 1 in paper, while n start from 0 in C++, but the
    // T_0 = ε_0 is useless in the propagation, so the index/notation used here
    // is same as in paper.
    TT = DM; // let T_n = ε_n firstly and T_1 ≡ ε_1
    for (int n = 2; n <= k_max; ++n) // loop from T_2, since T_1 ≡ ε_1
        for (int m = 1; m < n; ++m) {
            Complex_Matrix tmp;
            Matrix_Multiply(TT[n-m], DM[m], tmp); // tmp = T_{n-m}*ε_m
            Matrix_Add(TT[n], tmp, TT[n], -1.0); // T_n -= tmp
        }
    // Check the if TT[k_max] is close to 0, otherwise the k_max is too small
    std::cout << "The transfer tensor T[k_max] (a complex matrix) is:\n";
    for (int i = 0; i < DOFe*DOFe; ++i) { // index of inital state
        for (int f = 0; f < DOFe*DOFe; ++f) // index of final state
            std::cout << std::fixed << std::setprecision(8) << std::setw(26) << TT[k_max][i][f];
        std::cout << "\n";
    }
    std::cout << "Please check the above transfer tensor carefully.\n";
    std::cout << "If it is not close to 0, the results of TTM may not be resonable." << std::endl;
}

void DynamicsTTM::propagate() {
    // # Propagate the density matrix using transfer tensors and dynmacal maps
    // # if (n<k_max): ρ(n) = ε_n*ρ(0); else ρ(n) = Σ_{k=1,k_max}T_k*ρ(n-k)
    // Ref: J. Phys. Chem. Let. 2016, 7, 4809. Eqs. (1) and (3)
    // Note: here, ρ(n) is a DOFe^2-dimentional vector (DOFe*DOFe density matrix)
    for (int i = 0; i < init_state.size(); ++i) // init_state i = j*DOFe+k
        for (int n = 0; n <= nsteps; ++n)  // when n=0, ρ(n) is ρ(0)
            if (n < k_max) // ρ(n) = ε_n*ρ(0)
                // here, f=a*DOFe+b, ab is the subscripts of density matrix at time t
                // for ρ(0), only the element of init_state is 1, others are 0.
                // so matrix ε_n * vector ρ(0) is the cloumns of init_state in ε_n
                // i.e., ρ(n)[f] = ε_n[f][init_state], acutally it should be same as
                // input density matrix (learning period). And the ε_n is constructed
                // according to this relation.
                for (int f = 0; f < DOFe*DOFe; ++f)
                    RDM_output[i][n][f] = DM[n][f][init_state[i]];
            else         // ρ(n) = Σ_{k=1,k_max}T_k*ρ(n-k)
                for (int k = 1; k <= k_max; ++k) {
                    std::vector<Complex> tmp;
                    Matrix_Multiply(TT[k], RDM_output[i][n-k], tmp); // tmp = T_k*ρ_{n-k}
                    for (int f = 0; f < DOFe*DOFe; ++f)
                        RDM_output[i][n][f] += tmp[f]; // Σ{k=1,k_max}, accumalate each element
                }
}