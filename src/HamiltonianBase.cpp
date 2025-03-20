/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 3, 2021                                                 *
 * -------------------------------------------------------------------------- */

#include "HamiltonianBase.h"

void HamiltonianBase::init() {
    // * 0. Initialize data members and Hamiltonian.
    system_type          = param.getStr("system_type");
    model_type           = param.getStr("model_type");
    onthefly_type        = param.getStr("onthefly_type");
    allatom_type         = param.getStr("allatom_type");
    md_type              = param.getStr("md_type");
    PIMD_type            = param.getStr("PIMD_type");
    DOFe                 = param.getInt("DOFe");
    DOFn                 = param.getInt("DOFn");
    Condon_approximation = param.getBool("Condon_approximation");
    H.resize(DOFe, std::vector<double>(DOFe, 0));
    H_old.resize(DOFe, std::vector<double>(DOFe, 0));
    Heff.resize(DOFe, std::vector<double>(DOFe, 0));
    Heff_old.resize(DOFe, std::vector<double>(DOFe, 0));
    // Set Hamiltonian in diabatic or adiabatic representation
    representation = param.getStr("representation");
    if (representation.empty())  // If not specified, assume diabatic
        representation = "diabatic";
    else if (representation != "diabatic" && representation != "adiabatic" && representation != "quasi-diabatic")
        throw std::runtime_error("ERROR: Unsupported representation=" + representation +
            ". Only diabatic, adiabatic, and quasi-diabatic representations are supported now.");
    // Set unit conversion of energy to a.u., if real units are provided
    const std::string unit = param.getStr("energy_unit");
    unit2au = 1.0;
    if (unit == "cm-1")
        unit2au = cm2au;
    else if (unit == "eV")
        unit2au = eV2au;
    else if (unit == "kcal/mol")
        unit2au = kcal2au;
    else if (unit == "kj/mol")
        unit2au = kj2au;
    // if unit is au or empty, then no conversion
    // For SB/GOA model, if it is empty or "none", means to use the unitless model.
    else if (unit != "au" && unit != "" && unit != "none")
        throw std::runtime_error("ERROR: Unsupported energy_unit=" + unit);
    // Load Hamiltonian matrix from file, if an external file is provided, which
    // includes the off-diagonal (couplings) and diagonal (epsilons) elements.
    // The matrix should be a DOFe*DOFe real Hermitian matrix.
    // And the unit of value is same as energy_unit
    // In this case, the gamma_DA and epsilong in input will be ignored.
    const std::string H_load = param.getStr("H_load");
    if (!H_load.empty())
        loadHamiltonianMatrix(H_load);
    // Load value of gamma_DA (couplings) and epsilon with unit conversion from input.
    else {
        // Each value in input control file should be separated by ',' without space.
        // Note that if the number of input values is less than DOFe, the left values
        // will be set to 0, i.e., if the last few values are 0, they can be omitted.
        if (DOFe > 1) {
            // Read electronic couplings (gamma_DA) between each pair of states from
            // global paprameters. They are the off-diagonal elements of Hamiltonian.
            SplitString(diabatic_coupling, param.getStr("gamma_DA"));
            if (diabatic_coupling.size() > ((DOFe-1)*DOFe/2))
                throw std::runtime_error("ERROR: The number of values for electronic couplings is greater than (DOFe-1)*DOFe/2.");
            diabatic_coupling.resize((DOFe-1)*DOFe/2, 0);
            // Set the electronic couplings as off-diagonal element of Hamiltonian.
            // Note that the value of ij equals to the value of ji (real number).
            int index = 0;
            for (int i = 0; i < DOFe; ++i)
                for (int j = i+1; j < DOFe; ++j) {
                    diabatic_coupling[index] *= unit2au;
                    H[j][i] = H[i][j] = diabatic_coupling[index];
                    index++;
                }
        }
        // Get minimum energy of each state from input,
        SplitString(epsilon, param.getStr("epsilon"));
        if (epsilon.size() > DOFe)
            throw std::runtime_error("ERROR: The number of values for epsilon is greater than DOFe.");
        epsilon.resize(DOFe, 0);
        for (int i = 0; i < DOFe; ++i) {
            epsilon[i] *= unit2au;
            H[i][i] = epsilon[i];
        }
    }
    // Save Hamiltonian matrix from file (if required) in input energy unit
    // Do this in subclasses since the subclasses init() may modify it.
    //const std::string H_save= param.getStr("H_save");
    //if (!H_save.empty())
    //    saveHamiltonianMatrix(H_save);

    // * 1. Create an object of HamiltonianElec and initialize
    // Although, the electronic part is not always needed in all dynamics, for
    // example, the classical OpenMM MD doesn't require this. We create it here
    // for the consistense since it is reuqired in most of dynamics.
    Elec = std::make_shared<HamiltonianElec>(param);
    Elec->init();
}

double HamiltonianBase::getDiabaticCoupling(int i, int j) {
    if (i >= DOFe || i < 0 || j >= DOFe || j < 0 || i == j)
        throw std::runtime_error("ERROR: Called getDiabaticCoupling() with wrong state index.");
    // In the Condon case, the constant value is stored in diabatic_coupling
    // And only the value H_ij (i < j) is stored since H_ji = H_ij.
    if (Condon_approximation) {
        if (i < j)
            return diabatic_coupling[(2*DOFe-i-1)*i/2 + j-i - 1];
        else
            return diabatic_coupling[(2*DOFe-j-1)*j/2 + i-j - 1];
    }
    else // This should be defined in the subclasses if it supports non-Condon case
        getNonCondonCoupling(i, j);
}

void HamiltonianBase::loadHamiltonianMatrix(const std::string& file) {
    // Currently, only space separated file (.dat) is supported.
    std::ifstream inpfile(file);
    if (!inpfile)
        throw std::runtime_error("ERROR: Can't open Hamiltonin matrix file: " + file);
    int lineNumber = 0, count = 0;
    for (std::string line; getline(inpfile, line); ) {
        lineNumber++; // Let line number start from 1
        // A line starting with (#), or empty (including space, \t), will be ignored.
        if (line[0] == '#' || IsBlankString(line)) continue;
        // The elements of matrix should be divided with at least a space.
        std::vector<double> values;
        SplitLine(values, line);
        if (values.size() != DOFe)
            throw std::runtime_error("ERROR: Illegal row for Hamiltonian matrix at line" + std::to_string(lineNumber));
        for (int i = 0; i < DOFe; ++i)
            H[count][i] = values[i];
        count++;
        if (count == DOFe) break;
    }
    inpfile.close();
    if (count != DOFe)
        throw std::runtime_error("ERROR: Too few rows for Hamiltonian matrix in file:" + file);
    // Check if H is a real Hermitian matrix, H_ji = H_ij, i!=j;
    // Convert unit to au, and assigna values for epsilon and diabatic_coupling
    for (int i = 0; i < DOFe; ++i) {
        H[i][i] *= unit2au;
        epsilon.push_back(H[i][i]);
        for (int j = i+1; j < DOFe; ++j) {
            if (H[j][i] == H[i][j]) {
                H[i][j] *= unit2au;
                H[j][i] = H[i][j];
                diabatic_coupling.push_back(H[i][j]);
            }
            else
                throw std::runtime_error("ERROR: The input Hamiltonian matrix is not a real Hermitian matrix, i.e, H_jk != H_kj.");
        }
    }
}

void HamiltonianBase::saveHamiltonianMatrix(const std::string& file) {
    // Currently, only space separated file (.dat) is supported.
    FILE* outfile = CheckFile(file);
    // The unit of enengy in input, if not specified, assuming au.
    std::string energy_unit = param.getStr("energy_unit");
    if (energy_unit.empty())
        energy_unit = "au";
    fprintf(outfile, "# Hamiltonian matrix outputted by QCDyn in unit: %s\n", energy_unit.c_str());
    for (int i = 0; i < DOFe; ++i) {
        for (int j = 0; j < DOFe; ++j)
            fprintf(outfile, "%12.6g ", H[i][j]/unit2au);
        fprintf(outfile, "\n");
    }
    fclose(outfile);
}