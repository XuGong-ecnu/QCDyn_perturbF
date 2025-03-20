/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Nov. 23, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "HamiltonianModelBase.h"

/**
 * This class defines the Hamiltonian of the Multi-State (Global Bath) Harmonic (MSH)
 * model which is constructed from energy gap correlation functions obtained by
 * anharmonic all-atom MD simulation.
 *
 * And the minimum energy of each state (epsilon) and the couplings (gamma_DA)
 * between each pair of states should be provided in input control file (real unit).
 * And they will be converted to a.u. in the calculation.
 *
 * There are two models of MSH model, i.e., model I and model II. The
 * differnence is the size of S_j matrix: DOFe*(DOFe-1) (MSH-I); DOFe*(DOFe-2)
 * (MSH-II). However, we always build/load/save the full matrix of S_j (model I).
 * The total DOFn, calculations of energy and forces will depend on the model II
 * model I. When ground state is included, the model II can be used. The model II
 * and model can report identical results, and the advantage of model II is DOFn
 * decreasing and saving time.
 *
 * References:
 * [1] Zhubin Hu, Dominikus Brian, and Xiang Sun, J. Chem. Phys. 155, 124105 (2021)
 */
class HamiltonianMSH : public HamiltonianModelBase {
    friend class DynamicsBase;
    friend class DynamicsMQCBase;
    friend class DynamicsMF;
    friend class DynamicsLSC;
    friend class DynamicsSQC;
    friend class DynamicsSPM;
    friend class DynamicsCMM;
    friend class DynamicsECMM;
    friend class DynamicsFSSH;
    friend class DynamicsTBSH;
    friend class DynamicsMFRDM;
    friend class DynamicsECMMCV;

public:
    /**
     * Construct an HamiltonianMSH object.
     *
     * @param param   the global paramters
     */
    HamiltonianMSH(Parameters& param) : HamiltonianModelBase(param) {}
    ~HamiltonianMSH() {}
    /**
     * Initialize data members and model parameters.
     */
    void init();
    /**
     * Get potential energy.
     *
     * @param index  the index of state, start from 0
     * @return       the potential energy
     */
    double getPotentialEnergy(int index);
    /**
     * Compute and store diabatic forces into vector F. It is the negative
     * gradient of diabatic Hamiltonian matrix, i.e., F = - dH_ij/dR.
     *
     * If Condon_approximation is true, then the off-diagonal element of Hamiltonian
     * matrix (diabatic coupling) is a constant. So in this case (i !=j ), F is 0.
     *
     * @param i     the zero-based state index i
     * @param j     the zero-based state index j
     * @param F     the forces vector [out]
     */
    void getForces(int i, int j, std::vector<double>& F);

private:
    /**
     * Generate model paramters (frequencies: omega and equilbrium shifts: S)
     * based on input file containing energies of each state.
     * The input file contains DOFe columns corresponding to energies of each
     * state and no header, no index of row.
     * The unit of energy in  input file is specified by energy_unit in control
     * file.
     * The computed reorganization energy Er, reaction energy ΔE, average energy
     * gap <U>, and its variance σ^2 will be outputed to external file.
     */
    void buildModelParameters();
    /**
     * The function used in buildModelParameters().
     * To compute and return the value of f(ω_j) with given w and j.
     * f(ω_j) = 2Nω_j/(piCuu(0)) \int_0^∞ dt Cuu(t)/(ω_jt) sin(ω_jt) - j + 1/2
     * which is used to get frequency of discrete noromal modes from correlation
     * function (Cuu) of enenrgy gap.
     *
     * @param w    frequency in au
     * @param j    index, j = 1,...,N
     * @param dt   the time step in au for generating Cuu(t) from MD simulation
     * @param Cuu  time correlation function of energy gap (in au^2).
     */
    double fw(double w, int j, double dt, std::vector<double>& Cuu);
    /**
     * Similar as function fw() but compute and the value of J(ω) with given w.
     * J(ω) = βω/4 \int_0^∞ dt Cuu(t)/cos(ωt), β = 1/kT
     * which is used to get the continus spectral density from correlation
     * function (Cuu) of enenrgy gap.
     *
     * @param w    frequency in au
     * @param beta β = 1/kT, in au
     * @param dt   the time step in au for generating Cuu(t) from MD simulation
     * @param Cuu  time correlation function of energy gap (in au^2).
     */
    double Jw(double w, double beta, double dt, std::vector<double>& Cuu);
    /**
     * Load model paramters from external file.
     *
     * Requirments of file:
     * 1. First line is headers.
     * 2. Second - Last line  is the data, inlcude j, omega_j, S_j matrix in
     *    DOFe*(DOFe-1) base (but in the form of vector, i.e., vector[n*DOFe+m]
     *    = Matrix[n][m]).
     * 3. The number of cloumns is DOFe*(DOFe-1)+2
     * 4. The number of lines is N+1 (N is number of nromal modes).
     * 5. Each cloumn should be sepearted by one comma (,), i.e., csv file.
     * 6. The unit is au.
     */
    void loadModelParameters(const std::string& loadfile);
    /**
     * Save model paramters to external file.
     * The format of output is same as the descripition in loadModelParameters().
     */
    void saveModelParameters(const std::string& savefile);

private:
    // MSH model type: I (default) or II;
    std::string MSH_type;
};