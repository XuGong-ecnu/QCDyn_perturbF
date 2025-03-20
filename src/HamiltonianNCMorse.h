/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zengkui Liu @Sun Group @NYU-SH                                     *
 * First created: Dec. 15, 2021                                               *
 * Last updated:                                                              *
 * -------------------------------------------------------------------------- */

#pragma once
#include "HamiltonianModelBase.h"

/**
 * An HamiltonianNCMorse class defines the the non-Condon Morse potential model
 * Hamiltonian, which is a non-Condon model and used to study the gas-phase non-
 * adiabatic ultrafast photo-dissociation. 
 *
 * Here, the non-Condon Morse model is defined in a singular form for origin system.
 *
 * Reference:
 * [1] Chem. Phys. Lett. 349, 521 (2001).
 * [2] J. Phys. Chem. A. 125, 6845 (2021).
 */
class HamiltonianNCMorse : public HamiltonianModelBase {
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
     * Construct an HamiltonianNCMorse object.
     *
     * @param param   the global paramters
     */
    HamiltonianNCMorse(Parameters& param) : HamiltonianModelBase(param) {}
    ~HamiltonianNCMorse() {}
    /**
     * Initialize data members and model parameters.
     */
    void init();
    /**
     * Get total potential energy of specified state.
     *
     * @param index  the index of state
     * @return       the potential energy
     */
    double getPotentialEnergy(int index);
    /**
     * Compute and store diabatic forces into vector F. It is the negative
     * gradient of diabatic Hamiltonian matrix, i.e., F = - dH_ij/dR.
     *
     * Since LVC model is a non-Condon model, the off-diagonal element of
     * Hamiltonian matrix (diabatic coupling) is R-dependent. So in the case
     * of i !=j, F is not 0.
     *
     * @param i     the zero-based state index i
     * @param j     the zero-based state index j
     * @param F     the forces vector [out]
     */
    void getForces(int i, int j, std::vector<double>& F);

private:
    /**
     * Compute and return electronic non-Condon diabatic coupling i.e., the
     * off-diagonal element of Hamiltonian matrix, H_ij (i != j). Here, assuming
     * the Hamiltonian is a real symmetric matrix, and thus H_ij = H_ji (i != j).
     *
     * For LVC model, the non-Condon R-dependent diabatic coupling is defined as
     * H_ji = H_ij = Σγ_ij * R_k, where γ_ij is the linear coupling coefficient.
     *
     * @param i  the zero-based state index i
     * @param j  the zero-based state index j, i != j
     * @return   the non-Condon diabatic coupling
     */
    double getNonCondonCoupling(int i, int j);
    /**
     * Compute and return the reaction free energy ΔE for the LVC model.
     *
     * @return   the reaction free energy ΔE
     */
    double getReactionFreeEnergy();
    /**
     * Compute and return the reorganization energy Er for the LVC model.
     *
     * @return   the reorganization energy Er
     */
    double getReorganizationEnergy();
    /**
     * Generate the equilbrium shifts of nulcear which will be used in initial
     * nuclear sampling. The values of them will be stored in the base class
     * data member: shifts.
     *
     * For LVC model, it is decided from the system-bath coupling c_j of
     * specified sample_state.
     */
    void getSamplingShifts();
    /**
     * Generate model paramters based on spectral density, which is not
     * implemented yet.
     */
    void buildModelParameters();
    /**
     * Load model paramters from external file.
     * The format of file can be csv file or dat file with fixed width.
     *
     * Requirments of file (F-state LVC model):
     * 1. First line is headers: index,omgea,c0,c1,...,,c(F-1),gamma01,gamma02,
     *    gamma0(F-1),gamma12,gamma13,...,gamma1(F-1),...,gamma(F-2)(F-1). The
     *    number of c columns is F, and number of gamma columns is (F-1)*F/2.
     * 2. The number of lines is N+1 (N is number of nromal modes).
     * 3. The unit is au.
     */
    void loadModelParameters(const std::string& loadfile);
    /**
     * Save model paramters to external file.
     * The format of output is same as the descripition in loadModelParameters().
     */
    void saveModelParameters(const std::string& savefile);
    /*
    */
    std::string Morse_preset;

private:
    /**
     * Here we have the Morse potential given by
     * 
     * V_i (x) = D_i * [1 - exp(- b_i *(x - Req_i))]^2 + C_i, i = 1,2,3
     * 
     * where D_i serves as the depth of the Morse potential energy,
     * b_i  controls the width of the Morse potential well,
     * Req_i gives the equilbrium position of the PES,
     * and C_i serves as the vertical shift of the Morse potential.
     * 
     * The coupling term in this model is given in the Gaussian term given by
     * 
     * V_ij(x) = A_ij * exp{-a_ij * (x - Req_ij)^2 }, i,j = 1,2,3 and (i != j)
     * 
     * where the amplitudes of the coupling is given by A_ij, a_ij controls 
     * width of the peak, Req_ij serves the central position of the coupling 
     * distribution. 
     * 
     */
    std::vector<double> D, b, Req, C, A, Rcp, a;
    // Nuclaer mass 
    const double M = 20000.0;
};
