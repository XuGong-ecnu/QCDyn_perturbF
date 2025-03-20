/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * First created: Nov. 23, 2021                                               *
 * Last updated: Nov. 24, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "HamiltonianModelBase.h"

/**
 * An HamiltonianLVC class defines the the linear vibronic coupling (LVC) model
 * Hamiltonian, which is a non-Condon model and used to study the electronic
 * transitions through a conical intersection in gas phase molecules.
 *
 * Here, the LVC model is defined in a general form for multi-state system.
 *
 * Reference:
 * [1] J. Chem. Phys. 135, 234106 (2011).
 * [2] J. Chem. Phys. 144, 244105 (2016).
 * [3] J. Chem. Theory Comput. 16, 4479−4488 (2020).
 */
class HamiltonianLVC : public HamiltonianModelBase {
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
     * Construct an HamiltonianLVC object.
     *
     * @param param   the global paramters
     */
    HamiltonianLVC(Parameters& param) : HamiltonianModelBase(param) {}
    ~HamiltonianLVC() {}
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
     * Compute and return the donor-to-acceptor reaction free energy ΔE for the LVC model.
     *
     * @return   the donor-to-acceptor reaction free energy ΔE
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

private:
    // The c, gamma are special data member for LVC model, other model paramters
    // such as omega, shift, epsilon, and so on will be inherited from ModelBase.
    // It should be noted that the meaning of epsilon won't be the minima energies
    // since there will be extra terms for completing squares (in literature,
    // usually the "Δ" is used for it). And the c reprensents a N-dimentional
    // system-bath coupling in LVC model (in literature, the "d_j" called linear
    // shift is used). In other words, the "epsilon" and "c" used here are the
    // same as the "Δ" and "d" used in literature in LVC model. We use these
    // notations are consistent with the ones used in Spin-Boson model.
    // Each element of c represents a N-dimentional c_j in LVC model.
    // And the size of c is number of states, i.e., DOFe.
    // Each element of gamma represents a N-dimentional linear coupling coefficient
    // γ_ij in LVC model. And the size of gamma equals to the half number of
    // diabatic couplings, i.e., (DOFe-1)*DOFe/2. Only γ_ij (i < j) is stored
    // since since γ_ji = γ_ij. And the relation between the index of gamma and
    // the indices i, j is: index = (2*DOFe-i-1)*i/2 + j-i - 1, i < j. The order
    // γ_ij in gamma is 01, 02, 03, ... 0(F-1); 12, 13, ... 1(F-1); ...; (F-2)(F-1),
    // where F is number of states (DOFe).
    std::vector<std::vector<double>> c, gamma;
};