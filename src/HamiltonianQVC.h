/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * First created: Nov. 28, 2021                                               *
 * Last updated: Nov. 30, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "HamiltonianModelBase.h"

/**
 * An HamiltonianQVC class defines the the quadratic vibronic coupling (QVC) model
 * Hamiltonian, which is an extension of the LVC Hamiltonian, and more general
 * and allows the inclusion of effects such as the Duschinsky rotation of the
 * normal coordinates and changes in vibrational frequencies between two states.
 *
 * Here, the QVC model is defined for two-state system only.
 *
 * Reference:
 * [1] J. Chem. Phys. 141, 034104 (2014).
 */
class HamiltonianQVC : public HamiltonianModelBase {
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
    friend class DynamicsMixPES;

public:
    /**
     * Construct an HamiltonianQVC object.
     *
     * @param param   the global paramters
     */
    HamiltonianQVC(Parameters& param) : HamiltonianModelBase(param) {}
    ~HamiltonianQVC() {}
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
     * Since QVC model is a non-Condon model, the off-diagonal element of
     * Hamiltonian matrix (diabatic coupling) is R-dependent. So in the case
     * of i !=j, F is not 0.
     *
     * @param i     the zero-based state index i
     * @param j     the zero-based state index j
     * @param F     the forces vector [out]
     */
    void getForces(int i, int j, std::vector<double>& F);
    /**
     * Transform the normal modes or positions/momenta between the donor and
     * acceptor states using the donor-to-acceptor Duschinsky rotation matrix
     * and shift vector (if any).
     *
     * Transfrom donor to acceptor by using X_A,i = ΣJ_ij (X_D,j + shift_j)
     * Transfrom acceptor to donor by using X_D,i = (ΣJ_ji X_A,j) - shift_i
     *
     * @param X_old     the normal modes/positions to be transformed (read)
     * @param J         Duschinsky rotation matrix from donor to accrptor
     * @param shift     the shift vector between donor and acceptor state
     * @param X_new     the transformed normal modes/positions (write)
     * @param inverse   If true, do transfromation from acceptor to donor
     *                  Default is false and do transfromation from donor to acceptor
     */
    void transformNormalModes(const std::vector<double>& X_old, const Real_Matrix& J,
        const std::vector<double>& shift, std::vector<double>& X_new, bool inverse = false);

private:
    /**
     * Compute and return electronic non-Condon diabatic coupling i.e., the
     * off-diagonal element of Hamiltonian matrix, H_ij (i != j). Here, assuming
     * the Hamiltonian is a real symmetric matrix, and thus H_ij = H_ji (i != j).
     *
     * @param i  the zero-based state index i
     * @param j  the zero-based state index j, i != j
     * @return   the non-Condon diabatic coupling
     */
    double getNonCondonCoupling(int i, int j);
    /**
     * Generate model paramters based on spectral density, which is not
     * implemented yet.
     */
    void buildModelParameters();
    /**
     * Load model paramters from external file.
     * The format of file can be csv file or dat file with fixed width.
     *
     * Requirments of file (2-state QVC model with N normal modes):
     * 1. First line is headers: index,omgea_D,omgea_A,req_A,gamma,theta_0j,
     *    theta_1j,...,theta_(N-1)j,J_0j,J_1j,...,J_(N-1)j.
     * 2. The number of lines is N+1 (N is number of nromal modes).
     * 3. The unit is au.
     */
    void loadModelParameters(const std::string& loadfile);
    /**
     * Save model paramters to external file.
     * The format of output is same as the descripition in loadModelParameters().
     */
    void saveModelParameters(const std::string& savefile);
    /**
     * Check if the donor and acceptor states have same kinetic energy.
     * In principle, they should be identical.
     *
     * This function is used for debug only.
     *
     * @return true if the kinetic energies are same.
     */
    bool checkKineticEnergy();

private:
    // The gamma, J, and theta are special data member for QVC model, other model
    // paramters such as omega, shift, epsilon, req, and so on will be inherited
    // from ModelBase. It should be noted that epsilon[0] is zero and epsilon[1]
    // is the ΔE in literature. req is the shift vector of acceptor state (x_A
    // in the literature). And the constant (ΔDA) terms in the coupling in the
    // literature is the gamma_DA from input control file. Note that the omega
    // represents the N frequencies of donor state (omega_D) and the omgea_A
    // represents the N frequencies of acceptor state. They are different in QVC
    // model.
    std::vector<double> omega_A;
    // gamma represents a N-dimentional linear coupling coefficient
    // γ in QVC model, which is the same parameter in LVC model.
    std::vector<double> gamma;
    // N*N dimentional Duschinsky matrix of acceptor state.
    Real_Matrix J;
    // N*N dimentional quadratic coupling coefficient.
    Real_Matrix Theta;

};