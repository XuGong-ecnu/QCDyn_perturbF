/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Nov. 25, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "HamiltonianModelBase.h"

/**
 * An HamiltonianSpinBoson class defines the Hamiltonian of the well-known
 * Spin-Boson model (two-level system, DOFe=2).
 *
 * The SB, SB_req and GOA type Spin-Boson models with spec_density=Ohmic or Debye
 * (or load parameters from file) are supported now.
 *
 * Note that the SB_req is a equivalence of SB, but writting the Spin-Boson
 * model in terms of the relative shifts in the equilibrium positions (req), which
 * means the donor min is at 0, and the displacement between the donor and acceptor
 * equilibrium geometries along the j-th mode is req_j=2c_j/omega_j^2.
 *
 * In most papers, the Spin-Boson model means SB model, the spin-boson model
 * is symmetric when epsilon = 0 and asymmetric when epsilon ≠ 0.
 *
 * By defalut, the unit of SB model is not included explicitly, and the units of
 * parameters are based on the coupling (gamma_DA or Δ) or cutoff frequency
 * (hbar*omega_c). Then the relation between other parameters and Δ is as following:
 * ε = *Δ; β = /Δ; ω_c = *Δ; DT = /Δ; and η is dimensionless.
 *
 * But, you can set enegey_unit in input to use real unit, and in this case, the
 * temperatue should be in K (beta will be ignored), and DT should be in ps.
 *
 * The non-Condon case of SB/SB_req model are supported, which is a equivalence
 * of two-state LVC model. And the non-Condon GOA model is also supported, the
 * definition of it can be found: J. Chem. Phys. 144, 244105 (2016).
 */
class HamiltonianSpinBoson : public HamiltonianModelBase {
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
     * Construct an HamiltonianSpinBoson object.
     *
     * @param param   the global paramters
     */
    HamiltonianSpinBoson(Parameters& param) : HamiltonianModelBase(param) {}
    ~HamiltonianSpinBoson() {}
    /**
     * Initialize data members and model parameters.
     */
    void init();
    /**
     * Update Hamiltonian in adiabatic representation (potential energies, forces
     * of all states, and nonadiabatic coupling (NAC) vectors).
     *
     * For SB model, it is defined analyticallly and more efficent, which will
     * override the general verison in base class. While, for SB_req and GOA
     * model, I still use the general functions in base class.
     *
     * It will update data: H, H_old, T, T_dagger, F_all, NAC
     */
    void updateAdiabaticHamiltonian() override;
    /**
     * Get total potential energy of specified state.
     *
     * @param index  the index of state, 0 or 1
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
    /**
     * Get potential energy of specified normal mode of specified state.
     *
     * Note: the epsilon in Hamiltonian operator is not included here.
     *
     * @param state  the index of state, 0 or 1
     * @param mode   the index of normal mode, 0 ~ N-1
     * @return       the potential energy of this nomal mode of this state
     */
    double getNormalModePE(int state, int mode);
    /**
     * Get kinetic energy of specified normal mode.
     *
     * @param mode   the index of normal mode, 0 ~ N-1
     * @return       the kinetic energy of this nomal mode
     */
    double getNormalModeKE(int mode);

private:
    /**
     * Compute and return electronic non-Condon diabatic coupling i.e., the
     * off-diagonal element of Hamiltonian matrix, H_ij (i != j). Here, assuming
     * the Hamiltonian is a real symmetric matrix, and thus H_ij = H_ji (i != j).
     *
     * For SB/SB_req model, the non-Condon R-dependent diabatic coupling is defined
     * as H_ji = H_ij = Σγ * R_k, where γ is the linear coupling coefficient, which
     * is the same as LVC model.
     *
     * @param i  the zero-based state index i
     * @param j  the zero-based state index j, i != j
     * @return   the non-Condon diabatic coupling
     */
    double getNonCondonCoupling(int i, int j);
    /**
     * Generate the equilbrium shifts of nulcear which will be used in initial
     * nuclear sampling. The values of them will be stored in the base class
     * data member: shifts.
     */
    void getSamplingShifts();
    /**
     * Generate model paramters: omega, c, req, shift based on spectral density.
     * Ohmic and Debye spectral density are implemented now.
     */
    void buildModelParameters();
    /**
     * Load model paramters from external file.
     * The format of file can be csv file or dat file with fixed width.
     *
     * Requirments of file:
     * 1. First line is headers: index,omgea,c,req
     * 2. The number of lines is N+1 (N is number of nromal modes).
     * 3. If GOA model is used, the column of shift should be provided.
     * 4. If Condon_approximation is false (non-Condon), gamma should be provided.
     * 5. The unit is au.
     */
    void loadModelParameters(const std::string& loadfile);
    /**
     * Save model paramters to external file.
     * The format of output is same as the descripition in loadModelParameters().
     */
    void saveModelParameters(const std::string& savefile);
    /**
     * The analytical form for SB model to update adiabatic Hamiltionian and
     * forces and nonadiabatic coupling vectors.
     *
     * Note: this function is used for model_type = SB olu, not work for SB_req
     * or GOA model.
     */
    void updateAdiabaticHamiltonianForSB();

private:
    // gamma represents a N-dimentional linear coupling coefficients of diabatic
    // coupling, which is used in the case of non-Condon. H_ij = Σγ * R_k, i !=j
    // When non-Condon is used for SB/SB_req model, the model parameters including
    // gamma should be loaded from file. And in the non-Condon case, the SB/SB-req
    // model is an equivalent two-state. In the case of non-Condon GOA model,
    // the coupling coefficients can be dertermined from the input value of key
    // "GOA_gamma".
    std::vector<double> gamma;
};