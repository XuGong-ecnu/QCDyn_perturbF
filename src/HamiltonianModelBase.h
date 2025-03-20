/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 3, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "HamiltonianBase.h"
#include "HamiltonianElec.h"

/**
 * This is a base (abstract) class defines the Hamiltonian of model system.
 * Subclasses of it define particular definition of models, such as Spin-Boson
 * model, Multi-State Harmonic model, FMO model, LVC model, and so on.
 *
 * The common part of definitions for different models are defined here, such as
 * frequencies of normal modes, getKineticEnergy(), and so on.
 *
 * And the data (positions, velocities, and forces) that being propagated
 * by Dynamics object are stored in this class data members (in vector<double>
 * format, which are differnent with that in HamiltonianAllAtom).
 *
 * Note:
 * (1) For model system, all the units are a.u. in the calculation, and masses are 1.
 * (2) N is number of normal modes, and DOFn is totoal nulear DOF, they are different
 *     in MSH/FMO models.
 * (3) The model Hamiltonian is defined in the diabatic representation, but can
 *     be transformed to the adiabatic representation via rotation matrix, and
 *     the nonadiabatic coupling can also be computed.
 */
class HamiltonianModelBase : public HamiltonianBase {
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
    virtual ~HamiltonianModelBase() {}
    /**
     * Initialize data members of the HamiltonianModelBase object.
     */
    virtual void init();
    /**
     * Update Hamiltonian in diabatic representation (potential energies, forces
     * of all states, and averaged forces).
     *
     * It will update data: H, H_old, Heff, Heff_old, F_all, F_avg
     *
     * Note, Hamiltonian (H) denotes the standard one. If the averege potential
     * energy (H_avg) is removed from the diagonal elements of Hamiltonian (H),
     * then it is called effective Hamiltonian (Heff = H - H_avg), which is
     * commonly used in diabatic propagation of some dynamics methods.
     */
    virtual void updateDiabaticHamiltonian();
    /**
     * Update Hamiltonian in adiabatic representation (potential energies, forces
     * of all states, and nonadiabatic coupling (NAC) vectors).
     *
     * This will call internal functions updateHAndT() and updateFallAndNAC()
     * to finish this job via roation matrix, which are general for any model.
     * But, you can define analytical verisons (if any) to override it/them in
     * subclasses (such as Spin-Boson model).
     *
     * It will update data: H, H_old, T, T_dagger, F_all, NAC
     *
     * It is used for mapping dynamics with propagation in adiabatic
     * representation, which is proposed by Miller.
     *
     * References:
     * [1] J. Chem. Phys. 147, 064112 (2017)
     */
    virtual void updateAdiabaticHamiltonian();
    /**
     * Update Hamiltonian in quasi-diabatic representation.
     *
     * It will update data: H, H_old, F_all
     *
     * It is used for mapping dynamics with propagation in quasi-diabatic
     * representation, which is proposed by Pengfei Huo. And here, H_old the
     * current Hamiltonian in adiabatic basis, and H is the current Hamiltonian
     * in quasi-diabatic basis.
     *
     * References:
     * [1] J. Chem. Phys. 149, 044115 (2018)
     * [2] J. Chem. Theory Comput. 2018, 14, 1828−1840
     * [3] J. Phys. Chem. Lett. 2019, 10, 7062−7070
     * [4] J. Chem. Phys. 155, 084106 (2021)
     */
    virtual void updateQuasiDiabaticHamiltonian();
    /**
     * Compute and return diabatic potential energy.
     *
     * @param index  the index of state
     * @return       the potential energy
     */
    virtual double getPotentialEnergy(int index) = 0;
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
    virtual void getForces(int i, int j, std::vector<double>& F) = 0;
    /**
     * Get kinetic energy. All states share same velocities.
     *
     * @return       the kinetic energy
     */
    double getKineticEnergy();
    /**
     * Get the number of normal modes (N) of model.
     *
     * @return  the number of normal modes (N)
     */
    int getNumNormalModes() const {
        return N;
    }

protected:
    /**
     * Construct an HamiltonianModelBase object.
     *
     * @param param   the global paramters
     */
    HamiltonianModelBase(Parameters& param) : HamiltonianBase(param) {}
    /**
     * Initialize model paramters by loading from file or generating internally.
     * Resize the vector of parameters firstly before calling this.
     */
    void initializeModelParameters();
    /**
     * Generate model paramters internally.
     */
    virtual void buildModelParameters() = 0;
    /**
     * Load model paramters from external file.
     */
    virtual void loadModelParameters(const std::string& loadfile) = 0;
    /**
     * Save model paramters to external file.
     */
    virtual void saveModelParameters(const std::string& savefile) = 0;
    /**
     * Get frequencies (omega) and coupling coefficients (c) by discretize specific
     * spectral density function for the bath. Ohmic and Debye bath are supported.
     *
     * Various discretization schemes are implemented according to publications
     * so that can reproduce the result in publications.
     *
     * Note that the unit of omega is au.
     *
     * @param spec_density Ohmic or Debye bath
     * @param N            number of normal modes
     * @param omega_c      cutoff frequency or characteristic frequency for bath
     * @param eta          Kondo parameter for Ohmic bath
     * @param lambda       reorganization energy for Debye bath
     * @param omega        frequency of harmonic normal modes [out]
     * @param c            coupling coefficients [out]
     */
    void discretizeSpectralDensity(const std::string& spec_density, int N, double omega_c, double eta,
                                   double lambda, std::vector<double>& omega, std::vector<double>& c);
    /**
     * Internal function 1 used for updateAdiabaticHamiltonian().
     *
     * Update adiabatic Hamiltonian via the diagonalization of diabatic Hamiltonian.
     * This will also update the unitary rotation matrix, too, (T and T_dagger).
     *
     * Note, the indices of adiabatic states are in energy ascending order,
     * i.e., state 0 is ground state.
     *
     * It will update data: H, T, T_dagger.
     */
    virtual void updateHAndT();
    /**
     * Internal function 2 used for updateAdiabaticHamiltonian().
     *
     * Compute adiabatic forces and nonadiabatic coupling vectors from diabatic
     * forces via the rotation matrix.
     *
     * Before calling this, you should call updateHAndT() to update
     * the adiabatic Hamiltonian and unitary rotation matrix firstly.
     *
     * It will update data: Fall, NAC.
     */
    virtual void updateFallAndNAC();

protected:
    // the number of harmonic normal modes (oscillators)
    // for spin-boson model, DOFn = N, but not for MSH/FMO models
    int N;
    // model parameters: frequency of harmonic normal modes, size = N
    std::vector<double> omega;
    // model parameters: coupling coefficient between the electronic DOF and
    // nuclear normal modes
    std::vector<double> c;
    // model parameters: equilibrium shifts of each state along normal modes,
    // which are used in energy and forces calculations
    std::vector<double> req;
    // model parameters: shifts of initial state along normal modes (size=DOFn)
    // which are used in initial nuclear sampling in some models.
    std::vector<double> shift;
    // nuclear positions
    std::vector<double> R;
    // nuclear velocities (equal to momenta) since the all masess are 1.
    std::vector<double> V;
    // effective nuclear forces used for propagation
    std::vector<double> F;
    // gradient matrix in DOFe*DOFe base, i.e., -dH/dR
    // If Conda approximation/adiabatic is used, the off-diagonal is zero.
    // Here, use vector to store them: vector[n*DOFe+m] = Matrix[n][m]
    // Thus, the diabatic or adiabatic forces of state i is F_all[i*DOFe+i].
    Real_Matrix F_all;
    // average nuclear forces of each state, i.e. -dH_nn/dR (diagonal of H)
    // i.e., the average of diagonal elemnts of F_all. (diabatic only)
    std::vector<double> F_avg;
    // The nonadiabatic coupling (NAC) vectors between each pair of states, the order:
    // 01, 02, 03, ... 0i; 12, 13, ... 1i; ...; (i-1)i; in which i index of adiabatic
    // state. And, d_ij = -d_ji. The pairs of vectors is (DOFe-1)*DOFe/2 (DOFe > 1).
    // And the size of each nonadiabatic coupling vectors is DOFn. The index of
    // element (vector) in this matrix can be got by: index =
    // (2*DOFe-i-1)*i/2 + (j-i) -1, here, 0 =< i < j < DOFe, the ij element
    // is the index in full DOFe*DOFe matrix. e.g., DOFe=4, then the index of
    // d_12 in nonadiabatic_coupling is 3 (the 4th element).
    Real_Matrix NAC;
    // T is the unitary rotation matrix (transformation matrix), which can be used
    // to transofrm between diabatic and adiabatic density matrix, and can be used
    // to compute adiabatic forces from diabatic forces, also can be used to compute
    // the nonadiabatic coupling vectors between adiabatic states.
    // The adiabatic Hamiltonian (off-diagonal is 0) can be got by the diagonalization
    // of diabatic Hamiltonian. The diagonal element (adiabatic energies) of
    // adiabtic Hamiltonian are the eigenvalues, and the columns of rotation matrix
    // is are the coeffcients of eigenvectors. i.e. H_adiabtic = T_dagger*H_dibatic*T
    // T_dagger is conjugate transpose (Hermitian transpose) matrix of T, since T
    // is real matrix (when diabatic coupling is real number), then T_dagger =
    // transpose matrix of T. Note that the eigenvalues/vectors is in ascending order.
    Real_Matrix T, T_dagger;
};