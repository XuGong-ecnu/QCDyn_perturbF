/* --------------------------------------OpenMM------------------------------------ *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 3, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "Tools.h"
#include "OpenMM.h"
#include "Parameters.h"
#include "HamiltonianElec.h"
#include "Vec3.h"

/**
 * An Hamiltonian class defines the Hamiltonian of system including potential
 * energy in DOFe*DOFe electronic bases and diabatic electronic coupling. The
 * electronic part of Hamiltonian (such as electronic mapping variables) are
 * in the data member HamiltonianElec.
 *
 * The calculations of energy, force, and/or nonadiabatic coupling of system are
 * also defined here. The system can be a model or a realistic system.
 *
 * Note: the Hamiltionian can be on the diabatic or adiabatic representation,
 * which is depended on the system Hamiltonian or dynamics method.
 *
 * The data of system, such as nuclear positions, velocities, forces, and/or
 * electronic mapping variable positions, momenta, coefficients are stored in
 * Hamiltonian object, and will be propagated by particular Dynamics object.
 *
 * Each Hamiltonian object should be bound to a particular Dynamics object to do
 * dynamics. This connection is specified by passing the HamiltonianBase as
 * an argument to the constructor of the Dynamics.
 *
 * HamiltonianBase is an abstract class.
 * Subclasses define particular definition of Hamiltonian.
 */
class HamiltonianBase {
    friend class DynamicsBase;
    friend class DynamicsOpenMM;
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
    friend class DynamicsMD;
    friend class DynamicsReadTraj;
    friend class DynamicsRPMD;

public:
    virtual ~HamiltonianBase() {}
    /**
     * Initialize data members of HamiltonianBase class.
     */
    virtual void init();

    /**
     * Get smarter pointer to HamiltonianElec object.
     */
    std::shared_ptr<HamiltonianElec> getElec() {
        return Elec;
    }
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
    virtual void updateDiabaticHamiltonian() = 0;
    /**
     * Update Hamiltonian in adiabatic representation (potential energies, forces
     * of all states, and nonadiabatic coupling (NAC) vectors).
     *
     * It will update data: H, H_old, F_all, NAC
     */
    virtual void updateAdiabaticHamiltonian() = 0;
    /**
     * Update Hamiltonian in quasi-diabatic representation.
     *
     * It will update data: H, H_old, F_all
     */
    virtual void updateQuasiDiabaticHamiltonian() = 0;
    /**
     * Compute and return potential energy.
     *
     * @param index  the index of state
     * @return       the potential energy
     */
    virtual double getPotentialEnergy(int index) = 0;
    /**
     * Get kinetic energy. All states share same velocities.
     *
     * @return       the kinetic energy
     */
    virtual double getKineticEnergy() = 0;
    /**
     * Compute and return electronic diabatic coupling (in au) i.e., the off-diagonal
     * element of Hamiltonian matrix, H_ij (i != j). Here, assuming the
     * Hamiltonian is a real symmetric matrix, and thus H_ij = H_ji (i != j).
     *
     * If Condon_approximation is true, then the diabatic coupling is a constant.
     * Otherwise, it is R-dependent and thus should be computed based on current
     * positions.
     *
     * @param i  the zero-based state index i
     * @param j  the zero-based state index j, i != j
     * @return   the diabatic coupling in au
     */
    virtual double getDiabaticCoupling(int i, int j);
    /**
     * Get the number of electronic states (DOFe).
     *
     * @return    the number of electronic states (DOFe).
     */
    int getNumStates() const {
        return DOFe;
    }
    /**
     * Get the nuclear DOF (DOFn).
     *
     * @return    the nuclear DOF (DOFn).
     */
    int getNulearDOF() const {
        return DOFn;
    }
    /**
     * Load Hamiltonian matrix (DOFe*DOFe) from file.
     * Each value of column is separated by at least one space.
     * The unit of input value is specified by "energy_unit".
     * And the unit will be converted to au after loaded.
     *
     * @param  file  file name of input matrix file
     */
    void loadHamiltonianMatrix(const std::string& file);
    /**
     * Save current Hamiltonian matrix (DOFe*DOFe) to file.
     * The unit of output is au.
     *
     * @param  file  file name of output matrix file
     */
    void saveHamiltonianMatrix(const std::string& file);

protected:
    /**
     * Compute and return electronic non-Condon diabatic coupling (in au) i.e., the
     * off-diagonal element of Hamiltonian matrix, H_ij (i != j). Here, assuming
     * the Hamiltonian is a real symmetric matrix, and thus H_ij = H_ji (i != j).
     * The non-Condon diabatic coupling is R-dependent.
     *
     * This is a part of getDiabaticCoupling(int i, int j). Within HamiltonianBase
     * class,  this function will throw an Exception directly (assuming current
     * Hamiltonian is a Condon case). So if the system Hamiltonian defined in
     * subclass supports non-Condon case, you should define this function explictly
     * in the subclass to override this one.
     *
     * @param i  the zero-based state index i
     * @param j  the zero-based state index j, i != j
     * @return   the non-Condon diabatic coupling in au
     */
    virtual double getNonCondonCoupling(int i, int j) {
        // Please define it explictly in subclass, if it supports non-Condon case.
        throw std::runtime_error("ERROR: This Hamiltonian doesn't support non-Condon simulation.");
    }

protected:
    /**
     * This function is only used in the constructor of subclass, since you can
     * not construct an object for an abstract class.
     *
     * @param param   the global paramters
     */
    HamiltonianBase(Parameters& param) : param(param) {}
    // Parameters object controls the simulation
    Parameters& param;
    // The smarter pointer to HamiltonianElec object, which stores the electronic
    // part of Hamiltonian.
    std::shared_ptr<HamiltonianElec> Elec;
    // The system type: model, onthefly, allatom.
    // model: anlytical model Hamiltonian, SB, MSH, FMO and so on
    // onthefly: electronic structure interface,
    // allatom: calssical force filed MD interface, OpenMM
    std::string system_type, model_type, onthefly_type, allatom_type, md_type, PIMD_type;
    // Hamiltonian in diabatic or adiabatic or quasi-diabatic representation.
    // The model Hamiltonian supports these three representations.
    // The allatom (OpenMM) Hamiltonian supports diabatic representation only.
    std::string representation;
    // Whether Hamiltonian use Condon approximation. If it is true, the diabatic
    // couplings (off-diagonal element of Hamiltonian matrix) are real constants
    // which are specified by user with the key "gamma_DA" in the input control
    // file and stored in the data memeber "diabtaic_coupling". Otherwise, the
    // the diabatic coupling is R-dependent, and should be computed base on the
    // current positions at each step by getNonCondonCoupling().
    bool Condon_approximation;
    // Here, H is Hamiltonian of system in DOFe*DOFe electronic basis.
    // The diagonal elements are potential energies of each state (added energy
    // corrections (for all atom)). And the off-diagonal elements are electronic
    // couplings. When the Condon approximation is used, the diabatic couplings
    // are real constants. When removing average potential energy of all
    // states from diagonal elements of H, (Heff = H - H_avg), it is called
    // effective Hamiltion (Heff) matrix, which is used in some dymacis method.
    // H is always the Hamiltonian for current phase space, and H_old is the
    // previous Hamiltonian. The units of all elemets are a.u. (if any).
    // Note, H, H_old can be diabatic or adiabatic Hamiltonian depends on the
    // representation. Note, for adiabatic Hamiltonian, the index of
    // state is in energy ascending order, i.e., state 0 is ground state.
    Real_Matrix H, H_old, Heff, Heff_old;
    // DOFe: electronic DOF = number of states/topologies/surfaces.
    // DOFn: nuclear DOF, for model system, it is decided by user (read from input);
    // for allatom system, it is number of atoms (decided by structure file).
    int DOFe, DOFn;
    // The diabatic electronic couplings (gamma_DA) between each pair of states
    // from input, which should follow the order: 01, 02, 03, ... 0(F-1); 12, 13,
    // ... 1(F-1); ...; (F-2)(F-1), where F is DOFe. Here, gamma_DA = gamma_AD is
    // real number. Therefore, the size of vector is (DOFe-1)*DOFe/2 (DOFe > 1).
    // They are the off-diagonal elements of Hamiltonian, H_jk = H_jk, j != k.
    // The relation between the index of this vector and the matrix j, k is:
    // index = (2*DOFe-j-1)*j/2 + k-j - 1, j < k.
    // When using Condon approximation, they are constants. The unit is au.
    std::vector<double> diabatic_coupling;
    // For model, epsilon is minimum energy of each state (unit is au),
    // while the meaning of it may be different in different model.
    // For allatom, it is the energy correction of each state,
    // which is used for the potential energy correction.
    // If no vaules provided by user, then they are 0. (size = DOFe)
    std::vector<double> epsilon;
    // unit2au is the factor used to convert input unit of energy to a.u.
    // If energy_unit is provided in input, then convert it to a.u.
    // By defalut it is 1, which means no conversion.
    // All units used i Hamiltonian matrix during simulation is a.u. (if any)
    double unit2au;
};