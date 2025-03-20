/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "DynamicsElec.h"
#include "DynamicsBase.h"

/**
 * This class is a base (abstract) class for mixed quantum-classical (MQC)
 * dynamcis in diabatic or adiabatic representation. The subclasses of it define
 * dynamics methods, such as Ehrefest dynamics, surface hopping, mixed
 * quantum-classical Liouville (MQCL), and mapping dynamics (linearized
 * semiclassical (LSC), symmetrical quasi-classical (SQC), spin-mapping (SPM))....
 *
 * It defines the common behaviours (virtual functions) to complete a nonadiabatic
 * dynamimcs simulation. The explicit implementation of different methods should
 * be defined in the subclass.
 *
 * And some general base virtual functions are defined here, which can be directly
 * used for some methods. If them cannot applied to the subclass method, please
 * define a function with identical name in subclass to override the base one.
 */
class DynamicsMQCBase : public DynamicsBase {
public:
    virtual ~DynamicsMQCBase() {}
    /**
     * Initialize data members and check the legality of parameters.
     */
    virtual void init();
    /**
     * Make preparations for starting to propagate a trajectory.
     * You must call it before running a new trajectory.
     */
    virtual void beforeOneTraj();
    /**
     * Advance a simulation through time by taking a series of time steps.
     *
     * @param steps the number of time steps to take silently
     */
    virtual void dynamics(int steps);
    /**
     * Make conclusion for a trajectory.
     *
     * Pleases define it explicitly in subclass to override this one, if it is
     * required by the dynamcis method.
     */
    virtual void afterOneTraj() {}
    /**
     * Get averaged electronic reduced density matrix (RDM) or observables, which
     * is accumaleated and averaged by number of trajectories.
     *
     * This is a multi-frame reduced density matrix and each frame of it is
     * averaged by number of trajectories.
     *
     * It may inlcude sevral RDMs denpends on the methods, if they can computed
     * in one simulation.
     *
     * @return the constant reference to the RDM
     */
    virtual const std::vector<std::vector<Complex>>& getRDM() {
        return RDM;
    }
    virtual const std::vector<std::vector<Complex>>& getRDMLSC1() {
        return RDM_LSC1;
    }
    virtual const std::vector<std::vector<Complex>>& getRDMLSC2() {
        return RDM_LSC2;
    }
    virtual const std::vector<std::vector<Complex>>& getRDMRILSC1() {
        return RDM_RILSC1;
    }
    virtual const std::vector<std::vector<Complex>>& getRDMRILSC2() {
        return RDM_RILSC2;
    }
    virtual const std::vector<std::vector<Complex>>& getRDMRILSC3() {
        return RDM_RILSC3;
    }
    virtual const std::vector<std::vector<Complex>>& getAdiabaticRDM() {
        return AdiabaticRDM;
    }

protected:
    /**
     * Construct a DynamicsMQCBase object.
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     */
    DynamicsMQCBase(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsBase(param, Ha) {}

    // # The following three pure virtual functions must be defined explicitly
    // # in subclasses, since they are different in MQC dynamics methods.
    /**
     * Get initial sampling of electronic mapping variables.
     */
    virtual void samplingElec() = 0;
    /**
     * Get initial density (or window) at time zero of current trajectory.
     *
     * The results of it are stored in data member: init_density
     */
    virtual void getInitialDensity() = 0;
    /**
     * Update the electronic reduced density matrix at current time.
     *
     * Both the current (RDM_current) and average (RDM_average) ones are updated.
     */
    virtual void updateDensityMatrix() = 0;

    // # The following virtual functions are defined based on the Meyer-Miller
    // # mapping Hamiltonian, which can be inherited and used directly for
    // # general mapping dynamics, i.e., the subclasses: LSC, SQC, eCMM, SPM, MF.
    // # And for other methods they should be overrided by subclasses via
    // # defining identical functions explicitly (if necessary).
    // $ Model Hamiltonian with proprgation in diabatic basis
    /**
     * Propagate one step for a simulation of model using velocity Verlet integrator.
     *
     * This is original dynamics in diabatic representation.
     */
    virtual void oneStepForDiabaticModel();
    /**
     * Update the effective forces which will be used nuclear propagation.
     *
     * This is Ehrenfest-type mean-filed force derived from Hamiltonian (with
     * removing averaged potential energy), which works for most mapping dynamics
     * with two mapping variables.
     */
    virtual void updateDiabaticModelForces();
    // $ Model Hamiltonian with proprgation in quasi-diabatic basis
    /**
     * Propagate one step for a simulation of model.
     *
     * The EOMs used here are based on Meyer-Miller Hamiltonian in quasi-diabatic
     * (QD) representation by Pengfei Huo.
     *
     * References:
     * [1] J. Chem. Phys. 149, 044115 (2018)
     * [2] J. Chem. Theory Comput. 2018, 14, 1828−1840
     * [3] J. Phys. Chem. Lett. 2019, 10, 7062−7070
     * [4] J. Chem. Phys. 155, 084106 (2021)
     */
    virtual void oneStepForQuasiDiabaticModel();
    /**
     * Update the effective forces which will be used nuclear propagation.
     *
     * Do nothing here, since it has been included in oneStepForDiabaticModel().
     */
    virtual void updateQuasiDiabaticModelForces() {};
    // $ Model Hamiltonian with proprgation in adiabatic basis
    /**
     * Propagate one step for a simulation of model using velocity Verlet integrator.
     *
     * The EOMs used here are based on Meyer-Miller Hamiltonian in adiabatic
     * representation proposed by Miller [J. Chem. Phys. 147, 064112 (2017)].
     */
    virtual void oneStepForAdiabaticModel();
    /**
     * Update the effective forces which will be used nuclear propagation.
     */
    virtual void updateAdiabaticModelForces();
    /**
     * Transform the coefficient of electronic wavefunction between the diabatic
     * and adiabatic representations.
     *
     * Transfrom diabatic coefficient to adiabatic by using |i> = Σ_j T_ji |j>,
     * where i/j is the index of the adiabatic/diabatic state.
     *
     * Transfrom adiabatic coefficient to diabatic by using |i> = Σ_j T_ij |j>,
     * where i/j is the index of the diabatic/adiabatic state.
     *
     * The rotation matrix T is from Hamiltonian object.
     *
     * This function is used for a simulation of diabatic model with adiabatic
     * propgation. In the adiabatic propgation, all quantities are required to
     * be in the adiabatic representation. While the initial sampling of electronic
     * DOF is diabatic, and the calculations of reduced density matrix are based
     * on diabatic electronic DOF, too. Therefore, we need to transform the diabatic
     * coefficient to adiabatic before using them in the propagation. And after
     * propagation, before computing the reduced density matrix, we transfrom
     * the updated adiabatic electronic coefficient back to diabatic. Then, the
     * generated reduced density matrix is diabatic and can be compared directly
     * with the diabatic algorithm. In this case, we don't need to modify the
     * original sampling and calculations of observables, which is more convenient
     * than do the transformation in the initial sampling and observables.
     *
     * @param coeff_old the electronic coefficient to be transformed (read)
     * @param T         rotation/transformation matrix from Hamiltonian
     * @param coeff     the transformed electronic coefficient (write)
     * @param q         electronic mapping coordinates form coeff (write)
     * @param p         electronic mapping momenta form coeff (write)
     * @param inverse   If true, do transfromation from adiabatic to diabatic
     *                  Default is false and do transfromation from diabatic to adiabatic
     */
    virtual void transformCoefficient(const std::vector<Complex>& coeff_old, const Real_Matrix& T,
        std::vector<Complex>& coeff, std::vector<double>& q, std::vector<double>& p, bool inverse = false);
    // $ All-atom simulation with OpenMM in diabatic basis
    /**
     * Propagate one step for an all-atom simulation with OpenMM integrator using
     * external effective forces. [Plan B]
     */
    virtual void oneStepForAllAtom();
    /**
     * Propagate one step for an all-atom simulation with OpenMM integrator using
     * external effective forces. [Plan B]
     */
    virtual void oneStepForCustomAllAtom();
    /**
     * Update the effective forces (in kj/mol/nm) which will be used nuclear propagation.
     * This function is used for an all-atom simulation.
     */
    virtual void updateAllAtomForces();
    /**
     * Update the effective forces (in kj/mol/nm) which will be used nuclear propagation.
     * This function is used for an all-atom simulation.
     */
    virtual void updateCustomAllAtomForces();
    /**
     * Propagate one step for an all-atom simulation using leapfrog or velocity
     * Verlet integrator with the custom CPU code. [Plan A]
     * It should be equivalent of oneStepForAllAtom(), when constraints and
     * virtual sites are not used.
     * This function is used for reference and debug.
     */
    virtual void customCPUIntegrator();

protected:
    // The smarter pointer to HamiltonianElec object, which stores the electronic
    // mapping variables of Hamiltonian.
    std::shared_ptr<HamiltonianElec> Elec;
    // The smarter pointer to DynamicsElec object, which defines the electronic
    // propagation.
    std::shared_ptr<DynamicsElec> DyElec;
    // random number generator, which is used in electronic sampling.
    std::mt19937 elec_gen;
    // the frequency (in time steps) to update the electronic reduced density
    // matrix during the nuclear propagation, default is 1.
    int RDM_steps;
    // gamma is zero point energy (ZPE) parameter in the MMST mapping Hamiltonian.
    // It is used in mapping dynmaics, such as LSC, SQC, eCMM, SPM, Ehrenfest.
    // For LSC methods, gamma_MM is 1/2.
    // For SQC_square method, the optimal value is sqrt(3)-1)/2=0.366) and
    // For SQC_triangle method, the optimal value (1.0/3.0).
    // For Ehrenfest mean-filed method, it is zero.
    // For eCMM method, it can be specifed by user, in the range (-1/F, 1/2].
    // For SPM method, it is decided by the choice of phase space representation.
    // The specfic value of it will be defined in subclasses. If it is not
    // re-assigned in the subclasses, the value of it will be zero. We include it
    // here since it is commonly used in many mapping methods and used in general EOMs.
    // It will be used in electronic sampling and the calculation of observables
    // and equations of motions (if symmetrical Hamiltonian is used, the gamma
    // can be eliminated in the equations of motions).
    // Generally, the value of gamma should be between 0 and 1/2. When gamma=1/2,
    // it describes the whole ZPE effect. However, it has been observed that a
    // smaller value varying from 0 to 1/2, which indicates that only part of
    // ZPE is included, can lead to better results in practice[G. Stock, J. Chem.
    // Phys., 1995, 103, 2888–2902.]. Different choices of ZPE parameters as
    // well as several benchmarks have been discussed, and the optimal value of
    // gamma is found to be system-dependent. Refs: [S. J. Cotton and W. H. Miller,
    // J. Phys. Chem. A, 2013, 117, 7190–7194.] [S. J. Cotton and W. H. Miller,
    // J. Chem. Phys., 2016, 145, 144108.] [W. H. Miller and S. J. Cotton,
    // Faraday Discuss., 2016, 195, 9–30.]
    double gamma_MM;
    // The initial density matrix element indices |j><k|, e.g., "0,0"
    // Here, using single integer to represent two indices, i.e.,
    // init_state = j*DOFe+k, j=init_state/DOFe, k=init_state%DOFe.
    // Starting from all possible initial states is supported for some methods,
    // so here, defined it as a vector.
    std::vector<int> init_state;
    // Initial density cooresponding to init_state. (or window or A operator in
    // correlation function, the quantity of it is different in different methods).
    // Anyway, it is a qunatity that depends on the initial state and the phase
    // space and electronic mapping variables at time 0.
    // It keeps unchange in one trajectory, but changed in another traectory.
    // If the observables that starting from all initial states should be computed,
    // The init_density is a vector, and the index is the initial state .
    // Since some methods support all possible initial states in one simulation,
    // so here, defined it as a vector, in this case, the index of it represents
    // the init_state (using single integer j*DOFe+k to represent two indices).
    std::vector<Complex> init_density;
    // RDM stores the on-the-fly calculated reduced density matrix along the time (final results).
    // the density matrix of one step is stored as DOFe^2-dimentional vector.
    // RDM represents the RDM of current trajectory which is calculated on-the-fly
    // Since the density matrice that starting from differnent initial states
    // can be computed at one simulation for some methods, therefore it is a
    // vector of RDM. e.g., RDM[i][f], i is init_state, and f is the index of element in rho(t).
    // Here, i and f is, i (or f) = a*DOFe+b, where a, b is the subscripts of rho.
    std::vector<std::vector<Complex>> RDM, RDM_LSC1, RDM_LSC2, RDM_RILSC1, RDM_RILSC2, RDM_RILSC3;
    std::vector<std::vector<Complex>> AdiabaticRDM, AdiabaticRDM_LSC1, AdiabaticRDM_LSC2, AdiabaticRDM_RILSC1, AdiabaticRDM_RILSC2, AdiabaticRDM_RILSC3;
    // std::function is a polymorphic function wrapper. It stores a pointers to
    // member functions: oneStepForModel() or leapfrogVerletIntegrator() or
    // velocityVerletIntegrator() according to model_type and integrator parameters.
    // This binding is finished at init().
    std::function<void()> oneStep;
    // It stores a pointers to member functions: updateModelForces() or
    // updateAllAtomForces(), whose binding is finished at init().
    std::function<void()> updateF;
    // updateH() is to update Hamiltonian data that required
    // for dynamics, such as Hamiltonian matrix, forces, and/or nonadiabatic
    // coupling vectors, in diabtaic or adiabatic representation.
    // It stores a pointers to member functions in Hamiltonian object:
    // updateDiabaticHamiltonian() or updateAdiabaticHamiltonian() according to
    // representation. This binding is finished at init().
    std::function<void()> updateH;
};