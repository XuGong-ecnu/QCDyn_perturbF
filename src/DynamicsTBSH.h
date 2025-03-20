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
#include "DynamicsMQCBase.h"

/**
 * This class defines mixed quantum-classical Liouville dynamics with
 * Trotter-Based Surface-Hopping (TBSH)[1] algorithm.
 *
 * If the bath is harmonic and bilinearly coupled to the subsystem, then MQCL
 * dynamics entails no approximation and is equivalent to a full quantum
 * description of the dynamics of the entire system. And the momentum-jump
 * approximation is the only approximation used in simulations of the
 * quantum-classical dynamics with TBSH algorithms.
 *
 * In practice the applicability of this approach has been limited to short
 * times owing to unfavorable numerical scaling. There are some ways to combine
 * TBSH-MQCL method with other method to generate long time population dynamics,
 * such as generalized quantum master equation (GQME)[2], transfer tensor method
 * (TTM)[3].
 *
 * Note that the TBSH-MQCL method is defined in adiabatic base, i.e, the physical
 * quantities (including state/index/density) are in adiabatic base in the equation
 * 35 in orignal TBSH literature[1]. We need do the transformation for energy/force/denisty
 * between diabatic and adiabatic base via the unitary transformation matrix if
 * the Hamiltonian is defined in diabatic base.
 *
 * Unlike other methods, there is no explicit electron position/monmentum/coefficient,
 * since the adiabatic electronic density matrix element is propagated directly.
 *
 * Reference:
 * [1] D. Mac Kernan, G. Ciccotti, R. Kapral, J. Phys. Chem. B 2008, 112, 424.
 * [2] A. Kelly, T. E. Markland, J. Chem. Phys. 2013, 139, 014104.
 * [3] A. A. Kananenka, C.-Y. Hsieh, J. Cao, E. Geva, J. Phys. Chem. Let. 2016, 7, 4809.
 */
class DynamicsTBSH : public DynamicsMQCBase {
public:
    /**
     * Construct a DynamicsTBSH object.
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     */
    DynamicsTBSH(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsMQCBase(param, Ha) {}
    ~DynamicsTBSH() {}
    /**
     * Initialize data members and check the legality of parameters.
     */
    void init();
    /*
    * Here we set up a function to get RDM data member.
    * Notice here RDM here serves as a two-dimensional matrix, RDM
    * Here we take out one of the matrix at given time snapshot, RDM[i][j*DOFe+k],
    * where the first dimensionality gives the position in the initial state(s), 
    * the second dimensionality gives the order of the state which gets correlation with
    * the initial state[i], j*DOFe+k --> |j><k|.
    * 
    */
    const std::vector<std::vector<Complex>>& getRDM() {
        return RDM;
    }
    const std::vector<std::vector<Complex>>& getAdiabaticRDM() {
        return AdiabaticRDM;
    }

protected:
    /**
     * Generate uniform sampling of the initial index for pairs of quantum
     * states (s_init) and initialize s_prime_init, s_previous, s_current.
     * They are s_0, s'_0, s_j-1, s_j, respectively, in orignal TBSH literature
     * [JPCB,2008,112,424]
     */
    void samplingElec();
    /**
     * Get adiabatic electronic density matrix element of s_prime_init (s'_0)
     * via the transformation from diabatic initial denisty matrix to adiabatic.
     *
     * Specifically, it is the element with index s_prime_init in initial adiabatic
     * electronic denisty matrix (transform rho_init_diabatic to adiabtic).
     * Note that the rho(0) spcified by user is in diabatic base, so,
     * transform it to adiabatic via the unitary transformation matrix which is
     * the same matrix to transform diabatic Hamiltonian to adiabatic.
     * It keeps unchange in one trajectory, but changed in another traectory.
     * If the observables that starting from all initial states should be computed,
     * The init_density is a vector, and the index is the initial state.
     * The size of init_denisty depends on the initial state (one or DOFe*DOFe)
     *
     * The unitary transformation matrix is from Hamiltonian object, which is
     * used to transform diabatic Hamiltonian matrix to adiabatic Hamiltonian.
     * And it should be updated along the trajectory (time-dependent).
     *
     * The initial diabatic electronic denisty matrix is stored in data member:
     * rho_init_diabatic.
     * The value of result is stored in data member: init_density.
     */
    void getInitialDensity();
    /**
     * Update the adiabatic electronic reduced density matrix of current step.
     * And upadte the diabatic density matrix via the transformation from
     * adiabatic denisty matrix to diabatic.
     *
     * The unitary transformation matrix is from Hamiltonian object, which is
     * used to transform diabatic Hamiltonian matrix to adiabatic Hamiltonian.
     * And it should be updated along the trajectory (time-dependent).
     *
     * Both the current (RDM_current) and average (RDM_average) RDM within
     * diabatic and adiabatic representations are updated.
     *
     * How to compute obeservable or density matrix: rho(0)*rho(t)
     * Here, rho(0) is init_density, depends on the rho_init, s_prime_init
     * (s'_0) and R(0), P(0). rho(t) is density matrix of current step, depends
     * on the s_current and R(t), P(t). Both of them are in adiabatic base.
     *
     * Note that the elements is in adiabatic base in eq.35 in orignal TBSH
     * literature [JPCB,2008,112,424]. The unit operator of s_current element in
     * density matrix is 1, and the others are 0. The factor (complex) of it is
     * phase_factor*init_density. phase factor is a complex, and this behavoir
     * will update the real and imag part at the same time. So, for adiabatic
     * denisty matrix, just update and accumulate the single element of s_current.
     *
     * As for diabatic denisty matrix, we update it via the transformation from
     * adiabatic denisty matrix to diabatic. We transform the adiabatic unit matrix
     * firstly to diabatic (after transformation, all elemnets of diabtaic density
     * matrix have values), then multiply the phase_factor*init_density to all
     * elemnets of diabtic density matrix (the factor are same).
     */
    void updateDensityMatrix();
    /**
     * Update the effective forces which will be used nuclear propagation.
     * This function is used for a simulation of model.
     *
     * For TBSH-MQCL method, adiabatic forces are used.
     * If s_current is population, use forces of this state.
     * If s_current is coherence, use average forces of thes two states.
     */
    void updateAdiabaticModelForces() override;
    /**
     * Propagate one full step for a simulation of model, including:
     * (1) Adiabtic nuclear propagation with half-time step;
     * (2) Nonadiabatic propagation to decide s_j;
     * (3) Adiabtic nuclear propagation with another half-time step;
     * The phase factor will also be updated during propagation.
     */
    void oneStepForAdiabaticModel() override;
    /**
     * Adiabatic nuclear propagation with half-time step for model.
     * It is the first part and third part of in oneStepForAdiabaticModel().
     */
    void AdiabaticPropagateForModel();
    /**
     * Nonadiabatic propagation to decide s_j.
     * Filter and momentum shift (if any) will be applied here.
     * It is the second part in oneStepForAdiabaticModel().
     *
     * If no sufficent momentum, set probability to 0 directly before decide s_j
     * by stochastic probability algorithm.
     */
    void NonadiabaticPropagateForModel();

protected:
    // s_0, s'_0, s_j-1, s_j in orignal TBSH literature [JPCB,2008,112,424]
    // index for pairs of quantum states (index in DOFe^2-dimentional electornic
    // denisty vector), s= a*DOFe+b, s' = b*DOFe+a, a, b is the index of state,
    // DOFe is number of states. s' is the s with inerchanged indices.
    // Here, 0 <= a,b <= (DOFe-1), 0 <= s < (DOFe^2-1). a,b,s is interger.
    // s_init, s_previous, s_current, is used in time propagator
    // s_prime_init is used in initial electronic density matrix element.
    // Note that the state/iondex is in the adiabatic base.
    // s_init (random), s_prime_init keep unchange in one trajectory, but changed
    // in another traectory. s_previous, s_current will be updated along trajectory.
    int s_init, s_prime_init, s_previous, s_current;
    // Current indices of adiabatic states, a, b, from current pair of quantum
    // states (s_j). i.e., s_current = a*DOFe+b. They will be updated during
    // propagation as long as s_current changes.
    int a, b;
    // phase_factor is the factor of operator B^{s_N}_w(R,P) of s_current in
    // eq.35 in orignal TBSH literature [JPCB,2008,112,424], the normalize and average
    // fatcor N^2/M are also included, so the initial value of phase factor for each
    // trajectory is DOFe^2/ntraj. ntraj is total number of trajectories.
    // It will be updated along the trajectory.
    Complex phase_factor;
    // The initial diabatic density matrix (rho(0)) whose initial value are
    // specified by user via two indices, |j><k|, e.g., "0,0", which means rho_00
    // is 1, and other elements are 0.
    // It is constant in a simulation for all trajectory.
    // If init_state=all, then all possible initial states are considered, which
    // is useful when Transfer Tensor Method (TTM) will be applied.
    std::vector<Real_Matrix> rho_init_diabatic;
    // The RDM within diabatic and adiabatic representations.
    // Note, only the diabtaic one will be return when calling to get density
    // matrix.
    // std::vector<Complex_Matrix> RDM_current_diabatic, RDM_current_adiabatic,RDM_current;
    // std::vector<Complex_Matrix> RDM_average_diabatic, RDM_average_adiabatic,RDM_average;
    // If filter_bound > 0, then apply filter in the propagation.
    // By default, it is 0, which means no filter will use.
    // The filter should be used for medium and long time simulations.
    // And the value of bound should be tested until convergence is obtained.
    // (See discussion of filter in TBSH [JPCB,2008,112,424]).
    double filter_bound;
    // random number generator, which is used hopping probability.
    // Within TBSH method, the surface hopping probability at each step is
    // decided by a stochastic algorithm. So, a random number generator
    // is required.
    std::mt19937 TBSH_gen;
};