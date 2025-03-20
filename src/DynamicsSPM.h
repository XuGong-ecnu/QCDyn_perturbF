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
 * DynamicsSPM class define nonadiabatic dynamics with generalized spin-mapping
 * approach developed by Richardson[1-2].
 *
 * Since the correlation functions can be defined via the Q-, P-, and W- phase
 * space functions , there are 3 verisons of spin mapping methods (call them
 * the Q-, P-, and W-methods). All three methods are equivalent for an isolated
 * electornic system, but different when coupled to a nuclear environment.
 *
 * And there are 2 methods to generate initial electronic sampling, one is focused
 * method and another is full-sphere sampling method. The full-sphere method treats
 * initial conditions via a weighting procedure and sampling over the full sphere,
 * which is also used in Jian's CMM method (J. Chem. Phys. 2019, 151, 024105).
 * While, in the focused method, the initial distribution is limited to the points
 * in phase space which correspond directly to a pure initial state (population=1),
 * and let the phase space of other states is 0. The focused method is used in
 * Ehrenfest, LSC, SQC method. The advantage of full-sphere method is that it can
 * get any transition between n and m states in one simulation. One advantage of
 * using focused initial conditions is typically requires an order of magnitude
 * fewer trajectories to converge compared to full-sphere method, and it is possible
 * to define focused initial conditions also when starting from off-diagonal
 * elements of the density matrix (coherence).
 *
 * Therefore, there are 6 vesions of spin-mapping methods. For full-sphere initial
 * condition, the alternatives are to calculate the correlation function
 * <A_{s,0}B_{¯s,t}> in a symmetric way, meaning (s, ¯s)=(W, W) [W-method], or
 * in an asymmetric way, that is, (s, ¯s)=(Q, P) [Q-method] or (P, Q) [P-method].
 * The Jian's SPM can be described as (Q, Q), i.e., they used the Q-representation
 * at both initial and final times. Note that for focused methods the observable
 * must be calculated with the phase space function at initial and final time.
 * Among thess methods, the symmetric definition W-method is the most accurate
 * and robust. The (Q, P) sometimes can give same accuracy as W-method, but is
 * not always so reliable. While, (P, Q) was always less accurate. And these two
 * initial conditions give practically identical results for W-method. And as long
 * as the zero-point energy parameter is chosen to be γW, this asymmetric focusing
 * procedure is no worse than the symmetric full-sphere approach.
 *
 * The EOM of spin-mapping dynamics turns out to be equivalent to that of the
 * Meyer–Miller–Stock–Thoss (MMST) Hamiltonian, but with a new zero-point energy
 * parameter γ, for which is derived a closed formula as a function of DOFe. Note
 * that the γ in spin-mapping is different from the γ in MMST method (LSC/SQC).
 * The relation between them is γ_SPM = 2 * γ_MMST. But the value of γ has no
 * impact on the calculation of effective force.
 *
 * Note that the focused Q-method (γ=0) has the same initial distribution, dynamics,
 * and observables as the Ehrenfest method, the two methods are equivalent.
 *
 * The spin-mapping approach automatically treats the identity operator in the
 * same way as used in RI-LSC method, and Q-method can gives essentially the
 * same results as RI-LSC method. Like RI-LSC method, the total population is
 * always 1 in any trajectory. The methods are however not equivalent.
 *
 * The spin-mapping method can predict population dynamics in benchmark systems to
 * similar accuracy as other state-of-the-art mapping approaches such as SQC and
 * traceless MMST (RI-LSC). The advantage of it is no selection of γ, window, and
 * no projection operator used.
 *
 * Note that in this class, the hbar is 1 and not inlcuded explictly.
 *
 * References:
 * [1] J. E. Runeson and J. O. Richardson, J. Chem. Phys. 151, 044119 (2019)
 * [2] J. E. Runeson and J. O. Richardson, J. Chem. Phys. 152, 084110 (2020).
 */
class DynamicsSPM : public DynamicsMQCBase {
public:
    /**
     * Construct a DynamicsSPM object.
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     */
    DynamicsSPM(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsMQCBase(param, Ha) {}
    ~DynamicsSPM() {}
    /**
     * Initialize data members and check the legality of parameters.
     */
    void init();

private:
    /**
     * Get initial sampling of electronic mapping variables denpends on the
     * SPM_type (gamma_init) and/or init_state.
     */
    void samplingElec();
    /**
     * Get initial density at time zero (ρ(0)) according to initial state.
     *
     * The results of it are stored in data member: init_density
     */
    void getInitialDensity();
    /**
     * Update the electronic reduced density matrix at current time.
     * 
     * Both the current (RDM_current) and average (RDM_average) ones are updated.
     */
    void updateDensityMatrix();

private:
    // The method of electronic sampling (mapping variable positions and momenta).
    // The allowed methods are focused (only one initial population state is allowed)
    // and full (full-sphere) (starting from all possible initial states is allowed).
    // By defalut, full-sphere initial conditions is used.
    std::string elec_sample;
    // Here, SPM_type is Q-, P-, or W- phase space representation, which also
    // decides the type of spin-mapping method: Q-, P-, or W-method.
    // Note Q method with focused initial sampling is a equivalent as Ehrenfest.
    // By defalut, W-method will be used, which is recommanded.
    // In particular, window=MMST means standard MMST Hamiltonian, but with the
    // focused initial conditions method. Accutually. the SPM_type will decide
    // the choice of gamma_init, so only one window type can be used in one simulation.
    std::string SPM_type;
    // gamma is zero point energy parameter. The value of it is dependent on the
    // choiced representation. See: J. Chem. Phys. 152, 084110 (2020), Table I
    // Q: gamma_s = 0; W: gamma_s = 2/DOFe*(sqrt(DOFe+1)-1); P: gamma_s = 2
    // Here, gamma_init (gamma_s) is used in initial sampling/operator and Hamiltonian
    // And the gamma_final is used to calculate final operator (B(t))
    // For focused initial conditions, gamma_final = gamma_init
    // For full-sphere initial conditions, gamma_final is dual of gamma_init
    // i.e., Q_dual = P, P_dual = Q, W_dual = W.
    // In particular, gamma_s=1 is standard MMST, only for focused method in
    // Richardson's paper, but in Jian's J. Phys. Chem. Lett. 2021, 12, 2496−2501
    // paper, there is a general form for full-sphere (eCMM method), and here implemented.
    // Noth the gamma in spin mapping method is 2*gamma_MM (gamma used in LSC/SQC)
    // The values will affect the initial sampling and computation of RDM, but
    // has no impact on the effective force used for nulear propagation, since
    // the value is constant (you can inlcude it in your caluclation of nulcear
    // force like in the paper, but no influence, just like LSC/SQC method.).
    double gamma_init, gamma_final;
};