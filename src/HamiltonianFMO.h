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
 * This class defines the Hamiltonian of the multi-level Fenna–Matthews–Olson
 * (FMO) model, which has been used extensively as a benchmark.
 *
 * Here, the FMO model is implemented in general form (or called the general
 * site-exciton framework, see Reference [3]). The DOFe is not fixed,
 * and the reorganization energy (lambda) can be used different value for different
 * states to introduce the differnent strength of system–bath coupling (c). While,
 * frequency of harmonic normal modes (omega) is still one set for all states.
 *
 * Note that, the form of petential energy is defined in terms of the relative
 * shifts in the equilibrium positions according the equation (22) in Reference [1].
 *
 * The Debye spectral density with characteristic frequency (omega_c, real unit)
 * will be used to get omega and c.
 *
 * Here, the minimum energy of each state (epsilon), the couplings (gamma_DA),
 * the reorganization energy (lambda), and characteristic frequency (omega_c) should
 * be provided in input control file (real unit). And they will be converted to
 * a.u. to generate model parameters.
 *
 * Standard FMO model is 7-level model, and each diabatic state represents an
 * exciton localized on one of the sites. You can find the paprameters couplings
 * and epsilon (in cm-1) of this model in References [2] (spin-mapping) or [3]
 * (SQC-square). And one can specify model_type=FMO7 to use it directly without
 * input epsilon and gamma_DA. The reorganization energy (lambda) and characteristic
 * frequency (omega_c) still need to be provided in input control file (real unit),
 * since they can be changed for this 7-state FMO model.
 *
 * Addtionaly, the parameters for 4- and 5-state reduced-dimensional FMO
 * models can be found in References [4] (Appendix). And one can specify
 * model_type=FMO4 or FMO5 to use it directly.
 *
 * And the general FMO model also can be as 2-state model, parameters for 2-state
 * FMO model can be found in References [3].
 *
 * References:
 * [1] Xing Gao and Eitan Geva, J. Phys. Chem. A, 124, 11006−11016 (2020).
 * [2] J. E. Runeson and J. O. Richardson, J. Chem. Phys. 152, 084110 (2020).
 * [3] Stephen J. Cotton and William H. Miller, J. Chem. Theory Comput. 12, 983−991 (2016).
 * [4] Stephen J. Cotton and William H. Miller, J. Chem. Phys. 150, 104101 (2019)
 */
class HamiltonianFMO : public HamiltonianModelBase {
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
     * Construct an HamiltonianFMO object.
     *
     * @param param   the global paramters
     */
    HamiltonianFMO(Parameters& param) : HamiltonianModelBase(param) {}
    ~HamiltonianFMO() {}
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
     * Generate model paramters based on input paprameters and Debye spectral
     * denisty.
     */
    void buildModelParameters();
    /**
     * Load model paramters from external file.
     * The format of file can be csv file or dat file with fixed width.
     *
     * Requirments of file:
     * 1. First line is headers: index,omgea/au,c,req
     * 2. The number of lines is N+1 (N is number of nromal modes) for one set
     *    common c/req, or DOFn+1 (DOFn=DOFe*N) for DOFe sets c/req.
     * 3. The unit is au.
     */
    void loadModelParameters(const std::string& loadfile);
    /**
     * Save model paramters to external file.
     * The format of output is same as the descripition in loadModelParameters().
     */
    void saveModelParameters(const std::string& savefile);
};