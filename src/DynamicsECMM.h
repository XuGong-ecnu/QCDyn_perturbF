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
 * DynamicsECMM class defines nonadiabatic dynamics with extend classical
 * mapping models (eCMM) developed by Jian Liu [1].
 *
 * The feature of eCMM method is that it can run mapping dynamics with a specified
 * zero point energy (ZPE) parameter (γ: (-1/DOFe, ∞) in Meyer−Miller Hamiltonian.
 * A negative ZPE parameter is also supported and perform well in some case.
 *
 * Note that the CMM2 method [2-3], is also defined here, since it is the special
 * case of eCMM method with γ=0.
 *
 * The eCMM method use full-sphere electronic sampling, like spin-mapping and
 * CMM method.
 *
 * It should be note that the spin-mapping method (DynamicsSPM) with different
 * representation (Q, P, W, MMST) is special cases with specific value of γ.
 * (Note the γ defined in SPM is 2 times of standard ZPE parameter used here or
 * LSC/SQC method.) But, in SPM method, we don't need to specify the γ, sicne it
 * depends on the method.
 *
 *   The ZPE parameter used in mapping dynamics method:
 *              Method        ZPE parameter (γ)
 *                LSC           1/2
 *                SQC-square    (sqrt(3)-1)/2
 *                SQC-triangle  1/3
 *                SPM-Q         0
 *                SPM-P         1
 *                SPM-W         (sqrt(DOFe+1)-1)/DOFe
 *                SPM-MMST      1/2
 *                CMM2          0
 *                eCMM          specified by user
 *
 * Therefore, the CMM2 method (γ=0) and SPM-Q is equilvenlent.
 * And default γ value for eCMM method is 1/2, if γ is not specified by user.
 *
 * References:
 * [1] X. He, Z. Gong, B. Wu, J. Liu, J. Phys. Chem. Lett. 2021, 12, 2496.
 * [2] J. Liu, J. Chem. Phys. 2016, 145, 204105.
 * [3] X. He, J. Liu, J. Chem. Phys. 2019, 151, 024105.
 */
class DynamicsECMM : public DynamicsMQCBase {
public:
    /**
     * Construct a DynamicsECMM object.
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     */
    DynamicsECMM(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsMQCBase(param, Ha) {}
    ~DynamicsECMM() {}
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

private:
    /**
     * Get initial electronic mapping variables.
     * Same as full-sphere smapling in DynamcisSPM.
     * TODO: foucsed initial sampling
     */
    void samplingElec();
    /**
     * Get initial denisty (population or coherence operator) at time 0.
     * Same as initial denisty in DynamcisSPM.
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
};