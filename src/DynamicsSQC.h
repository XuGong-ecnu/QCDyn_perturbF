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
 * This class defines the symmetrical quasi-classical (SQC) mapping dynamics.
 *
 * For SQC mapping dynamics, there are 2 window types[1]: square and triangle,
 * and the latter may be a bit better (in most cases, they are very simillar).
 * The type of window will decide the initial electronic sampling and
 * calculation of observables.
 *
 * The results of SQC method for 7 state FMO model can be found in ref[2].
 *
 * Unlike LSC, for SQC method, only one reduced density matrix (RDM) of specified
 * initial state and window type can be computed in one simulation.
 *
 * It should be noted that the final RDM should be renormalized (including
 * coherence).
 *
 * Reference:
 * [1] S. J. Cotton, W. H. Miller, J. Chem. Phys. 2013, 139, 234112.
 * [1] S. J. Cotton, W. H. Miller, J. Chem. Phys. 2016, 145, 144108.
 * [2] S. J. Cotton, W. H. Miller, J. Chem. Phys. 2019, 150, 104101.
 */
class DynamicsSQC : public DynamicsMQCBase {
public:
    /**
     * Construct a DynamicsSQC object.
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     */
    DynamicsSQC(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsMQCBase(param, Ha) {}
    ~DynamicsSQC() {}
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
     * Get initial sampling of electronic mapping variables aacording to
     * window type and initial state.
     */
    void samplingElec();
    /**
     * For SQC method, the init_density is 1 when starting from population,
     * and it is the phase factor when starting from coherence.
     */
    void getInitialDensity();
    /**
     * Update the electronic reduced density matrix at current time.
     *
     * Both the current (RDM_current) and average (RDM_average) ones are updated.
     */
    void updateDensityMatrix();

private:
    // The type of window function, square or triangle.
    std::string window_type;
    // norm_pop and norm_coh are used to renormalize the population and coherence
    // in final averaged reduced denisty matrix, respectively. This normalization
    // makes total population of all states along the time is always 1.
    // This alleviates the issue of “wondering” trajectories that sample regions
    // of action space that do not lie within windows, in which the trajectories
    // lie outside the window has no contribution to final RDM.
    // Here, the value of norm_pop and norm_coh is the number of trajectories
    // lie within population or coherence window. Then, final RDM is normalized
    // by: RDM_norm(t) = RDM_raw(t) / norm(t).
    // Ref: J. Chem. Phys. 148, 181102 (2018) and the code in github:
    // https://github.com/jprov410/mqds,
    // file: mqds-master/mqds/src/general_src/windows.f90
    std::vector<double> norm_pop, norm_coh;
};