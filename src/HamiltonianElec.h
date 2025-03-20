/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Oct. 28, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "Tools.h"
#include "Parameters.h"

/**
 * An HamiltonianElec class defines the electronic mapping variables and coefficient
 * of electronic wavefunction used in different mapping methods. It is the electronic
 * part of Hamiltonian used in mapping dynamics (same for all-atom or model system).
 */
class HamiltonianElec {
    friend class DynamicsElec;
    friend class DynamicsBase;
    friend class DynamicsMQCBase;
    friend class DynamicsMF;
    friend class DynamicsLSC;
    friend class DynamicsSQC;
    friend class DynamicsSPM;
    friend class DynamicsCMM;
    friend class DynamicsECMM;
    friend class DynamicsFSSH;
    friend class DynamicsMFRDM;
    friend class DynamicsECMMCV;
    
public:
    /**
     * Construct an HamiltonianElec object
     *
     * @param param   the global paramters
     */
    HamiltonianElec(Parameters& param) : param(param) {}
    ~HamiltonianElec() {}
    /**
     * Initialize data members of HamiltonianElec class.
     */
    void init();
    /**
     * Get electronic mapping variable positions (read only).
     *
     * @return reference to the q in HamiltonianElec object
     */
    const std::vector<double>& get_q() const {
        return q;
    }
    /**
     * Get electronic mapping variable momenta (read only).
     *
     * @return reference to the p in HamiltonianElec object
     */
    const std::vector<double>& get_p() const {
        return p;
    }
    /**
     * Get coefficient of electronic wavefunction (read only).
     *
     * @return reference to the coeff in HamiltonianElec object
     */
    const std::vector<Complex>& get_coeff() const {
        return coeff;
    }

private:
    // Parameters object controls the simulation
    Parameters& param;
    // electronic DOF = number of states/topologies/surfaces.
    int DOFe;
    // electronic mapping variables: positions and momenta (LSC and SQC method)
    // Note that, for SQC method, the electronic mapping variables are action (n)
    // and angle (u), but they can be expressed to variables q and p that are used
    // in propagation. Therefore, we only store mapping variable q and p here.
    std::vector<double> q, p;
    // coefficient of electronic wavefunction, coeff[j]=(q[j]+I*p[j])/sqrt(2)
    // It is used drectly in Ehrenfest/surface hopping methods.
    // And it is can be propageted instead of q, p.
    std::vector<Complex> coeff;
    // electronic mapping variables used in Jian Liu's CMM methods, Refs:
    // [1] J. Liu, J. Chem. Phys. 2016, 145, 204105.
    // [2] X. He, J. Liu, J. Chem. Phys. 2019, 151, 024105.
    std::vector<double> x, y, px, py;
    // commutator variables used in Jian Liu's eCMMcv methods, Refs:
    // [1] X. He, B. Wu, Z. Gong, J. Liu, J. Phys. Chem. A. 2021, 125, 6845. 
    std::vector<double> xcv, pcv, scv;
};