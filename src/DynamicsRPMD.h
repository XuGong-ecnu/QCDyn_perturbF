/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "DynamicsBase.h"
#include "Vec3.h"

class DynamicsRPMD : public DynamicsBase {
public:
    /**
     * Construct a DynamicsElec object.
     *
     * @param param   the global paramters
     */
    DynamicsRPMD(Parameters& param, std::shared_ptr<HamiltonianBase> Ha) : DynamicsBase(param, Ha) {}
    ~DynamicsRPMD() {}
    /**
     * Initialize data members and check the legality of parameters.
     */
    void init();
     /**
     * Do nothing for a general MD simulation with ForceFieldBase.
     */
    void beforeOneTraj();
    /**
     * Advance a simulation through time by taking a series of time steps.
     *
     * @param steps the number of time steps to take silently
     */
    void dynamics(int steps);
    /**
     * Do nothing for a general MD simulation.
     */
    void afterOneTraj() {}
private:
    /**
     * Do nuclear propagation with one step with providing the external forces
     * (instead of the forces which is always computed from force fileds based
     * on current positions) at each step.
     */
    void myOneStep();
    void RingPolymer_init(const std::vector<double>& masses, const double& omega_n);
    void RingPolymer_nm_propagate(const std::vector<double>& masses, std::vector<std::vector<Vec3>>& R_RP, std::vector<std::vector<Vec3>>& V_RP);
    void RingPolymer_Langevin(std::vector<std::vector<Vec3>>& V_RP, const std::vector<double>& masses, const double& beta_n);
private:
    std::vector<double> C_mat, Evo_mat, c_lngv1, c_lngv2;
    // The smarter pointer to HamiltonianForceFieldBase object.
    std::shared_ptr<HamiltonianRingPolymer> ha;
    // Integrator implements an algorithm for advancing the simulation through time, 
    // such as leapfrog Verlet, velocity Verlet, Langvin integrators.
    std::string integrator;
    std::string PIMD_thermostat;
    int ilngv_flag;
    std::vector<double> CMD_sigma;
};