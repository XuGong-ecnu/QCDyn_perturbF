/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: Juen. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsRPMD.h"

void DynamicsRPMD::init() {
    std::cout.precision(12);
    DynamicsBase::init();
    if (system_type != "allatom")
        throw std::runtime_error("ERROR: Unsupported system_type=" + system_type + " for DynamicsRPMD.");
    if (allatom_type != "CustomForceField")
        throw std::runtime_error("ERROR: Unsupported allatom_type=" + allatom_type + " for DynamicsRPMD.");
     if (PIMD_type != "RPMD" || PIMD_type != "PIMD" || PIMD_type != "CMD")
        throw std::runtime_error("ERROR: Unsupported PIMD_type=" + PIMD_type + " for DynamicsRPMD.");   
    if (dyn_type != "MD")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for DynamicsRPMD.");   
    integrator = param.getStr("integrator");
    // Create smarter pointer to RingPolymer intergrator
    ha = std::static_pointer_cast<HamiltonianRingPolymer>(Ha);
    PIMD_thermostat = param.getStr("PIMD_thermostat");
    ilngv_flag = param.getInt("ilngv_flag");
    C_mat.resize(nbeads * nbeads);
    Evo_mat.resize(DOFn * nbeads * 2 * 2);
    c_lngv1.resize(nbeads);
    c_lngv2.resize(nbeads);
    CMD_sigma.resize(nbeads);
    RingPolymer_init(ha->masses, ha->omega_n);
}

void DynamicsRPMD::beforeOneTraj() {
    step = 0; // reset current step to 0.
    ha->updateAllForces(); 
}

void DynamicsRPMD::dynamics(int steps) {
    for (int i = 0; i < steps; ++i) {
        myOneStep();  
    }
    this->step += steps;
}

void DynamicsRPMD::myOneStep() {
    if (integrator == "velocityVerlet") {
        // TO DO: for half step Langevin thermostat 
        if (PIMD_thermostat == "Langevin")
            RingPolymer_Langevin(ha->V_RP, ha->masses, ha->beta_n);
        // Update velocities with a half step using effective force.
        RP_MOVE_Half_V(ha->masses, ha->F_RP, ha->V_RP);
        //propagate RP in normal modes
        RingPolymer_nm_propagate(ha->masses, ha->R_RP, ha->V_RP);
        // Update Forces with new positions
        ha->updateAllForces();
        // Update momenta with a half step using new effective forces
        RP_MOVE_Half_V(ha->masses, ha->F_RP, ha->V_RP);// masses: cmd and RPMD different
    }
    else if (integrator == "leapfrog") {
        // Update velocities with a full step using effective force.
        RP_MOVE_V(ha->masses, ha->F_RP, ha->V_RP);
        //propagate RP in normal modes
        RingPolymer_nm_propagate(ha->masses, ha->R_RP, ha->V_RP);
        // At the start of each time step, we need to update force
        ha->updateAllForces(); 
    }
    else
        throw std::runtime_error("ERROR: Unsupported integrator=" + integrator + " for CustomForceField.");
}

void DynamicsRPMD::RingPolymer_nm_propagate(const std::vector<double>& masses, 
                std::vector<std::vector<Vec3>>& R_RP, std::vector<std::vector<Vec3>>& V_RP) {
    
    double p_nm[nbeads][DOFn][3];
    double q_nm[nbeads][DOFn][3];
    double p, q;
    //Transform to normal mode basis
    for (int i = 0; i < DOFn; i++)
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < nbeads; k++) {
                p_nm[k][i][a] = 0;
                q_nm[k][i][a] = 0;
            }
    
    for (int i = 0; i < DOFn; i++)
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < nbeads; k++) {
                for (int j = 0; j < nbeads; j++) {
                    p_nm[k][i][a] += V_RP[j][i][a] * C_mat[j * nbeads + k];
                    q_nm[k][i][a] += R_RP[j][i][a] * C_mat[j * nbeads + k];
                }
                p_nm[k][i][a] *= masses[i] * CMD_sigma[k] * CMD_sigma[k];;
            }
    //Exact evolution of RP normal modes for one time step
    for (int i = 0; i < DOFn; i++)
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < nbeads; k++) {
                p  =  Evo_mat[nbeads * 4 * i + 4 * k + 0] * p_nm[k][i][a]
                    + Evo_mat[nbeads * 4 * i + 4 * k + 1] * q_nm[k][i][a];
                q  =  Evo_mat[nbeads * 4 * i + 4 * k + 2] * p_nm[k][i][a]
                    + Evo_mat[nbeads * 4 * i + 4 * k + 3] * q_nm[k][i][a];
                p_nm[k][i][a] = p;
                q_nm[k][i][a] = q;
            }
    //Transform back to Cartesian coordinate
    for (int i = 0; i < DOFn; i++)
        for (int a = 0; a < 3; a++)
            for (int j = 0; j < nbeads; j++) {
                V_RP[j][i][a] = 0;
                R_RP[j][i][a] = 0;
            }
    for (int i = 0; i < DOFn; i++)
        for (int a = 0; a < 3; a++)
            for (int j = 0; j < nbeads; j++)  {
                for (int k = 0; k < nbeads; k++) {
                    V_RP[j][i][a] += C_mat[j * nbeads + k] * p_nm[k][i][a] / CMD_sigma[k] / CMD_sigma[k];
                    R_RP[j][i][a] += C_mat[j * nbeads + k] * q_nm[k][i][a];
                }
                V_RP[j][i][a] /= masses[i];
            }
}
 
void DynamicsRPMD::RingPolymer_init(const std::vector<double>& masses, const double& omega_n) {

    //Set harmonic constants of normal modes of ring polymer
    std::vector<double> omega_nm;
    
    omega_nm.resize(nbeads);
    for (int k = 0; k < nbeads; k++) {
        omega_nm[k] = 2 * omega_n * sin( k * pi / nbeads);
    }

    if (PIMD_type == "RPMD") {
        for (int k = 0; k < nbeads; k++) CMD_sigma[k] = 1.0 ;
    }
    else if (PIMD_type == "CMD") {
        CMD_sigma[0] = 1.0 ;
        double mmm;
        mmm = pow(nbeads, nbeads/(nbeads-1));
        for (int k = 1; k < nbeads; k++) { 
            double beta1 = ha->beta_n;
            // This unit need to check!
            CMD_sigma[k] = omega_nm[k] * hbar_SI * beta1 /mmm ; // / nbeads**(nbeads/(nbeads-1));
        }
    }
    else {
        // The PIMD_type == PIMD is not supported in this part.
        throw std::runtime_error("ERROR: Unsupported PIMD_type=" + PIMD_type + " for DynamicsRPMD.");  
    }

    //Transformation matrix between original beads to normal mode beads
    //k=0,...,n-1, but j=1,...,n. in Eq.(17) in Ceriotti, JCP 133, 124104 (2010)
    //but in this code, j=0,...,n-1
    
    for (int j = 0; j < nbeads; j++) {
        C_mat[j * nbeads + 0] = sqrt(1.0 / nbeads);//k=0
        C_mat[j * nbeads + nbeads / 2] = sqrt(1.0 / nbeads) * pow(-1, j + 1);//k=n/2
        for (int k = 1; k < nbeads / 2; k++)
            C_mat[j * nbeads + k] = sqrt(2.0 / nbeads) * cos(2 * pi * (j + 1) * k / nbeads);
        for (int k = nbeads / 2 + 1; k < nbeads ; k++)
            C_mat[j * nbeads + k] = sqrt(2.0 / nbeads) * sin(2 * pi * (j + 1) * k / nbeads);
    }
    //Evolution matrix for normal modes p and q
    for (int i = 0; i < DOFn; i++) {
        for (int k = 0; k < nbeads; k++) {
            if (omega_nm[k] == 0) {
                Evo_mat[nbeads * i * 4 + k * 4 + 0] = 1.0;
                Evo_mat[nbeads * i * 4 + k * 4 + 1] = 0;
                Evo_mat[nbeads * i * 4 + k * 4 + 2] = DT / masses[i];
                Evo_mat[nbeads * i * 4 + k * 4 + 3] = 1.0;
            } else {
                Evo_mat[nbeads * i * 4 + k * 4+ 0] = cos(omega_nm[k] * DT);
                Evo_mat[nbeads * i * 4 + k * 4 + 1] = -1 * masses[i] * omega_nm[k] * sin(omega_nm[k] * DT);
                Evo_mat[nbeads * i * 4 + k * 4 + 2] = 1.0/(masses[i] * omega_nm[k]) * sin(omega_nm[k] * DT);
                Evo_mat[nbeads * i * 4 + k * 4 + 3] = cos(omega_nm[k] * DT);
            }
        }
    }
    //Langevin thermostat coefficients
    double gamma_0, gamma_k;
    // This is for Thermostated RPMD to all normal modes.
    double tau_lngv_ps;
    tau_lngv_ps = param.getDouble("tau_lngv_ps");
    // No thermostat.
    if (ilngv_flag == 0) gamma_0 = 0;
    // Thermostated RPMD to all normal modes.
    if (ilngv_flag == 1) gamma_0 = 1.0 / tau_lngv_ps ;
    // TRPMD, thermostat to higher freq normal modes, not to centroid.
    if (ilngv_flag == 2) gamma_0 = 0; 

    if (nbeads > 1) {
        for (int k = 0; k < nbeads; k++) {
            if (k == 0) {
                c_lngv1[0] = exp(-DT * 0.5 * gamma_0);
                c_lngv2[0] = sqrt(1 - c_lngv1[0] * c_lngv1[0]);
            }
            else {
                gamma_k = 2 * omega_nm[k];
                c_lngv1[k] = exp(-DT * 0.5 * gamma_k);
                c_lngv2[k] = sqrt(1 - c_lngv1[k] * c_lngv1[k]);
            }
        }
    }
}

void DynamicsRPMD::RingPolymer_Langevin(std::vector<std::vector<Vec3>>& V_RP, 
                    const std::vector<double>& masses, const double& beta_n) {
    double fct;
    std::normal_distribution<double> normal_dist_xi(0.0,1.0);
    std::mt19937 gen;
    gen.seed(step * 273 + 1453);//1894
    double p_nm[nbeads][DOFn][3];
    //Transform to normal mode basis
    for (int i = 0; i < DOFn; i++)
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < nbeads; k++) {
                p_nm[k][i][a] = 0;
            }
    for (int i = 0; i < DOFn; i++)
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < nbeads; k++) {
                for (int j = 0; j < nbeads; j++) {
                    p_nm[k][i][a] += V_RP[j][i][a] * C_mat[j * nbeads + k];
                }
                //p_nm[k][i][a] *= masses[i]; //Note: we are using NM velocity for thermostat
            }
    for (int i = 0; i < DOFn; i++) {
        fct = 1.0 / sqrt(masses[i] * beta_n);
        //we use velocity here, factor = sqrt(m/beta_n)-> sqrt(1/m beta_n) in Eq. 28 in JCP 133, 124104 (2010).
        for (int a = 0; a < 3; a++)
            for (int k = 0; k < nbeads; k++) {
                p_nm[k][i][a] = c_lngv1[k] * p_nm[k][i][a] + fct * c_lngv2[k] * normal_dist_xi(gen);
            }
    }
    //Transform back to Cartesian coordinate
    for (int i = 0; i < DOFn; i++)
        for (int a = 0; a < 3; a++)
            for (int j = 0; j < nbeads; j++) {
                V_RP[j][i][a] = 0;
            }
    for (int i = 0; i < DOFn; i++)
        for (int a = 0; a <3; a++)
            for (int j = 0; j < nbeads; j++)  {
                for (int k = 0; k < nbeads; k++) {
                    V_RP[j][i][a] += C_mat[j * nbeads + k] * p_nm[k][i][a];
                }
            }
}
