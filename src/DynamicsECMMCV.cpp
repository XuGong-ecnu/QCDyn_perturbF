/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zengkui Liu @Sun Group @NYU-SH                                     *
 * First created: Dec. 14, 2021                                               *
 * Last updated: Dec. 27, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsECMMCV.h"

void DynamicsECMMCV::init(){
    int DOFe = param.getInt("DOFe");
    double dofe = static_cast<double>(DOFe);
    // 0.Initialization of data members of DynamicsMQCBase
    DynamicsMQCBase::init();
    if (dyn_type != "eCMMcv")
        throw std::runtime_error("ERROR: Unsupported dyn_type" + dyn_type + "for DynamicsECMMCV");
    // 1. Read and check parameters for eCMMcv mapping dynamics
    //  Here we should specify the model type and the dynamics type
    if (representation == "quasi-diabatic")
        throw std::runtime_error("ERROR: Only diabatic representation is supported for dyn_type=" + dyn_type);
    if (system_type != "model")
        throw std::runtime_error("ERROR: Realistic system is not supported for dyn_type=" + dyn_type);
    // Bind functions of electronic propagation according to dyn_type
    switch (representation[0]) {
        case 'd': //diabatic
            MOVE_elec = std::bind(&DynamicsECMMCV::MOVE_elec_eCMMcv, this); break;
        // case 'a': //adiabatic
        //     MOVE_elec = std::bind(&DynamicsECMMCV::MOVE_elec_eCMMcv1, this); break;

    }
    // If γ is not specified, we employ the default values of the Eq 22 
    // in the reference: J. Phys. Chem. A. 2021, 125, 6845.
    // std::cout << "dofe " << dofe << std::endl;
    if (param.getStr("gamma_MM").empty())
        gamma_MM = (sqrt(dofe + 1)  - 1) / dofe;  // if ZPE parameters, γ_MM, is not specified, we take the default value
    else
        gamma_MM = param.getDouble("gamma_MM");
    if (gamma_MM <= (-1.0/DOFe))
        throw std::runtime_error("ERROR: Illegal value for gmmma_MM (requires > -1/DOFe) for eCMMcv method");      
    // std::cout << "gamma_MM" << gamma_MM << std::endl;
    // std::cout << "Das ist eCMMcv Ver Dec 28, 2021" << std::endl;
    traj_debug = 0;
    // std::cout << "inin good" << std::endl;
}


void DynamicsECMMCV::samplingElec(){
    // * Same as full-sphere sampling in DynamicsSPM
    // Note the gamma_MM used here is 1/2 of that used in DynamicsSPM
    std::vector<double>& q = Elec->q;
    std::vector<double>& p = Elec->p;
    std::vector<Complex>& coeff = Elec->coeff;
    std::vector<double>& xcv = Elec->xcv;
    std::vector<double>& pcv = Elec->pcv;
    std::vector<double>& scv = Elec->scv;
    std::normal_distribution<double> normal_dist(0.0, 1.0);
    //normalization factors for spherical sampling
    double NF = 0.0;
    double cv2 = 0.0;
    double xp2 = 0.0;
 
    //1. electronic sampling
    //yes it is exactly the same as ECMM methods
    double sum = 0.0;
    const double rad2 = 2.0 + DOFe*gamma_MM*2.0; // the squared {q,p}-radius

    for (int i = 0; i < DOFe; i++) {
        q[i] = normal_dist(elec_gen);
        p[i] = normal_dist(elec_gen);
        sum += q[i]*q[i] + p[i]*p[i];
    }
    for (int i = 0; i < DOFe; i++) {
        q[i] *= sqrt(rad2/sum);
        p[i] *= sqrt(rad2/sum);
    }
    for (int i = 0; i < DOFe; i++)
        coeff[i] = (q[i] + I*p[i])/sqrt(2.0);
    // debug codes by Cesare
    // std::cout << "p0, p1 = " << p[0] << "," << p[1] << std::endl;
    // std::cout << "q0, q1 = " << q[0] << "," << q[1] << std::endl;

    // test codes by Cesare to make sure the initial electronic DOF focusd

    // given one sampling of ppqq
    /*
    for (int i = 0; i < DOFe; i++) {
        if (i == 0 ) {
            q[i] = sqrt(1+gamma_MM);
            p[i] = sqrt(1+gamma_MM);
        }
        else {
            q[i] = sqrt(gamma_MM);
            p[i] = sqrt(gamma_MM);
        }
        // std::cout << q[i] << ", " << p[i] << ", ";
    }
    */
    /*
     for (int i = 0; i < DOFe; i++) {
        if (i == 0) {
            q[i] = sqrt(1+gamma_MM);
            p[i] = sqrt(1+gamma_MM);
        }
        else {
            q[i] = sqrt(gamma_MM);
            p[i] = sqrt(gamma_MM);
        }
        // std::cout << (q[i]*q[i] + p[i]*q[i]) << ", ";
    } 
    */ 
   // Focus sampling
   /*
   for (int i = 0; i < DOFe; i++) {
       xp2 = 0.5 * (q[i]*q[i] + p[i]*p[i]);
       if (i == 0) {
           NF = sqrt((1.0 +gamma_MM)/xp2);
           q[i] *= NF;
           p[i] *= NF;
       }
       else {
           NF = sqrt((gamma_MM)/xp2);
           q[i] *= NF;
           p[i] *= NF;
       }
   }
   */
    // std::cout << std::endl;


    // 2. commutator variable sampling
    // This section samples commutator variables following 
    // So that the initial nuclear sampling are taken on the donor state
    // or decleared ininital state (which is a pure state)
    // Eq. (41) (42) (43) in the reference J. Phys. Chem. A. 2021, 125, 6845.
    for (int k = 0; k < init_state.size(); ++k){
        for (int i = 0; i < DOFe; i++){
            for (int j = 0; j < DOFe; j++) {
                cv2 = 0.0;
                xp2 = q[i]*q[j] + p[i]*p[j];
                if (i != j) {
                    xcv[i*DOFe+j] = 0;
                    pcv[i*DOFe+j] = 0;
                }
                else{ //i == j condition (diagonal elements of Γ matrix)
                    int n = init_state[k] / DOFe; // initial state i = n*DOFe+m, |n><m|
                    int m = init_state[k] % DOFe; 
                    // So the initial state is |n><n| or |m><m|, which is the |j_occ> in the Eq 
                    // 40, 41, 42 in JPCA, 2021, 125, 6845
                    if (n != m) {
                        // std::cout << "initial state" << init_state[k] << std::endl;
                        throw std::runtime_error("ERROR: coherent state as initial state is not supported by eCMMcv method at present.");
                    }
                    else {  // n == m case 
                        xcv[i*DOFe+j] = normal_dist(elec_gen);
                        pcv[i*DOFe+j] = normal_dist(elec_gen);
                        // xcv[i*DOFe+j] = -q[i];
                        // pcv[i*DOFe+j] = p[i];
                        // std::cout << "xcv = " << xcv[i*DOFe+j] << ", pcv = " << pcv[i*DOFe+j] << ", ";
                        // std::cout << i << ", " << j << ", ";
                        cv2 = xcv[i*DOFe+j]*xcv[i*DOFe+j] + pcv[i*DOFe+j]*pcv[i*DOFe+j];
                        // std::cout << ", cv2 " << cv2 << ", ";
                        if (i == n) {
                            scv[i] = abs(xp2-2.0)/(xp2-2.0);
                            NF = sqrt(abs(xp2-2.0)/cv2);
                        }
                        else {
                            scv[i] = 1.0;
                            // std::cout << "xp2 " << xp2 << ", cv2 " << cv2 << ", ";
                            NF = sqrt(xp2/cv2);
                        }
                        // std::cout << "NF " << NF << std::endl;
                        xcv[i*DOFe+j] *= NF;
                        pcv[i*DOFe+j] *= NF;
                        // try codes by Cesare
                        // make it positive 
                        // xcv[i*DOFe+j] = abs(xcv[i*DOFe+j]);
                        // pcv[i*DOFe+j] = abs(pcv[i*DOFe+j]);
                    }
                }
            }
        }
    }
    // debug code by Cesare
    // std::cout << "NF 4cv = " << NF << std::endl;
    // std::cout << "ini_st = " << init_state[0] << std::endl;
    // std::cout << "int_size " << init_state.size() << std::endl;
    // std::cout << "xcv[0] = " << xcv[0] << ", xcv[1] = " << xcv[1] << std::endl;
    // std::cout << "pcv[0] = " << pcv[0] << ", pcv[1] = " << pcv[1] << std::endl;
    // for (int i = 0; i < DOFe; i++) {
    //    std::cout << " xcv[" << i << "] = " << xcv[i*DOFe+i] << ", x[" << i << "] = " << q[i] << std::endl;
    //    std::cout << " pcv[" << i << "] = " << pcv[i*DOFe+i] << ", p[" << i << "] = " << p[i] << std::endl;
    //    // std::cout << (scv[i]*(pcv[i*DOFe+i]*pcv[i*DOFe+i]+xcv[i*DOFe+i]*xcv[i*DOFe+i]) - (q[i]*q[i]+p[i]*p[i])) << std::endl;
    //    std::cout << "γ = " << (scv[i]*(pcv[i*DOFe+i]*pcv[i*DOFe+i]+xcv[i*DOFe+i]*xcv[i*DOFe+i])) << std::endl;    
    // }
    // for (int i = 0; i < DOFe; i++) {
    //    for (int j = 0; j < DOFe; j++) {
    //        std::cout << (scv[i]*(pcv[i*DOFe+j]*pcv[i*DOFe+j]+xcv[i*DOFe+j]*xcv[i*DOFe+j])) << ", ";
    //    }
    //    std::cout << std::endl;
    // }
    // double checkxp = 0.0;
    // double checkxpcv= 0.0;
    // double gamma_ij = 0.0;
    // for (int i = 0; i < DOFe; i++) {
    //     checkxp += 0.5 * (q[i]*q[i] + p[i]*p[i]);
    //     gamma_ij = 0.0 ;
    //     for (int j = 0; j < DOFe; j++) {
    //        gamma_ij  += 0.5 * (scv[i]*(pcv[i*DOFe+j]*pcv[i*DOFe+j]+xcv[i*DOFe+j]*xcv[i*DOFe+j]));
    //        checkxpcv += 0.5 * (scv[i]*(pcv[i*DOFe+j]*pcv[i*DOFe+j]+xcv[i*DOFe+j]*xcv[i*DOFe+j]));
    //     }
    //     std::cout << "eq 43 check i = " << i << ", γ = " << gamma_ij << 
    //     ", check δ = " << (gamma_ij -  (0.5 * (q[i]*q[i] + p[i]*p[i]))) << std::endl;
    // }
    // std::cout << "initial sampling check " << std::endl;
    // std::cout << "mv " << checkxp << ", cv " << checkxpcv << std::endl;
    // debug codes for initialization of the method 
    // double xxpp = 0.0;
    // for (int i = 0; i < DOFe; i++) {
    //     xxpp += (q[i] * q[i] + p[i] * p[i]);
    // }
    // xxpp *= 0.5; 
    // xxpp -= 1.0;
    // xxpp /= static_cast<double>  (DOFe);
    // double cvc = 0.0; //cv check
    // for (int i = 0; i < DOFe; i++) {
    //     for (int j = 0; j < DOFe; j++) {
    //        cvc += (scv[i] * (xcv[i*DOFe+j]*xcv[i*DOFe+j] + pcv[i*DOFe+j]*pcv[i*DOFe+j]));
    //     }
    // }
    // cvc *= 0.5;
    // cvc /= static_cast<double> (DOFe);
    // std::cout << "qqpp " << xxpp << ", cvc " << cvc << ", gamma " << gamma_MM << std::endl;
    // std::cout << "qqpp " << xxpp << ", gamma " << gamma_MM << std::endl;

    // double cvtemp = 0.0;
    // double xptemp = 0.0;
    // for (int j = 0; j < DOFe; j++) {
    //     cvtemp = 0.0;
    //     xptemp = 0.0;
    //     for (int i = 0; i < DOFe; i++) {
    //          cvtemp += (scv[i] * (xcv[i*DOFe+j]*xcv[i*DOFe+j]+ pcv[i*DOFe+j]*pcv[i*DOFe+j]));
    //      }
    //      xptemp = (q[j]*q[j] + p[j]*p[j]);
    //      std::cout << "ratio " << (xptemp / cvtemp) << std::endl;
    // }
    traj_debug += 1;
    // std::cout << "sampling good" << std::endl;
}

void DynamicsECMMCV::getInitialDensity() {
    // For eCMMcv method, the init_density is the A(0) operator in correlation function.
    // Same as getInitialDensity() in DynamicsECMM
    const std::vector<double>& q = Elec->q;
    const std::vector<double>& p = Elec->p;
    const std::vector<double>& xcv = Elec->xcv;
    const std::vector<double>& pcv = Elec->pcv;
    const std::vector<double>& scv = Elec->scv;

    for (int i = 0; i < init_state.size(); i++) {
        int j = init_state[i] / DOFe; // initial state i = j*DOFe+k, |j><k|
        int k = init_state[i] % DOFe;
        if (j == k)  // population operator
            init_density[i] = 0.5 * (q[j]*q[j] + p[k]*p[k]) - gamma_MM;
        else
            init_density[i] = 0.5 * (q[j] - I*p[j]) * (q[k] + I*p[k]);
    }
    // double gamma_jk = 0.0; 
    // for (int i = 0; i < init_state.size(); i++) {
    //     int j = init_state[i] / DOFe; // initial state i = j*DOFe+k, |j><k|
    //     int k = init_state[i] % DOFe;
    //     for (int n = 0; n < DOFe; n++) {
    //         gamma_jk += 0.5*scv[n]*(xcv[n*DOFe+j]*xcv[n*DOFe+k]+pcv[n*DOFe+j]*pcv[n*DOFe+k]);
    //     }
    // }

    // for (int i = 0; i < init_state.size(); ++i) {
    //    int j = init_state[i] / DOFe; // initial state i = j*DOFe+k, |j><k|
    //    int k = init_state[i] % DOFe;
    //    if (j == k) // population operator
    //        init_density[i] = 0.5 * (q[j]*q[j] + p[k]*p[k]) - gamma_jk;
    //    else        // coherence operator (jk)
    //        init_density[i] = 0.5 * (q[j] - I*p[j]) * (q[k] + I*p[k]) - gamma_jk;
    // }
    // std::cout << "init_den = " << init_density[0] << ", gamma_jk " << gamma_jk << std::endl;
    // for debug only
    // double sum = 0.0; 
    // for (int i = 0; i < DOFe; i++) 
    //     for (int j = 0; j < DOFe; j++)
    //        sum += 0.5 * scv[i] * (xcv[i*DOFe+j]*xcv[i*DOFe+j]+pcv[i*DOFe+j]*pcv[i*DOFe+j])/DOFe;
    // std::cout << "sum over CV " << sum << std::endl;
    // std::cout << "gamma_jk    " << gamma_jk << std::endl;
    // std::cout.precision(16);
    // std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    // const Real_Matrix& H_new = ha->Heff;
    // std::cout << q[0] << ", " << p[0] << ", " << q[1] << ", " << p[1] << ", " << q[2] << ", " << p[2] << ", " 
    //     << H_new[0][0] << ", " << H_new[1][1] << ", " << H_new[2][2] << ", " << H_new[1] << ", " << H_new[1] << ", " << H_new[2] << ", " << H_new[2] << ", " 
    // double initial_potential = 0;
    // commutator variables 
    // double cv = 0.0;
    // Mapping variables 
    // double mv = 0.0;
    // for (int j = 0; j < DOFn; j++) {
    //     for (int n = 0; n < DOFe; n++){
    //         for (int m = 0; m < DOFe; m++) {
    //             cv = 0.0;
    //             mv = 0.5 * (p[n]*p[m] + q[n]*q[m]); 
    //             for (int k = 0; k < DOFe; k++) {
    //                 cv += (0.5 * scv[k] * (xcv[k*DOFe+n]*xcv[k*DOFe+m] + pcv[k*DOFe+n]*pcv[k*DOFe+m]));
    //             }
    //             if (n == m) // from diagnoal of H (population)
    //                 initial_potential += H_new[n][m] * (mv - cv);
    //                 // ha->F[j] += ha->F_all[n*DOFe+m][j]  * (mv - cv);
    //             else if (n < m) // from off-diagnoal of H (coherence)
    //                 // multiply 2 since nm is same as mn (real part)
    //                 // ha->F[j] += ha->F_all[n*DOFe+m][j] * 2.0 * (mv - cv);
    //                 initial_potential += H_new[n][m] * 2.0 * (mv - cv);
    //         }
    //     }
    //     //ha->F[j] /= hbar;
    // }
    // std::cout << ha->R[0] << ", " << ha->V[0] << ", " << initial_potential << std::endl;
    // std::cout << "init den good" << std::endl;
}   

void DynamicsECMMCV::updateDensityMatrix() {
    const std::vector<double>& q = Elec->q;
    const std::vector<double>& p = Elec->p;
    const std::vector<double>& xcv = Elec->xcv;
    const std::vector<double>& pcv = Elec->pcv;
    const std::vector<double>& scv = Elec->scv;
    // The formula can be found in Jian's JPCL paper (eCMM):
    // Ref2: X. He, J. Liu, J Phys. Chem. Lett. 2021, 12, 2496. Eq.20
    // However, they can be used for population-population operator only.
    // The coherence operator can be found in DynamcisSPM.
    // Actually, the formula used in SPM and eCMM are equivalent.
    // γ=0 -> SPM-Q; γ=1 -> SPM-P; γ=1/2 -> SPM-MMST.
    // The population/coherence at t has a scale factor, and for population,
    // there is a shift value. they depend on the DOFe and gamma_MM.
    // And the normalized factor (DOFe) arised from full-sphere sampling will be
    // multiplied in getDensityMatrix() to get final RDMs.
    const double factor = (1.0+DOFe)/(1.0+DOFe*gamma_MM)/(1.0+DOFe*gamma_MM);
    const double shift = (1.0-gamma_MM)/(1.0+DOFe*gamma_MM);
    for (int i = 0; i < init_density.size(); ++i)
        for (int j = 0; j < DOFe; ++j)
            for (int k = 0; k < DOFe; ++k) {
                Complex B = 0 ; // B(q_k, q_j)(t)
                if (j == k)   // population operator
                    B = 0.5 * factor * (q[j]*q[j] + p[k]*p[k]) - shift;
                else // coherence operator (kj) (same as DynamicsSPM)
                    B = 0.5 * factor * (q[k] - I*p[k]) * (q[j] + I*p[j]);
                // Here, RDM is stored as DOFe^2-dimentional vector (index j*DOFe+k).
                // Note that the value of each element is accumaleated (+=).
                Complex sigma = (double) (DOFe) * init_density[i] * B;
                // RDM_current[i][step/RDM_steps][j*DOFe+k] =  sigma;
                // RDM_average[i][step/RDM_steps][j*DOFe+k] += sigma / (double)(ntraj);
                RDM[i][j*DOFe+k] = sigma;
            }
    // Previous BUG codes by Cesare
    // The old version of observables (δnm -> Γnm / γ_MM)
    // double gamma_kj = 0;
    // double Fgamma = 0;
    // double shift = 0;
    // double oneFgamma = 0;
    // The formula can be found in Jian's JPCA paper (eCMMcv):
    // Ref:  X. He, B. Wu, Z. Gong, J. Liu, J. Phys. Chem. A. 2021, 125, 6845.
    // Eq. 19, 20, 40  
    // However, they can be used for population-population operator only.
    // for (int n = 0; n < DOFe; n++) {
    //     oneFgamma += 0.5 * (q[n]*q[n] + p[n]*p[n]);
    //    for (int k = 0; k < DOFe; k++) {
    //         Fgamma += 0.5*scv[k]*(xcv[k*DOFe+n]*xcv[k*DOFe+n]+pcv[k*DOFe+n]*pcv[k*DOFe+n]);
    //     }
    // }
    // std::cout << "Fgamma  " << Fgamma << " one and Fgamma " << oneFgamma << std::endl; 
    // The formula can not be found in Jian's JPCA paper (eCMMcv):
    // Instead the formula used in spin-mapping to compute σ_jk(t).
    // Ref:  J. Chem. Phys. 152, 084110 (2020) Eq. 40-42
    // double factor = (1.0+DOFe)/(1.0+Fgamma)/(1+Fgamma);
    // The coherence operator can be found in DynamcisSPM.
    // Actually, the formula used in SPM and eCMM are equivalent.
    // γ=0 -> SPM-Q; γ=1 -> SPM-P; γ=1/2 -> SPM-MMST.
    // The population/coherence at t has a scale factor, and for population,
    // there is a shift value. they depend on the DOFe and gamma_MM.
    // And the normalized factor (DOFe) arised from full-sphere sampling will be
    // multiplied in getDensityMatrix() to get final RDMs.
    // for (int i = 0; i < init_density.size(); ++i)
    //    for (int j = 0; j < DOFe; ++j)
    //         for (int k = 0; k < DOFe; ++k) {
    //             Complex B = 0 ; // B(q_k, q_j)(t)
    //             gamma_kj = 0 ;
    //             for (int n = 0; n < DOFe; ++n) {
    //                 gamma_kj += scv[n]*(xcv[n*DOFe+j]*xcv[n*DOFe+k]+pcv[n*DOFe+j]*pcv[n*DOFe+k]);
    //             }
    //             // shift = (1.0 - gamma_kj)/(1.0 - DOFe*gamma_kj);
    //             shift = (1.0-gamma_MM)/(1.0+DOFe*gamma_MM);
    //             if (j == k)  { // population operator 
    //                 B = 0.5 * factor * (q[j]*q[j] + p[k]*p[k]) - shift*gamma_kj/gamma_MM;
    //             }
    //             else { // coherence operator (kj) (same as DynamicsSPM)
    //                 B = 0.5 * factor * (q[k] - I*p[k]) * (q[j] + I*p[j])- shift*gamma_kj/gamma_MM;
    //             }
    //             // Here, RDM is stored as DOFe^2-dimentional vector (index j*DOFe+k).
    //             // Note that the value of each element is accumaleated (+=).
    //             rdm[i][step/RDM_steps][j*DOFe+k] += init_density[i] * B;
    //        }
    // for debug only
    // double sum = 0.0; 
    // for (int i = 0; i < DOFe; i++) 
    //    for (int j = 0; j < DOFe; j++)
    //         sum += 0.5 * scv[i] * (xcv[i*DOFe+j]*xcv[i*DOFe+j]+pcv[i*DOFe+j]*pcv[i*DOFe+j])/DOFe;
    // std::cout << "sum over CV " << sum << std::endl;
    // double sum = 0.0;
    // for (int i = 0; i < DOFe; i++) 
    //    sum += 0.5 * (q[i]*q[i]+p[i]*p[i]);
    // sum /= (1.0 + DOFe * gamma_MM * 1.0);
    // std::cout << "sum over MV " << sum << std::endl;
    // std::cout << "traj " << traj_debug << std::endl;
    // if (traj_debug ==  17620) {
    //     std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    //     std::cout.precision(16);
    //     std::cout << ha->R[0] << ", " << ha->V[0] << ", " << ha->F[0] << ", " <<
    //         q[0] << ", " << p[0] << ", " << q[1] << ", " << p[1] << ", " <<
    //         q[2] << ", " << p[2] << std::endl;
    // }
    // if (traj_debug ==  17620) {
    //     std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    //     const Real_Matrix& H_new = ha->Heff;
    //     std::cout.precision(16);
    //     std::cout <<  H_new[0][0] << ", " <<  H_new[1][1] << ", " <<  H_new[2][2] << ", "<< H_new[0][1] << ", " <<  H_new[1][0] << ", " <<  H_new[0][2]   << ", " <<  H_new[2][0]  << ", " <<  H_new[1][2]  << ", " <<  H_new[2][1]  <<  std::endl;
    // }
}

void DynamicsECMMCV::oneStepForDiabaticModel(){

    //This is velocity-Verlet integrator
    std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    // debug codes of Cesare
    // std::cout.precision(16);
    // std::cout << " R " << ha->R[0] << ", V " << ha->V[0] << ", F " << ha->F[0] ;//<< std::endl;
    // Update velocities/momenta with a half step using effective force.
    MOVE_Half_V(ha->F, ha->V);
    // Update positions with a full step
    MOVE_R(ha->V, ha->R);
    // std::cout << "R = " << ha -> R[0];
    // std::cout << "R = " << ha -> R[0] << ", F = " << ha -> F[0] 
    // << ", V = " << ha -> V[0] << ", " << std::endl;
    // Update H and F_all since R has been changed.
    // std::cout << ", pot " << (ha->H[0][0] + ha-> H[1][1] + ha -> H[2][2]) << std::endl;
    updateH();
    // Update electronic mapping variables.
    // std::cout << ", F = " << ha->F[0] << std:: endl;
    MOVE_elec();
    // Update effective forces
    updateF();
    // Update velocities/momenta with a half step using new effective forces
    MOVE_Half_V(ha->F, ha->V);
    // debug code by Cesare
    // std::cout << " R " << ha->R[0] << ", V " << ha->V[0] << ", F " << ha->F[0] << std::endl;
    // const std::vector<double>& q = Elec->q;
    // const std::vector<double>& p = Elec->p;
    // double xxpp = 0.0;
    // for (int i = 0; i < DOFe; i++) {
    //     xxpp += (q[i] * q[i] + p[i] * p[i]);
    // }
    // xxpp *= 0.5; 
    // xxpp -= 1.0;
    // xxpp /= static_cast<double>  (DOFe);
    // std::cout <<"ppqq " << xxpp << std::endl;
    // CoutCode by Cesare
    // debug nuclear seed 42 
    // debug elec seed 73
    // dt = 1e-5ps
    // nsteps=8000


}

void DynamicsECMMCV::updateDiabaticModelForces() {
    std::vector<double>& p = Elec->p;
    std::vector<double>& q = Elec->q;
    std::vector<double>& xcv = Elec->xcv;
    std::vector<double>& pcv = Elec->pcv;
    std::vector<double>& scv = Elec->scv;
    // This section implemented the Eq. (37) of JPCA. 2021. 125, 6845
    std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    // Reset all elements to 0 with keeping same size
    std::fill(ha->F.begin(), ha->F.end(), 0.0);
    // debug codes by Cesare
    // std:: cout << "Dynamics ECMMCV DOFn check " << DOFn << std::endl;
    // commutator variables 
    double cv = 0.0;
    // Mapping variables 
    double mv = 0.0;
    for (int j = 0; j < DOFn; j++) {
        for (int n = 0; n < DOFe; n++){
            for (int m = 0; m < DOFe; m++) {
                cv = 0.0;
                mv = 0.5 * (p[n]*p[m] + q[n]*q[m]); 
                for (int k = 0; k < DOFe; k++) {
                    cv += (0.5 * scv[k] * (xcv[k*DOFe+n]*xcv[k*DOFe+m] + pcv[k*DOFe+n]*pcv[k*DOFe+m]));
                }
                if (n == m) // from diagnoal of H (population)
                    ha->F[j] += ha->F_all[n*DOFe+m][j]  * (mv - cv);
                else if (n < m) // from off-diagnoal of H (coherence)
                    // multiply 2 since nm is same as mn (real part)
                    // If Condon approximation is used (coupling is contant), it is 0.
                    ha->F[j] += ha->F_all[n*DOFe+m][j] * 2.0 * (mv - cv);
            }
        }
        ha->F[j] /= hbar;
        // ha->F[j] += ha->F_avg[j];
    }
}

Complex DynamicsECMMCV::getOperator(int j, int k) {
    const std::vector<double>& q = Elec->q;
    const std::vector<double>& p = Elec->p;
    if (j == k)   // population operator
        return 0.5 * (q[j]*q[j] + p[k]*p[k]) - gamma_MM;
    else
        return 0.5 * (q[j]-I*p[j]) * (q[k]+I*p[k]);
    // const std::vector<double>& xcv = Elec->xcv;
    // const std::vector<double>& pcv = Elec->pcv;
    // const std::vector<double>& scv = Elec->scv;
    // Complex gamma = 0.0;
    // for (int i = 0; i < DOFe; i++) {
    //     gamma += 0.5 * scv[i] * (xcv[i*DOFe+j]*xcv[i*DOFe+k]+pcv[i*DOFe+j]*pcv[i*DOFe+k]); 
    // }
    //  if (j == k)   // population operator   
    //     return 0.5 * (q[j]*q[j] + p[k]*p[k]) - gamma;
    // else        // coherence operator 
    //     return 0.5 * (q[j]-I*p[j]) * (q[k]+I*p[k]) - gamma;
}

void DynamicsECMMCV::MOVE_elec_eCMMcv(){
    std::vector<double>& p = Elec->p;
    std::vector<double>& q = Elec->q;
    std::vector<double>& xcv = Elec->xcv;
    std::vector<double>& pcv = Elec->pcv;
    std::vector<double>& scv = Elec->scv;
    // const Real_Matrix& H_old = Ha->Heff_old;
    const Real_Matrix& H_new = Ha->Heff;
    // For all CMM family methods, electronic steps per nuclear step is 1
    const double dt = DT;
    // debug codes by Cesare
    // std::cout << "dt  "  << dt << std::endl;
    // electronic sampling variables updated according to Eq (38) 
    // in J. Phys. Chem. A. 2021, 125, 6845.
    for (int n = 0; n < DOFe; n++) {
        for (int m = 0; m < DOFe; m++) {
            q[n] += ( dt * (H_new[n][m] * p[m]));
            p[n] += (-dt * (H_new[n][m] * q[m]));
        }
    }
    // commutator variables update in Eq (38)
    //  J. Phys. Chem. A. 2021, 125, 6845.
    for (int k = 0; k < DOFe; k++) {
        for (int n = 0; n < DOFe; n++) {
            for (int m = 0; m < DOFe; m++) {
                xcv[k*DOFe+n] += (-scv[k] * dt * (H_new[n][m] * pcv[k*DOFe+m]));
                pcv[k*DOFe+n] += ( scv[k] * dt * (H_new[n][m] * xcv[k*DOFe+m]));
            }
        }
    }
}


