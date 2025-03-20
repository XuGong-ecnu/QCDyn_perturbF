/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                       *
 * Last updated: Jan. 6, 2022                                                 *
 * -------------------------------------------------------------------------- */

#include "ObservableTCF.h"

void ObservableTCF::init() {
    LEN_TRAJ  = param.getInt("LEN_TRAJ");
    LEN_CORR  = param.getInt("LEN_CORR");
    DT    = param.getDouble("DT");
    LEN_TCF  = param.getInt("LEN_TCF");
    LEN_SKIP = param.getInt("LEN_SKIP");
    nsteps  = param.getInt("nsteps");
    steps_for_config = param.getInt("steps_for_config");
    TCF.resize(9*LEN_TRAJ);
    TCF_accum.resize(9*LEN_TRAJ);
    TCF_fly.resize(9*LEN_TRAJ);
    TCF_fly_accum.resize(9*LEN_TRAJ);
    A.resize(LEN_CORR);
    A_fly.resize(LEN_CORR);
    B.resize(LEN_CORR);
    B_fly.resize(LEN_CORR);
    Avec.resize(LEN_CORR);
    Avec_fly.resize(LEN_CORR);
    Bvec.resize(LEN_CORR);
    Bvec_fly.resize(LEN_CORR);
}

void ObservableTCF::setEachStepPiTensor(int step, std::vector<double> PiCurrentTensor) {
        for (int a = 0; a < 9; a++){
     
    //    TCF[step * 9 - 9 + a] = PiCurrentTensor[a];
        int index = step * 9 + a - 9;
        if(index>=9*LEN_TRAJ){break;}
        TCF[index] = PiCurrentTensor[a];
    }
}

void ObservableTCF::setEachPositiveStepPiTensor(int step, std::vector<double> PiCurrentTensor) {
    for (int a = 0; a < LEN_TCF; a++)
        TCF[step * 9 + a] = PiCurrentTensor[a];// TCF means positive perturb
}

void ObservableTCF::setEachNegtiveStepPiTensor(int step, std::vector<double> PiCurrentTensor) {
    for (int a = 0; a < LEN_TCF; a++)
        TCF_fly[step * 9 + a] = PiCurrentTensor[a];// TCF_fly means negtive perturb
                
    // This is for Rxzxz result. xz = 0 * 3 + 2; 
    std::vector<int> perturb_parameters;
    SplitString(perturb_parameters, param.getStr("perturb_parameters"));
    int mu, nu;
    if (perturb_parameters[2] == 0) {
        mu = 1;
        nu = 2;
    }
    else if (perturb_parameters[3] == 0) {
        mu = 0;
        nu = 2;
    }
    else if (perturb_parameters[4] == 0) {
        mu = 0;
        nu = 1;
    }
    double munu = 3 * mu + nu;
    
    B[step] = TCF[step * 9 + munu] - TCF_fly[step * 9 + munu];

    if (step == LEN_CORR) {
        for (int i = 0; i < LEN_CORR; i++) {
            A[i] += B[step];
        }
    }
}

void ObservableTCF::writeEqRxzxzResult(std::string& file) {
    // Save Rxzxz to file.
    FILE* outfile = CheckFile(file);
    // calculate the Rxzxz result
    Rxzxz();
    // save the Rxzxz response function out
    for (int i = 0; i < LEN_CORR; i++) {
        fprintf(outfile, "%.12g     %.12g", i * DT * LEN_SKIP, A[i]);
        fprintf(outfile, "\n");
    }
    fclose(outfile);
}

void ObservableTCF::writeNoneqRxzxzResult(std::string& file) {
    // Save Rxzxz to file.
    FILE* outfile = CheckFile(file);
    // calculate the Rxzxz result
    double count = nsteps/steps_for_config;
    for(int i = 0; i < LEN_CORR; i++) {
        A[i] /= count;
        A[i] *= 1;// This is the unit for the noneq part
    }
    // save the Rxzxz response function out
    for (int i = 0; i < LEN_CORR; i++) {
        fprintf(outfile, "%.12g     %.12g", i * DT * LEN_SKIP, A[i]);
        fprintf(outfile, "\n");
    }
    fclose(outfile);
}

void ObservableTCF::Rxzxz() { 
//calculate OKE anisotropic/depolarized response fn Rxzxz(t): t=0, 1, ... tcor. Where trun is actually input polarizability tensors. 
//tgap=1, with time unit = tgap * STEPS(=10) dt
    double norm[LEN_CORR/LEN_SKIP];
    double A0[9], Att0[9];
    double h;
    int n = LEN_CORR/LEN_SKIP;
    for (int t = 0; t < LEN_CORR; ++t){
        A[t] = 0;
        norm[t] = 0;
    }
    for (int t0 = 0; t0 < LEN_TRAJ; t0++) {
        for (int i = 0; i < 9; ++i) A0[i] = TCF[t0 * 9 + i];    
        int tt0max = min(LEN_TRAJ, t0 + LEN_CORR);      
        for (int tt0 = t0; tt0 < tt0max; tt0 += LEN_SKIP) {
            for (int i = 0; i < 9; ++i) Att0[i] = TCF[tt0 * 9 + i];
            int t = (tt0 - t0) / LEN_SKIP;
            A[t] += (1/10.0)* Pairwise_product(A0, Att0) - (1/30.0) * Trace(Att0)*Trace(A0);
            norm[t]++;
        }
    }
    for (int t = 0; t < n; t++){
        norm[t] = A[t]/norm[t]; //norm[tcor/tgap] for temporary storage for Corr[tcor/tgap]
    }
    // apply -d/dt Corr[t], notice: need to times beta=1/kT
    h = 2 * LEN_SKIP * DT; // h=2dt
    A[0] = -(-3 * norm[0] + 4 * norm[1] - norm[2]) / h;
    A[n-1] = -(norm[n-3] - 4 * norm[n-2] + 3 * norm[n-1]) / h;
    for( int t = 1; t < n-1; t++) {
        A[t] = -(norm[t+1] - norm[t-1])/h;
    }
}

double ObservableTCF::Pairwise_product(double AA[9], double BB[9]) {
    double pp(0);
    int i;
    for (i=0; i<9; i++)
        pp += AA[i]*BB[i];
    return pp;
}

double ObservableTCF::Trace(double Aa[9]) {
    double tr(0);
    tr += Aa[0]+Aa[4]+Aa[8];
    return tr;
}

int ObservableTCF::min(int a, int b) {
    if (a <= b) {
        return a;
    }
    else {
        return b;
    }
}
