/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: Apr. 11, 2022                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "vector"
#include "Tools.h"
#include "Parameters.h"

class ObservableTCF {
    
public:
    /**
     * Construct an object of ObservableTCF
     *
     * @param param   the global paramters
     * @param Ha      Hamiltonian object
     * @param Dy      Dynamics object
     */
    ObservableTCF(Parameters& param) : param(param) {}
    ~ObservableTCF() {}
        /** 
     * In order to construct the new operator to run solve the resize problem
     */
//    ObservableTCF(const ObservableTCF& observableTCF): param(observableTCF.param) {}
//    ObservableTCF& operator=(const ObservableTCF& other)
//    {
//        param = other.param;
//        return *this;
//    }


ObservableTCF(const ObservableTCF& other) : param(other.param) {
        TCF = other.TCF;
        TCF_accum = other.TCF_accum;
        TCF_fly = other.TCF_fly;
        TCF_fly_accum = other.TCF_fly_accum;
        A = other.A;
        A_fly = other.A_fly;
        B = other.B;
        B_fly = other.B_fly;
        Avec = other.Avec;
        Avec_fly = other.Avec_fly;
        Bvec = other.Bvec;
        Bvec_fly = other.Bvec_fly;
        
    }
    
    ObservableTCF& operator=(const ObservableTCF& other) {
        if (this != &other) { 
            
            TCF = other.TCF;
            TCF_accum = other.TCF_accum;
            TCF_fly = other.TCF_fly;
            TCF_fly_accum = other.TCF_fly_accum;
            A = other.A;
            A_fly = other.A_fly;
            B = other.B;
            B_fly = other.B_fly;
            Avec = other.Avec;
            Avec_fly = other.Avec_fly;
            Bvec = other.Bvec;
            Bvec_fly = other.Bvec_fly;
            
        }
        return *this;
    }

    
    /**
     * Initialize data members and check the legality of parameters.
     */
    void init();
    /**
     * Set one step Pi tensor to the TCF[].
     */
    void setEachStepPiTensor(int step, std::vector<double> PiCurrentTensor);
    /**
     * Set one step Pi tensor to the TCF[] in positive perturb MD.
     */
    void setEachPositiveStepPiTensor(int step, std::vector<double> PiCurrentTensor);
    /**
     * Set one step Pi tensor to the TCF[] in negtive perturb MD.
     */
    void setEachNegtiveStepPiTensor(int step, std::vector<double> PiCurrentTensor);
    /**
     * Uing the eq method to get the Rxzxz result file.
     */
    void writeEqRxzxzResult(std::string& file);
    /**
     * Uing the noeq method to get the Rxzxz result file.
     */
    void writeNoneqRxzxzResult(std::string& file);
    /**
     * This is the eq method to calculate the Rxzxz.
     */
    void Rxzxz();
protected:
    /**
     * Calculate the matrix.
     */
    double Pairwise_product(double AA[9], double BB[9]);
    /**
     * Get the trace of the matrix.
     */
    double Trace(double Aa[9]);
    /**
     * Judge the min parameter.
     */
    int min(int a, int b);
    
private:
    Parameters& param;
    int ntraj;
    int LEN_TRAJ;
    int LEN_CORR;
    int LEN_TCF;
    int LEN_SKIP;
    int steps_for_config;
    int nsteps;
    double DT;
    std::vector<double> TCF;
    std::vector<double> TCF_accum;
    std::vector<double> TCF_fly;
    std::vector<double> TCF_fly_accum;
    std::vector<double> A;
    std::vector<double> A_fly;
    std::vector<double> B;
    std::vector<double> B_fly;
    std::vector<double> Avec;
    std::vector<double> Avec_fly;
    std::vector<double> Bvec;
    std::vector<double> Bvec_fly;
};
