/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: June. 10, 2022                                                *
 * -------------------------------------------------------------------------- */
#pragma once
#include "Tools.h"
#include "Parameters.h"
#include "StructureVec3.h"
#include "Vec3.h"
#include "ForceFieldNeighborList.h"
#include "ForceFieldBase.h"
//#include "gamma.h"
#include "pme.h" 
#include "algorithm"

class ForceFieldPolar : public ForceFieldBase{

public: 
    /**
     * Construct a ForceFieldPolar object.
     *
     * @param param   the global paramters
     */    
    ForceFieldPolar(Parameters& param) : ForceFieldBase(param) {}
    ~ForceFieldPolar() {}
    /**
     * Make preparations for calculating the forces and energy.
     * 1. Initialize force_contribution to control the calculation of force part 
     * 2. Initialize data members such as R, V, F, forceParts and energyParts.
     */
    void init();
    /**
     * Fill the polar part with the force field parameters from topology.
     */
    void setForceFieldPolarParameters();
    /**
     * Calculete the polarizability tensor
     */
    void calculatePiTensor();
    /**
     * Calculete the dipole.R left and right direction
     */
    void calculateDipoleRleftAright();
    /**
     * Calculete the induced dipole tensor
     */
    void calculateInducedDipole();
    /**
     * Compute and return potential energy of given state. If includeForces = true,
     * the forces will be updated, and you can use getForces() to retrieve the 
     * forces that were calculated after calling this.
     * @param includeForces  true if forces should be calculated
     * @return               the potential energy
     */
    double getPotentialEnergy(bool includeForces); 
    /**
     * Compute and return polar force field energy of given state. If includeForces=true,
     * the forces will be updated, and you can use getForces() to
     * retrieve the forces that were calculated after calling this.
     * @param includeForces  true if forces should be calculated
     * @return               the Polar energy
     */
    double getPolarPotentialEnergy(bool includeForces);
    /**
     * We should divide the polarizable force and energy into Two different methods.
     * First is the point dipole model for polarizable energy calculation Ref[1]
     * [1]. Sala J, Guardia E, Masia M. The polarizable point dipoles method with electrostatic 
     * damping: Implementation on a model system[J]. The Journal of chemical physics,
     * 2010, 133(23): 234101.
     * Using the Ewald method.
     * The dammping has the ES screening functions.
     * Compute the Polar force and energy cosmponents
     * @param forces           total force that to compute
     * @param polarforce        angleforce that to compute
     * @param polarenergy      angleenergy that to compute
     * @param includeForces    true if forces should be calculated
     */
    // TO DO
    //void calculateESPolarForceAndEnergy(std::vector<Vec3>& forces, std::vector<Vec3>& polarforce, 
    //                            double& polarenergy, bool includeForces); 
    /**
     *  Second is the  for polarizable energy calculation Ref[2]
     * [2]. Toukmaji A, Sagui C, Board J, et al. Efficient particle-mesh Ewald based approach to 
     * fixed and induced dipolar interactions[J]. The Journal of chemical physics, 2000,
     *  113(24): 10913-10927. 
     * We will compare the energy and force with the Tinker software using the Amoeba force field.
     * Focus: Ref: Ponder, J. W.; Wu, C.; Ren, P.; Pande, V. S.; Chodera, J. D.; Schnieders, M. J.; 
     * Haque, I.; Mobley, D. L.; Lambrecht, D. S.; DiStasio, R. A., Jr.; Head-Gordon, M.; Clark, 
     * G. N. I.; Johnson, M. E.; Head-Gordon, T. Current Status of the AMOEBA Polarizable Force Field. 
     * J. Phys. Chem. B 2010, 114, 2549−2564.
     * In the AMOEBA polarization model, the damping factor provides another control over the ability 
     *  of an atom to polarize; the universal damping factor adopted by AMOEBA is a = 0.39, which 
     * effectively leads to a stronger damping and less shortrange polarization than the original 
     *  value of 0.572 suggested by Thole. We have kept the same atomic polarizabilities (Å3) given by 
     * Thole, eg: 1.334 for carbon, 0.496 for hydrogen, 1.073 for nitrogen, and 0.837 for oxygen. 
     * The only exception is for carbon and hydrogen in aromatic rings, where we found the use of somewhat 
     * larger values greatly improves the molecular polarizability tensor of benzene and polycyclic aromatics.
     * Compute the Polar force and energy cosmponents
     * @param forces           total force that to compute
     * @param polarforce        angleforce that to compute
     * @param polarenergy      angleenergy that to compute
     * @param includeForces    true if forces should be calculated
     */
    void calculateTTM3PolarForceandEnergy(std::vector<Vec3>& forces, std::vector<Vec3>& polarforce, 
                                double& totalenergy, double& polarenergy, bool includeForces);   
    /**
     *  Third is the  for polarizable energy calculation Ref[3]
     * The dammping has the TTM-4 model damping functions.
     * [3]. Ren P, Ponder J W. Polarizable atomic multipole water model for molecular 
     * mechanics simulation[J].The Journal of Physical Chemistry B, 2003, 107(24): 5933-5947.
     * The dammping has the TTM-3 model damping functions.
     * We will compare the energy and force with the MBX software.
     * Compute the Polar force and energy cosmponents
     * @param forces           total force that to compute
     * @param polarforce        angleforce that to compute
     * @param polarenergy      angleenergy that to compute
     * @param includeForces    true if forces should be calculated
     */
    void calculateTTM4PolarForceandEnergy(std::vector<Vec3>& forces, std::vector<Vec3>& polarforce, 
                                double& totalenergy, double& polarenergy, bool includeForces);  
    
    /**
     * Get the polarizable tensor from the ForceFieldPolar
     */
    void getPiTensor(std::vector<double>& Pi_tensor); 
    /**
     * Get the dipole tensor from the ForceFieldPolar
     */
    void getDipoleRleftAright(double& DEScosleft, double& DEScosright, int& DESnumLeft, 
                            int& DESnumRight, double& DESdirLx, double& DESdirRx);    
    /**
     * Get the Induced dipole tensor from the ForceFieldPolar
     */
    void getMuTensor(Vec3& Mu_tensor);  
    /**
     * This is the most simplest model for the interation
     * s[1] = 1;
     * s[1] = 1;
     */
    void TTM0Model(double& S0r, double& S1r, double& S2r, double& S3r);
    /**
     * This is the piecewise function model
     * 
     */
    void TTM1Model(double& S0r, double& S1r, double& S2r, double& S3r);
    /**
     * This is the Gaussian ES damping model for the polarizable force field in the Ref[1]
     * s[0] = erfc(alpha_Ew * r),
     * s[1] = s[0] + 2 * alpha_Ew * r * exp(-(alpha_Ew * r)^2)/sqrt(π),
     * s[2] = s[1] + 4 * (alpha_Ew * r)^3 * exp(-(alpha_Ew * r)^2)/(3 * sqrt(π)),
     * s[3] = s[2] + 8 * (alpha_Ew * r)^5 * exp(-(alpha_Ew * r)^2)/(15 * sqrt(π)),
     * need modified
     */ 
    void TTM2Model(double alphaEwald, double deltaR, double& S0r, 
                        double& S1r, double& S2r, double& S3r);
    /**
     * This is the Exponential damping model(TTM3) for the polarizable force field, 
     * s[0] = 1 - exp(-a * (r/A)^3) + a^(1/3) * (r/A) * γ[2/3, a * (r/A)^3],
     * s[1] = 1 - exp(-a * (r/A)^3),
     * s[2] = 1 - [1 + a * (r/A)^3] * exp(-a * (r/A)^3),
     * s[3] = 1 - [1 + a * (r/A)^3 + 3 * a^2 * (r/A)^6 / 5] * exp(-a * (r/A)^3),
     * where A = (alpha_Ew * alpha_Ew)^1/6
     * a is Thole damping parameter for chosen pair of atoms
     * a is a scaling parameter equal to 0.39 for all of the simulated systems 
     * but CS2, where it is set to 0.77.
     * This isn't surport the CS2 sysytem.
     */ 
    void TTM3Model(double dampingParam, double alphaij, double& S0r, double& S1r, double& S2r, double& S3r);   
    /**
     * This is the order 4 exponential damping (TTM-4) model for the polarizable force field
     * For MB-pol model, the order 4 exponential damping (TTM-4) is used
     * s[0] = 1 - exp[-a * (r/A)^4] + a^1/4 * (r/A) * γ[3/4, a * (r/A)^4],
     * s[1] = 1 - exp[-a * (r/A)^4],
     * s[2] = s[1] - 4 * a * (r/A)^4 * exp[-a * (r/A)^4]/3,
     * s[3] = s[2] - 4 * a(r/A)^4 * exp[-a * (r/A)^4] * (4 * a * (r/A)^4 - 1),
     * where A = (alpha_Ew * alpha_Ew)^1/6
     * a is Thole damping parameter for chosen pair of atoms
     * If the system don' t have water molecule, so the bond 1-2 and 1-3 pair atoms 
     * will chose dampingParam2 = 0.3, and so the bond 1-4 and no-boned pair atoms 
     * will chose dampingParam3 = 0.055, else if the system have water molecule, 
     * bond 1-2 will chose dampingParam1 = 0.626, and so the bond 1-4 and no-boned 
     * pair atoms will chose dampingParam3 = 0.055.
     * Only surport the no water sysytem.
     */ 
    void TTM4Model(double dampingParam, double alphaij, double& S0r, double& S1r, double& S2r, double& S3r);  
    /**
     * This is the Gaussian damping model for the polarizable force field Ref[1].APPENDIX: 
     * SCREENING FUNCTIONS
     * s[0] = erf(r/a),
     * s[1] = s[0] - (2/sqrt(π))* (r/a) * exp(-(r/a)^2),
     * s[2] = s[1] - (4/3 *sqrt(π))* (r/a)^3 * exp(-(r/a)^2),
     * s[3] = s[2] - (8/15 *sqrt(π))* (r/a)^5 * exp(-(r/a)^2).
     * need modified
     */      
    void Gaussian_model(double alphaij, double& S0r, double& S1r, double& S2r, double& S3r);
    /**
     * This is the judge the Thole model pair.
     */ 
    void addTholeExclusionsToSet(std::vector<std::set<int> >& bonded12, std::vector<std::set<int> >& bonded13);
    void calcualtePerturbPolarForce(std::vector<Vec3>& forces);
    void calcualtePerturbMuForce(std::vector<Vec3>& forces);
    void setTrueMdPerturb();
    void setFalseMdPerturb();
    void setPositiveDirection();
    void setNegativeDirection();
    void Pi_tensor(std::vector<double>& iPi_all, std::vector<double>& iPi);
    void Mu_tensor(std::vector<double>& iMu_all, std::vector<double>& iMu);
    void Dk_Pi_tensor(int k, std::vector<double>& iPi_all, std::vector<double>& Tij_all, int conf_diff, std::vector<double>& DkPi);
    void Dk_Mu_tensor(int k, std::vector<double>& iMu_all, std::vector<double>& Tij_all, int conf_diff, std::vector<double>& DkMu);
    void Dk_Tij_tensor(int k, int i, int j, double DkTij[27]);
    void Perturb_Polar(std::vector<double>& DPi, std::vector<Vec3>& perturbforce);
    void Perturb_Dipole(std::vector<double>& DMu, std::vector<Vec3>& perturbforce);
    void Dk_E_perm(int k, std::vector<double>& DkE_perm, std::vector<double>& DkE);
protected:
    /** 
     * This is to get the atom or molecule polarizability
     */
    struct FFPolars{int atomIndex; double alpha;};     
    /**
     * For thole damping a pair for 12 
     */
    struct Exclusions12{int atomIndex[2];};
    /**
     * For thole damping a pair for 13
     */
    struct Exclusions13{int atomIndex[2];}; 
      
private:
    // force field parameters of the system
    std::vector<FFPolars>           ffPolars;
    std::vector<Exclusions12>   exclusions12;
    std::vector<Exclusions13>   exclusions13;
    // The damping_type : TTM0, TTM1, TTM2, TTM3 and TTM4;
    std::string damping_type;
    // This is the the number of iterations
    int iter_max;
    // This is the TOLERANCE of the iterations
    double TOLERANCE;
    // The induced dipole of all atoms
    std::vector<double> Mu_all;
    Vec3 Mu;
    //The polarizability tensor of all atoms
    std::vector<double> Pi_all;
    double Pi[3][3];
    //std::vector<std::vector<double>> Pi;
    std::vector<double> Dipole_mol;
    // This is for the DES dipole . R calcualtion part.
    double cosleft;
    double cosright;
    double Rdipoleleft;
    double Rdipoleright;
    double dirLx;
    double dirRx;
    Vec3 dipole_right;
    Vec3 dipole_left;
    int numRight;
    int numLeft;
    std::vector<double> Perturbforce;
    // Control the force field polar and nonpolar
    std::string forcefield_type;
    // choose to calculate the perturb force;
    bool md_perturb;
    // choose the direction of the perturb
    bool direction;

};
