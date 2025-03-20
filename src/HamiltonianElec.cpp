/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Mar. 02, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "HamiltonianElec.h"

void HamiltonianElec::init() {
    DOFe = param.getInt("DOFe");
    // Mapping variable position and momenta (LSC/SQC method)
    q.resize(DOFe, 0);
    p.resize(DOFe, 0);
    // Mapping variables (JianLiu's CMM method)
    x.resize(DOFe, 0);
    y.resize(DOFe, 0);
    px.resize(DOFe, 0);
    py.resize(DOFe, 0);
    // Complex coefficient of wavefunction, coeff[j] = (q[j] + I * p[j])/sqrt(2)
    // It is used drectly in Ehrenfest/surface hopping methods.
    coeff.resize(DOFe, 0);
    // commutator variables for eCMMcv method
    scv.resize(DOFe, 0);
    xcv.resize(DOFe*DOFe, 0);
    pcv.resize(DOFe*DOFe, 0);
}