/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 13, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "ObservableBase.h"

void ObservableBase::init() {
    // * Get parameters from Hamilatonian object, the legality of them should
    // be checked by Hamilatonian.
    DOFe     = param.getDouble("DOFe");
    dyn_type = param.getStr("dyn_type");
    ntraj    = param.getInt("ntraj");
    nsteps   = param.getInt("nsteps");
    step     = 0;
    DeltaT   = 0; // initialize it in subclass
    nucl_start = param.getInt("nucl_start");
    nucl_skip  = param.getInt("nucl_skip");
    nucl_end   = param.getInt("nucl_end");
}