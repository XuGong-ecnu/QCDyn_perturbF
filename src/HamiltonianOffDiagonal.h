/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: Jan. 14, 2021                                                *
 * -------------------------------------------------------------------------- */
#pragma once
#include "Tools.h"
#include "Parameters.h"
#include "Topology.h"
#include "Vec3.h"

class HamiltonianOffDiagonal {

public: 
    /**
     * Construct a ForceFieldBase object.
     *
     * @param param   the global paramters
     */    
    HamiltonianOffDiagonal(Parameters& param) : param(param) {}
    ~HamiltonianOffDiagonal() {}
    /**
     * Make preparations for calculating the forces and energy.
     * 1. Initialize force_contribution to control the calculation of force part 
     * 2. Initialize data members such as R, V, F, forceParts and energyParts.
     */
    void init();
    /**
     * Fill the each part with thenewvdw parameters from topology.
     * @param topology     the Topology object loaded from topology file
     * @param param        the global simulation parameters
     */
    void setNewVdwParameters();
    double getDiabaticCoupling(int i, int j);
    std::vector<Vec3> getDiabaticCouplingforce(int i, int j);
    /**
     * Get the dR in PBC system
     */ 
    void GetDeltaRPeriodic(const Vec3& atomCoordinatesI, const Vec3& atomCoordinatesJ,
                        const Vec3 (&periodicBoxVectors)[3], double deltaR[5]);
    /**
     * set the posiation from ForceFieldBase 
     */    
    void setPositions(const std::vector<Vec3>& positions);
    /**
     * Set the default values of the vectors defining the axes of the periodic 
     * box (measured in nm).  Any newly created system will have its box vectors 
     * set to these.  They will affect any Force added to the System that uses 
     * periodic boundary conditions.
     *
     * Triclinic boxes are supported, but the vectors must satisfy certain 
     * requirements.  In particular, a must point in the x direction, b must 
     * point "mostly" in the y direction, and c must point "mostly"
     * in the z direction.  See the documentation for details.
     *
     * @param a      the vector defining the first edge of the periodic box
     * @param b      the vector defining the second edge of the periodic box
     * @param c      the vector defining the third edge of the periodic box
     */
    void setPeriodicBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c);
    /**
     * set the Topologies from ForceFieldBase
     */
    void setTopologies(const Topology& fftopologies);
    /** 
     * This is to get the atom or molecule polarizability
     */
    struct FFNewVDW{ int atomElec; double vdwparam0; double vdwparam1; double vdwparam2; double vdwparam3; 
                    double vdwparam4; double vdwparam5; double vdwparam6; double vdwparam7; double vdwparam8;};    

protected:
    // Parameters object controls the simulation
    Parameters& param;
    // force field parameters of the system
    Topology topologies; 
    std::vector<FFNewVDW>         ffnewvdw;
    // box size
    Vec3 periodicBoxVectors[3];
    // this is the total atom in the sysytem
    int DOFn;
    int DOFe;
    std::vector<Vec3> R;
};