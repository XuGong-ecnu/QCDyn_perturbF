#ifndef OPENMM_VERLETINTEGRATOR_H_
#define OPENMM_VERLETINTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2020 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "Integrator.h"
#include "openmm/Kernel.h"
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This is an Integrator which simulates a System using the leap-frog Verlet algorithm.
 */

class OPENMM_EXPORT VerletIntegrator : public Integrator {
public:
    /**
     * Create a VerletIntegrator.
     *
     * @param stepSize the step size with which to integrate the system (in picoseconds)
     */
    explicit VerletIntegrator(double stepSize);
   /**
     * Advance a simulation through time by taking a series of time steps.
     *
     * @param steps   the number of time steps to take
     */
    void step(int steps);
    /**
     * Advance a simulation through time by taking one step with the providing
     * external forces instead of the forces which is computed from force fileds
     * based on current positions.
     *
     * Note that, you should call updateContextState() firstly before calling
     * this function. (Added by zhubin, Nov. 14, 2021)
     *
     * @param forces the vector contains forces
     */
    void oneStep(const std::vector<Vec3>& forces);

protected:
    /**
     * This will be called by the Context when it is created.  It informs the Integrator
     * of what context it will be integrating, and gives it a chance to do any necessary initialization.
     * It will also get called again if the application calls reinitialize() on the Context.
     */
    void initialize(ContextImpl& context);
    /**
     * This will be called by the Context when it is destroyed to let the Integrator do any necessary
     * cleanup.  It will also get called again if the application calls reinitialize() on the Context.
     */
    void cleanup();
    /**
     * Get the names of all Kernels used by this Integrator.
     */
    std::vector<std::string> getKernelNames();
    /**
     * Compute the kinetic energy of the system at the current time.
     */
    double computeKineticEnergy();
    /**
     * Get the time interval by which velocities are offset from positions.  This is used to
     * adjust velocities when setVelocitiesToTemperature() is called on a Context.
     */
    double getVelocityTimeOffset() const {
        return getStepSize()/2;
    }
private:
    Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_VERLETINTEGRATOR_H_*/
