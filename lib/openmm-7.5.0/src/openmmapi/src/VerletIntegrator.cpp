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

#include "openmm/VerletIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <string>

using namespace OpenMM;
using std::string;
using std::vector;

VerletIntegrator::VerletIntegrator(double stepSize) {
    setStepSize(stepSize);
    setConstraintTolerance(1e-5);
}

void VerletIntegrator::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    context = &contextRef;
    owner = &contextRef.getOwner();
    kernel = context->getPlatform().createKernel(IntegrateVerletStepKernel::Name(), contextRef);
    kernel.getAs<IntegrateVerletStepKernel>().initialize(contextRef.getSystem(), *this);
}

void VerletIntegrator::cleanup() {
    kernel = Kernel();
}

vector<string> VerletIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateVerletStepKernel::Name());
    return names;
}

double VerletIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateVerletStepKernel>().computeKineticEnergy(*context, *this);
}

void VerletIntegrator::step(int steps) {
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");
    for (int i = 0; i < steps; ++i) {
        context->updateContextState();
        context->calcForcesAndEnergy(true, false, getIntegrationForceGroups());
        kernel.getAs<IntegrateVerletStepKernel>().execute(*context, *this);
    }
}

// Advance a simulation through time by taking one step with the providing
// external forces instead of the forces which is computed from force fileds
// based on current positions.
// Note that, you should call updateContextState() firstly before calling
// this function. (Added by zhubin, Nov. 14, 2021)
void VerletIntegrator::oneStep(const std::vector<Vec3>& forces) {
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");
    context->setForces(forces);
    kernel.getAs<IntegrateVerletStepKernel>().execute(*context, *this);
}