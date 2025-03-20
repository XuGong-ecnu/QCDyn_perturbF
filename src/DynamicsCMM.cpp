/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 16, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "DynamicsCMM.h"

void DynamicsCMM::init() {
    // 0. Initialize date members of DynamicsMQCBase.
    DynamicsMQCBase::init();
    if (dyn_type.substr(0,3) != "CMM")
        throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for DynamicsCMM.");
    if (representation != "diabatic")
            throw std::runtime_error("ERROR: Only diabatic representation is supported for dyn_type=" + dyn_type);
    if (system_type != "model")
        throw std::runtime_error("ERROR: Realistic system is not supported for dyn_type=" + dyn_type);
    // Bind functions of electronic propagation according to dyn_type
    switch (dyn_type[3]) {
        case '1':
            MOVE_elec = std::bind(&DynamicsCMM::MOVE_elec_CMM1, this); break;
        case '3':
            MOVE_elec = std::bind(&DynamicsCMM::MOVE_elec_CMM3, this); break;
        case '4':
            MOVE_elec = std::bind(&DynamicsCMM::MOVE_elec_CMM4, this); break;
        case '5':
            MOVE_elec = std::bind(&DynamicsCMM::MOVE_elec_CMM5, this); break;
        case '6':
            MOVE_elec = std::bind(&DynamicsCMM::MOVE_elec_CMM6, this); break;
        default:
            throw std::runtime_error("ERROR: Unsupported dyn_type=" + dyn_type + " for DynamicsCMM.");
    }
}

void DynamicsCMM::samplingElec() {
    // * CMM1 and CMM3-6 use same initial sampling (4 mapping variables).
    // Ref: J. Chem. Phys. 151, 024105 (2019) Eq.28
    // For initial conditions, requires total population Σ_n{P_n}=1,
    // in which, 0 << P_n << 1 , P_n = x_n*py_n-y_n*px_n.
    // which means all {x,y,px,py} are uniformly sampled on a surface of sphere.
    // Note that it’s a difference in sampling of the initial state from the in
    // ref: J. Liu, J. Chem. Phys. 2016, 145, 204105. eq. C3
    // A standard method is to generate standard normal (Gaussian) distribution
    // and construct a unit vector from them. That is, when Xi∼N(0,1) and
    // λ=srqt(Σ_i{Xi^2}), then Xi/λ is uniformly distributed on the sphere.
    // This method works well for d-dimensional spheres, too.
    // We can call this method as normalised Gaussian sampling.
    // Note, this method is similar/same as full-sphere sampling in spin-mapping.
    std::vector<double>& x = Elec->x;
    std::vector<double>& y = Elec->y;
    std::vector<double>& px = Elec->px;
    std::vector<double>& py = Elec->py;
    std::normal_distribution<double> normal_dist(0.0, 1.0);
    double sum = 0.0;
    for (int i = 0; i < DOFe; i++) {
        x[i] = normal_dist(elec_gen);
        y[i] = normal_dist(elec_gen);
        px[i] = -y[i]; // only x, y is independent, So, 4 mapping variables (CMM1)
        py[i] = x[i];  // is very simillar as 2 mapping variables (CMM2)
        sum += x[i]*x[i] + y[i]*y[i] + px[i]*px[i] + py[i]*py[i];
    }
    for (int i = 0; i < DOFe; i++) {
        x[i] *= sqrt(2.0/sum);
        y[i] *= sqrt(2.0/sum);
        px[i] *= sqrt(2.0/sum);
        py[i] *= sqrt(2.0/sum);
    }
}

void DynamicsCMM::getInitialDensity() {
    // For CMM method, the init_density is the A(0) operator in correlation function.
    // The defination of population or coherence operator see the Hamiltonian:
    // Ref: J. Liu, J. Chem. Phys. 2016, 145, 204105.
    for (int i = 0; i < init_state.size(); ++i)  // initial state i = j*DOFe+k, |j><k|
        init_density[i] = getOperator(init_state[i]/DOFe, init_state[i]%DOFe);
}

void DynamicsCMM::updateDensityMatrix() {
    // For CMM method, correlation function (observable) is defined as in Eq.32:
    // Ref: X. He, J. Liu, J. Chem. Phys. 2019, 151, 024105. Eq.32
    // σ_jk(t) = DOFe*(DOFe+1)* <A_nm(0)B_kj(t)> - 1. <...> is averaged by ntraj.
    // A_nm(0) has been computed and stored in init_density. nm is initial indices.
    // Similar formula can be found in Jian's JPCL paper (eCMM):
    // Ref2: X. He, J. Liu, J Phys. Chem. Lett. 2021, 12, 2496. Eq.20
    // However, they can be used for population-population operator only.
    // Here, we use the formula used in spin-mapping to compute σ_jk(t).
    // Ref: J. Chem. Phys. 152, 084110 (2020) Eq. 40-42
    // Actually, they are equivalent. e.g., when γ=0,(CMM2) is the SPM-Q method.
    // Note, in this formula, the population or coherence operator at time
    // t has a factor. Since the γ in all CMMs method is zero (like SPM-Q method),
    // therefore, then factor is DOFe+1, and for population, minus 1 is needed.
    // The details, please refer to spin-mapping code/paper.
    for (int i = 0; i < init_density.size(); ++i)
        for (int j = 0; j < DOFe; ++j)
            for (int k = 0; k < DOFe; ++k) {
                Complex B = (DOFe+1.0) * getOperator(k, j); // B(q_k, q_j)(t)
                if (j == k)   // population operator, minus 1 is needed. (γ=0)
                    B -= 1.0; // The factor and -1 also can be found in JPCL,Eq.20
                // jk element of RDM matrix at time t: σ_jk(t)
                // DOFe is normalized factor arised from full-sphere sampling.
                Complex sigma = (double)(DOFe) * init_density[i] * B;
                                // Here, RDM is stored as DOFe^2-dimentional vector (index j*DOFe+k).
                // RDM_current[i][step/RDM_steps][j*DOFe+k] = sigma;
                // For RDM_average, the value of each element is accumaleated (+=)
                // and averaged by number of trajectories.
                // RDM_average[i][step/RDM_steps][j*DOFe+k] += sigma / (double)(ntraj);

                //calculate RDM of the snapshot
                RDM[i][j*DOFe+k] = sigma;
            }
}

// TODO: perhaps can use the general one by taking MOVE_elec as parameter
void DynamicsCMM::oneStepForDiabaticModel() {
    // This is Velocity Verlet integrator
    std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    // Update velocities/momenta with a half step using effective force.
    MOVE_Half_V(ha->F, ha->V);
    // Update positions with a full step
    MOVE_R(ha->V, ha->R);
    // Update H and F_all since R has been changed.
    updateH();
    // Update electronic mapping variables.
    MOVE_elec();
    // Update effective forces
    updateF();
    // Update velocities/momenta with a half step using new effective forces
    MOVE_Half_V(ha->F, ha->V);
}

void DynamicsCMM::updateDiabaticModelForces() {
    // * Same as in DynamicsMQCBase() but for different operators
    std::shared_ptr<HamiltonianModelBase> ha = std::static_pointer_cast<HamiltonianModelBase>(Ha);
    // Reset all elements to 0 with keeping same size
    std::fill(ha->F.begin(), ha->F.end(), 0.0);
    for (int j = 0; j < DOFn; j++) {
        for (int n = 0; n < DOFe; n++)
            for (int m = 0; m < DOFe; m++)
                if (n == m) // from diagnoal of H (population)
                    ha->F[j] += (ha->F_all[n*DOFe+m][j] - ha->F_avg[j]) * getOperator(n, m).real();
                else if (n < m) // from off-diagnoal of H (coherence)
                    // multiply 2 since nm is same as mn (real part)
                    // If Condon approximation is used (coupling is contant), it is 0.
                    ha->F[j] += ha->F_all[n*DOFe+m][j] * 2.0 * getOperator(n, m).real();
        ha->F[j] /= hbar;
        ha->F[j] += ha->F_avg[j];
    }
}

Complex DynamicsCMM::getOperator(int j, int k) {
    const std::vector<double>& x = Elec->x;
    const std::vector<double>& y = Elec->y;
    const std::vector<double>& px = Elec->px;
    const std::vector<double>& py = Elec->py;
    // The formula of population or coherence operator can be found in:
    // Ref: J. Liu, J. Chem. Phys. 2016, 145, 204105. (Eq.29,49,54,56,59)
    if (j == k) // population
        if (dyn_type == "CMM1" || dyn_type == "CMM6")      // CMM1 and CMM6 have same population
            return (x[j]*py[j] - y[k]*px[k]);
        else                                               // CMM3-5 have same population
            return 0.25*((x[j]+py[j])*(x[j]+py[j]) + (y[k]-px[k])*(y[k]-px[k]));
    else // coherence
        if (dyn_type == "CMM1" || dyn_type == "CMM3")      // CMM1 and CMM3 have same coherence
            return ((x[j]*py[k] - y[j]*px[k]) + I*(x[j]*px[k] + y[j]*py[k]));
        else if (dyn_type == "CMM4" || dyn_type == "CMM6") // CMM4 and CMM6 have same coherence
            return 0.5*((x[j]*x[k] + px[j]*px[k]) + (y[j]*y[k] + py[j]*py[k])
                + I*((x[j]*px[k] - px[j]*x[k]) + (y[j]*py[k] - py[j]*y[k])));
        else                                               // CMM5 coherence
            return 0.25*((x[j]+py[j])*(x[k]+py[k]) + (y[j]-px[j])*(y[k]-px[k])
                + I*((x[j]+py[j])*(y[k]-px[k]) - (y[j]-px[j])*(x[k]+py[k])));
}

// TODO multi-steps for electronic propagation
void DynamicsCMM::MOVE_elec_CMM1() {
    std::vector<double>& x = Elec->x;
    std::vector<double>& y = Elec->y;
    std::vector<double>& px = Elec->px;
    std::vector<double>& py = Elec->py;
    const Real_Matrix& H_old = Ha->Heff_old;
    const Real_Matrix& H_new = Ha->Heff;
    // For all CMM methods, electronic steps per nuclear step is 1
    const double dt = DT;
    // CMM1 propagation Ref: J. Chem. Phys. 145, 204105 (2016) eq.C1, eq.30
    // Velocity Verlet integrtation ref: J. Chem. Phys. 138, 144106 (2013):
    // 1. p(t+0.5dt) = p(t) - 0.5dt*h(t)*q(t) [h(t) = H_old]
    // 2. q(t+dt) = q(t)+ dt*h[t+dt]*p(t+0.5dt) [h(t) = H_new]
    // 3. p(t+0.5dt) = p(t+0.5dt) - 0.5*dt*h(t+dt)*q(t+dt) [h(t) = H_new]
    // here, {x, px}, like momenta p; {y, py} like position q
    // Step 1. propagate {x, px} for a half time interval 0.5dt using H(t) and {y, py}(t)
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe; m++) {
            x[n] += -0.5 * dt * (H_old[n][m] * y[m]);
            px[n] += -0.5 * dt * (H_old[n][m] * py[m]);
        }
    // Step 2. propagate {y, py} for a time interval dt using H(t+dt) and {x, px}(t+0.5dt)
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe; m++) {
            y[n] += dt * (H_new[n][m] * x[m]);
            py[n] += dt * (H_new[n][m] * px[m]);
        }
    // Step 3. propagate {x, px} for a half time interval 0.5dt using H(t+dt) and {y, py}(t+dt)
    // same as step 1 but using H(t+dt) (H_new)
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe; m++) {
            x[n] += -0.5 * dt * (H_new[n][m] * y[m]);
            px[n] += -0.5 * dt * (H_new[n][m] * py[m]);
        }
}

void DynamicsCMM::MOVE_elec_CMM2() {
    // This original numerical algorithm to do electronic propagation for CMM2
    // method that defined in JianLiu's paper [J. Chem. Phys. 2016, 145, 204105.].
    // However, we don't use it here, we use RK4 method to do this like other
    // mapping dynmaics (LSC/SQC/SPM).
    std::vector<double>& q = Elec->q;
    std::vector<double>& p = Elec->p;
    const Real_Matrix& H_old = Ha->Heff_old;
    const Real_Matrix& H_new = Ha->Heff;
    // For all CMM methods, electronic steps per nuclear step is 1
    const double dt = DT;
    // CMM2 propagation Ref: J. Chem. Phys. 145, 204105 (2016) eq.C4, eq.41
    // Velocity Verlet integrtation ref: J. Chem. Phys. 138, 144106 (2013):
    // 1. p(t+0.5dt) = p(t) - 0.5dt*h(t)*q(t) [h(t) = H_old]
    // 2. q(t+dt) = q(t)+ dt*h[t+dt]*p(t+0.5dt) [h(t) = H_new]
    // 3. p(t+0.5dt) = p(t+0.5dt) - 0.5*dt*h(t+dt)*q(t+dt) [h(t) = H_new]
    // Step 1. propagate {p} for a half time interval 0.5dt using H(t) and {q}(t)
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe; m++)
            p[n] += -0.5 * dt * (H_old[n][m] * q[m]);
    // Step 2. propagate {q} for a time interval dt using H(t+dt) and {p}(t+0.5dt)
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe; m++)
            q[n] += dt * (H_new[n][m] * p[m]);
    // Step 3. propagate {p} for a half time interval 0.5dt using H(t+dt) and {q}(t+dt)
    // same as step 1 but using H(t+dt) and {q}(t+dt)
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe; m++)
            p[n] += -0.5 * dt * (H_new[n][m] * q[m]);
}

void DynamicsCMM::MOVE_elec_CMM3() {
    std::vector<double>& x = Elec->x;
    std::vector<double>& y = Elec->y;
    std::vector<double>& px = Elec->px;
    std::vector<double>& py = Elec->py;
    const Real_Matrix& H_old = Ha->Heff_old;
    const Real_Matrix& H_new = Ha->Heff;
    // For all CMM methods, electronic steps per nuclear step is 1
    const double dt = DT;
    // CMM3 propagation Ref: J. Chem. Phys. 145, 204105 (2016) eq.C7 - C10, eq.50
    // Velocity Verlet integrtation ref: J. Chem. Phys. 138, 144106 (2013):
    // 1. p(t+0.5dt) = p(t) - 0.5dt*h(t)*q(t) [h(t) = H_old]
    // 2. q(t+dt) = q(t)+ dt*h[t+dt]*p(t+0.5dt) [h(t) = H_new]
    // 3. p(t+0.5dt) = p(t+0.5dt) - 0.5*dt*h(t+dt)*q(t+dt) [h(t) = H_new]
    // here, {x, px}, like momenta p; {y, py} like position q
    // Step 1. propagate {x, px} for a half time interval 0.5dt using H(t) and {y, py}(t)
    for (int n = 0; n < DOFe; n++) {
        // (1) get px[t+0.25dt] using x[t] and py[t] (half step: 0.25dt)
        // Note: px[n] <- px[n] - 0.25*dt* (0.5*H[n][n]py[n] + ΣmH[n][m]py[m](m!=n))
        // ======> px[n] <- px[n] + 0.125*dt*(H[n][n]py[n]) -0.25dt*ΣmH[n][m]y[m]
        px[n] += -0.125 * dt * H_old[n][n] * (x[n]-py[n]);
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.25 * dt * (H_old[n][m] * py[m]);
        // (2) get x[t+0.5dt] using y[t] and px[t+0.25dt] (full step: 0.5dt)
        x[n] += 0.25 * dt * H_old[n][n] * (px[n]+y[n]);
        for (int m = 0; m < DOFe; m++)
            x[n] += -0.5 * dt * (H_old[n][m] * y[m]);
        // (3) get px[t+0.5dt] using x[t+0.5dt] and py[t] (half step: 0.25dt) same as (1)
        px[n] += -0.125 * dt * H_old[n][n] * (x[n]-py[n]);
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.25 * dt * (H_old[n][m] * py[m]);
    }
    // Step 2. propagate {y, py} for a time interval dt using H(t+dt) and {x, px}(t+0.5dt)
    for (int n = 0; n < DOFe; n++) {
        // (1) get py[t+0.5dt] using px[t+0.5dt] and y[t] (half step: 0.5dt)
        py[n] += -0.25 * dt * H_new[n][n] * (px[n]+y[n]);
        for (int m = 0; m < DOFe; m++)
            py[n] += 0.5 * dt * (H_new[n][m] * px[m]);
        // (2) get y[t+dt] using x[t+0.5dt] and py[t+0.5dt] (full step: dt)
        y[n] += 0.5 * dt * H_new[n][n] * (py[n]-x[n]);
        for (int m = 0; m < DOFe; m++)
            y[n] += dt * (H_new[n][m] * x[m]);
        // (3) get py[t+dt] using px[t+0.5dt] and y[t+dt] (half step: 0.5dt) same as (1)
        py[n] += -0.25 * dt * H_new[n][n] * (px[n]+y[n]);
        for (int m = 0; m < DOFe; m++)
            py[n] += 0.5 * dt * (H_new[n][m] * px[m]);
    }
    // Step 3. propagate {x, px} for a half time interval 0.5dt using H(t+dt) and {y, py}(t+dt)
    // same as step 1 but using H(t+dt) (H_new)
    for (int n = 0; n < DOFe; n++) {
        px[n] += -0.125 * dt * H_new[n][n] * (x[n]-py[n]);
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.25 * dt * (H_new[n][m] * py[m]);
        x[n] += 0.25 * dt * H_new[n][n] * (px[n]+y[n]);
        for (int m = 0; m < DOFe; m++)
            x[n] += -0.5 * dt * (H_new[n][m] * y[m]);
        px[n] += -0.125 * dt * H_new[n][n] * (x[n]-py[n]);
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.25 * dt * (H_new[n][m] * py[m]);
    }
}

void DynamicsCMM::MOVE_elec_CMM4() {
    std::vector<double>& x = Elec->x;
    std::vector<double>& y = Elec->y;
    std::vector<double>& px = Elec->px;
    std::vector<double>& py = Elec->py;
    const Real_Matrix& H_old = Ha->Heff_old;
    const Real_Matrix& H_new = Ha->Heff;
    // For all CMM methods, electronic steps per nuclear step is 1
    const double dt = DT;
    // * The EOM of CMM4 is very similar as CMM3, but exchange x <-> py
    // CMM4 propagation Ref: J. Chem. Phys. 145, 204105 (2016) eq.C11 - C13, eq.55
    // Velocity Verlet integrtation ref: J. Chem. Phys. 138, 144106 (2013):
    // 1. p(t+0.5dt) = p(t) - 0.5dt*h(t)*q(t) [h(t) = H_old]
    // 2. q(t+dt) = q(t)+ dt*h[t+dt]*p(t+0.5dt) [h(t) = H_new]
    // 3. p(t+0.5dt) = p(t+0.5dt) - 0.5*dt*h(t+dt)*q(t+dt) [h(t) = H_new]
    // here, {px, py}, like momenta p; {x, y} like position q
    // Step 1. propagate {px, py} for a half time interval 0.5dt using H(t) and {x, y}(t)
    for (int n = 0; n < DOFe; n++) {
        // (1) get px[t+0.25dt] using py[t] and x[t] (half step: 0.25dt)
        px[n] += -0.125 * dt * H_old[n][n] * (py[n]-x[n]);
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.25 * dt * (H_old[n][m] * x[m]);
        // (2) get py[t+0.5dt] using y[t] and px[t+0.25dt] (full step: 0.5dt)
        py[n] += 0.25 * dt * H_old[n][n] * (px[n]+y[n]);
        for (int m = 0; m < DOFe; m++)
            py[n] += -0.5 * dt * (H_old[n][m] * y[m]);
        // (3) get px[t+0.5dt] using py[t+0.5dt] and x[t] (half step: 0.25dt) same as (1)
        px[n] += -0.125 * dt * H_old[n][n] * (py[n]-x[n]);
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.25 * dt * (H_old[n][m] * x[m]);
    }
    // Step 2. propagate {x, y} for a time interval dt using H(t+dt) and {px, py}(t+0.5dt)
    for (int n = 0; n < DOFe; n++) {
        // (1) get x[t+0.5dt] using px[t+0.5dt] and y[t] (half step: 0.5dt)
        x[n] += -0.25 * dt * H_new[n][n] * (px[n]+y[n]);
        for (int m = 0; m < DOFe; m++)
            x[n] += 0.5 * dt * (H_new[n][m] * px[m]);
        // (2) get y[t+dt] using x[t+0.5dt] and py[t+0.5dt] (full step: dt)
        y[n] += 0.5 * dt * H_new[n][n] * (x[n]-py[n]);
        for (int m = 0; m < DOFe; m++)
            y[n] += dt * (H_new[n][m] * py[m]);
        // (3) get x[t+dt] using px[t+0.5dt] and y[t+dt] (half step: 0.5dt) same as (1)
        x[n] += -0.25 * dt * H_new[n][n] * (px[n]+y[n]);
        for (int m = 0; m < DOFe; m++)
            x[n] += 0.5 * dt * (H_new[n][m] * px[m]);
    }
    // Step 3. propagate {px, py} for a half time interval 0.5dt using H(t+dt) and {x, y}(t+dt)
    // same as step 1 but using H(t+dt) (H_new)
    for (int n = 0; n < DOFe; n++) {
        px[n] += -0.125 * dt * H_new[n][n] * (py[n]-x[n]);
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.25 * dt * (H_new[n][m] * x[m]);
        py[n] += 0.25 * dt * H_new[n][n] * (px[n]+y[n]);
        for (int m = 0; m < DOFe; m++)
            py[n] += -0.5 * dt * (H_new[n][m] * y[m]);
        px[n] += -0.125 * dt * H_new[n][n] * (py[n]-x[n]);
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.25 * dt * (H_new[n][m] * x[m]);
    }
}

void DynamicsCMM::MOVE_elec_CMM5() {
    std::vector<double>& x = Elec->x;
    std::vector<double>& y = Elec->y;
    std::vector<double>& px = Elec->px;
    std::vector<double>& py = Elec->py;
    const Real_Matrix& H_old = Ha->Heff_old;
    const Real_Matrix& H_new = Ha->Heff;
    // For all CMM methods, electronic steps per nuclear step is 1
    const double dt = DT;
    // CMM5 propagation derived from: J. Chem. Phys. 145, 204105 (2016) eq.57,
    // which is similar as CMM3 method.
    // Velocity Verlet integrtation ref: J. Chem. Phys. 138, 144106 (2013):
    // 1. p(t+0.5dt) = p(t) - 0.5dt*h(t)*q(t) [h(t) = H_old]
    // 2. q(t+dt) = q(t)+ dt*h[t+dt]*p(t+0.5dt) [h(t) = H_new]
    // 3. p(t+0.5dt) = p(t+0.5dt) - 0.5*dt*h(t+dt)*q(t+dt) [h(t) = H_new]
    // here, {x, px}, like momenta p; {y, py} like position q
    // Step 1. propagate {x, px} for a half time interval 0.5dt using H(t) and {y, py}(t)
    // (1) get px[t+0.25dt] using py[t] and x[t] (half step: 0.25dt)
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.125 * dt * H_old[n][m] * (py[m]+x[m]);
    // (2) get x[t+0.5dt] using px[t+0.25dt] and y[t] (full step: 0.5dt)
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe; m++)
            x[n] += 0.25 * dt * H_old[n][m] * (px[m]-y[m]);
    // (3) get px[t+0.5dt] using py[t] and x[t+0.5dt] (same as (1)) (half step: 0.25dt)
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.125 * dt * H_old[n][m] * (py[m]+x[m]);
    // Step 2. propagate {y, py} for a time interval dt using H(t+dt) and {x, px}(t+0.5dt)
    // (1) get py[t+0.5dt] using px[t+0.5dt] and y[t] (half step: 0.5dt)
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe; m++)
            py[n] += 0.25 * dt * H_new[n][m] * (px[m]-y[m]);
    // (2) get y[t+dt] using x[t+0.5dt] and py[t+0.5dt] (full step: dt)
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe; m++)
            y[n] += 0.5 * dt * H_new[n][m] * (x[m]+py[m]);
    // (3) get py[t+dt] using px[t+0.5dt] and y[t+dt] (same as (1)) (half step: 0.5dt)
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe; m++)
            py[n] += 0.25 * dt * H_new[n][m] * (px[m]-y[m]);
    // Step 3. propagate {x, px} for a half time interval 0.5dt using H(t+dt) and {y, py}(t+dt)
    // same as step 1 but using H(t+dt) (H_new)
    // (1) get px[t+0.75dt] using py[t+dt] and x[t+0.5dt] (half step: 0.25dt)
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.125 * dt * H_old[n][m] * (py[m]+x[m]);
    // (2) get x[t+dt] using px[t+0.75dt] and y[t+dt] (full step: 0.5dt)
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe; m++)
            x[n] += 0.25 * dt * H_old[n][m] * (px[m]-y[m]);
    // (3) get px[t+dt] using py[t+dt] and x[t+dt] (same as (1)) (half step: 0.25dt)
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.125 * dt * H_old[n][m] * (py[m]+x[m]);
}

void DynamicsCMM::MOVE_elec_CMM6() {
    std::vector<double>& x = Elec->x;
    std::vector<double>& y = Elec->y;
    std::vector<double>& px = Elec->px;
    std::vector<double>& py = Elec->py;
    const Real_Matrix& H_old = Ha->Heff_old;
    const Real_Matrix& H_new = Ha->Heff;
    // For all CMM methods, electronic steps per nuclear step is 1
    const double dt = DT;
    // CMM6 propagation derived from: J. Chem. Phys. 145, 204105 (2016) eq.60
    // which is similar as CMM4 method (the difference is factor).
    // Velocity Verlet integrtation ref: J. Chem. Phys. 138, 144106 (2013):
    // 1. p(t+0.5dt) = p(t) - 0.5dt*h(t)*q(t) [h(t) = H_old]
    // 2. q(t+dt) = q(t)+ dt*h[t+dt]*p(t+0.5dt) [h(t) = H_new]
    // 3. p(t+0.5dt) = p(t+0.5dt) - 0.5*dt*h(t+dt)*q(t+dt) [h(t) = H_new]
    // here, {px, py}, like momenta p; {x, y} like position q
    // Step 1. propagate {px, py} for a half time interval 0.5dt using H(t) and {x, y}(t)
    for (int n = 0; n < DOFe; n++) {
        // (1) get px[t+0.25dt] using py[t] and x[t] (half step: 0.25dt)
        px[n] += -0.25 * dt * H_old[n][n] * (py[n]-x[n]);
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.25 * dt * (H_old[n][m] * x[m]);
        // (2) get py[t+0.5dt] using y[t] and px[t+0.25dt] (full step: 0.5dt)
        py[n] += 0.5 * dt * H_old[n][n] * (px[n]+y[n]);
        for (int m = 0; m < DOFe; m++)
            py[n] += -0.5 * dt * (H_old[n][m] * y[m]);
        // (3) get px[t+0.5dt] using py[t+0.5dt] and x[t] (half step: 0.25dt) same as (1)
        px[n] += -0.25 * dt * H_old[n][n] * (py[n]-x[n]);
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.25 * dt * (H_old[n][m] * x[m]);
    }
    // Step 2. propagate {x, y} for a time interval dt using H(t+dt) and {px, py}(t+0.5dt)
    for (int n = 0; n < DOFe; n++) {
        // (1) get x[t+0.5dt] using px[t+0.5dt] and y[t] (half step: 0.5dt)
        x[n] += -0.5 * dt * H_new[n][n] * (px[n]+y[n]);
        for (int m = 0; m < DOFe; m++)
            x[n] += 0.5 * dt * (H_new[n][m] * px[m]);
        // (2) get y[t+dt] using x[t+0.5dt] and py[t+0.5dt] (full step: dt)
        y[n] += dt * H_new[n][n] * (x[n]-py[n]);
        for (int m = 0; m < DOFe; m++)
            y[n] += dt * (H_new[n][m] * py[m]);
        // (3) get x[t+dt] using px[t+0.5dt] and y[t+dt] (half step: 0.5dt) same as (1)
        x[n] += -0.5 * dt * H_new[n][n] * (px[n]+y[n]);
        for (int m = 0; m < DOFe; m++)
            x[n] += 0.5 * dt * (H_new[n][m] * px[m]);
    }
    // Step 3. propagate {px, py} for a half time interval 0.5dt using H(t+dt) and {x, y}(t+dt)
    // same as step 1 but using H(t+dt) (H_new)
    for (int n = 0; n < DOFe; n++) {
        px[n] += -0.25 * dt * H_old[n][n] * (py[n]-x[n]);
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.25 * dt * (H_old[n][m] * x[m]);
        py[n] += 0.5 * dt * H_old[n][n] * (px[n]+y[n]);
        for (int m = 0; m < DOFe; m++)
            py[n] += -0.5 * dt * (H_old[n][m] * y[m]);
        px[n] += -0.25 * dt * H_old[n][n] * (py[n]-x[n]);
        for (int m = 0; m < DOFe; m++)
            px[n] += -0.25 * dt * (H_old[n][m] * x[m]);
    }
}