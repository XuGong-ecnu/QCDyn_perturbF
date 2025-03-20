/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 20, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "HamiltonianMSH.h"

void HamiltonianMSH::init() {
    // * 0. Initialize data members and Hamiltonian.
    HamiltonianModelBase::init();
    MSH_type = param.getStr("MSH_type");
    if (model_type != "MSH")
        throw std::runtime_error("ERROR: Undefined model_type=" + model_type + " for multi-state Harmonic model.");
    if (MSH_type != "I" && MSH_type != "II")
        throw std::runtime_error("ERROR: Undefined MSH_type=" + MSH_type + " for multi-state Harmonic model.");
    if (MSH_type == "II" && DOFe == 2)
        throw std::runtime_error("ERROR: MSH model II can be used only for DOFe > 2.");
    // For MSH model I, number of normal modes, N = DOFn/(DOFe-1).
    // For MSH model II, number of normal modes, N = DOFn/(DOFe-2).
    // If groud state is included, then model I can be reduced to model II,
    // which can save time and get identical result.
    const int length = MSH_type == "I" ? (DOFe-1) : (DOFe-2);
    if (DOFn % length != 0)
        throw std::runtime_error("ERROR: For MSH model, total nuclear DOFn should "
            "be divided exactly by DOFe-1 (MSH model I) or DOFe-2 (MSH model II).");
    N = DOFn/length;

    // * 1. Initialize and get model parameters by building or loading
    omega.resize(N, 0);
    c.resize(N, 0);
    // For MSH model, req is the equilibrium shifts between each pair of states
    // (off-diagonal) along each mode.  And for each harmonic nomral mode, the shift
    // is a DOFe*(DOFe-1) matrix. So, the size of req (1D vector) is N*[DOFe*(DOFe-1)].
    req.resize(N*DOFe*(DOFe-1), 0);
    shift.resize(DOFn, 0);
    HamiltonianModelBase::initializeModelParameters();

    // * 2. Extend parameters for nulcear sampling and energy/forces calculations
    // Note that after initializeModelParameters(), the model parameters have been
    // created or loaded, but for convenience and general, we extend the length
    // of omega, req, and set shift accordingly.
    // Let omega size = DOFn, the R of multi-state share same set omega j.
    // DOFn/N = DOFe-1 (model I), or DOFe-2 (model II), the number of columns of
    // S matrix used in nulear sampling and calculations.
    omega.resize(DOFn, 0);
    for (int i = 1; i < DOFn/N; ++i)
        for (int j = 0; j < N; ++j)
            omega[i*N+j] = omega[j];
    // Set equilibrium shifts of each R according the state index used to sample.
    // Generally, the last state is ground state and as sample state. Shift of
    // R_j is -S_j_(F-1)F (Note negative sign). j*DOFe*(DOFe-1) is the starting
    // index of S_j, and sample_state*(DOFe-1)+i is the index of S_(F-1)F in S
    // matrix. Here, S_(F-1)F is the values of Fth row (F-1)th col of S matrix.
    // Note that the state index in input file should start from 0.
    // Note that, here req is always a full matrix in DOFe*(DOFe-1) base,
    // regardless of using Model 3 or Model 2.
    const int sample_state = param.getInt("sample_state");
    if (sample_state >= DOFe || sample_state < 0)
        throw std::runtime_error("ERROR: input sample_state=" + std::to_string(sample_state) + " is out of scope (starting from 0).");
    for (int i = 0; i < DOFn/N; ++i)
        for (int j = 0; j < N; ++j)
            shift[i*N+j] = -req[j*DOFe*(DOFe-1)+sample_state*(DOFe-1)+i];

    // * 3. Save Hamiltonian matrix from file (if required) in input energy unit
    // Do this in subclasses since the subclasses init() may modify it.
    const std::string H_save= param.getStr("H_save");
    if (!H_save.empty())
        saveHamiltonianMatrix(H_save);
}

void HamiltonianMSH::buildModelParameters() {
    // Ref: Multistate-Harmonic-Model_report
    // * 1. Get file name of energy file from global paramters (param)
    // Energy input file: potential energies of each state.
    // File requirment: one column for one state in order and
    // number of cols = DOFe, no headers, no index of row.
    std::string energyfile = param.getStr("energy_load");
    if (energyfile.empty())
        throw std::runtime_error("ERROR: Please provide an input energy file.");
    // Get energy coorection (if any) of each state from input.
    // Each value should be separated by one comma "," with the energy_unit
    const std::string energy_correction_str = param.getStr("energy_correction");
    std::vector<double> energy_correction;
    if (!energy_correction_str.empty()) {
        SplitString(energy_correction, energy_correction_str);
        if (energy_correction.size() != DOFe)
            throw std::runtime_error("ERROR: The number of values for energy_correction should be equal to DOFe.");
    }
    else // By defalut, it is zero.
        energy_correction.resize(DOFe, 0);
    // * 2. Load energy data from file
    std::vector<std::vector<double>> V(DOFe); // to store energies of each state
    std::ifstream inp(energyfile);
    if (!inp)
        throw std::runtime_error("ERROR: Can't open input energy file: " + energyfile);
    // According to the file extension name (csv or dat) to decide the
    // function used to process each line of file
    std::function<void(std::vector<double>&, const std::string&)> split;
    std::string type = GetFileSuffix(energyfile);
    if (type == "dat") // fields are separated by at least one space.
        // For overloading functions, we must specify the type of function pointer
        // e.g., int(*)(int) is the type of function pointer, which can be as parameter
        // for function, and int(int) is the type of function.
        split = (void(*)(std::vector<double>&, const std::string&))SplitLine;
    else if (type == "csv") // fields are separated by one comma (,).
        split = std::bind((void(*)(std::vector<double>&, const std::string&, char))
            SplitString, std::placeholders::_1, std::placeholders::_2, ',');
    else
        throw std::runtime_error("ERROR: Unknown file type for input: " + type);
    // loop for loading values of each line until end of file (unknown length)
    int lineNumber = 0;
    for (std::string line; getline(inp, line); lineNumber++) {
        std::vector<double> fields;
        split(fields, line);
        if (fields.size() < DOFe) // number of cols must be DOFe
            throw std::runtime_error("ERROR: Too few columns at line " + std::to_string(lineNumber+1));
        for (int i = 0; i < DOFe; ++i)
            V[i].emplace_back(fields[i]+energy_correction[i]); // unit same as input
    }
    inp.close();
    const int V_size = V[0].size(); // the number/length of energy data (V)
    // * 3. Compute energy gap: U_XY = V_X - V_Y
    // Note, the number of sets of U, Er, S_j is DOFe*(DOFe-1)/2.
    // The order of U is U_12, U_13, U_23,...U_(DOFe-1)(DOFe),..., (order by rows)
    // Note that the index of state (XY) starts from 1 in paper
    // Using the order number (single index) represents a pair of states (XY):
    // U_index -> U_xy, 1 <= x < y <= DOFe, index = (y-1)(y-2)/2 + x - 1
    // Note that first index is 0 in code (thus -1)
    // If we know index, then (index+1) = (y-1)(y-2)/2 + x
    // We find the largest y makes (y-1)(y-2)/2 < (index+1)
    // y < [3 + sqrt(1+8(index+1))] = b, then y is the interger of b-1,
    // e.g., b = 7.8, then y = 7; b = 6 y = 5. Then, x can be decided.
    const int U_count =  DOFe*(DOFe-1)/2; // number of sets of U_xy: DOFe*(DOFe-1)/2
    // U vector: U_12, U_13, U_23,...U_(DOFe-1)(DOFe). Each element is a vector.
    // Note that the unit of it is au at last.
    std::vector<std::vector<double>> U(U_count, std::vector<double>(V_size, 0));
    for (int y = 2; y <= DOFe; ++y)
        for (int x = 1; x < y; ++x) {
            // The index for U_XY in U vector is: index = (y-1)(y-2)/2 + x-1
            int index = (y-1)*(y-2)/2 + x-1;
            for (int i = 0; i < V_size; ++i)
                // state index from 0 in code thus x-1, y-1
                // and here we convert input energy unit to au.
                U[index][i] = (V[x-1][i] - V[y-1][i]) * unit2au;
        }
    // * 4. Compute energy-gap variance: σ^2 = C_UU(0) = <U^2>-<U>^2
    std::vector<double> U_avg(U_count, 0); // average of U_XY: <U>  in au
    std::vector<double> sigma2(U_count, 0); // variance of U_XY: σ^2 in au
    for (int n = 0; n < U_count; ++n) {
        double sum = 0, sum2 = 0; // sum = ΣU; sum2 = Σ(U^2)
        for (int i = 0; i < V_size; ++i) {
            sum += U[n][i];
            sum2 += U[n][i]*U[n][i];
        }
        U_avg[n] = sum/V_size;
        sigma2[n] = sum2/V_size - U_avg[n]*U_avg[n];
    }
    // * 5. Compute reorganization energy: Er = σ^2/2kT (k is k_B constant)
    const double temperature = param.getDouble("temperature"); // real unit in K
    if (temperature <= 0)
        throw std::runtime_error("ERROR: Illegal value for temperature (requires > 0).");
    // Note the unit of U has been au, then convert kT to au
    const double kT = temperature * kT2au; // kT in unit of au
    std::vector<double> Er(U_count, 0); // Er in au
    for (int n = 0; n < U_count; ++n)
        Er[n] = 0.5 * sigma2[n] / kT; // sigma2 in au^2
    // * 6. Compute reaction free energy: ΔE = - Er - <U>
    std::vector<double> deltaE(U_count, 0); // ΔE in au
    for (int n = 0; n < U_count; ++n)
        deltaE[n] = -Er[n] - U_avg[n];
    // Output average of U, reorganization energy and free energy to file
    // TCF_name denotes the prefix of the files of intermediate results, by
    // default it is "MSH".
    const std::string TCF_prefix = param.getStr("TCF_prefix");
    std::string outfile = TCF_prefix + "_U_sigma_Er_deltaE.dat";
    FILE* out = CheckFile(outfile);
    // The unit of enengy in input, if not specified, assuming au.
    std::string energy_unit = param.getStr("energy_unit");
    if (energy_unit.empty())
        energy_unit = "au";
    fprintf(out, "%s%s\n%s%s\n","Parameters computed from input energy file: ", energyfile.c_str(),
        "The units of them are the same as input: ", energy_unit.c_str());
    fprintf(out, "XY denotes the pair of states; <U_XY> is the mean of energy gap;\n"
        "sigma^2 is the variance of U; Er_XY is the reorganization energy;\n"
        "deltaE_XY = -Er_XY - <U_XY> is the reaction free energy.\n");
    fprintf(out, "----------------------------------------------------------------------\n");
    fprintf(out, "%6s%12s %12s %12s %12s\n", "XY  ", "<U_XY>", "sigma^2", "Er_XY", "deltaE_XY");
    for (int y = 2; y <= DOFe; ++y)
        for (int x = 1; x < y; ++x) {
            // The index for U_XY in U vector is: index = (y-1)(y-2)/2 + x-1
            // And Er, sigma2, deltaE, U_avg are the same as U
            // Here, we use same energy unit as the input.
            int index = (y-1)*(y-2)/2 + x-1;
            fprintf(out, "%3d%-3d%12.6g %12.6g %12.6g %12.6g\n", x, y, U_avg[index]/unit2au,
                sigma2[index]/unit2au/unit2au, Er[index]/unit2au, deltaE[index]/unit2au);
        }
    fclose(out);
    // * 7. Compute spectral density to get frequency: omega
    std::string spec_density = param.getStr("spec_density");
    if (spec_density.empty())
        spec_density = "TCF";  // By default, using spectral density from TCF
    // (1) Get frequency from Debye spectral density
    if (spec_density == "Debye") {
        // Get characteristic frequency (omega_c) and reorganization energy
        // (lambda) with unit converison used for Debye sepectral density.
        // Remember to do the unit converison.
        const double omega_c = param.getDouble("omega_c") * unit2au;
        const double lambda = param.getDouble("lambda") * unit2au;
        // Get omega by discretization of Debye spectral density.
        // This will give you one set omega and c (c is uesless here).
        discretizeSpectralDensity("Debye", N, omega_c, 1, lambda, omega, c);
    }
    // (2) Compute spectral density from time correlation function (TCF) of the
    // energy gap between X and Y PESs: Cuu(t)_XY = <U(t)_XY U(0)_XY> - <U_XY>^2
    // from MD simulations. Note Cuu(0) = <U(0)^2> - <U>^2 = σ^2.
    else if (spec_density == "TCF") {
        // 1. Compute Cuu(t) from U_XY
        // TCF_ntraj is the numer of trajectories in input energy file.
        // By default (not specified in input control file), it is 1.
        const int TCF_ntraj = param.getInt("TCF_ntraj");
        if (TCF_ntraj < 0 || V_size % TCF_ntraj != 0)
            throw std::runtime_error("ERROR: Illegal value for TCF_ntraj, it should be larger than 0 and "
                "the total length of energy data can be divided by it exactly.");
        // TCF_totsteps is the total steps of each traj in input energy file.
        // Here, it is computed by (length of enenrgy data in file)/ntraj
        const int TCF_totsteps = V_size / TCF_ntraj;
        // TCF_corsteps is the number of correlation stpes = correlation time/DT
        // By default (not specified in input control file), it is 800.
        const int TCF_corsteps = param.getInt("TCF_corsteps");
        if (TCF_corsteps <= 0 || TCF_corsteps > TCF_totsteps)
            throw std::runtime_error("ERROR: Illegal value for TCF_corsteps, it should be larger than 0 and "
                " less than the total steps of MD simulation.");
        // TCF_XY denotes which pair of U to be used for frequency. the input is
        // two indices separeted by one comma ",", e.g., 1,2 represnets U_12
        // If not specified, U_12 will be used by default. (starting from 1)
        const std::string TCF_XY = param.getStr("TCF_XY");
        std::vector<int> index; // index[0] = X, index[1] = Y
        SplitString(index, TCF_XY);
        if (index.size() != 2)
            throw std::runtime_error("ERROR: Illegal TCF_XY=" + TCF_XY +
                ". You must provide 2 indices (starting from 1) to represent which "
                "pair of energy gap to be use for spectral density.");
        if (index[0] > DOFe || index[0] < 1 || index[1] > DOFe || index[1] < 1 || index[0] >= index[1])
                throw std::runtime_error("ERROR: Illegal TCF_XY, require 1 <= X < Y <= DOFe.");
        // The index for U_XY in U vector is: index = (y-1)(y-2)/2 + x-1
        const int U_XY = (index[1]-1)*(index[1]-2)/2 + index[0]-1;
        // CUU will store the final CUU(t) (averaged by ntraj) with unit of au^2
        // Here, use the equivalent form of Cuu(t) = <dU(t)dU(0)> = <(Ut-<U>)*(U0-<U>)>
        // = < UtU0 + <U>^2 - <U>U0 - <U>Ut > = <UtU0> + <U>^2 - <U><U0> - <U><Ut>
        // = <UtU0> + <U>^2 - <U>^2 - <U>^2 (since each point of U can be as U0 and
        // Ut, then <U0> = <Ut> = <U>) = <UtU0> - <U>^2
        // Here, dU = U - <U> is the fluctuations of U.
        std::vector<double> CUU(TCF_corsteps, 0);
        for (int n = 0; n < TCF_ntraj; ++n) {
            // Get U_avg of current traj firstly
            double U_avg = 0;
            for (int i = 0; i < TCF_totsteps; ++i)
                U_avg += U[U_XY][n*TCF_totsteps + i];
            U_avg /= TCF_totsteps;
            // Compute dU = U - <U> for each step of current trajectory
            std::vector<double> dU(TCF_totsteps, 0);
            for (int i = 0; i < TCF_totsteps; ++i)
                dU[i] = U[U_XY][n*TCF_totsteps + i] - U_avg;
            // Compute Cuu(t) = <dU(0)dU(t)> = <U(t)U(0)>-<U>^2
            // Here, use step instead of time, and step = t/DT
            std::vector<double> CUU_temp(TCF_corsteps, 0); // Cuu(t) of current traj
            // count is the number of time for Cuu(t) be computed and used to
            // get mean of Cuu(t). One long trajectory, can be splitted a lot of
            // windows with lenth of correlation time/step. The dU(0) is shiftted
            // step by step. Note, at the tail of trajectory, the window is not
            // complete, but also can be used to compute Cuu(t).
            std::vector<int> count(TCF_corsteps, 0);
            for (int s0 = 0; s0 < TCF_totsteps; ++s0) { // step for dU(0) in traj
                int st_max = std::min(TCF_totsteps, s0+TCF_corsteps);
                for (int st = s0; st < st_max; ++st) { // step for dU(t) in traj
                    int s = st - s0; // step at time t in a window
                    CUU_temp[s] += dU[s0]*dU[st];
                    count[s]++;
                }
            }
            // Get Cuu(t) of current traj and acuumalate (average) it to final Cuu
            // Note that the unit of Cuu is au^2.
            for (int s = 0; s < TCF_corsteps; ++s) {
                CUU_temp[s] /= count[s];
                CUU[s] += CUU_temp[s]/TCF_ntraj;
            }
        }
        // Additional output the time correlation function Cuu(t)
        const double TCF_DT =  param.getDouble("TCF_DT"); // defalut 0.005 (ps)
        outfile = TCF_prefix + "_CUU_" + std::to_string(index[0]) + std::to_string(index[1]) + ".dat";
        out = CheckFile(outfile);
        fprintf(out,"%8s%24s\n","time(ps)","CUU(eV^2)"); // time in ps is step * TCF_DT
        for (int t = 0; t < CUU.size(); ++t)
            fprintf(out, "%8.3f%24.16g\n", (t+1)*TCF_DT, CUU[t]*au2eV*au2eV);
        fclose(out);
        // 2. Get frequency (omega) of discrete normal modes from the Cuu(t)
        // by solving the following equation numerically (j = 1,...,N):
        // 2Nω_j/(piCuu(0)) \int_0^∞ dt Cuu(t)/(ω_jt) sin(ω_jt) = j - 1/2,
        // which must be treated via a root solving algorithm, such as the
        // secant method (https://en.wikipedia.org/wiki/Secant_method).
        // The secant method requires two initial values, x0 and x1,
        // the tolerance of convergence, max iteration.
        // In our case, the root we need find is ω_j of function f(ω_j)
        const double DT_au = param.getDouble("TCF_DT") * ps2au; // defalut 0.005 (ps)
        const double tol = param.getDouble("TCF_tol"); // default 1e-8
        const int max_iter = param.getDouble("TCF_maxiter"); // default 1000
        if (DT_au < 0 || tol < 0 || max_iter <= 0)
        throw std::runtime_error("ERROR: Illegal value for TCF_DT (requires > 0), "
            "or TCF_tol (requires >= 1), or TCF_maxiter (requires > 0).");
        for (int j = 0; j < N; ++j) {
            // Here, we always use same initial values. If they are not close to
            // root, it may take more iterations to converge or not converge.
            double x0 = 0.1*cm2au, x1 = 0.5*cm2au; // initial 0.1 and 0.5 cm-1
            for (int i = 0; abs(x1-x0) > tol; ++i) {
                // Compute f(x0) f(x1)
                // Note the units for all parameters should be au or au^2 (Cuu)
                // and index j is 1,...,N (thus using j+1)
                double f0 = fw(x0, j+1, DT_au, CUU);
                double f1 = fw(x1, j+1, DT_au, CUU);
                // Get new value of x
                double x = (x0*f1 - x1*f0) / (f1 - f0);
                // using x1 and x instead of x0 and x1 to repeat the process
                // until we reach a sufficiently high level of precision
                // (a sufficiently small difference (tol) between x and x1)
                // then the value of x1 is the ω_j we want.
                x0 = x1;
                x1 = x;
                // Not converge within the max_iter
                if (i ==  max_iter-1) {
                    std::cout << "Try to get frequency (omega) from TCF with following conditions:\n";
                    std::cout << "TCF_tol     = " << tol << "\n";
                    std::cout << "TCF_maxiter = " << max_iter << "\n";
                    std::cout << "Current j   = " << j+1 << std::endl;
                    throw std::runtime_error("ERROR: Failed to converge for getting this frequency.");
                }
            }
            // Get the root ω_j in au.
            omega[j] = x1;
        }
        // 3. Compute the continus spectral density J(ω)
        // J(ω) = βω/4 \int_0^∞ dt Cuu(t)/cos(ωt), β = 1/kT
        // Here, we take 1.2*w_N as w_max in Jw, w_N is the largest discrete w
        const double w_max = omega[N-1] * 1.2; // in au
        const double dw = 0.2*cm2au; // frequency interval for J(ω)
        const int nw = w_max / dw; // number of w in J(ω)
        const double beta = 1.0 / kT; // β = 1/kT in au
        std::vector<double> Jw(nw, 0); // in au
        for (int j = 1; j < nw; ++j)  // strating from 1, since J(0) = 0
            Jw[j] = HamiltonianMSH::Jw(dw*j, beta, DT_au, CUU); // all unit in au
        std::string outfile = TCF_prefix + "_TCF_spectral_density_Jw.csv";
        FILE* out = CheckFile(outfile);
        fprintf(out, "%s%s%s\n", "w/au",",","J(w)/au");
        for (int j = 0; j < nw; ++j)
            fprintf(out, "%.16g%s%.16g\n", dw*j, ",", Jw[j]);
        fclose(out);
    }
    else
       throw std::runtime_error("ERROR: Unsupported spectral density: " + spec_density);
    // * 8. Determine equilibrium shifts S_j
    // * 8.1 Compute a_j: a_j = sqrt(2.0/N)*(1.0/w_j)
    std::vector<double> a(N, 0); // N is number of normal modes
    for (int j = 0; j < N; ++j)
        a[j] = sqrt(2.0/N)/omega[j];
    // * 8.2 Compute angle cos(θ_jk) = (v_1j^2 + v_1k^2 - v_jk^2) / 2|v_1j||v_1k|,
    // where |v_jk| = v_jk = sqrt(Er_jk). For a F-state system (F > 2), we need
    // θ_23, θ_24, θ_34, ... , θ_(F-1)F, in total (F-1)(F-2)/2
    // Using the order number (single index) represents a pair of states (jk):
    // θ_index -> θ_jk, 2 <= x < y <= DOFe, index = [(k-2)(k-3)/2 + j-1] - 1
    // Note that first index is 0 in code (thus -1)
    const int theta_count = (DOFe-1)*(DOFe-2)/2; // number of θ_jk
    std::vector<double> cos_theta; // cosθ_jk
    std::vector<double> sin_theta; // sinθ_jk
    // Note a vector cannot be resized to 0 element (when DOFe=2, no θ)
    if (theta_count > 0) {
        cos_theta.resize(theta_count, 0);
        sin_theta.resize(theta_count, 0);
        for (int k = 3; k <= DOFe; ++k)
            for (int j = 2; j < k; ++j) {
                // The index of Er_jk in Er vector is (k-1)(k-2)/2 + j - 1
                int Er_1j = (j-1)*(j-2)/2 + 1 - 1;
                int Er_1k = (k-1)*(k-2)/2 + 1 - 1;
                int Er_jk = (k-1)*(k-2)/2 + j - 1;
                // index of θ_jk in cos/sin_θ_jk vector is (k-2)(k-3)/2 + j-2
                int theta_jk = (k-2)*(k-3)/2 + j - 2;
                // cos(θ_jk) = (v_1j^2 + v_1k^2 - v_jk^2) / 2|v_1j||v_1k|
                // v_jk^2 = (sqrt(Er_jk))^2 = abs(Er_jk)
                cos_theta[theta_jk] = 0.5 * (abs(Er[Er_1j]) + abs(Er[Er_1k]) - abs(Er[Er_jk]))
                    / sqrt(abs(Er[Er_1j]*Er[Er_1k]));
                // sinθ_jk = sqrt(1-cosθ_jk^2) and > 0 since 0 < θ_jk < pi
                sin_theta[theta_jk] = sqrt(1-cos_theta[theta_jk]*cos_theta[theta_jk]);
            }
        // Additional output the value of cosθ, sinθ and angles in degree.
        std::string outfile = TCF_prefix + "_theta.dat";
        FILE* out = CheckFile(outfile);
        fprintf(out, "%s\n","The angles theta in degree are:");
        fprintf(out, "%6s%12s %12s %12s\n", "jk  ", "theta", "cos_theta", "sin_theta");
        for (int k = 3; k <= DOFe; ++k)
            for (int j = 2; j < k; ++j) {
                // index of θ_jk in cos/sin_θ_jk vector is (k-2)(k-3)/2 + j-2
                int theta_jk = (k-2)*(k-3)/2 + j - 2;
                double angle = acos(cos_theta[theta_jk])/pi*180.0; // in degree
                fprintf(out, "%3d%-3d%12.6g %12.6g %12.6g\n", j, k, angle, cos_theta[theta_jk], sin_theta[theta_jk]);
            }
        fclose(out);
    }
    // * 8.3 Compute angle θ_jk'
    // For a F-state system (F > 3), we need θ_34', θ_35', θ_45', ... ,
    // θ_(F-1)F', in total (F-2)(F-3)/2;
    // Using the order number (single index) represents a pair of states (jk):
    // θ_index -> θ_jk, 3 <= x < y <= DOFe, index = [(k-3)(k-4)/2 + j-2] - 1
    // Note that first index is 0 in code (thus -1)
    const int theta_prime_count = (DOFe-2)*(DOFe-3)/2; // number of θ_jk'
    std::vector<double> cos_theta_prime; // cosθ_jk'
    std::vector<double> sin_theta_prime; // sinθ_jk'
    // Note a vector cannot be resized to 0 element (when DOFe<=3, no θ')
    if (theta_prime_count > 0) {
        cos_theta_prime.resize(theta_prime_count, 0);
        sin_theta_prime.resize(theta_prime_count, 0);
        for (int k = 4; k <= DOFe; ++k)
            for (int j = 3; j < k; ++j) {
                double numerator = 0, denumerator = 1;
                // The index of θ_jk in cos/sin_θ_jk vector is (k-2)(k-3)/2 + j-2
                int theta_2j = (j-2)*(j-3)/2 + 2-2;
                int theta_2k = (k-2)*(k-3)/2 + 2-2;
                int theta_jk = (k-2)*(k-3)/2 + j-2;
                denumerator = sin_theta[theta_2j] * sin_theta[theta_2k];
                for (int i = 3; i < j; ++i) {
                    // The index of θ_ik' in cos/sin_θ_ik' vector is (k-3)(k-4)/2 + i-3
                    int theta_prime_ij = (j-3)*(j-4)/2 + i-3;
                    int theta_prime_ik = (k-3)*(k-4)/2 + i-3;
                    denumerator *= sin_theta_prime[theta_prime_ij] * sin_theta_prime[theta_prime_ik];
                }
                numerator = cos_theta[theta_jk] - cos_theta[theta_2j]*cos_theta[theta_2k];
                // numerator has the following temp_sum term when j > 3.
                double temp_sum = 0;
                for (int i = 3; i < j; ++i) {
                    double temp_prod = 1;
                    for (int l = 3; l <= i; ++l) {
                        // The index of θ_lk' in cos/sin_θ_lk' vector is (k-3)(k-4)/2 + l-3
                        int theta_prime_lj = (j-3)*(j-4)/2 + l-3;
                        int theta_prime_lk = (k-3)*(k-4)/2 + l-3;
                        if (l == i)
                            temp_prod *= cos_theta_prime[theta_prime_lj] * cos_theta_prime[theta_prime_lk];
                        else
                            temp_prod *= sin_theta_prime[theta_prime_lj] * sin_theta_prime[theta_prime_lk];
                    }
                    temp_sum += temp_prod;
                }
                numerator -=  sin_theta[theta_2j] * sin_theta[theta_2k] * temp_sum;
                // The index of θ_jk' in cos/sin_θ_jk' vector is (k-3)(k-4)/2 + j-3
                int theta_prime_jk = (k-3)*(k-4)/2 + j-3;
                cos_theta_prime[theta_prime_jk] = numerator / denumerator;
                // sinθ_jk' = sqrt(1-cosθ_jk'^2) and > 0 since 0 < θ_jk' < pi
                sin_theta_prime[theta_prime_jk] = sqrt(1-cos_theta_prime[theta_prime_jk]*cos_theta_prime[theta_prime_jk]);
            }
        // Additional output the value of cosθ', sinθ' and angles in degree.
        std::string outfile = TCF_prefix + "_theta_prime.dat";
        FILE* out = CheckFile(outfile);
        fprintf(out, "%s\n","The angles theta_prime in degree are:");
        fprintf(out, "%6s%12s %12s %12s\n", "jk  ", "theta'", "cos_theta'", "sin_theta'");
        for (int k = 4; k <= DOFe; ++k)
            for (int j = 3; j < k; ++j) {
                // The index of θ_jk' in cos/sin_θ_jk' vector is (k-3)(k-4)/2 + j-3
                int theta_prime_jk = (k-3)*(k-4)/2 + j-3;
                double angle = acos(cos_theta_prime[theta_prime_jk])/pi*180.0; // in degree
                fprintf(out, "%3d%-3d%12.6g %12.6g %12.6g\n", j, k, angle,
                    cos_theta_prime[theta_prime_jk], sin_theta_prime[theta_prime_jk]);
            }
        fclose(out);
    }
    // * 8.4 Compute equilibrium shifts S_XY
    // They are S_12, S_13, S_23,...S_(DOFe-1)(DOFe),..., (order by rows)
    // Using the order number (single index) represents a pair of states (XY):
    // S_index -> S_xy, 1 <= x < y <= DOFe, index = (y-1)(y-2)/2 + x - 1
    // The index of S_XY is the same as U_XY Er_XY. in total (DOFe-1)(DOFe)/2
    // For each S_XY, there are N elements, N is number of normal modes.
    std::vector<std::vector<double>> S(U_count, std::vector<double>(N, 0));
    // Here, position is the non-zero equilibrium position component of PESs,
    // such as Ax, Bx, By, Cx, Cy, Cz, in our paper.
    std::vector<double> position(U_count, 0);
    for (int y = 2; y <= DOFe; ++y)
        for (int x = 1; x < y; ++x) {
            // factor is the factor from theta/theta' in S formula
            // For S_12 only, no factor from theta/theta' and thus it is 1.
            double factor = 1;
            if (y > 2) { // S_13, S_23, S_14, ..., S_(DOFe-1)(DOFe)
                // the first item comes from θ_2y, for S_1Y is cosθ_2y, others are sinθ_2y
                // The index of θ_jk in cos/sin_θ_jk vector is (k-2)(k-3)/2 + j-2
                int theta_2y = (y-2)*(y-3)/2 + 2-2;
                if (x == 1)
                    factor *= cos_theta[theta_2y];
                else
                    factor *= sin_theta[theta_2y];
                // the second item comes from product of sin θ_3y', θ_4y', ..., θ_xy'
                for (int i = 3; i <= x; ++i) {
                    // The index of θ_iy' in cos/sin_θ_iy' vector is (y-3)(y-4)/2 + i-3
                    int theta_prime_iy = (y-3)*(y-4)/2 + i-3;
                    factor *= sin_theta_prime[theta_prime_iy];
                }
                // then multiply cos θ_(x+1)y' only when x < y-1 and x > 1
                // Note for the last elment in a row, i.e., S_(y-1)y, doesn't
                // multiply this cos θ_(x+1)y' item.
                if (x > 1 && x < y-1) {
                    // The index of θ_(x+1)y' in cos/sin_θ_xy' vector is (y-3)(y-4)/2 + x+1-3
                    int theta_prime_x1y = (y-3)*(y-4)/2 + x+1-3;
                    factor *= cos_theta_prime[theta_prime_x1y];
                }
            }
            // Compute S_j^XY based on the a_j, Er_1y, and factor
            // The index for S_XY in S vector is: index = (y-1)(y-2)/2 + x-1
            // which is the same as Er_XY, U_XY.
            int S_xy = (y-1)*(y-2)/2 + x-1;
            int Er_1y = (y-1)*(y-2)/2 + 1-1;
            // position is the Ax, Bx, By, Cx, Cy, Cz,..., in our paper.
            position[S_xy] = sqrt(Er[Er_1y]) * factor;
            for (int j = 0; j < N; ++j)
                S[S_xy][j] = a[j] * position[S_xy];
        }
    // Additional output the equilibrium positions in au of the PESs.
    outfile = TCF_prefix + "_equilibrium_positions.dat";
    out = CheckFile(outfile);
    fprintf(out, "%s\n","The equilibrium positions in a.u. of the PESs are:");
    char c = 'A';
    for (int y = 2; y <= DOFe; ++y) {
        fprintf(out, "%3c  ", c);
        for (int x = 1; x < y; ++x) {
            // The index for S_XY in S vector is: index = (y-1)(y-2)/2 + x-1
            int S_xy = (y-1)*(y-2)/2 + x-1;
            fprintf(out, "%12.6f", position[S_xy]);
            //if (x < y-1)
            //    fprintf(out, ","); // comma separated
        }
        fprintf(out, "\n");
        c++;
    }
    fclose(out);
    // * 9. Assign values of shifts to data member: req from S matrix
    // The meanning of req is the same as S (equilibrium shifts), but req is
    // a full DOFe*(DOFe-1) matrix stored as (1D vector) N*[DOFe*(DOFe-1)].
    // The req matrix is S_XY, 1 <= X <= DOFe-1 and 1 <= Y <= DOFe-1,
    // where, XY is the state index in paper starting from 1. Each element of
    // matrix, S_XY, has N values, i.e., N matrices of req. Then each matrix in
    // req is stored as a vector DOFe*(DOFe-1), that is, (S_11[1], S_21[1],...
    // S_(DOFe-1)DOFe[1]), ..., (S_11[j], ..., S_(DOFe-1)DOFe[j]), j is from 1 to N.
    // The rest values of req that doesn't exist in S, are 0. The rest elements
    // are S_11, S_12, S_1(DOFe-1),  S_22,..., S_2(DOFe-1),.... (the row with
    // Y=1 and the up diagonal of matrix, Y <= X). S is the down diagonal of
    // matrix without the first row (X < Y <= DOFe and Y != 1).
    // The req vector has been resized and initialized to 0.
    for (int y = 2; y <= DOFe; ++y)
        for (int x = 1; x < y; ++x) {
            // The index for S_XY in S vector is: S_xy = (y-1)(y-2)/2 + x -1
            // while, in the req full matrix is: req_xy = (y-1)(DOFe-1)+ x -1
            // (-1 makes req_xy starting from 0).
            int S_xy = (y-1)*(y-2)/2 + x-1;
            int req_xy = (y-1)*(DOFe-1) + x -1;
            for (int j = 0; j < N; ++j)
                req[j*DOFe*(DOFe-1) + req_xy] = S[S_xy][j];
        }
}

double HamiltonianMSH::fw(double w, int j, double dt, std::vector<double>& Cuu) {
    // Compute f(ω_j) = [2Nω_j/(piCuu(0))]*[\int_0^∞ dt Cuu(t)/(ω_jt) sin(ω_jt)] - j + 1/2
    // Note that the unit of them is au or au^2.
    double factor = 2*N*w/pi/Cuu[0]; // factor = 2Nω_j/(piCuu(0))
    // Get f(t) = Cuu(t)/(ω_jt) sin(ω_jt)
    const int size = Cuu.size();
    std::vector<double> ft(size, 0);
    // To avoid dividing by zero, the integrand at t = 0 is replaced by its
    // analytical limit, i.e., ft[0] = Cuu[0]
    ft[0] = Cuu[0];
    for (int s = 1; s < size; ++s) {  // s is step at time t, and t = s*DT
        double wt = w * s * dt;       // Note starting from 1.
        ft[s] = Cuu[s] * sin(wt) / wt;
    }
    // Compute integral \int f(t)dt
    double integral = Integrate(ft, dt);
    // Final result of f(w_j)
    return factor*integral + 0.5 - j;
}

double HamiltonianMSH::Jw(double w, double beta, double dt, std::vector<double>& Cuu) {
    // Compute the continus spectral density J(ω)
    // J(ω) = βω/4 \int_0^∞ dt Cuu(t)/cos(ωt), β = 1/kT
    // Note that the unit of them is au or au^2.
    double factor = 0.25 * beta * w;
    // Compute f(t) = Cuu(t)/cos(ωt)
    const int size = Cuu.size();
    std::vector<double> ft(size, 0);
    for (int s = 0; s < size; ++s) {  // s is step at time t, and t = s*DT
        double wt = w * s * dt;
        ft[s] = Cuu[s] * cos(wt);
    }
    // Compute integral \int f(t)dt
    double integral = Integrate(ft, dt);
    // Final result of J(w)
    return factor*integral;
}

void HamiltonianMSH::loadModelParameters(const std::string& loadfile) {
    std::ifstream inpfile(loadfile);
    if (!inpfile)
        throw std::runtime_error("ERROR: Can't open MSH model parameters file: " + loadfile);
    const int cols = DOFe*(DOFe-1) + 2;
    int lineNumber = 0;
    for (std::string line; getline(inpfile, line); ) {
        lineNumber++; // Let line number start from 1
        if (lineNumber == 1) continue; // First line is headers
        std::vector<std::string> values;
        SplitString(values, line, ','); // csv file
        if (values.size() < cols)
            throw std::runtime_error("ERROR: Too few columns at line " + std::to_string(lineNumber));
        // index (j) is lineNumber-2, starting form 0.
        omega[lineNumber-2] = std::stod(values[1]);
        // Here, the req is 1D vector: S[j][n][m] = req[j*n*m + n*DOFe+m]
        for (int i = 0; i < (cols-2); ++i)
            req[(lineNumber-2)*(cols-2) + i] = std::stod(values[i+2]);
    }
    if (lineNumber != N+1) // The required number of lines is N+1
        throw std::runtime_error("ERROR: Too few lines of data in " + loadfile);
    inpfile.close();
}

void HamiltonianMSH::saveModelParameters(const std::string& savefile) {
    FILE* outfile = CheckFile(savefile);
    // CSV file, fields are separated by one comma (,).
    fprintf(outfile, "%s,%s", "j", "omega_j/au");
    for (int n = 0; n < DOFe; n++)
        for (int m = 0; m < DOFe-1; m++)
                fprintf(outfile, ",%s%d%d", "Sj_", m+1, n+1);
    fprintf(outfile, "\n");
    for (int j = 0; j < N; ++j) {
        fprintf(outfile, "%d,%.16g", j+1, omega[j]);
        for (int i = 0; i < DOFe*(DOFe-1); ++i)
            fprintf(outfile, ",%.16g", req[j*DOFe*(DOFe-1) + i]);
        fprintf(outfile, "\n");
    }
    fclose(outfile);
}

double HamiltonianMSH::getPotentialEnergy(int index) {
    if (index >= DOFe || index < 0)
        throw std::runtime_error("ERROR: State index out of scope.");
    double PE = 0.0;
    // Ref: J. Chem. Phys. 155, 124105 (2021) Eq. 18
    // Here, N is the number of harmonic modes. size is the size of S matrix.
    // row is the index of row of S matrix should be used. DOFn/N = DOFe-1 (model I)
    // or DOFe-2 (model II), is the number of columns of S matrix used for calculation.
    const int size = DOFe*(DOFe-1);
    const int row  = index*(DOFe-1);
    for (int i = 0; i < DOFn/N; ++i)
        for (int j = 0; j < N; ++j) {
            double X = R[i*N+j] - req[j*size + row + i];
            PE += omega[j] * omega[j] * X * X;
        }
    return 0.5 * PE + epsilon[index];
}

void HamiltonianMSH::getForces(int i, int j, std::vector<double>& F) {
    if (i >= DOFe || i < 0 || j >= DOFe || j < 0)
        throw std::runtime_error("ERROR: Called getForces() with wrong state index.");
    F.resize(DOFn, 0);
    if (i == j) { // froces of each state from diagonal Hamiltonian
        const int size = DOFe*(DOFe-1);
        const int row  = i*(DOFe-1);
        for (int s = 0; s < DOFn/N; ++s)
            for (int k = 0; k < N; ++k)
                F[s*N+k] = -omega[k] * omega[k] * (R[s*N+k] - req[k*size + row + s]);
    }
    // forces from off-diagonal Hamiltonian (diabatic coupling)
    else if (!Condon_approximation)
        throw std::runtime_error("ERROR: The force from non-Condon diabatic coupling for MSH model is undefined.");
    else // In the case of Condon_approximation, forces are zero
        std::fill(F.begin(), F.end(), 0);
}