/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is main function of the Quantum Classical Dynamics (QCDyn) program.   *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Nov. 09, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "Simulation.h"

int main(int argc, char* argv[]) {
    // Allways enclose all calls with a try/catch block to make sure that
    // usage and runtime errors are caught and reported.
    try {
        // Start to run main program and record current time.
        std::cout << "Start to run QCDyn program at " << CurrentTime() << ".\n" << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        // Create a static Parameters object.
        Parameters& param = Parameters::get();
        // Create a Simulation object and run.
        Simulation simulation(param);
        simulation.init(argc, argv);
        simulation.run();
        // Termination of program.
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsedTime = end - start;
        std::cout << "The total elapsed time is " << std::fixed << std::setprecision(3) <<
            elapsedTime.count()  << " seconds.\n";
        std::cout << "Normal termination of QCDyn program at " << CurrentTime() << "." << std::endl;
        return 0;
    }
    // Catch and report usage and runtime errors detected by program and fail.
    catch(const std::exception& e) {
        std::cout << "EXCEPTION: " << e.what() << std::endl;
        return 1;
    }
}