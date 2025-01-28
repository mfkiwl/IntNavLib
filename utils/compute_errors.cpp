#include <iostream>
#include <fstream>
#include <stdexcept>
#include "intnavlib/intnavlib.h"

using namespace intnavlib;

int main(int argc, char** argv) {
    if (argc != 3) {
        throw std::runtime_error("Usage: ./compute_errors <ground_truth_profile> <estimated_profile>");
    }

    // Input filenames
    std::string ground_truth_file = argv[1];
    std::string estimated_file = argv[2];

    // Initialize readers for ground truth and estimated profiles
    MotionProfileReader ground_truth_reader(ground_truth_file);
    MotionProfileReader estimated_reader(estimated_file);

    // Output errors filename
    std::string errors_output_file = estimated_file + "_errors.csv";
    ErrorsWriter errors_writer(errors_output_file);

    // Navigation solutions
    NavSolutionNed ground_truth_nav;
    NavSolutionNed estimated_nav;
    
    // Throw first timestamp
    ground_truth_reader.readNextRow(ground_truth_nav);

    // Loop through the profiles and compute errors
    while (ground_truth_reader.readNextRow(ground_truth_nav) && 
           estimated_reader.readNextRow(estimated_nav)) {
        // Compute errors in NED
        // std::cout << ground_truth_nav.time << " " << estimated_nav.time << std::endl;
        ErrorsNed errors = calculateErrorsNed(ground_truth_nav, estimated_nav);
        // Write errors to the output file
        errors_writer.writeNextRow(errors);
    }

    std::cout << "Errors computation completed. Results saved to " << errors_output_file << std::endl;

    return 0;
}
