#include "Lattice.hpp"

#include <filesystem>
#include <stdexcept>

// Created by Steven Cotterill
// Evolves either global field theories or gauged field theories on the lattice.
// Numerical and system parameters are defined in external config files.

int main()
{

    try 
    {
        // Create the standard output directory and warn before existing files
        // with matching names can be overwritten.
        const std::filesystem::path data_directory(DATA_DIR);

        if (std::filesystem::exists(data_directory))
        {
            if (!std::filesystem::is_empty(data_directory))
            {
                std::cerr
                    << "WARNING: The output directory "
                    << data_directory
                    << " is not empty. Existing output files with matching names may be overwritten."
                    << std::endl;
            }
        }

        std::filesystem::create_directories(data_directory);
        std::filesystem::create_directories(data_directory / "Continual");




        // Construct the lattice and load in from config files.
        Lattice lattice;

        // Run any analyses of the initial state that have been assigned.
        lattice.initialAnalysis();

        // Main part of the simulation, evolves the system until the end of the sim and handles any assigned continual analyses.
        lattice.evolve();

        // Run any analyses of the final state that have been assigned.
        lattice.finalAnalysis();

        return EXIT_SUCCESS;
    }
    catch (const std::exception& exception)
    {
        std::cerr << "Error in " << exception.what() << std::endl;

        return EXIT_FAILURE;
    }
    catch (...)
    {
        std::cerr << "Unknown error occurred" << std::endl;

        return EXIT_FAILURE;
    }

}



