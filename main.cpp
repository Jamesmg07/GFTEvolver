#include "Lattice.hpp"

#include <filesystem>
#include <stdexcept>
#include <iostream>

#ifdef GFT_ENABLE_MPI
#include <mpi.h>
#endif

// Created by Steven Cotterill
// Evolves either global field theories or gauged field theories on the lattice.
// Numerical and system parameters are defined in external config files.

int main(int argc, char* argv[])
{

    int rank = 0;
    int numRanks = 1;

    #ifdef GFT_ENABLE_MPI
    bool mpiInitialised = false;

    if (MPI_Init(&argc, &argv) != MPI_SUCCESS)
    {
        std::cerr << "Failed to initialise MPI." << std::endl;
        return EXIT_FAILURE;
    }
    mpiInitialised = true;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

    #else
        (void)argc;
        (void)argv;
    #endif



    try 
    {
        if (rank == 0){
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
        }

        #ifdef GFT_ENABLE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
        #endif

        {
            // Existing constructors print their configuration through std::cout.
            // Suppress this duplicate start-up output on worker ranks.
            #ifdef GFT_ENABLE_MPI
            if (rank != 0)
                std::cout.setstate(std::ios_base::failbit);
            #endif

            // Construct the lattice and load in from config files.
            Lattice lattice(rank, numRanks);

            #ifdef GFT_ENABLE_MPI
            if (rank != 0)
                std::cout.clear();
            #endif

            // Run any analyses of the initial state that have been assigned.
            lattice.initialAnalysis();

            // Main part of the simulation, evolves the system until the end of the sim and handles any assigned continual analyses.
            lattice.evolve();

            // Run any analyses of the final state that have been assigned.
            lattice.finalAnalysis();
        }

        #ifdef GFT_ENABLE_MPI
        MPI_Finalize();
        mpiInitialised = false;
        #endif

        return EXIT_SUCCESS;
    }
    catch (const std::exception& exception)
    {
        if (numRanks > 1)
        {
            std::cerr
                << "Error on MPI rank "
                << rank
                << " in "
                << exception.what()
                << std::endl;
        }
        else
        {
            std::cerr
                << "Error in "
                << exception.what()
                << std::endl;
        }

        #ifdef GFT_ENABLE_MPI
        if (mpiInitialised)
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        #endif

        return EXIT_FAILURE;
    }
    catch (...)
    {
        if (numRanks > 1)
        {
            std::cerr
                << "Unknown error occurred on MPI rank "
                << rank
                << std::endl;
        }
        else
        {
            std::cerr
                << "Unknown error occurred"
                << std::endl;
        }

        #ifdef GFT_ENABLE_MPI
        if (mpiInitialised)
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        #endif

        return EXIT_FAILURE;
    }

}



