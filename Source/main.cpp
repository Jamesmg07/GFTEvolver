#include "Lattice.hpp"

// Created by Steven Cotterill
// Evolves either global field theories or gauged field theories on the lattice.
// Numerical and system parameters are defined in external config files.

int main()
{
    // Construct the lattice and load in from config files.
    Lattice lattice;

    // Run any analyses of the initial state that have been assigned.
    lattice.initialAnalysis();

    // Main part of the simulation, evolves the system until the end of the sim and handles any assigned continual analyses.
    lattice.evolve();

    // Run any analyses of the final state that have been assigned.
    lattice.finalAnalysis();

    return 0;

}



