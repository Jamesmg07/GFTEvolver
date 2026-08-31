#!/bin/bash

#SBATCH --job-name=GFTEvolver
#SBATCH --output=/home/jmg/fork_rep/GFTEvolver/GFTEvolver_%j.out
#SBATCH --error=/home/jmg/fork_rep/GFTEvolver/GFTEvolver_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=400:00:00
#SBATCH --mem=120000

# Load required software
module purge
module load gcc/12.5.0
module load openmpi4.1.8-gcc12.5.0
module load cmake/4.1.2

# Add my user-installed Ninja to PATH
export PATH="$HOME/.local/bin:$PATH"

# Go to GFTEvolver
cd /home/jmg/fork_rep/GFTEvolver

# Configure
cmake -S . -B build-mpi -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DGFT_ENABLE_MPI=ON \
    -DCMAKE_C_COMPILER="$(which gcc)" \
    -DCMAKE_CXX_COMPILER="$(which g++)"

# Build
cmake --build build-mpi

# Run with the number of MPI tasks requested from SLURM
mpirun -np "$SLURM_NTASKS" ./build-mpi/GFTEvolver