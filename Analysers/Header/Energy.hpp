#pragma once

#include "Analyser.hpp"
#include "Model.hpp"
#include <array>

class Energy
    : public Analyser
{
private:

    /////////////////////////////////////////////////////////  Variables  //////////////////////////////////////////////////////////////////////

    const Model &model;

    // These need to match the number of options provided.
    // Effectively they are hard-coded but in such a way that changing it later is slightly less error prone.
    // They are usually equal because of the set-up of the energy calculations, but global must strictly be less than local.
    static constexpr unsigned numGlobalOptions = 6U;
    static constexpr unsigned numLocalOptions = 6U;

    float energy, potential, kinetic, gradient, magnetic, electric;
    float* const energyPointers[numGlobalOptions];

    std::vector<float> energyDensity, potentialDensity, kineticDensity, gradientDensity, magneticDensity, electricDensity;
    std::vector<float>* const densityPointers[numLocalOptions];

    const double &dx, &dy, &dz;
    int rank, numRanks;
    unsigned nx, ny, nz;      // GLOBAL lattice dimensions
    unsigned globalXStart;    // this rank's offset into the global x-axis

    bool globalOptions[numGlobalOptions], localOptions[numLocalOptions], anyLocalOptions, anyGlobalOptions;
    unsigned globalFrequency, localFrequency, counter;
    bool globalOutput, localOutput;
    long long unsigned ownedSiteBegin, ownedSiteEnd;

    std::string globalQuantitiesPath, localQuantitiesPath;

    bool monoSeparationEnabled;
    std::string monoSeparationPath;
    unsigned monoSeparationFrequency;
    double monoMinimumSeparation;
    bool monoSeparationOutput;

    

    ///////////////////////////////////////////////////////  Initialisers  /////////////////////////////////////////////////////////////////////

    /*
     * Load in chosen options and parameters from the associated config file.
     *
     * @param        string path                Path to the config file.
     * @param        bool debug                 Outputs loaded choices and parameters if true.
     */
    void configure(const std::string path, const bool debug = false);

    /*
     * Initialise any member variables not dealt with in configure.
     *
     * @param        long long unsigned grid_size                Total spatial size of the lattice.
     */
    void initVariables(const long long unsigned grid_size);

    /*
    * Add the MPI rank suffix before the filename extension.
    * For a one-rank run the original path is returned unchanged.
    */
    std::string rankLocalPath(const std::string &path, const int file_rank) const;

    /*
    * Reassemble one local-energy output file from the rank-local pieces.
    */
    void mergeRankLocalOutput(const std::string &path, const std::vector<unsigned long long> &owned_site_counts) const;

public:

    //////////////////////////////////////////////////  Constructors/Destructors  //////////////////////////////////////////////////////////////

    Energy(const Model &model,
    const double &dx,
    const double &dy,
    const double &dz,
    const long long unsigned grid_size,
    const long long unsigned owned_site_begin,
    const long long unsigned owned_site_end,
    const int rank,
    const int num_ranks,
    const unsigned nx,
    const unsigned ny,
    const unsigned nz, 
    const unsigned global_x_start);
    virtual ~Energy();

    //////////////////////////////////////////////////////  Public Functions  //////////////////////////////////////////////////////////////////

    /*
     * Not required.
     */
    void initialAnalysis();

    /*
     * Calculates energy densities at each visited site in the
     * evolved-site diagnostic domain. Fixed support sites are not
     * visited; their full-grid local-array entries are placeholders.
     * 
     * @param        long long unsigned index                              Index for the density arrays
     * @param        float* local_scalar_pointers[2]                       Pointers to the scalar field at current location for both timesteps.
     * @param        vector<vector<float*>> scalar_pointers                Pointers to the scalar fields at grid locations required by the stencil
     *                                                                     (assumed the same as the 2nd derivative locations)
     * @param        float* local_vector_pointers[2]                       Pointers to the vector field at current location for both timesteps.
     * @param        vector<vector<float*>> vector_pointers                Pointers to the vector fields at grid locations required.
     */
    void preEvolveLocationAnalysis(const long long unsigned index,
                                   const float* const local_scalar_pointers[2], const std::vector<std::vector<const float*>> &scalar_pointers,
                                   const float* const local_vector_pointers[2], const std::vector<std::vector<const float*>> &vector_pointers);

    /*
     * Not required.
     * 
     * @param        unsigned t_now                                        Index to determine locations in array that correspond to "now" (other is future)
     * @param        long long unsigned index                              Index for the density arrays
     * @param        float* local_scalar_pointers[2]                       Pointers to the scalar field at current location for both timesteps.
     * @param        vector<vector<float*>> scalar_pointers                Pointers to the scalar fields at grid locations required by the stencil
     *                                                                     (assumed the same as the 2nd derivative locations)
     * @param        float* local_vector_pointers[2]                       Pointers to the vector field at current location for both timesteps.
     * @param        vector<vector<float*>> vector_pointers                Pointers to the vector fields at grid locations required.
     */
    void postEvolveLocationAnalysis(const unsigned &t_now, const long long unsigned index,
                                    const float* const local_scalar_pointers[2], const std::vector<std::vector<const float*>> &scalar_pointers,
                                    const float* const local_vector_pointers[2], const std::vector<std::vector<const float*>> &vector_pointers);

    /*
     * Outputs energy integrated over the evolved-site diagnostic
     * domain and resets the energy accumulators.
     *
     * @param        unsigned time_step                                    The current timestep.
     */
    void timestepAnalysis(const unsigned &time_step);

    /*
     * Not required.
     */
    void finalAnalysis();
};