#pragma once

#include "Analyser.hpp"
#include "Model.hpp"

class GaugeCondition
    : public Analyser
{
private:

    /////////////////////////////////////////////////////////  Variables  //////////////////////////////////////////////////////////////////////

    const Model &model;

    // These need to match the number of options provided.
    // Effectively they are hard-coded but in such a way that changing it later is slightly less error prone.
    static constexpr unsigned numGlobalOptions = 2U;
    static constexpr unsigned numLocalOptions = 1U;

    unsigned numEquations;

    std::vector<float> integratedAbsViolation, maxAbsViolation;
    std::vector<float>* const globalPointers[numGlobalOptions];

    std::vector<std::vector<float>> localViolation;
    std::vector<std::vector<float>>* const localPointers[numLocalOptions];

    const double &dx, &dy, &dz;
    const long long unsigned gridSize;
    const unsigned &numVectorComponents;
    int rank, numRanks;

    bool globalOptions[numGlobalOptions], localOptions[numLocalOptions], anyLocalOptions, anyGlobalOptions;
    unsigned globalFrequency, localFrequency, counter;
    bool globalOutput, localOutput;
    long long unsigned ownedSiteBegin, ownedSiteEnd;

    std::string globalQuantitiesPath, localQuantitiesPath;

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
    void initVariables(const long long unsigned &grid_size);

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

    GaugeCondition(const Model &model,
               const double &dx,
               const double &dy,
               const double &dz,
               const long long unsigned grid_size,
               const unsigned &num_vector_components,
               const long long unsigned owned_site_begin,
               const long long unsigned owned_site_end,
               const int rank,
               const int num_ranks);
    virtual ~GaugeCondition();

    //////////////////////////////////////////////////////  Public Functions  //////////////////////////////////////////////////////////////////

    /*
     * Not required.
     */
    void initialAnalysis();

    /*
     * Not required.
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
     * Calculates the gauge-constraint violation at each visited site
     * in the evolved-site diagnostic domain. Fixed support sites are
     * excluded; their full-grid local-array entries are placeholders.
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
    * Outputs integrated and maximum absolute violations over the
    * evolved-site diagnostic domain and resets the accumulators.
    *
    * @param        unsigned time_step       The current timestep.
    */
    void timestepAnalysis(const unsigned &time_step);

    /*
     * Not required.
     */
    void finalAnalysis();
};