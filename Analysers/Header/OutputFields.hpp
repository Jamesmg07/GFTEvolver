#pragma once

#include "Analyser.hpp"

class OutputFields:
    public Analyser
{
private:

    ////////////////////////////////////////////////////  Variables  ////////////////////////////////////////////////////////

    bool outputBothTimesteps, outputInitial, outputFinal, outputContinual;
    unsigned outputFrequency, completedTimesteps;

    std::string initialAnalysisPath, continualAnalysisPath, finalAnalysisPath;

    const std::vector<float> &scalarFields, &vectorFields;
    const unsigned &numScalarComponents, &numVectorComponents;

    // Rank-local storage geometry. The owned x-slab is a contiguous range in
    // flattened [x][y][z] site ordering; halos lie outside this range.
    long long unsigned storageVolume, ownedSiteBegin, ownedSiteEnd;
    int rank, numRanks;

    //////////////////////////////////////////////////  Initialisers  ///////////////////////////////////////////////////////

    /*
     * Loads in chosen options and parameters from the associated config file.
     * 
     * @param        string path                Path to the config file.
     * @param        bool debug                 Outputs loaded choices and parameters if true.
     */
    void configure(const std::string path, const bool debug = false);

    /*
    * Add the requested rank suffix before the filename extension for
    * multi-rank output. The one-rank filename is left unchanged.
    */
    std::string rankLocalPath(const std::string &path, const int file_rank) const;

    /*
    * Reassemble one logical field output from its rank-local pieces.
    * Rank files are removed only after the merged file is complete.
    */
    void mergeRankOutput(
        const std::string &path,
        const std::vector<unsigned long long> &owned_site_counts) const;

    /*
    * Outputs the requested field configuration.
    *
    * @param        string path                Output file path.
    * @param        unsigned time_step         Number of completed field updates.
    */
    void outputFields(const std::string &path,
                    const unsigned &time_step) const;


public:

    ////////////////////////////////////////////  Constructors/Destructors  /////////////////////////////////////////////////

    OutputFields(const std::vector<float> &scalar_fields,
             const unsigned &num_scalar_components,
             const std::vector<float> &vector_fields,
             const unsigned &num_vector_components,
             const long long unsigned storage_volume,
             const long long unsigned owned_site_begin,
             const long long unsigned owned_site_end,
             const int rank,
             const int num_ranks);
    virtual ~OutputFields();

    ////////////////////////////////////////////////  Public Functions  /////////////////////////////////////////////////////

    /*
     * Outputs the initial field configuration.
     */
    void initialAnalysis();

    /*
     * Not required.
     *
     * @param        long long unsigned index                              Index for the density arrays
     * @param        float* local_scalar_pointers[2]                       Pointers to the scalar field at current location for both timesteps.
     * @param        vector<vector<float*>> scalar_pointers                Pointers to the scalar fields at the grid locations required by the stencil
     *                                                                     (assumed the same as the 2nd derivative locations)
     * @param        float* local_vector_pointers[2]                       Pointers to the vector fields at current location for both timesteps.
     * @param        vector<vector<float*>> vector_pointers                Pointers to the vector fields at the grid locations required.
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
     * @param        vector<vector<float*>> scalar_pointers                Pointers to the scalar fields at the grid locations required by the stencil
     *                                                                     (assumed the same as the 2nd derivative locations)
     * @param        float* local_vector_pointers[2]                       Pointers to the vector fields at current location for both timesteps.
     * @param        vector<vector<float*>> vector_pointers                Pointers to the vector fields at the grid locations required.
     */
    void postEvolveLocationAnalysis(const unsigned &t_now, const long long unsigned index,
                                    const float* const local_scalar_pointers[2], const std::vector<std::vector<const float*>> &scalar_pointers,
                                    const float* const local_vector_pointers[2], const std::vector<std::vector<const float*>> &vector_pointers);

    /*
     * Outputs the field every timestep.
     *
     * @param        unsigned &time_step                The current timestep
     */
    void timestepAnalysis(const unsigned &time_step);

    /*
     * Outputs the final field configuration.
     */
    void finalAnalysis();

};