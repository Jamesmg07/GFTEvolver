#pragma once

#include <vector>
#include <fstream>
#include <iostream>

class Analyser
{
private:

    ////////////////////////////////////////////  Initialisers  ///////////////////////////////////////////////////

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to load in any option choices or parameters from the config file
     * associated with the child class.
     * 
     * @param        string path                Path to the config file.
     * @param        bool debug                 Outputs loaded choices and parameters if true.
     */
    virtual void configure(const std::string path, const bool debug = false) = 0;

public:

    //////////////////////////////////////  Constructors/Destructors  /////////////////////////////////////////////

    Analyser();
    virtual ~Analyser();

    //////////////////////////////////////////  Public Functions  /////////////////////////////////////////////////

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to do any processing of data that is required before the evolution of the fields (after initial conditions).
     * If an Analyser does not need to do that, the function must still be defined but it can be an empty
     * function.
     */
    virtual void initialAnalysis() = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to do processing of data that is required at every location in the grid, 
     * and occurs just before the fields are evolved.
     * If an Analyser does not need to do that, the function must still be defined but it can be an empty
     * function.
     * 
     * @param        long long unsigned index                              Index for the density arrays
     * @param        float* local_scalar_pointers[2]                       Pointers to the scalar field at current location for both timesteps.
     * @param        vector<vector<float*>> scalar_pointers                Pointers to the scalar field at grid locations required by the stencil.
     *                                                                     (assumed the same as the 2nd derivative locations)
     * @param        float* local_vector_pointers[2]                       Pointers to the vector field at current location for both timesteps.
     * @param        vector<vector<float*>> vector_pointers                Pointers to the vector field at grid locations required by the stencil.
     */
    virtual void preEvolveLocationAnalysis(const long long unsigned index, 
                                           const float* const local_scalar_pointers[2], const std::vector<std::vector<const float*>> &scalar_pointers,
                                           const float* const local_vector_pointers[2], const std::vector<std::vector<const float*>> &vector_pointers  ) = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to do processing of data that is required at every location in the grid, 
     * and occurs just after the fields are evolved.
     * If an Analyser does not need to do that, the function must still be defined but it can be an empty
     * function.
     * 
     * @param        unsigned t_now                                        Index to determine locations in array that correspond to "now" (other is future)
     * @param        long long unsigned index                              Index for the density arrays
     * @param        float* local_scalar_pointers[2]                       Pointers to the scalar field at current location for both timesteps.
     * @param        vector<vector<float*>> scalar_pointers                Pointers to the scalar field at grid locations required by the stencil.
     *                                                                     (assumed the same as the 2nd derivative locations)
     * @param        float* local_vector_pointers[2]                       Pointers to the vector field at current location for both timesteps.
     * @param        vector<vector<float*>> vector_pointers                Pointers to the vector field at grid locations required by the stencil.
     */
    virtual void postEvolveLocationAnalysis(const unsigned &t_now, const long long unsigned index, 
                                            const float* const local_scalar_pointers[2], const std::vector<std::vector<const float*>> &scalar_pointers,
                                            const float* const local_vector_pointers[2], const std::vector<std::vector<const float*>> &vector_pointers  ) = 0;

    /*
     * Pure virtual function that must defined in each child class.
     * Intended use is to do processing of data that is required during the evolution of the fields,
     * but only once - not at every location in the grid.
     * Runs after the fields at each grid site have been evolved forward one timestep.
     * If an Analyser does not need to do that, the function must still be defined but it can be an empty
     * function.
     * 
     * @param        unsigned &time_step                The current timestep.
     */
    virtual void timestepAnalysis(const unsigned &time_step) = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to do any processing of data that is required after the evolution has completed.
     * If an Analyser does not need to do that, the function must still be defined but it can be an empty
     * function.
     */
    virtual void finalAnalysis() = 0;
};