#pragma once

#include "BoundaryCondition.hpp"
#include "Model.hpp"

class Fixed:
    public BoundaryCondition
{
private:

    /////////////////////////////////////////////////////  Initialisers  ////////////////////////////////////////////////////////////

    /*
     * Load in any parameters relating to the global derivatives from the associated config file.
     * Currently there are none for this case, but an empty function and config file are still provided
     * in case they are needed in the future.
     * 
     * @param        string path                Path to the config file.
     * @param        bool debug                 Outputs loaded parameters if true.
     */
    void configure(const std::string path, const bool debug = false);

    ////////////////////////////////////////////////////  Private Functions  ////////////////////////////////////////////////////////

    /*
     * Add on the contributions from the x position (and neighbours) to the running indices array
     * 
     * @param        vector<vector<long long unsigned>> running_indices                 Copy of the starting running indices.
     * @param        unsigned &loc                                                      Location along the specified axis.
     * @param        unsigned axis                                                      The specified axis.
     * 
     * @return       vector<vector<long long unsigned>>                                 The next level of running indices.
     */
    std::vector<std::vector<long long unsigned>> updateRunningIndices(std::vector<std::vector<long long unsigned>> running_indices, 
                                                                      const unsigned &loc, const unsigned axis) const;

    /*
     * Creates an array of pointers to the locations determined by the running_indices.
     * 
     * @param        vector<vector<long long unsigned>> running_indices                 The running indices to be used to create the array.
     * 
     * @return       vector<vector<float*>>                                             Array of pointers to locations in the lattice.
     */
    std::vector<std::vector<const float*>> generateScalarPointers(const std::vector<std::vector<long long unsigned>> &running_indices) const;

    /*
     * Same as above but also performs a final modification of the running indices in the process of creating the pointers.
     * 
     * @param        vector<vector<long long unsigned>> running_indices                 The running indices to be used to create the array.
     * @param        unsigned &loc                                                      Location along the specified axis.
     * @param        unsigned axis                                                      The specified axis.
     * 
     * @return       vector<vector<float*>>                                             Array of pointers to locations in the lattice.
     */
    std::vector<std::vector<const float*>> generateScalarPointers(const std::vector<std::vector<long long unsigned>> &running_indices,
                                                                   const unsigned &loc, const unsigned axis) const;

    /*
     * Same as above but this time for the vector fields.
     * WARNING: THIS PROBABLY WON'T ADAPT VERY WELL TO DIFFERENT SIZED STENCILS FOR THE WILSON LOOPS
     * 
     * @param        vector<vector<long long unsigned>> running_indices                 The running indices to be used to create the array.
     * 
     * @return       vector<vector<float*>>                                             Array of pointers to locations in the lattice.
     */
    std::vector<std::vector<const float*>> generateVectorPointers(const std::vector<std::vector<long long unsigned>> &running_indices,
                                                                  const unsigned &x_loc, const unsigned &y_loc, const unsigned &z_loc) const;

    /*
     * Returns a pointer to the first scalar field component at the specified location.
     *
     * @param        unsigned &t_step               Specifies which timestep is current.
     * @param        unsigned &x_loc                Specifies x location.
     * @param        unsigned &y_loc                Specifies y location.
     * @param        unsigned &z_loc                Specifies z location.
     * 
     * @return       float*                         Pointer to the first field component at that location.
     */
    float* getScalarFieldPointer(const unsigned &t_step, const unsigned &x_loc, const unsigned &y_loc, const unsigned &z_loc);

    /*
     * Returns a pointer to the first vector field component at the specified location.
     *
     * @param        unsigned &t_step               Specifies which timestep is current.
     * @param        unsigned &x_loc                Specifies x location.
     * @param        unsigned &y_loc                Specifies y location.
     * @param        unsigned &z_loc                Specifies z location.
     * 
     * @return       float*                         Pointer to the first field component at that location.
     */
    float* getVectorFieldPointer(const unsigned &t_step, const unsigned &x_loc, const unsigned &y_loc, const unsigned &z_loc);

public:

    ///////////////////////////////////////////////  Constructors/Destructors  //////////////////////////////////////////////////////

    Fixed(std::vector<float> &scalar_fields, std::vector<float> &vector_fields, std::vector<Analyser*> &analysers,
             const unsigned &num_scalar_components, const unsigned &num_vector_components,
             Model &model, const unsigned &nx, const unsigned &ny, const unsigned &nz, const double &dt,
             const std::vector<int> bound_vector);
    virtual ~Fixed();

    ///////////////////////////////////////////////////  Public Functions  /////////////////////////////////////////////////////////

    /*
    * Leaves sites assigned to this fixed boundary unchanged. These
    * allocated sites form a frozen support shell for neighbouring
    * evolved stencils.
    *
    * Both initial time buffers must contain identical values on this shell.
    *
    * @param        unsigned t_now          Specifies the current buffer.
    * @param        unsigned stencil_size   Size of the largest default stencil.
    */
    void evolve(const unsigned &t_now, const unsigned &stencil_size);

    /*
    * Performs no post-evolution analysis because fixed support sites
    * lie outside the evolved-site diagnostic domain.
    *
    * @param        unsigned t_now          Specifies the current buffer.
    * @param        unsigned stencil_size   Size of the largest default stencil.
    */
    void postEvolveAnalysis(const unsigned &t_now, const unsigned &stencil_size);

};