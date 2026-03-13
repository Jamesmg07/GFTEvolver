#pragma once

#include "Analyser.hpp"
#include "Model.hpp"

class Energy
    : public Analyser
{
private:

    /////////////////////////////////////////////////////////  Variables  //////////////////////////////////////////////////////////////////////

    const Model &model;

    float energy, potential, kinetic, gradient, magnetic, electric;
    float* const energy_pointers[6];

    std::vector<float> energy_density, potential_density, kinetic_density, gradient_density, magnetic_density, electric_density;
    std::vector<float>* const density_pointers[6];

    const double &dx, &dy, &dz;

    bool globalOptions[6], localOptions[6], anyLocalOptions, anyGlobalOptions;
    unsigned globalFrequency, localFrequency, counter;
    bool globalOutput, localOutput;

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
    void initVariables(const long long unsigned grid_size);

public:

    //////////////////////////////////////////////////  Constructors/Destructors  //////////////////////////////////////////////////////////////

    Energy(const Model &model, const double &dx, const double &dy, const double &dz, const long long unsigned grid_size);
    virtual ~Energy();

    //////////////////////////////////////////////////////  Public Functions  //////////////////////////////////////////////////////////////////

    /*
     * Not required.
     */
    void initialAnalysis();

    /*
     * Calculates the energy density at each location.
     * For now it just adds this to a total but could add the ability to regularly output the energy density.
     * 
     * @param        long long unsigned index                              Index for the density arrays
     * @param        float* local_scalar_pointers[2]                       Pointers to the scalar field at current location for both timesteps.
     * @param        vector<vector<float*>> scalar_pointers                Pointers to the scalar fields at grid locations required by the stencil
     *                                                                      (assumed the same as the 2nd derivative locations)
     * @param        float* local_vector_pointers[2]                       Pointers to the vector field at current location for both timesteps.
     * @param        vector<vector<float*>> vector_pointers                Pointers to the vector fields at grid locations required.
     */
    void locationAnalysis(const long long unsigned index,
                          const float* const local_scalar_pointers[2], const std::vector<std::vector<const float*>> &scalar_pointers,
                          const float* const local_vector_pointers[2], const std::vector<std::vector<const float*>> &vector_pointers);

    /*
     * Output the total energy to file and reset the total energy.
     */
    void timestepAnalysis(const unsigned &time_step);

    /*
     * Not required.
     */
    void finalAnalysis();
};