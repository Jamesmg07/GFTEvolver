#pragma once

#include "Model.hpp"

#include "InitialCondition.hpp"
#include "RandomUniform.hpp"

#include "BoundaryCondition.hpp"
#include "Periodic.hpp"
#include "Fixed.hpp"

#include "Analyser.hpp"
#include "OutputFields.hpp"
#include "Energy.hpp"
#include "GaugeCondition.hpp"

#include <chrono>

enum Initial_Condition_Types {ZERO_IC = 0, RANDOM_UNIFORM_IC};
enum Boundary_Condition_Types {UNASSIGNED_BOUNDARY_CONDITION = 0, FIXED, NEUMANN, PERIODIC};

class Lattice
{
private:

    ///////////////////////////////////  Variables  ////////////////////////////////////////

    unsigned nx, ny, nz, nt; // Lattice size
    double dx, dy, dz, dt; // Lattice spacings

    unsigned numScalarFieldComponents, numVectorFieldComponents; // Number of field components
    std::vector<float> scalarFields, vectorFields; // Field arrays

    unsigned stencilSize; // The largest default stencil size.

    // Contains the potential, gradients and wilson loops.
    Model model;

    // Classes to deal with the initial conditions. Child classes will be instantiated based on lattice configuration.
    InitialCondition *scalarInitialConditions, *vectorInitialConditions;

    // Classes to deal with the boundary conditions. Child classes will be instantiated based on lattice configuration.
    // 6 face boundaries, 12 edge boundaries and 8 corner boundaries, so 26 in total.
    std::vector<BoundaryCondition*> boundaryConditions;

    // Classes relating to all of the analyses that should be performed.
    std::vector<Analyser*> analysers;

    // Optional reporting of progression throughout the simulation.
    bool progressReport;
    unsigned reportFrequency; // How often to output the number of timesteps that have finsihed.

    //////////////////////////////////  Initialisers  //////////////////////////////////////

    /*
     * Initialises all variables to placeholder or default values.
     * Pointers are set to nullptr.
     */
    void initVariables();

    /*
     * Loads in the numerical parameters of the lattice (nx, ny, nz, nt, dx, dy, dz, dt) from a config file
     * Also defines the number of scalar and vector (gauge) field components.
     * 
     * @param        string path             The file path to the lattice config file.
     * @param        bool debug              Outputs the parameters read in from the config file if true.
     */
    void configure(const std::string path, const bool debug = false);

    /*
     * Initialise InitialCondition classes based on chosen child class in config file.
     * The chosen child class will be configured based on its own associated config file.
     * 
     * @param        vector<int> initial_condition_type                Contains all IC choices.
     */
    void initInitialCondition(const std::vector<int> &initial_condition_types);

    /*
     * Instantiates the existing 26 boundary-region objects from the six
     * user-assigned face boundary conditions.
     *
     * @param        vector<int> face_boundary_condition_types
     *               Choices for (-x), (+x), (-y), (+y), (-z), (+z).
     */
    void initBoundaryCondition(const std::vector<int> &face_boundary_condition_types);

    /*
     * Allocate the appropriate amount of memory and run a function to generate the initial conditions
     */
    void initFields();

    ////////////////////////////////  Private Functions  /////////////////////////////////////

    /*
     * Check the size of the largest default stencil being used by the gradients and wilson loops.
     */
    unsigned getDefaultStencilSize() const;

    /*
     * Use the stencil size to determine how much of the evolution will be performed in the default manner.
     * The rest need adjustments which will be handled by the boundary conditions classes.
     * 
     * @param         const unsigned &stencil_size            The largest default stencil size.
     * 
     * @return        vector<vector<unsigned>>                Limits of the interiod grid loops.
     */
    std::vector< std::vector<unsigned> > determineResponsibilities(const unsigned &stencil_size) const;


    /*
    * Runs post-evolution location analyses after the new timestep
    * has been completed across the full dynamic grid.
    *
    * @param        unsigned t_now
    * @param        vector<vector<unsigned>> loop_limits
    */
    void postEvolveAnalysis(const unsigned &t_now, const std::vector<std::vector<unsigned>> &loop_limits);


    /*
     * Add on the contributions from the x position (and neighbours) to the running indices array
     * 
     * @param        vector<vector<long long unsigned>> running_indices                 Copy of the running indices to be updated.
     * @param        unsigned &loc                                                      Location along the specified axis.
     * @param        unsigned axis                                                      The specified axis.
     * 
     * @return       vector<vector<long long unsigned>>                                 Next version of running indices.
     */
    std::vector<std::vector<long long unsigned>> updateRunningIndices(std::vector<std::vector<long long unsigned>> running_indices, 
                                                                      const unsigned &loc, const unsigned axis) const;

    /*
     * Creates an array of pointers for the scalar fields at the locations determined by the running_indices.
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
    std::vector<std::vector<const float*>> generateVectorPointers(const std::vector<std::vector<long long unsigned>> &running_indices) const;

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

    /////////////////////////////  Constructors/Destructors  ////////////////////////////////

    Lattice();
    virtual ~Lattice();

    ////////////////////////////////  Public Functions  /////////////////////////////////////

    /*
     * For all assigned analysers, this will run their own initialAnalysis functions.
     */
    void initialAnalysis() const;

    /*
     * Loops through nt timesteps and evolves the fields according to the equations of motion.
     * Optionally outputs simulation data throughout the evolution.
     */
    void evolve();

    /*
     * For all assigned analysers, this will run their own finalAnalysis functions.
     */
    void finalAnalysis() const;

};