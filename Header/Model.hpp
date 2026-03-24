#pragma once

#include "Potential.hpp"
#include "SO_N.hpp"

#include "Gradient.hpp"
#include "GlobalGradients.hpp"
#include "StandardModelGradients.hpp"

#include "WilsonLoop.hpp"
#include "StandardModelLoops.hpp"

enum Potential_Types {UNASSIGNED_POTENTIAL = 0, SO_N_POTENTIAL, TWO_HDM_POTENTIAL};
enum Gradient_Types {UNASSIGNED_GRADIENT = 0, GLOBAL_GRADIENT, SM_GRADIENT};
enum Wilson_Loop_Types {UNASSIGNED_WILSON_LOOP = 0, SM_WILSON_LOOP};

class Model
{
private:

    ///////////////////////////////  Variables  //////////////////////////////////////

    // Allows a choice of predefined model types.
    // I would prefer to move away from this method and allow a more general approach
    // but that is yet to be implemented and needs some thought as to how it would work.

    Potential* potential;
    Gradient* gradient;
    WilsonLoop* wilsonLoop;

    // Damping parameters
    unsigned ntDamped;
    float dampingFactor;
    float currentDamping;
    bool evolveGauge;

    // Contributions to the equations of motion (at one position) from the potential, gradients and wilson loops.
    // Scalar equations of motion
    std::vector<double> potentialContributions;
    std::vector<double> derivativeContributions;
    // Vector equations of motion
    std::vector<double> currentContributions;
    std::vector<double> wilsonLoopMagneticContributions;
    std::vector<double> wilsonLoopElectricContributions;

    /////////////////////////////  Initialisers  /////////////////////////////////////

    /*
     * Initialise member variables to default values.
     * It's expected that these values will be overwritten.
     */
    void initVariables();

    /*
     * Initialise potential based on chosen child class.
     * Will configure the chosen potential based on the config file associated with the child class.
     * 
     * @param        int potential_type                     Corresponds to a chosen potential read from Model.cfg.
     * @param        unsigned &num_scalar_components        Reference to the number of scalar field components.
     */
    void initPotential(const int& potential_type, const unsigned &num_scalar_components);

    /*
     * Initialise gradient based on chosen child class.
     * Will configure the chosen gradient based on the config file associated with the child class.
     * 
     * @param        int gradient_type                    Corresponds to a chosen gradient read from Model.cfg
     * @param        unsigned &num_scalar_components      Reference to the number of scalar field components.
     * @param        unsigned &num_vector_components      Reference to the number of vector field components.
     * @param        unsigned &nx, &ny, &nz               References to the number of grid points in each direction.
     * @param        double &dt, &dx, &dy, &dz            References to the timestep size and lattice spacing in each direction.
     */
    void initGradient(
        const int& gradient_type,
        const unsigned &num_scalar_components, const unsigned &num_vector_components,
        const unsigned &nx, const unsigned &ny, const unsigned &nz,
        const double &dt, const double &dx, const double &dy, const double &dz
    );

    /*
     * Initialise wilsonLoop based on chosen child class.
     * Will configure the chosen wilsonLoop based on the config file associated with the child class.
     * 
     * @param        int wilson_loop_type                Corresponds to a chosen potential read from Model.cfg
     * @param        double &dt, &dx, &dy, &dz           References to the timestep size and lattice spacing in each direction.
     */
    void initWilsonLoop(const int& wilson_loop_type, const double &dt, const double &dx, const double &dy, const double &dz);


public:

    /////////////////////////  Constructors/Destructors  /////////////////////////////

    Model();
    virtual ~Model();

    ///////////////////////////  Public Functions  ///////////////////////////////////  

    /*
     * Loads in the model type (from predefined sets) and some parameters relating to the model
     * such as damping, expansion etc.
     * Any parameters of the potential, gradients or wilson loops are left to be loaded by the
     * respective classes.
     * 
     * @param        string path                          The file path to the model config file.
     * @param        unsigned &num_scalar_components      Reference to the number of scalar field components.
     * @param        unsigned &num_vector_components      Reference to the number of vector field components.
     * @param        unsigned &nx, &ny, &nz               References to the number of grid points in each direction.
     * @param        double &dx, &dy, &dz, &dt            References to the lattice spacing in each direction and timestep.
     * @param        bool debug                           Outputs the parameters read in from the config file if true.
     */
    void configure(
        const std::string path, 
        const unsigned &num_scalar_components, const unsigned &num_vector_components,
        const unsigned &nx, const unsigned &ny, const unsigned &nz,
        const double &dx, const double &dy, const double &dz, const double &dt,
        const bool debug = false);

    /*
     * Set up for the energy calculation. Only runs if energy analyser is being used.
     *
     * @param        bool store_energy                Boolean that determines whether energy calculations should be performed during the evolution loop.
     */
    void energyPreparation(const bool store_energy) const;

    /*
     * Updates internal parameters relating to the evolution that affect how evolution works, e.g damping.
     * 
     * @param        unsigned time_iter                Current timestep.
     */
    void update(unsigned time_iter);

    /*
     * Get the largest size of the default stencil (assumed the same in each direction) from the gradient and wilson loop classes.
     * It is assumed that the stencil is symmetric, so ,for example, 1 means that the stencil looks at neighbours
     * that are one lattice site away on both sides.
     * 
     * @return        vector<unsigned>                The stencil size.
     */
    unsigned getDefaultStencilSize() const;

    /*
     * Returns the number of constraint equations that should be satisfied throughout the evolution.
     * This is used to check that the equivalent of Gauss's law remains satisfied.
     * 
     * @return        unsigned                The number of constraint equations.
     *
     */
    unsigned getNumberOfConstraintEquations() const;

    /*
     * Calculate the contribution to the equations of motion from the potential.
     * 
     * @param        float* local_scalar_field                A pointer to the first scalar field component at this position.
     */
    void calcPotentialContributions(const float* const local_scalar_field);

    /*
     * Calculate the potential energy associated with the fields at this position.
     *
     * @param        float* local_scalar_field                A pointer to the first scalar field component at this position.
     * 
     * @return       float                                    The potential energy (density) associated with this position.
     */
    float calcPotentialEnergy(const float* const local_scalar_field) const;

    /*
     * Calculate the contribution to the equations of motion from the gradient term.
     * It contributes to both the scalar and vector equations of motion so this will determine
     * both the derivativeContribution (scalar part) and the currentContribution (vector part).
     * 
     * @param        vector<vector<float*>> &scalar_pointers                Array of pointers to the scalar fields at the required grid positions.
     * @param        vector<vector<float*>> &vector_pointers                Array of pointers to the vector fields at the required grid positions.
     */
    void calcGradientContributions(const std::vector<std::vector<const float*>> &scalar_pointers, const std::vector<std::vector<const float*>> &vector_pointers);

    /*
     * Calculate the gradient energy associated with the fields at this position.
     *
     * @param        vector<vector<float*>> &scalar_pointers                 Array of pointers to the scalar fields at the required grid positions.
     * @param        vector<vector<float*>> &vector_pointers                 Array of pointers to the vector fields at the required grid positions.
     * 
     * @return       float                                                   The gradient energy (density) associated with this position.
     */
    float calcGradientEnergy(const std::vector<std::vector<const float*>> &scalar_pointers, const std::vector<std::vector<const float*>> &vector_pointers) const;

    /*
     * Calculates the kinetic energy associated with the fields at this position.
     * 
     * @param        float* local_scalar_fields[2]                           Pointers to scalar fields at this location, for both timesteps.
     * 
     * @return       float                                                   The kinetic energy (density) at this position.
     */
    float calcKineticEnergy(const float* const local_scalar_fields[2]) const;

    /*
     * Calculate the contribution to the equations of motion from the Yang-Mills terms.
     * 
     * @param        vector<vector<float*>> &vector_pointers                Array of pointers to the vector fields at the required grid positions.
     * @param        float* local_vector_fields[2]                          Pointers to vector fields at this location, for both timesteps.
     */
    void calcYangMillsContributions(const std::vector<std::vector<const float*>> &vector_pointers, const float* const local_vector_fields[2]);

    /*
     * Calculates the magnetic energy associated with the fields at this position.
     *
     * @param        vector<vector<float*>> &vector_pointers                Array of pointers to the vector fields at the required grid positions.
     * 
     * @return       float                                                  The magnetic energy (density) at this position.
     */
    float calcMagneticEnergy(const std::vector<std::vector<const float*>> &vector_pointers) const;

    /*
     * Calculates the electric energy associated with the fields at this position.
     *
     * @param        float* local_vector_fields[2]                           Pointers to the vector fields at this location, for both timesteps.
     * 
     * @return       float                                                   The electric energy (density) at this position.
     */
    float calcElectricEnergy(const float* const local_vector_fields[2]) const;

    /*
     * Calculates the violation of the constraint equations associated with the temporal gauge fixing.
     *
     * @param        unsigned num_equations                                  The number of constraint equations.
     * @param        long long int t_future_index                            Index to add (may be negative) that swaps from present indices to future.
     * @param        float* local_scalar_fields[2]                           Pointers to scalar fields at this location, for both timesteps.
     * @param        vector<vector<float*>> &vector_pointers                 Array of pointers to the vector fields at the required grid positions.
     * 
     * @return       vector<float>                                           Violation of the constraint equations.
     */
    std::vector<float> calcConstraintViolation(const unsigned &num_equations, const long long int &t_future_index,
                                               const float* const local_scalar_fields[2], 
                                               const std::vector<std::vector<const float*>> &vector_pointers) const;

    /*
     * Responsible for calculating the value of the field at the next timestep.
     *
     * @param        float* local_scalar_fields[2]                           Pointers to scalar fields at this location, for both timesteps.
     * @param        float* local_vector_fields[2]                           Pointers to vector fields at this location, for both timesteps.
     * @param        float &dt                                               Reference to the timestep size.
     * @param        unsigned &num_scalar_components                         Reference to the number of scalar field components.
     * @param        unsigned &num_vector_components                         Reference to the number of vector field components.
     */
    void evolve(float* const local_scalar_fields[2], float* const local_vector_fields[2],
                const double &dt, const unsigned &num_scalar_components, const unsigned &num_vector_components);

    
};