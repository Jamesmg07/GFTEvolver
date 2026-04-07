#pragma once

#include <vector>
#include <fstream>
#include <iostream>

class Gradient
{
private:

protected:

    ////////////////////////////////////////////  Variables  ///////////////////////////////////////////////

    // Outermost vector contains the different sizes of stencils which may be used near the boundaries.
    // Different stencils for each dimension are stored in the interior 2D array, so that they may be
    // different in each direction (e.g different grid spacings).
    std::vector< std::vector< std::vector<double> > > stencils;

    // Same as above but this time for calculating the gradient energy density.
    // Assumes that the stencils are the same size.
    std::vector< std::vector< std::vector<double> > > energy_stencils;

    // Store the sizes of each stencil type so that they don't have to be computed during evolution.
    std::vector<unsigned> stencil_sizes;

public:

    //////////////////////////////////////  Constructors/Destructors  ///////////////////////////////////////

    Gradient();
    virtual ~Gradient();

    /////////////////////////////////////////  Public Functions  ////////////////////////////////////////////

    /*
     * Pure virtual function that must be defined in each child class.
     * Intention is to get the size of the default stencil (assumed the same in each direction).
     * It is assumed that the stencil is symmetric, so ,for example, 1 means that the stencil looks at neighbours
     * that are one lattice site away on both sides.
     * 
     * @return    unsigned                The stencil size.
     */
    virtual unsigned getDefaultStencilSize() const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to return the number of constraint equations (due to the fixing of the temporal gauge) that should be satisfied throughout.
     * 
     * @return        unsigned                                             The number of constraint equations.
     */
    virtual unsigned getNumberOfConstraintEquations() const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to calculate the gradient energy term.
     * 
     * @param        vector<vector<float*>> &scalar_pointers                Array of pointers to the scalar fields at the required grid positions.
     * @param        vector<vector<float*>> &vector_pointers                Array of pointers to the vector fields at the required grid positions.
     * 
     * @return       float                                                  The gradient energy (density) at this position.
     */
    virtual float calcGradientEnergy(const std::vector<std::vector<const float*>> &scalar_pointers, 
                                     const std::vector<std::vector<const float*>> &vector_pointers) const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to calculate the kinetic energy term.
     * 
     * @param        float* local_scalar_fields[2]                           Pointers to scalar fields at this location, for both timesteps.
     * 
     * @return       float                                                   The kinetic energy (density) at this position.
     */
    virtual float calcKineticEnergy(const float* const local_scalar_fields[2]) const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to calculate the contributions from the gradients to the constraint equations associated with the temporal gauge fixing.
     * 
     * @param        float* local_scalar_fields[2]                           Pointers to scalar fields at this location, for both timesteps.
     * 
     * @return       vector<float>                                           Contribution to the constraint equations.
     */
    virtual std::vector<float> calcConstraintContributions(const float* const local_scalar_fields[2]) const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to calculate the contribution from the gradient
     * energy term to the scalar equations of motion.
     * 
     * @param        vector<vector<float*>> &scalar_pointers                 Array of pointers to the scalar fields at the required grid positions.
     * @param        vector<vector<float*>> &vector_pointers                 Array of pointers to the vector fields at the required grid positions.
     * 
     * @return       vector<double>                                          Contribution to the scalar equations of motion.
     */
    virtual std::vector<double> calcDerivatives(const std::vector<std::vector<const float*>> &scalar_pointers, 
                                                const std::vector<std::vector<const float*>> &vector_pointers) const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to calculate the contribution from the gradient
     * energy term to the vector (gauge) equations of motion.
     * 
     * @param        vector<vector<float*>> &scalar_pointers                 Array of pointers to the scalar fields at the required grid positions.
     * @param        vector<vector<float*>> &vector_pointers                 Array of pointers to the vector fields at the required grid positions.
     * 
     * @return       vector<double>                                          Contribution to the vector equations of motion.
     */
    virtual std::vector<double> calcCurrents(const std::vector<std::vector<const float*>> &scalar_pointers, 
                                             const std::vector<std::vector<const float*>> &vector_pointers) const = 0;

};



//////////////////////////////////////////////  Null Gradient  //////////////////////////////////////////////////
// Trivial version of Gradient, where all functions either do nothing or return zero(s).

class NullGradient:
    public Gradient
{
private:

    ///////////////////////////////////////////  Variables  /////////////////////////////////////////////////////

    const unsigned &numScalarComponents, &numVectorComponents;

    //////////////////////////////////////////  Initialisers  ///////////////////////////////////////////////////

    void configure(const std::string path, const bool debug = false);

public:

    /////////////////////////////////////  Constructors/Destructors  ////////////////////////////////////////////

    NullGradient(const unsigned &num_scalar_components, const unsigned &num_vector_components);
    virtual ~NullGradient();

    /////////////////////////////////////////  Public Functions  ////////////////////////////////////////////////

    unsigned getDefaultStencilSize() const;
    
    unsigned getNumberOfConstraintEquations() const;

    float calcGradientEnergy(const std::vector<std::vector<const float*>> &scalar_pointers, 
                             const std::vector<std::vector<const float*>> &vector_pointers) const;

    float calcKineticEnergy(const float* const local_scalar_fields[2]) const;

    std::vector<float> calcConstraintContributions(const float* const local_scalar_fields[2]) const;

    std::vector<double> calcDerivatives(const std::vector<std::vector<const float*>> &scalar_pointers, 
                                        const std::vector<std::vector<const float*>> &vector_pointers) const;

    std::vector<double> calcCurrents(const std::vector<std::vector<const float*>> &scalar_pointers, 
                                     const std::vector<std::vector<const float*>> &vector_pointers) const;

};