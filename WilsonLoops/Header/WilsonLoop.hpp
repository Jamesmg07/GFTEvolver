#pragma once

#include <vector>
#include <fstream>
#include <iostream>

class WilsonLoop
{
private:

public:

    //////////////////////////////////////////  Constructors/Destructors  /////////////////////////////////////////

    WilsonLoop();
    virtual ~WilsonLoop();

    //////////////////////////////////////////////  Public Functions  /////////////////////////////////////////////

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to return whether the generator representation is being used for the gauge fields.
     * This is written with SU(2) in mind, where it can be useful to use the quaternion representation instead.
     * May need to be expanded later if there end up being more than two convenient representations for some group elements.
     * 
     * @return        bool                                                  True if using the generator representation.
     */
    virtual bool isUsingGeneratorRepresentation() const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to return the relevant squares of the gauge couplings, for each component.
     * 
     * @param        unsigned comp_iter                                     Integer that decides which component of the gauge field.
     * 
     * @return       float                                                  Square of the relevant gauge coupling.
     */
    virtual float getSqrCouplings(const unsigned comp_iter) const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to return the number of evolution equations for the vector fields.
     * This can differ from the number of vector field components if there is some redundancy in the description of the gauge fields.
     * It would use less memory to eliminate this redundancy but this may come at the cost of numerical accuracy so that may not always
     * be the preferred approach.
     * 
     * @return        unsigned                                              The number of evolution equations.
     */
    virtual unsigned getNumberOfEvolutionEquations() const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to return the number of constraint equations (due to the fixing of the temporal gauge) that should be satisfied throughout.
     * 
     * @return        unsigned                                              The number of constraint equations.
     */
    virtual unsigned getNumberOfConstraintEquations() const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is just to set up for the energy calculation. Only runs if energy analyser is being used.
     * 
     * @param        bool store_energy                                      Boolean that determines if energy calculations are performed during evolution.
     */
    virtual void energyPreparation(const bool store_energy) = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to calculate the energy density coming from the spatial part of the Yang-Mills term(s).
     * The contribution from the time components will be treated separetely.
     * 
     * @param        vector<vector<float*>> &vector_pointers                Array of pointers to the vector fields at the required grid positions.
     * 
     * @return       float                                                  Energy density at this position from the (spatial part of the) Yang-Mills term.
     */
    virtual float calcMagneticEnergy(const std::vector<std::vector<const float*>> &vector_pointers) const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to calculate the energy density coming the temporal part of the Yang-Mills term(s).
     * 
     * @param       float* local_vector_fields[2]                           Pointers to the vector fields at this location, for both timesteps.
     * 
     * @return      float                                                   Energy density at this position from the (temporal part of the) Yang-Mills term.
     */
    virtual float calcElectricEnergy(const float* const local_vector_fields[2]) const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to calculate the contributions from the wilson loops to the constraint equations associated with the temporal gauge fixing.
     * 
     * @param        long long int t_future_index                            Index to add (may be negative) that swaps from present indices to future.
     * @param        vector<vector<float*>> &vector_pointers                 Array of pointers to the vector fields at the required grid positions.
     * 
     * @return       vector<float>                                           Contribution to the constraint equations.
     */
    virtual std::vector<float> calcConstraintContributions(const long long int &t_future_index, 
                                                           const std::vector<std::vector<const float*>> &vector_pointers) const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to calculate the contribution to the equations of motion coming from the spatial part of the Yang-Mills term.
     * The contribution from the time components will be treated separately.
     * 
     * @param        vector<vector<float*>> &vector_pointers                Array of pointers to the vector fields at the required grid positions.
     * 
     * @return       vector<double>                                         Contribution to the vector equations of motion.
     */
    virtual std::vector<double> calcMagneticContributions(const std::vector<std::vector<const float*>> &vector_pointers) const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to calculate the contribution to the equations of motion coming from half of the temporal part of the Yang-Mills term,
     * the half that depends upon the current timestep and the previous one.
     * The other half will be used to evolve the fields and dealt with separately.
     * 
     * @param        float* local_vector_fields[2]                          Pointers to the vector fields at this location, for both timesteps.
     * 
     * @return       vector<double>                                         Contribution to the vector equations of motion.
     */
    virtual std::vector<double> calcElectricContributions(const float* const local_vector_fields[2]) const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to determine the group element at the next timestep, given that the Im part of the trace of \sigma^a U_i(t+dt)U_i^\dagger(t) = RHS_i^a.
     * It is assumed that, for sufficiently small dt, U_i(t+dt)U_i^\dagger(t) is close to the identity which resolves the remaining ambiguity.
     * Finally, multiply on the right by U_i(t) to get the group element at the next timestep.
     * 
     * @param        float* local_vector_fields[2]                          Pointers to the vector fields at this location, for both timesteps.
     * @param        vector<float> equation_RHS                             Array containing the right-hand side of the equation, for all components.
     */
    virtual void evolve(float* const local_vector_fields[2], std::vector<double> equation_RHS) = 0;

    /*
    * Pure virtual function that must be defined in each child class.
    * Intention is to get the largest spatial distance from the local lattice site required by the wilson loops.
    * A value of 1 means that only the local site and nearest-neighbour sites are required.
    *
    * @return        unsigned                                              The required spatial stencil size.
    */
    virtual unsigned getDefaultStencilSize() const = 0;

};



/////////////////////////////////////////////  Null Wilson Loop  /////////////////////////////////////////////////////////
// Trivial version of WilsonLoop, where all functions either do nothing or return zero(s).

class NullWilsonLoop:
    public WilsonLoop
{
private:

    ////////////////////////////////////////////  Variables  ///////////////////////////////////////////////

    const unsigned &numVectorComponents;

public:

    ////////////////////////////////////  Constructors/Destructors  ////////////////////////////////////////

    NullWilsonLoop(const unsigned &num_vector_components);
    virtual ~NullWilsonLoop();

    ////////////////////////////////////////  Public Functions  ////////////////////////////////////////////

    bool isUsingGeneratorRepresentation() const;

    float getSqrCouplings(const unsigned comp_iter) const;

    unsigned getNumberOfEvolutionEquations() const;

    unsigned getNumberOfConstraintEquations() const;

    void energyPreparation(const bool store_energy);

    float calcMagneticEnergy(const std::vector<std::vector<const float*>> &vector_pointers) const;

    float calcElectricEnergy(const float* const local_vector_fields[2]) const;

    std::vector<float> calcConstraintContributions(const long long int &t_future_index, 
                                                           const std::vector<std::vector<const float*>> &vector_pointers) const;

    std::vector<double> calcMagneticContributions(const std::vector<std::vector<const float*>> &vector_pointers) const;

    std::vector<double> calcElectricContributions(const float* const local_vector_fields[2]) const;

    void evolve(float* const local_vector_fields[2], std::vector<double> equation_RHS);

    unsigned getDefaultStencilSize() const;

};