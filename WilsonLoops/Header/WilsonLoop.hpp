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
     * Intended use it to return the relevant squares of the gauge couplings, for each component.
     * 
     * @param        unsigned comp_iter                                     Integer that decides which component of the gauge field.
     * 
     * @return       float                                                  Square of the relevant gauge coupling.
     */
    virtual float getSqrCouplings(const unsigned comp_iter) const = 0;

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
     * Intended use is to calculate the contribution to the equations of motion coming from the spatial part of the Yang-Mills term.
     * The contribution from the time components will be treated separately.
     * 
     * @param        vector<vector<float*>> &vector_pointers                Array of pointers to the vector fields at the required grid positions.
     * 
     * @return       vector<double>                                         Contribution to the vector equations of motion.
     */
    virtual std::vector<double> calcMagneticContributions(const std::vector<std::vector<const float*>> &vector_pointers, long long unsigned index) const = 0;

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

};