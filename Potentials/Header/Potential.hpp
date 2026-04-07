#pragma once

#include <vector>
#include <fstream>
#include <iostream>

class Potential
{
private:

public:

    ///////////////////////////////////////  Constructors/Destructors  ///////////////////////////////////////

    Potential();
    virtual ~Potential();

    //////////////////////////////////////////  Public Functions  ////////////////////////////////////////////

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use to to calculate the potential energy density at this location.
     * Results is to be fed into output for storage until the calculation has been
     * performed for the whole grid.
     * 
     * @param        float* field                Pointer to the scalar fields at this position.
     */
    virtual float calcPotentialEnergy(const float* field) const = 0;

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use to to calculate the contribution from the potential to each equation of motion.
     * 
     * @param        float* field                Pointer to the scalar fields at this position.
     * 
     * @return       vector<double>              Array containing the contribution to each equation of motion.  
     */
    virtual std::vector<double> calcPotentialDerivatives(const float* field) = 0;


};



///////////////////////////////////////////  Null Potential  ////////////////////////////////////////////////////
// Trivial version of Potential, where all functions either return zero(s) or do nothing.

class NullPotential:
    public Potential
{
private:

    //////////////////////////////////////////  Variables  ///////////////////////////////////////////////

    const unsigned &numScalarComponents;

public:

    ////////////////////////////////////  Constructors/Destructors  //////////////////////////////////////

    NullPotential(const unsigned &num_scalar_components);
    virtual ~NullPotential();

    ///////////////////////////////////////  Public Functions  ///////////////////////////////////////////

    float calcPotentialEnergy(const float* field) const;

    std::vector<double> calcPotentialDerivatives(const float* field);
};