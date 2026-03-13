#pragma once

#include <vector>
#include <fstream>
#include <iostream>

class Potential
{
private:

    /////////////////////////////////////////////  Initialisers  /////////////////////////////////////////////

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to load in the parameters of the potential, with the possibility to 
     * include overrides for alternative parameterisations.
     * 
     * @param        string path                Path to the config file.
     * @param        bool debug                 Outputs loaded parameters if true.
     */
    virtual void configure(const std::string path, const bool debug = false) = 0;

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