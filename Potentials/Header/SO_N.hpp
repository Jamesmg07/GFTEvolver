#pragma once

#include "Potential.hpp"

#include <cmath>

class SO_N:
    public Potential
{
private:

    //////////////////////////////////////////////  Variables  //////////////////////////////////////////////////

    const unsigned &numComponents;
    double lambda;
    double etaSqr;
    double fieldSqrMagnitude;

    /////////////////////////////////////////////  Initialisers  ////////////////////////////////////////////////

    /*
     * Initialises potential parameters to zero.
     */
    void initVariables();

    /*
     * Loads in the parameters of the potential from a config file.
     * 
     * @param        string path                Path to the config file.
     * @param        bool debug                 Outputs loaded parameters if true.
     */
    void configure(const std::string path, const bool debug = false);

public:

    ////////////////////////////////////////  Constructors/Destructors  /////////////////////////////////////////

    SO_N(const unsigned &num_components);
    virtual ~SO_N();

    ///////////////////////////////////////////  Public Functions  //////////////////////////////////////////////

    /*
     * Calculates the contribution from the potential to each equation of motion.
     * 
     * @param        float* field                Pointer to the scalar fields at this position.
     * 
     * @return       vector<double>              Array containing the contribution to each equation of motion.  
     */
    std::vector<double> calcPotentialDerivatives(const float* field);

    /*
     * Calculates the potential energy density at this location in the lattice.
     * Assumes that the square magnitude of the field has already been calculated by the
     * calcPotentialDerivatives function.
     * 
     * @param        float* field                Pointer to the scalar fields at this position.
     * 
     * @return       float                       The potential energy (density) at this location.
     */
    float calcPotentialEnergy(const float* field) const;

};