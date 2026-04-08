#pragma once

#include "Potential.hpp"

#include <cmath>

class Double_SO_N:
    public Potential
{
private:

    //////////////////////////////////////////////  Variables  //////////////////////////////////////////////////

    const unsigned &totalNumComponents;
    unsigned numComponents[2];
    double m1, m2, m12r, m12i;
    double l1, l2, l3, l4p5, l4m5, l5i, l6r, l6i, l7r, l7i;

    bool breakToSON, breakToUHalfN;

    // Quantities that regularly appear in the calculations. 
    double fieldSqrMagnitude[2], fieldsDot, fieldsAntisym;

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

    // Possibly add a function to determine the ground state of the potential.
    // This may be quite computationally taxing but only needs to run once, then store potential energy of ground state to be subtracted.

public:

    ////////////////////////////////////////  Constructors/Destructors  /////////////////////////////////////////

    Double_SO_N(const unsigned &num_components);
    virtual ~Double_SO_N();

    ///////////////////////////////////////////  Public Functions  //////////////////////////////////////////////

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

    /*
     * Calculates the contribution from the potential to each equation of motion.
     * 
     * @param        float* field                Pointer to the scalar fields at this position.
     * 
     * @return       vector<double>              Array containing the contribution to each equation of motion.  
     */
    std::vector<double> calcPotentialDerivatives(const float* field);
};