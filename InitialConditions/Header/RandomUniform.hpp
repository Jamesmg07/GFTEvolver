#pragma once

#include "InitialCondition.hpp"

#include <fstream>
#include <iostream>
#include <random>

class RandomUniform:
    public InitialCondition
{
private:

    ////////////////////////////////////////////////  Variables  /////////////////////////////////////////////////////

    int seed;
    float minRange, maxRange;
    bool setTimestepsEqual;

    ///////////////////////////////////////////////  Initialisers  ///////////////////////////////////////////////////

    /*
     * Initialises member variables to zero.
     */
    void initVariables();

    /*
     * Loads in parameters for the initial conditions from a config file.
     *
     * @param        string path                Path to the config file.
     * @param        bool debug                 Outputs loaded parameters if true.
     */
    void configure(const std::string path, const bool debug = false);

public:

    //////////////////////////////////////////  Constructors/Destructors  ////////////////////////////////////////////

    RandomUniform();
    virtual ~RandomUniform();

    /////////////////////////////////////////////  Public Functions  /////////////////////////////////////////////////

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to set the values of the fields at the initial two timesteps.
     * 
     * @param        vector<float>& scalarFields                Reference to scalarFields. To be assigned values.
     * @param        vector<float>& vectorFields                Reference to vectorFields. To be assigned values.
     */
    void setInitialFields(std::vector<float>& field) const;
};