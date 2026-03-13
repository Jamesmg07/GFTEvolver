#pragma once

#include <vector>

class InitialCondition
{
private:

public:

    /////////////////////////////////////////  Constructors/Destructors  ///////////////////////////////////////////

    InitialCondition();
    virtual ~InitialCondition();

    /////////////////////////////////////////////  Public Functions  ///////////////////////////////////////////////

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to set the values of the fields at the initial two timesteps.
     * 
     * @param        vector<float>& field                Reference to a field array. To be assigned values.
     */
    virtual void setInitialFields(std::vector<float>& field) const = 0;
};