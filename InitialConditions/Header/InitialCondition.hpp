#pragma once

#include <stdexcept>
#include <vector>


/*
 * Read-only description of the global lattice and this rank's owned/storage
 * geometry. Initial-condition classes use this mapping to evaluate the same
 * physical prescription independently of the MPI decomposition.
 */
struct InitialConditionGeometry
{
    unsigned globalNx, globalNy, globalNz;
    unsigned localNx, storageNx;
    unsigned globalXStart;
    unsigned ownedXBegin, ownedXEnd;
    double dx, dy, dz, dt;

    unsigned localToGlobalX(const unsigned local_x) const
    {
        if (local_x < this->ownedXBegin || local_x >= this->ownedXEnd)
            throw std::out_of_range(
                "INITIALCONDITION:: Global x-index requested for a non-owned site."
            );

        return this->globalXStart
             + (local_x - this->ownedXBegin);
    }
};


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
    virtual void setInitialFields(std::vector<float>& field, const InitialConditionGeometry& geometry, const unsigned num_components) const = 0;
};