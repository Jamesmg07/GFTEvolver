#pragma once

#include "Model.hpp"
#include "Analyser.hpp"

class BoundaryCondition
{
private:

    ////////////////////////////////////////////////////  Initialisers  ///////////////////////////////////////////////////////////

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to load in any parameters relating to the boundary condition
     * 
     * @param        string path                Path to the config file.
     * @param        bool debug                 Outputs loaded parameters if true.
     */
    virtual void configure(const std::string path, const bool debug = false) = 0;

protected:

    std::vector<float> &scalarFields, &vectorFields;
    std::vector<Analyser*> &analysers;
    std::vector< std::vector<unsigned> > loopLimits;
    const std::vector<int> boundVector;
    Model &model;
    const unsigned &numScalarComponents, &numVectorComponents;
    const unsigned &nx, &ny, &nz;
    const double &dt;
    unsigned boundaryDepth;



public:

    ///////////////////////////////////////////////  Constructors/Destructors  ////////////////////////////////////////////////////

    BoundaryCondition(std::vector<float> &scalar_fields, std::vector<float> &vector_fields, std::vector<Analyser*> &analysers,
                      const unsigned &num_scalar_components, const unsigned &num_vector_components,
                      Model &model, const unsigned &nx, const unsigned &ny, const unsigned &nz, const double &dt,
                      const std::vector<int> bound_vector);
    virtual ~BoundaryCondition();

    ///////////////////////////////////////////////////  Public Functions  ////////////////////////////////////////////////////////

    /*
     * Determine the grid positions that this class will be responsible for evolving.
     * 
     * @param        unsigned stencil_size                Size of the largest default stencil being used.
     * 
     * @return       bool                                 Returns true if this class is not responsible for any grid points.
     */
    bool determineResponsibilities(const unsigned &stencil_size);

    /*
     * Pure virtual function that must be defined in each child class.
     * Intended use is to evolve the fields that are the responsibility of the
     * boundary condition classes. This may be more than just the edges of the lattice,
     * it depends on the size of the largest default stencil.
     * 
     * @param        unsigned t_now                       Specify the current timestep
     * @param        usigned stencil_size                 Size of the largest default stencil being used.
     */
    virtual void evolve(const unsigned &t_now, const unsigned &stencil_size) = 0;

};