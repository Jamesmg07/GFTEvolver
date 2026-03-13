#include "BoundaryCondition.hpp"

////////////////////////////////////////////////  Constructors/Destructors  /////////////////////////////////////////////////

BoundaryCondition::BoundaryCondition(std::vector<float> &scalar_fields, std::vector<float> &vector_fields, std::vector<Analyser*> &analysers,
                                     const unsigned &num_scalar_components, const unsigned &num_vector_components,
                                     Model &model, const unsigned &nx, const unsigned &ny, const unsigned &nz, const double &dt,
                                     const std::vector<int> bound_vector)
    : scalarFields(scalar_fields), vectorFields(vector_fields), analysers(analysers), 
      numScalarComponents(num_scalar_components), numVectorComponents(num_vector_components), 
      model(model), nx(nx), ny(ny), nz(nz), dt(dt), 
      boundVector(bound_vector)
{
    this->loopLimits.resize(3, std::vector<unsigned>(2, 0) );
}

BoundaryCondition::~BoundaryCondition()
{
}

bool BoundaryCondition::determineResponsibilities(const unsigned &stencil_size)
{
    // Useful to define an array of the number of grid points in each direction
    std::vector<unsigned> lattice_sizes{this->nx, this->ny, this->nz};

    // Check whether this boundary has no responsibilities
    bool is_empty = false;

    for (unsigned iter = 0; iter < 3; iter++)
    {
        switch (this->boundVector[iter])
        {
        case -1:
            
            this->loopLimits[iter][0] = 0;
            this->loopLimits[iter][1] = stencil_size;
            break;

        case 0:

            this->loopLimits[iter][0] = stencil_size;

            if (lattice_sizes[iter] < 2*stencil_size)
                this->loopLimits[iter][1] = stencil_size;
            else
                this->loopLimits[iter][1] = lattice_sizes[iter] - stencil_size;
            break;

        case 1:
            
            if (lattice_sizes[iter] < 2*stencil_size)
                this->loopLimits[iter][0] = stencil_size;
            else
                this->loopLimits[iter][0] = lattice_sizes[iter] - stencil_size;

            this->loopLimits[iter][1] = lattice_sizes[iter];
            break;
        
        default:

            std::cout << "BOUNDARYCONDITION::ERROR: Invalid boundVector component: " << this->boundVector[iter] << std::endl;
            break;
        }

        if (this->loopLimits[iter][0] >= this->loopLimits[iter][1])
            is_empty = true;
    }

    return is_empty;
}
