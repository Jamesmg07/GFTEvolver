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

bool BoundaryCondition::determineResponsibilities(
    const unsigned &stencil_size,
    const unsigned &owned_x_begin,
    const unsigned &owned_x_end)
{
    // nx is the allocated local x-width including halos.
    // Responsibilities are restricted to the owned physical interval.
    std::vector<unsigned> lattice_begins{
        owned_x_begin, 0U, 0U
    };

    std::vector<unsigned> lattice_ends{
        owned_x_end, this->ny, this->nz
    };

    std::vector<unsigned> lattice_sizes{
        owned_x_end - owned_x_begin,
        this->ny,
        this->nz
    };

    bool is_empty = false;

    for (unsigned iter = 0; iter < 3; iter++)
    {
        switch (this->boundVector[iter])
        {
        case -1:

            this->loopLimits[iter][0] = lattice_begins[iter];
            this->loopLimits[iter][1] =
                lattice_begins[iter] + stencil_size;
            break;

        case 0:

            this->loopLimits[iter][0] =
                lattice_begins[iter] + stencil_size;

            if (lattice_sizes[iter] < 2*stencil_size)
                this->loopLimits[iter][1] =
                    lattice_begins[iter] + stencil_size;
            else
                this->loopLimits[iter][1] =
                    lattice_ends[iter] - stencil_size;
            break;

        case 1:

            if (lattice_sizes[iter] < 2*stencil_size)
                this->loopLimits[iter][0] =
                    lattice_begins[iter] + stencil_size;
            else
                this->loopLimits[iter][0] =
                    lattice_ends[iter] - stencil_size;

            this->loopLimits[iter][1] = lattice_ends[iter];
            break;

        default:

            std::cout
                << "BOUNDARYCONDITION::ERROR: Invalid boundVector component: "
                << this->boundVector[iter]
                << std::endl;
            break;
        }

        if (this->loopLimits[iter][0] >=
            this->loopLimits[iter][1])
        {
            is_empty = true;
        }
    }

    return is_empty;
}
