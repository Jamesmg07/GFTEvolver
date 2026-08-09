#pragma once

#include "Periodic.hpp"

/*
 * Boundary region containing an internal MPI process interface.
 *
 * Interior currently reuses Periodic's evolution and pointer machinery.
 * The distinction is semantic: x-neighbours are supplied through halo
 * storage rather than by wrapping within this rank's owned slab.
 */
class Interior:
    public Periodic
{
public:

    Interior(std::vector<float> &scalar_fields,
             std::vector<float> &vector_fields,
             std::vector<Analyser*> &analysers,
             const unsigned &num_scalar_components,
             const unsigned &num_vector_components,
             Model &model,
             const unsigned &nx,
             const unsigned &ny,
             const unsigned &nz,
             const double &dt,
             const std::vector<int> bound_vector);

    virtual ~Interior();
};