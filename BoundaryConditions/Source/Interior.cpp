#include "Interior.hpp"

Interior::Interior(std::vector<float> &scalar_fields,
                   std::vector<float> &vector_fields,
                   std::vector<Analyser*> &analysers,
                   const unsigned &num_scalar_components,
                   const unsigned &num_vector_components,
                   Model &model,
                   const unsigned &nx,
                   const unsigned &ny,
                   const unsigned &nz,
                   const double &dt,
                   const std::vector<int> bound_vector)
    : Periodic(scalar_fields, vector_fields, analysers,
               num_scalar_components, num_vector_components,
               model, nx, ny, nz, dt, bound_vector)
{
}

Interior::~Interior()
{
}