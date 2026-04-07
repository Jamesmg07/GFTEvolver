#include "Gradient.hpp"

//////////////////////////////////////////  Constructors/Destructors  ///////////////////////////////////////////

Gradient::Gradient()
{
}

Gradient::~Gradient()
{
}




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                              Null Gradient                                                  //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////  Constructors/Destructors  /////////////////////////////////////////////

NullGradient::NullGradient(const unsigned &num_scalar_components, const unsigned &num_vector_components)
    : numScalarComponents(num_scalar_components), numVectorComponents(num_vector_components)
{
}

NullGradient::~NullGradient()
{
}

/////////////////////////////////////////////  Public Functions  ////////////////////////////////////////////////

unsigned NullGradient::getDefaultStencilSize() const
{
    return 0;
}

unsigned NullGradient::getNumberOfConstraintEquations() const
{
    return 0;
}

float NullGradient::calcGradientEnergy(const std::vector<std::vector<const float *>> &scalar_pointers, 
                                       const std::vector<std::vector<const float *>> &vector_pointers) const
{
    return 0.0f;
}

float NullGradient::calcKineticEnergy(const float *const local_scalar_fields[2]) const
{
    return 0.0f;
}

std::vector<float> NullGradient::calcConstraintContributions(const float *const local_scalar_fields[2]) const
{
    return std::vector<float>();
}

std::vector<double> NullGradient::calcDerivatives(const std::vector<std::vector<const float *>> &scalar_pointers, 
                                                  const std::vector<std::vector<const float *>> &vector_pointers) const
{
    return std::vector<double>(this->numScalarComponents, 0.0);
}

std::vector<double> NullGradient::calcCurrents(const std::vector<std::vector<const float *>> &scalar_pointers, const std::vector<std::vector<const float *>> &vector_pointers) const
{
    return std::vector<double>(this->numVectorComponents, 0.0);
}
