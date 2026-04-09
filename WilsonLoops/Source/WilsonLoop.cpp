#include "WilsonLoop.hpp"

/////////////////////////////////////////////  Constructors/Destructors  /////////////////////////////////////////

WilsonLoop::WilsonLoop()
{
}

WilsonLoop::~WilsonLoop()
{
}




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                NullWilsonLoop                                                 //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////  Constructors/Destructors  ///////////////////////////////////////////

NullWilsonLoop::NullWilsonLoop(const unsigned &num_vector_components)
    : numVectorComponents(num_vector_components)
{
}

NullWilsonLoop::~NullWilsonLoop()
{
}

////////////////////////////////////////////////  Public Functions  ///////////////////////////////////////////////

bool NullWilsonLoop::isUsingGeneratorRepresentation() const
{
    return false;
}

float NullWilsonLoop::getSqrCouplings(const unsigned comp_iter) const
{
    return 0.0f;
}

unsigned NullWilsonLoop::getNumberOfEvolutionEquations() const
{
    return 0;
}

unsigned NullWilsonLoop::getNumberOfConstraintEquations() const
{
    return 0;
}

void NullWilsonLoop::energyPreparation(const bool store_energy)
{
}

float NullWilsonLoop::calcMagneticEnergy(const std::vector<std::vector<const float *>> &vector_pointers) const
{
    return 0.0f;
}

float NullWilsonLoop::calcElectricEnergy(const float *const local_vector_fields[2]) const
{
    return 0.0f;
}

std::vector<float> NullWilsonLoop::calcConstraintContributions(const long long int &t_future_index, const std::vector<std::vector<const float *>> &vector_pointers) const
{
    return std::vector<float>();
}

std::vector<double> NullWilsonLoop::calcMagneticContributions(const std::vector<std::vector<const float *>> &vector_pointers) const
{
    return std::vector<double>(this->numVectorComponents, 0.0);
}

std::vector<double> NullWilsonLoop::calcElectricContributions(const float *const local_vector_fields[2]) const
{
    return std::vector<double>(this->numVectorComponents, 0.0);
}

void NullWilsonLoop::evolve(float *const local_vector_fields[2], std::vector<double> equation_RHS)
{
}
