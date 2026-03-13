#include "GlobalGradients.hpp"

//////////////////////////////////////////////////  Initialisers  /////////////////////////////////////////////////////

void Gradients::Global::configure(const std::string path, const bool debug)
{
    // Nothing to configure. Add different orders of derivatives to config at some point.
    // For now there is just 2nd order.

    // Set-up the stencil for the equations of motion.
    this->stencils.resize(1, std::vector< std::vector<double> >());
    this->stencils[0].resize(3, std::vector<double>(3, 0.0));

    this->stencils[0][0][0] = 1.0/(this->dx*this->dx);
    this->stencils[0][0][1] = -2.0/(this->dx*this->dx);
    this->stencils[0][0][2] = 1.0/(this->dx*this->dx);

    this->stencils[0][1][0] = 1.0/(this->dy*this->dy);
    this->stencils[0][1][1] = -2.0/(this->dy*this->dy);
    this->stencils[0][1][2] = 1.0/(this->dy*this->dy);

    this->stencils[0][2][0] = 1.0/(this->dz*this->dz);
    this->stencils[0][2][1] = -2.0/(this->dz*this->dz);
    this->stencils[0][2][2] = 1.0/(this->dz*this->dz);

    // Allocate size of stencil_sizes based on number of stencil options
    stencil_sizes.resize(1, 0);
    stencil_sizes[0] = 1;


    // Set up the stencil for calculating the gradient energy.
    this->energy_stencils.resize(1, std::vector< std::vector<double> >());
    this->energy_stencils[0].resize(3, std::vector<double>(3, 0.f));

    this->energy_stencils[0][0][1] = -1.0/this->dx;
    this->energy_stencils[0][0][2] = 1.0/this->dx;

    this->energy_stencils[0][1][1] = -1.0/this->dy;
    this->energy_stencils[0][1][2] = 1.0/this->dy;

    this->energy_stencils[0][2][1] = -1.0/this->dz;
    this->energy_stencils[0][2][2] = 1.0/this->dz;

    if (debug)
    {
        std::cout << "GRADIENTS::GLOBAL::\n"
                  << "No need to load anything from the config file for this option."
                  << "\n" << std::endl;
    }
}

////////////////////////////////////////////  Constructors/Destructors  ///////////////////////////////////////////////

Gradients::Global::Global(
    const unsigned &num_scalar_components, const unsigned &num_vector_components, 
    const unsigned &nx, const unsigned &ny, const unsigned &nz, 
    const double &dt, const double &dx, const double &dy, const double &dz
) : numScalarComponents(num_scalar_components), numVectorComponents(num_vector_components),
    dt(dt), dx(dx), dy(dy), dz(dz)
{
    this->configure(std::string(SOURCE_DIR) + "/Config/Global.cfg", true);
}

Gradients::Global::~Global()
{
}

//////////////////////////////////////////////// Public Functions  ///////////////////////////////////////////////////

unsigned Gradients::Global::getDefaultStencilSize() const
{
    // Assumption is that the stencil has an odd number of elements, so an integer divide by 2 gives the size.
    // Other assumption is that the default stencil has the same size in each direction.
    return this->stencils[0][0].size()/2;
}

float Gradients::Global::calcGradientEnergy(const std::vector<std::vector<const float*>> &scalar_pointers,
                                            const std::vector<std::vector<const float*>> &vector_pointers) const
{
    float gradient_energy = 0.f;

    for (unsigned comp_iter = 0; comp_iter < this->numScalarComponents; comp_iter++)
    {
        for (unsigned axis_iter = 0; axis_iter < 3; axis_iter++)
        {
            float derivative = 0.f;

            for (unsigned stencil_iter = 0; stencil_iter < 2*this->stencil_sizes[0] + 1; stencil_iter++)
            {
                derivative += this->energy_stencils[0][axis_iter][stencil_iter]*scalar_pointers[axis_iter][stencil_iter][comp_iter];
            }

            gradient_energy += powf(derivative, 2);
        }
    }

    return 0.5f*gradient_energy;
}

float Gradients::Global::calcKineticEnergy(const float* const local_scalar_fields[2]) const
{
    float kinetic_energy = 0.f;

    for (unsigned comp_iter = 0; comp_iter < this->numScalarComponents; comp_iter++)
    {
        kinetic_energy += powf((local_scalar_fields[1][comp_iter] - local_scalar_fields[0][comp_iter])/this->dt, 2);
    }

    return 0.5f*kinetic_energy;
}

std::vector<double> Gradients::Global::calcDerivatives(const std::vector<std::vector<const float*>> &scalar_pointers,
                                                       const std::vector<std::vector<const float*>> &vector_pointers) const
{
    std::vector<double> derivativeContributions(this->numScalarComponents, 0.f);

    for (unsigned comp_iter = 0; comp_iter < this->numScalarComponents; comp_iter++)
    {
        for (unsigned axis_iter = 0; axis_iter < 3; axis_iter++)
        {
            for (unsigned stencil_iter = 0; stencil_iter < 2*this->stencil_sizes[0] + 1; stencil_iter++)
            {
                derivativeContributions[comp_iter] += this->stencils[0][axis_iter][stencil_iter]*static_cast<double>(scalar_pointers[axis_iter][stencil_iter][comp_iter]);
            }
        }
    }

    return derivativeContributions;
}

std::vector<double> Gradients::Global::calcCurrents(const std::vector<std::vector<const float*>> &scalar_pointers, 
                                                   const std::vector<std::vector<const float*>> &vector_pointers) const
{
    return std::vector<double>();
}
