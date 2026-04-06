#include "StandardModelGradients.hpp"

/////////////////////////////////////////  Initialisers  ////////////////////////////////////////////////////

void Gradients::StandardModel::configure(const std::string path, const bool debug)
{
    // Add different orders of derivatives to config at some point.
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
    this->energy_stencils[0].resize(3, std::vector<double>(3, 0.0));

    this->energy_stencils[0][0][1] = -1.0/this->dx;
    this->energy_stencils[0][0][2] = 1.0/this->dx;

    this->energy_stencils[0][1][1] = -1.0/this->dy;
    this->energy_stencils[0][1][2] = 1.0/this->dy;

    this->energy_stencils[0][2][1] = -1.0/this->dz;
    this->energy_stencils[0][2][2] = 1.0/this->dz;



    // Load in from file
    std::ifstream ifs(path);

    // Read in the parameter values from the configuration file
    std::string description; // Descriptions in the config file (to be dumped)
    if (ifs.is_open())
    {
        // Skip the initial configuration file description
        for(int iter = 0; iter < 12; iter++) std::getline(ifs, description);

        std::getline(ifs, description, ':');
        ifs >> this->divisionByZeroTolerance; // For calculating unit vector.

    } 

    ifs.close();

    if (debug)
    {
        std::cout << "GRADIENTS::STANDARDMODEL::\n"
                  << "Division by zero tolerance: " << this->divisionByZeroTolerance
                  << "\n" << std::endl;
    }
}


std::vector<double> Gradients::StandardModel::transform(const std::vector<double> &scalar_fields, const float *const vector_pointer,
                                                        const unsigned dir_index, const bool conjugate) const
{
    std::vector<double> transformed_field(4, 0.0), hyper_field(4, 0.0);

    int conj_fac = 1;
    if (conjugate)
        conj_fac = -1;

    // Firstly deal with the hypercharge part which is much easier.

    for (int comp_iter = 0; comp_iter < 4; comp_iter++)
    {
        int pm = -1 + 2*(comp_iter%2);
        hyper_field[comp_iter] = std::cos(static_cast<double>(vector_pointer[dir_index]))*static_cast<double>(scalar_fields[comp_iter])
                               + pm*std::sin(conj_fac*static_cast<double>(vector_pointer[dir_index]))*static_cast<double>(scalar_fields[comp_iter-pm]);
    }

    // Second step is to deal with the isospin transformation. First need to compute the matrix (without doing an infinite sum!).
    // Use U_L = e^{iw_i^a\sigma^a} = \cos|w_i^a|\sigma^0 + i\sin|w_i^a|\hat{w}_i^a\sigma^a, where |w_i^a| is the sqrt of quadrature sum over a (not i!).

    // double w_mag = 0.0;
    // for (unsigned comp_iter = 1; comp_iter < 4; comp_iter++)
    //     w_mag += std::pow(static_cast<double>(vector_pointer[dir_index + comp_iter]), 2);

    // w_mag = std::sqrt(w_mag);

    // double w_unit[3] = {0.0, 0.0, 0.0};
    // if (w_mag > this->divisionByZeroTolerance)
    // {
    //     for (unsigned comp_iter = 0; comp_iter < 3; comp_iter++)
    //         w_unit[comp_iter] = conj_fac*static_cast<double>(vector_pointer[dir_index + comp_iter + 1])/w_mag;
    // }
    // else
    // {
    //     w_unit[0] = 1.0; // Arbitrarily choose this direction.
    // }

    // // Perform the SU(2) transformation.
    // transformed_field[0] = std::cos(w_mag)*hyper_field[0] + std::sin(w_mag)*(-w_unit[0]*hyper_field[3] + w_unit[1]*hyper_field[2] - w_unit[2]*hyper_field[1]);
    // transformed_field[1] = std::cos(w_mag)*hyper_field[1] + std::sin(w_mag)*(w_unit[0]*hyper_field[2] + w_unit[1]*hyper_field[3] + w_unit[2]*hyper_field[0]);
    // transformed_field[2] = std::cos(w_mag)*hyper_field[2] + std::sin(w_mag)*(-w_unit[0]*hyper_field[1] - w_unit[1]*hyper_field[0] + w_unit[2]*hyper_field[3]);
    // transformed_field[3] = std::cos(w_mag)*hyper_field[3] + std::sin(w_mag)*(w_unit[0]*hyper_field[0] - w_unit[1]*hyper_field[1] - w_unit[2]*hyper_field[2]);

    // Alternative method that directly maps to the c_representation
    std::vector<double> U_rep = this->getSU2Representation(vector_pointer, dir_index, conjugate);

    //std::cout << std::setprecision(12) << U_rep[0] << " " << U_rep[1] << " " << U_rep[2] << " " << U_rep[3] << std::endl;

    transformed_field[0] =  U_rep[0]*hyper_field[0] - U_rep[3]*hyper_field[1] + U_rep[2]*hyper_field[2] - U_rep[1]*hyper_field[3];
    transformed_field[1] =  U_rep[3]*hyper_field[0] + U_rep[0]*hyper_field[1] + U_rep[1]*hyper_field[2] + U_rep[2]*hyper_field[3];
    transformed_field[2] = -U_rep[2]*hyper_field[0] - U_rep[1]*hyper_field[1] + U_rep[0]*hyper_field[2] + U_rep[3]*hyper_field[3];
    transformed_field[3] =  U_rep[1]*hyper_field[0] - U_rep[2]*hyper_field[1] - U_rep[3]*hyper_field[2] + U_rep[0]*hyper_field[3];

    return transformed_field;
}


std::vector<double> Gradients::StandardModel::getSU2Representation(const float *const vector_pointer, const unsigned dir_index, const bool conjugate) const
{
    std::vector<double> c_representation(4, 0.f);

    int conj_fac = 1;
    if (conjugate)
        conj_fac = -1;

    // Approach when 3 dofs are stored:

    // double w_mag = 0.0;
    // for (unsigned comp_iter = 1; comp_iter < 4; comp_iter++) // comp_iter starts at 1 because 0 is the hypercharge, which I don't care about here.
    //     w_mag += static_cast<double>(vector_pointer[dir_index + comp_iter])*static_cast<double>(vector_pointer[dir_index + comp_iter]);

    // w_mag = std::sqrt(w_mag);
    // c_representation[0] = std::cos(w_mag);

    // if (w_mag > this->divisionByZeroTolerance)
    // {
    //     for (unsigned comp_iter = 1; comp_iter < 4; comp_iter++)
    //         c_representation[comp_iter] = conj_fac*std::sin(w_mag)*static_cast<double>(vector_pointer[dir_index + comp_iter])/w_mag;
    // }
    // else
    // {
    //     c_representation[1] = std::sin(w_mag); // Arbitrarily choose this direction.
    // }

    // Approach when 4 dofs are stored:

    c_representation[0] = static_cast<double>(vector_pointer[dir_index + 1]);
    for (unsigned comp_iter = 2; comp_iter < 5; comp_iter++) // comp_iter starts at 2 because 0 is the hypercharge and 1 the c0 component.
        c_representation[comp_iter-1] = static_cast<double>(conj_fac*vector_pointer[dir_index + comp_iter]);

    return c_representation;
}


std::vector<double> Gradients::StandardModel::getScalarField(const float *const scalar_pointer) const
{
    std::vector<double> scalar_field = {static_cast<double>(scalar_pointer[0]), static_cast<double>(scalar_pointer[1]), 
                                        static_cast<double>(scalar_pointer[2]), static_cast<double>(scalar_pointer[3])};
    return scalar_field;
}


/////////////////////////////////////  Constructors/Destructors  ////////////////////////////////////////////


Gradients::StandardModel::StandardModel(
    const unsigned &num_scalar_components, const unsigned &num_vector_components,
    const unsigned &nx, const unsigned &ny, const unsigned &nz,
    const double &dt, const double &dx, const double &dy, const double &dz) : 
    numScalarComponents(num_scalar_components), numVectorComponents(num_vector_components),
    dt(dt), dx(dx), dy(dy), dz(dz)
{
    this->configure(std::string(SOURCE_DIR) + "/Config/StandardModel.cfg", true);

    if (this->numScalarComponents != 4)
        throw std::runtime_error("GRADIENTS::STANDARDMODEL:: 4 scalar components are required but " + std::to_string(this->numScalarComponents) 
                                + " have been assigned.");
    if (this->numVectorComponents != 12 && this->numVectorComponents != 15)
        throw std::runtime_error("GRADIENTS::STANDARDMODEL:: Either 12 or 15 vector components are required (depending on whether the generator "
                                "or quaternion representation is used) but " + std::to_string(this->numVectorComponents) + " have been assigned.\n"
                                "Note that 3 spatial dimensions are assumed so requested number of components is multiplied by 3 internally.");
}

Gradients::StandardModel::~StandardModel()
{
}

////////////////////////////////////////  Public Functions  /////////////////////////////////////////////////


unsigned Gradients::StandardModel::getDefaultStencilSize() const
{
    // Assumption is that the stencil has an odd number of elements, so an integer divide by 2 gives the size.
    // Other assumption is that the default stencil has the same size in each direction.
    return this->stencils[0][0].size()/2;
}


float Gradients::StandardModel::calcGradientEnergy(const std::vector<std::vector<const float *>> &scalar_pointers,
                                                   const std::vector<std::vector<const float *>> &vector_pointers) const
{
    float gradient_energy = 0.f;

    // Assumes that this is a purely forward derivative (i.e no negative side contributes)

    std::vector<float> derivatives(3*this->numScalarComponents, 0.f);

    for (unsigned stencil_iter = this->stencil_sizes[0]; stencil_iter < 2*this->stencil_sizes[0] + 1; stencil_iter++)
    {
        unsigned stencil_dif = 2*this->stencil_sizes[0] - stencil_iter;
        
        for (unsigned axis_iter = 0; axis_iter < 3; axis_iter++)
        {
            // Approach when using 3 dofs:
            //unsigned dir_index = axis_iter*4; // To index the correct spatial component of the gauge fields

            // Approach when using 4 dofs:
            unsigned dir_index = axis_iter*5;

            unsigned dir_scalar_index = axis_iter*this->numScalarComponents; // To store scalar fields in derivative before squaring.

            // Transform the fields with the gauge fields.

            std::vector<double> transformed_field = this->getScalarField(scalar_pointers[axis_iter][stencil_iter]);

            for (int transform_iter = 0; transform_iter < stencil_dif; transform_iter++)
            {
                transformed_field = this->transform(transformed_field, vector_pointers[axis_iter][stencil_iter + transform_iter], dir_index, false);
            }

            // Now can sum over the stencil
            for (unsigned comp_iter = 0; comp_iter < 4; comp_iter++)
            {
                derivatives[dir_scalar_index + comp_iter] += this->energy_stencils[0][axis_iter][stencil_iter]*transformed_field[comp_iter];
            }


            
        }
    }

    for (float derivative : derivatives)
        gradient_energy += std::pow(derivative, 2);

    return 0.5f*gradient_energy;
    
}


float Gradients::StandardModel::calcKineticEnergy(const float *const local_scalar_fields[2]) const
{
    float kinetic_energy = 0.f;

    for (unsigned comp_iter = 0; comp_iter < this->numScalarComponents; comp_iter++)
    {
        kinetic_energy += powf((local_scalar_fields[1][comp_iter] - local_scalar_fields[0][comp_iter])/static_cast<float>(this->dt), 2);
    }

    return 0.5f*kinetic_energy;
}

std::vector<float> Gradients::StandardModel::calcConstraintContributions(const float *const local_scalar_fields[2]) const
{
    std::vector<float> contribution(4, 0.f);

    // Hypercharge current
    contribution[0] = 0.25*( local_scalar_fields[1][0]*local_scalar_fields[0][1] - local_scalar_fields[1][1]*local_scalar_fields[0][0]
                           + local_scalar_fields[1][2]*local_scalar_fields[0][3] - local_scalar_fields[1][3]*local_scalar_fields[0][2] );


    // Isospin currents
    contribution[1] = 0.25*( local_scalar_fields[1][0]*local_scalar_fields[0][3] - local_scalar_fields[1][1]*local_scalar_fields[0][2]
                           + local_scalar_fields[1][2]*local_scalar_fields[0][1] - local_scalar_fields[1][3]*local_scalar_fields[0][0]  );

    contribution[2] = 0.25*( -local_scalar_fields[1][0]*local_scalar_fields[0][2] - local_scalar_fields[1][1]*local_scalar_fields[0][3]
                            + local_scalar_fields[1][2]*local_scalar_fields[0][0] + local_scalar_fields[1][3]*local_scalar_fields[0][1] );

    contribution[3] = 0.25*( local_scalar_fields[1][0]*local_scalar_fields[0][1] - local_scalar_fields[1][1]*local_scalar_fields[0][0]
                           - local_scalar_fields[1][2]*local_scalar_fields[0][3] + local_scalar_fields[1][3]*local_scalar_fields[0][2] );

    return contribution;
}

std::vector<double> Gradients::StandardModel::calcDerivatives(const std::vector<std::vector<const float *>> &scalar_pointers,
                                                              const std::vector<std::vector<const float *>> &vector_pointers) const
{
    std::vector<double> derivativeContributions(this->numScalarComponents, 0.0);

    // In this case, I will need stencil pointers for the scalar and the vector fields.
    // Probably the global case will need to be given pointers that do nothing too...

    // First thing to do is to define a function that transforms the scalar fields under the action of a given link variable.

    for (unsigned stencil_iter = 0; stencil_iter < 2*this->stencil_sizes[0] + 1; stencil_iter++)
    {
        int stencil_dif = static_cast<int>(stencil_iter)- static_cast<int>(this->stencil_sizes[0]) ;
        bool is_negative_side = stencil_dif < 0;

        for (unsigned axis_iter = 0; axis_iter < 3; axis_iter++)
        {
            // Approach when using 3 dofs:
            //unsigned dir_index = axis_iter*4; // For indexing the correct spatial component of the vector fields

            // Approach when using 4 dofs:
            unsigned dir_index = axis_iter*5;

            // Transform the field
            std::vector<double> transformed_field = this->getScalarField(scalar_pointers[axis_iter][stencil_iter]);

            // Central point with stencil_dif=0 will not be transformed.
            if (is_negative_side)       // Transform for larger stencils would go like U(x-dx)U(x-2dx)f(x-2dx) etc...
            {
                for (int transform_iter = 0; transform_iter < -stencil_dif; transform_iter++)
                {
                    transformed_field = this->transform(transformed_field, vector_pointers[axis_iter][transform_iter], dir_index, false);
                }
            }
            else                        // Transform for larger stencils would go like U^\dagger(x)U^\dagger(x+dx)f(x+2dx) etc...
            {
                for (int transform_iter = 0; transform_iter < stencil_dif; transform_iter++)
                {
                    transformed_field = this->transform(transformed_field, vector_pointers[axis_iter][stencil_iter - 1 - transform_iter], dir_index, true);
                }
            }           
            
            // Now that the fields are transformed properly, the sum over the stencil proceeds as in the global case.

            for (unsigned comp_iter = 0; comp_iter < this->numScalarComponents; comp_iter++)
            {
                derivativeContributions[comp_iter] += this->stencils[0][axis_iter][stencil_iter]*transformed_field[comp_iter];
            }
        }
    }

    return derivativeContributions;
}


std::vector<double> Gradients::StandardModel::calcCurrents(const std::vector<std::vector<const float*>> &scalar_pointers,
                                                           const std::vector<std::vector<const float*>> &vector_pointers) const
{
    // THIS IS ANOTHER FUNCTION THAT CURRENTLY ASSUMES A STENCIL SIZE OF 1
    // NEEDS GENERALISING FOR LARGER STENCILS

    std::vector<double> currentContribution(12, 0.0);

    for (unsigned axis_iter = 0; axis_iter < 3; axis_iter++)
    {
        // Approach when using 3 dofs:
        //unsigned eq_dir_index = axis_iter*4;
        //unsigned dir_index = axis_iter*4; // For indexing the correct spatial component of the vector fields

        // Approach when using 4 dofs:
        unsigned eq_dir_index = axis_iter*4;
        unsigned dir_index = axis_iter*5;

        // First step is to transform phi(x) with U(x) to get tilde{phi(x)}
        std::vector<double> transformed_field = this->getScalarField(scalar_pointers[axis_iter][1]);
        transformed_field = this->transform(transformed_field, vector_pointers[axis_iter][1], dir_index, false);

        std::vector<double> neighbour_field = this->getScalarField(scalar_pointers[axis_iter][2]);

        // Next step is to calculate the imaginary part of tilde{phi(x)}^dagger sigma^mu phi(x+s_i).
        // Also needs multiplying by -1/4 and the gauge coupling^2

        // Hypercharge current
        currentContribution[eq_dir_index] = 0.25*( transformed_field[0]*neighbour_field[1] - transformed_field[1]*neighbour_field[0]
                                                 + transformed_field[2]*neighbour_field[3] - transformed_field[3]*neighbour_field[2] );


        // Isospin currents
        currentContribution[eq_dir_index + 1] = 0.25*( transformed_field[0]*neighbour_field[3] - transformed_field[1]*neighbour_field[2]
                                                     + transformed_field[2]*neighbour_field[1] - transformed_field[3]*neighbour_field[0]  );

        currentContribution[eq_dir_index + 2] = 0.25*( -transformed_field[0]*neighbour_field[2] - transformed_field[1]*neighbour_field[3]
                                                      + transformed_field[2]*neighbour_field[0] + transformed_field[3]*neighbour_field[1] );

        currentContribution[eq_dir_index + 3] = 0.25*( transformed_field[0]*neighbour_field[1] - transformed_field[1]*neighbour_field[0]
                                                     - transformed_field[2]*neighbour_field[3] + transformed_field[3]*neighbour_field[2] );

    }

    // Returns components of 0.25*(g or gp)^2\Phi^\dagger(x)U_i^\dagger\sigma^\mu\Phi(x+s_i)
    return currentContribution;
}
