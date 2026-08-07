#include "TwoHDMGradients.hpp"

/////////////////////////////////////////  Initialisers  ////////////////////////////////////////////////////

void Gradients::TwoHDM::configure(const std::string path, const bool debug)
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
        std::cout << "GRADIENTS::TWOHDM::\n"
                  << "Division by zero tolerance: " << this->divisionByZeroTolerance
                  << "\n" << std::endl;
    }
}


std::vector<double> Gradients::TwoHDM::transform(const std::vector<double> &scalar_fields, const float *const vector_pointer,
                                                 const unsigned dir_index, const bool conjugate) const
{
    std::vector<double> transformed_field(8, 0.0), hyper_field(8, 0.0), SM_field(8, 0.0);

    int conj_fac = 1;
    if (conjugate)
        conj_fac = -1;

    // Firstly deal with the hypercharge part which is much easier.

    for (int comp_iter = 0; comp_iter < 8; comp_iter++)
    {
        int pm = -1 + 2*(comp_iter%2); // Alternates between -1 and 1
        hyper_field[comp_iter] = std::cos(static_cast<double>(vector_pointer[dir_index]))*static_cast<double>(scalar_fields[comp_iter])
                               + pm*std::sin(conj_fac*static_cast<double>(vector_pointer[dir_index]))*static_cast<double>(scalar_fields[comp_iter - pm]);
    }

    // Second step is to deal with the isospin transformation. First need to compute the matrix (without doing an infinite sum!).
    // Call function that transforms e^{iw^a\sigma^a} into c_0\sigma^0 + ic_a\sigma^a since it is easier to work with c_0 and c_a.
    std::vector<double> U_rep = this->getSU2Representation(vector_pointer, dir_index, 1U, conjugate);

    for (unsigned doublet_iter = 0; doublet_iter < 2; doublet_iter++) // Loop over both doublets
    {
        unsigned doublet_index = doublet_iter*4U;
        SM_field[doublet_index]     =  U_rep[0]*hyper_field[doublet_index] - U_rep[3]*hyper_field[doublet_index + 1] 
                                    +  U_rep[2]*hyper_field[doublet_index + 2] - U_rep[1]*hyper_field[doublet_index + 3];
        SM_field[doublet_index + 1] =  U_rep[3]*hyper_field[doublet_index] + U_rep[0]*hyper_field[doublet_index + 1] 
                                    +  U_rep[1]*hyper_field[doublet_index + 2] + U_rep[2]*hyper_field[doublet_index + 3];
        SM_field[doublet_index + 2] = -U_rep[2]*hyper_field[doublet_index] - U_rep[1]*hyper_field[doublet_index + 1] 
                                    +  U_rep[0]*hyper_field[doublet_index + 2] + U_rep[3]*hyper_field[doublet_index + 3];
        SM_field[doublet_index + 3] =  U_rep[1]*hyper_field[doublet_index] - U_rep[2]*hyper_field[doublet_index + 1] 
                                    -  U_rep[3]*hyper_field[doublet_index + 2] + U_rep[0]*hyper_field[doublet_index + 3];
    }

    unsigned offset;
    if (this->gaugeNum > 4U)
    {
        if (this->usingGeneratorRepresentation)
            offset = 4U;
        else
            offset = 5U;
    }
    else
        return SM_field;

    // Only proceeds if there is an extra gauged U(1) or SU(2) (gaugeNum > 4)

    if (this->gaugeNum == 5U) // Also transform under the extra U(1)
    {
        for (int comp_iter = 0; comp_iter < 8; comp_iter++)
        {
            int pm = -1 + 2*(comp_iter%2); // Alternates between -1 and 1 (period is 2)
            int pm4 = 1 - 2*(comp_iter/4); // Alternates between 1 and -1 (period is 4 as long as max comp_iter is 7)
            transformed_field[comp_iter] = std::cos(static_cast<double>(vector_pointer[dir_index + offset]))*SM_field[comp_iter]
                                         + pm*std::sin(pm4*conj_fac*static_cast<double>(vector_pointer[dir_index + offset]))*SM_field[comp_iter - pm];
        }
    }
    else if (this->gaugeNum == 7U) // Also transform under the extra SU(2)
    {
        U_rep = this->getSU2Representation(vector_pointer, dir_index, offset, conjugate);

        for (unsigned inner_iter = 0; inner_iter < 2; inner_iter++)
        {
            unsigned inner_index = inner_iter*2;
            transformed_field[inner_index]     =  U_rep[0]*SM_field[inner_index] - U_rep[3]*SM_field[inner_index + 1] 
                                               +  U_rep[2]*SM_field[inner_index + 4] - U_rep[1]*SM_field[inner_index + 5];
            transformed_field[inner_index + 1] =  U_rep[3]*SM_field[inner_index] + U_rep[0]*SM_field[inner_index + 1] 
                                               +  U_rep[1]*SM_field[inner_index + 4] + U_rep[2]*SM_field[inner_index + 5];
            transformed_field[inner_index + 4] = -U_rep[2]*SM_field[inner_index] - U_rep[1]*SM_field[inner_index + 1] 
                                               +  U_rep[0]*SM_field[inner_index + 4] + U_rep[3]*SM_field[inner_index + 5];
            transformed_field[inner_index + 5] =  U_rep[1]*SM_field[inner_index] - U_rep[2]*SM_field[inner_index + 1] 
                                               -  U_rep[3]*SM_field[inner_index + 4] + U_rep[0]*SM_field[inner_index + 5];
        }
    }

    return transformed_field;
}


std::vector<double> Gradients::TwoHDM::getSU2Representation(const float *const vector_pointer, const unsigned dir_index, const unsigned offset,
                                                            const bool conjugate) const
{
    std::vector<double> c_representation(4, 0.f);

    int conj_fac = 1;
    if (conjugate)
        conj_fac = -1;

    if (usingGeneratorRepresentation)    // Approach when 3 dofs are stored
    {
        double w_mag = 0.0;
        for (unsigned comp_iter = offset; comp_iter < offset + 3; comp_iter++)
            w_mag += std::pow(static_cast<double>(vector_pointer[dir_index + comp_iter]), 2);

        w_mag = std::sqrt(w_mag);
        c_representation[0] = std::cos(w_mag);

        if (w_mag > this->divisionByZeroTolerance)
        {
            for (unsigned comp_iter = 1; comp_iter < 4; comp_iter++)
                c_representation[comp_iter] = conj_fac*std::sin(w_mag)*static_cast<double>(vector_pointer[dir_index + comp_iter + offset - 1])/w_mag;
        }
        else
        {
            c_representation[1] = std::sin(w_mag); // Arbitrarily choose this direction.
        }
    }
    else                                // Approach when 4 dofs are stored
    {
        c_representation[0] = static_cast<double>(vector_pointer[dir_index + offset]);
        for (unsigned comp_iter = 1; comp_iter < 4; comp_iter++)
            c_representation[comp_iter] = static_cast<double>(conj_fac*vector_pointer[dir_index + comp_iter + offset]);
    }

    return c_representation;
}


std::vector<double> Gradients::TwoHDM::getScalarField(const float *const scalar_pointer) const
{
    std::vector<double> scalar_field = {static_cast<double>(scalar_pointer[0]), static_cast<double>(scalar_pointer[1]), 
                                        static_cast<double>(scalar_pointer[2]), static_cast<double>(scalar_pointer[3]),
                                        static_cast<double>(scalar_pointer[4]), static_cast<double>(scalar_pointer[5]),
                                        static_cast<double>(scalar_pointer[6]), static_cast<double>(scalar_pointer[7])};
    return scalar_field;
}


/////////////////////////////////////  Constructors/Destructors  ////////////////////////////////////////////


Gradients::TwoHDM::TwoHDM(
    const unsigned &num_scalar_components, const unsigned &num_vector_components,
    const unsigned &nx, const unsigned &ny, const unsigned &nz,
    const double &dt, const double &dx, const double &dy, const double &dz,
    const bool using_generator_representation) : 
    numScalarComponents(num_scalar_components), numVectorComponents(num_vector_components),
    dt(dt), dx(dx), dy(dy), dz(dz),
    usingGeneratorRepresentation(using_generator_representation)
{
    this->configure(std::string(SOURCE_DIR) + "/Config/TwoHDM.cfg", true);

    if (this->numScalarComponents != 8)
        throw std::runtime_error("GRADIENTS::TWOHDM:: 8 scalar components are required but " + std::to_string(this->numScalarComponents) 
                                + " have been assigned.");
    
    if (this->usingGeneratorRepresentation)
    {
        if (this->numVectorComponents == 12)
            this->gaugeNum = 4U; // Just SM
        else if (this->numVectorComponents == 15)
            this->gaugeNum = 5U; // SM x U(1)
        else if (this->numVectorComponents == 21)
            this->gaugeNum = 7U; // SM x SU(2)
        else
            throw std::runtime_error("GRADIENTS::TWOHDM:: When using the generator representation, either 12 (just SM), 15 (SM x U(1)) "
                                     "or 21 (SM x SU(2)) vector components are required but " + std::to_string(this->numVectorComponents) +
                                     " have been assigned.\n"
                                     "Note that 3 spatial dimensions are assumed so requested number of components is multiplied by 3 internally.");
    }
    else
    {
        if (this->numVectorComponents == 15)
            this->gaugeNum = 4U; // Just SM
        else if (this->numVectorComponents == 18)
            this->gaugeNum = 5U; // SM x U(1)
        else if (this->numVectorComponents == 27) 
            this->gaugeNum = 7U; // SM x SU(2)
        else
            throw std::runtime_error("GRADIENTS::TWOHDM:: When using the quaternion representation, either 15 (just SM), 18 (SM x U(1)) "
                                     "or 27 (SM x SU(2)) vector components are required but " + std::to_string(this->numVectorComponents) +
                                     " have been assigned.\n"
                                     "Note that 3 spatial dimensions are assumed so requested number of components is multiplied by 3 internally.");
    }

    std::cout << "GRADIENTS::TWOHDM:: " << this->gaugeNum << std::endl;
}

Gradients::TwoHDM::~TwoHDM()
{
}

////////////////////////////////////////  Public Functions  /////////////////////////////////////////////////


unsigned Gradients::TwoHDM::getDefaultStencilSize() const
{
    // Assumption is that the stencil has an odd number of elements, so an integer divide by 2 gives the size.
    // Other assumption is that the default stencil has the same size in each direction.
    return this->stencils[0][0].size()/2;
}

unsigned Gradients::TwoHDM::getNumberOfConstraintEquations() const
{
    return this->gaugeNum;
}


float Gradients::TwoHDM::calcGradientEnergy(const std::vector<std::vector<const float *>> &scalar_pointers,
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
            unsigned dir_index = axis_iter*this->numVectorComponents/3U;
            unsigned dir_scalar_index = axis_iter*this->numScalarComponents; // To store scalar fields in derivative before squaring.

            // Transform the fields with the gauge fields.

            std::vector<double> transformed_field = this->getScalarField(scalar_pointers[axis_iter][stencil_iter]);

            for (int transform_iter = 0; transform_iter < stencil_dif; transform_iter++)
            {
                transformed_field = this->transform(transformed_field, vector_pointers[axis_iter][stencil_iter + transform_iter], dir_index, false);
            }

            // Now can sum over the stencil
            for (unsigned comp_iter = 0; comp_iter < 8; comp_iter++)
            {
                derivatives[dir_scalar_index + comp_iter] += this->energy_stencils[0][axis_iter][stencil_iter]*transformed_field[comp_iter];
            }


            
        }
    }

    for (float derivative : derivatives)
        gradient_energy += std::pow(derivative, 2);

    return 0.5f*gradient_energy;
    
}


float Gradients::TwoHDM::calcKineticEnergy(const float *const local_scalar_fields[2]) const
{
    float kinetic_energy = 0.f;

    for (unsigned comp_iter = 0; comp_iter < this->numScalarComponents; comp_iter++)
    {
        kinetic_energy += powf((local_scalar_fields[1][comp_iter] - local_scalar_fields[0][comp_iter])/static_cast<float>(this->dt), 2);
    }

    return 0.5f*kinetic_energy;
}

std::vector<float> Gradients::TwoHDM::calcConstraintContributions(const float *const local_scalar_fields[2]) const
{
    std::vector<float> contribution(this->gaugeNum, 0.f);

    for (unsigned iter = 0; iter < 2; iter++)
    {
        unsigned d_ind = iter*4U; // Swaps between \Phi_1 and \Phi_2
        unsigned i_ind = iter*2U; // Swaps between upper and lower parts of \Phi_a

        // Hypercharge current
        contribution[0] += 0.25*( local_scalar_fields[1][d_ind]*local_scalar_fields[0][d_ind + 1] - local_scalar_fields[1][d_ind + 1]*local_scalar_fields[0][d_ind]
                                + local_scalar_fields[1][d_ind + 2]*local_scalar_fields[0][d_ind + 3] - local_scalar_fields[1][d_ind + 3]*local_scalar_fields[0][d_ind + 2] );


        // Isospin currents
        contribution[1] += 0.25*( local_scalar_fields[1][d_ind]*local_scalar_fields[0][d_ind + 3] - local_scalar_fields[1][d_ind + 1]*local_scalar_fields[0][d_ind + 2]
                                + local_scalar_fields[1][d_ind + 2]*local_scalar_fields[0][d_ind + 1] - local_scalar_fields[1][d_ind + 3]*local_scalar_fields[0][d_ind]  );

        contribution[2] += 0.25*( -local_scalar_fields[1][d_ind]*local_scalar_fields[0][d_ind + 2] - local_scalar_fields[1][d_ind + 1]*local_scalar_fields[0][d_ind + 3]
                                 + local_scalar_fields[1][d_ind + 2]*local_scalar_fields[0][d_ind] + local_scalar_fields[1][d_ind + 3]*local_scalar_fields[0][d_ind + 1] );

        contribution[3] += 0.25*( local_scalar_fields[1][d_ind]*local_scalar_fields[0][d_ind + 1] - local_scalar_fields[1][d_ind + 1]*local_scalar_fields[0][d_ind]
                                - local_scalar_fields[1][d_ind + 2]*local_scalar_fields[0][d_ind + 3] + local_scalar_fields[1][d_ind + 3]*local_scalar_fields[0][d_ind + 2] );

        // Higgs family currents
        if (this->gaugeNum == 5U)
        {
            contribution[4] += 0.25*( local_scalar_fields[1][i_ind]*local_scalar_fields[0][i_ind + 1] - local_scalar_fields[1][i_ind + 1]*local_scalar_fields[0][i_ind]
                                    - local_scalar_fields[1][i_ind + 4]*local_scalar_fields[0][i_ind + 5] + local_scalar_fields[1][i_ind + 5]*local_scalar_fields[0][i_ind + 4] );
        }
        else if (this->gaugeNum == 7U)
        {
            contribution[4] += 0.25*( local_scalar_fields[1][i_ind]*local_scalar_fields[0][i_ind + 5] - local_scalar_fields[1][i_ind + 1]*local_scalar_fields[0][i_ind + 4]
                                    + local_scalar_fields[1][i_ind + 4]*local_scalar_fields[0][i_ind + 1] - local_scalar_fields[1][i_ind + 5]*local_scalar_fields[0][i_ind] );

            contribution[5] += 0.25*( -local_scalar_fields[1][i_ind]*local_scalar_fields[0][i_ind + 4] - local_scalar_fields[1][i_ind + 1]*local_scalar_fields[0][i_ind + 5]
                                     + local_scalar_fields[1][i_ind + 4]*local_scalar_fields[0][i_ind] + local_scalar_fields[1][i_ind + 5]*local_scalar_fields[0][i_ind + 1] );

            contribution[6] += 0.25*( local_scalar_fields[1][i_ind]*local_scalar_fields[0][i_ind + 1] - local_scalar_fields[1][i_ind + 1]*local_scalar_fields[0][i_ind]
                                    - local_scalar_fields[1][i_ind + 4]*local_scalar_fields[0][i_ind + 5] + local_scalar_fields[1][i_ind + 5]*local_scalar_fields[0][i_ind + 4] );
        }
    }

    return contribution;
}

std::vector<double> Gradients::TwoHDM::calcDerivatives(const std::vector<std::vector<const float *>> &scalar_pointers,
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
            unsigned dir_index = axis_iter*this->numVectorComponents/3U;

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


std::vector<double> Gradients::TwoHDM::calcCurrents(const std::vector<std::vector<const float*>> &scalar_pointers,
                                                    const std::vector<std::vector<const float*>> &vector_pointers) const
{
    // THIS IS ANOTHER FUNCTION THAT CURRENTLY ASSUMES A STENCIL SIZE OF 1
    // NEEDS GENERALISING FOR LARGER STENCILS

    std::vector<double> currentContribution(this->gaugeNum*3U, 0.0);

    for (unsigned axis_iter = 0; axis_iter < 3; axis_iter++)
    {
        unsigned eq_dir_index = axis_iter*this->gaugeNum;
        unsigned dir_index = axis_iter*this->numVectorComponents/3U;

        // First step is to transform phi(x) with U(x) to get tilde{phi(x)}
        std::vector<double> transformed_field = this->getScalarField(scalar_pointers[axis_iter][1]);
        transformed_field = this->transform(transformed_field, vector_pointers[axis_iter][1], dir_index, false);

        std::vector<double> neighbour_field = this->getScalarField(scalar_pointers[axis_iter][2]);

        // Next step is to calculate the imaginary part of tilde{phi(x)}^dagger sigma^mu phi(x+s_i).
        // Also needs multiplying by 1/4 and the gauge coupling^2 (done outside class)

        for (unsigned iter = 0; iter < 2; iter++)
        {
            unsigned d_ind = iter*4U; // Swaps between \Phi_1 and \Phi_2
            unsigned i_ind = iter*2U; // Swaps between upper and lower parts of \Phi_a

            // Hypercharge current
            currentContribution[eq_dir_index] += 0.25*( transformed_field[d_ind]*neighbour_field[d_ind + 1] - transformed_field[d_ind + 1]*neighbour_field[d_ind]
                                                      + transformed_field[d_ind + 2]*neighbour_field[d_ind + 3] - transformed_field[d_ind + 3]*neighbour_field[d_ind + 2] );


            // Isospin currents
            currentContribution[eq_dir_index + 1] += 0.25*( transformed_field[d_ind]*neighbour_field[d_ind + 3] - transformed_field[d_ind + 1]*neighbour_field[d_ind + 2]
                                                          + transformed_field[d_ind + 2]*neighbour_field[d_ind + 1] - transformed_field[d_ind + 3]*neighbour_field[d_ind]  );

            currentContribution[eq_dir_index + 2] += 0.25*( -transformed_field[d_ind]*neighbour_field[d_ind + 2] - transformed_field[d_ind + 1]*neighbour_field[d_ind + 3]
                                                          + transformed_field[d_ind + 2]*neighbour_field[d_ind] + transformed_field[d_ind + 3]*neighbour_field[d_ind + 1] );

            currentContribution[eq_dir_index + 3] += 0.25*( transformed_field[d_ind]*neighbour_field[d_ind + 1] - transformed_field[d_ind + 1]*neighbour_field[d_ind]
                                                          - transformed_field[d_ind + 2]*neighbour_field[d_ind + 3] + transformed_field[d_ind + 3]*neighbour_field[d_ind + 2] );

            // Higgs family currents
            if (this->gaugeNum == 5U)
            {
                currentContribution[eq_dir_index + 4] += 0.25*( transformed_field[i_ind]*neighbour_field[i_ind + 1] - transformed_field[i_ind + 1]*neighbour_field[i_ind]
                                                              - transformed_field[i_ind + 4]*neighbour_field[i_ind + 5] + transformed_field[i_ind + 5]*neighbour_field[i_ind + 4] );
            }
            else if (this->gaugeNum == 7U)
            {
                currentContribution[eq_dir_index + 4] += 0.25*( transformed_field[i_ind]*neighbour_field[i_ind + 5] - transformed_field[i_ind + 1]*neighbour_field[i_ind + 4]
                                                              + transformed_field[i_ind + 4]*neighbour_field[i_ind + 1] - transformed_field[i_ind + 5]*neighbour_field[i_ind] );

                currentContribution[eq_dir_index + 5] += 0.25*( -transformed_field[i_ind]*neighbour_field[i_ind + 4] - transformed_field[i_ind + 1]*neighbour_field[i_ind + 5]
                                                              + transformed_field[i_ind + 4]*neighbour_field[i_ind] + transformed_field[i_ind + 5]*neighbour_field[i_ind + 1] );

                currentContribution[eq_dir_index + 6] += 0.25*( transformed_field[i_ind]*neighbour_field[i_ind + 1] - transformed_field[i_ind + 1]*neighbour_field[i_ind]
                                                              - transformed_field[i_ind + 4]*neighbour_field[i_ind + 5] + transformed_field[i_ind + 5]*neighbour_field[i_ind + 4] );
            }
        }
    }

    return currentContribution;
}
