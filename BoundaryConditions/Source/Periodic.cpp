#include "Periodic.hpp"

/////////////////////////////////////////////////////  Initialisers  /////////////////////////////////////////////////////////////

void Periodic::configure(const std::string path, const bool debug)
{
    if (debug)
    {
        std::cout << "BOUNDARYCONDITIONS::PERIODIC::\n"
                  << "No need to load anything from the config file for this option."
                  << "\n" << std::endl;
    }
}

///////////////////////////////////////////////////  Private Functions  //////////////////////////////////////////////////////////


std::vector<std::vector<long long unsigned>> Periodic::updateRunningIndices(std::vector<std::vector<long long unsigned>> running_indices, 
                                                                            const unsigned &loc, const unsigned axis) const
{

    // Perhaps some of this can be moved so that it is only done once in the constructor or configuration function.

    unsigned stencil_size = (running_indices[0].size()-1)/2;
    long long unsigned sub_array_size;
    unsigned modulus;
    if (axis == 0)
    {
        sub_array_size = 1ULL*this->ny*this->nz;
        modulus = this->nx;
    }
    else if (axis == 1)
    {
        sub_array_size = 1ULL*this->nz;
        modulus = this->ny;
    }
    else if (axis == 2)
    {
        sub_array_size = 1ULL;
        modulus = this->nz;
    }
    else
        std::cout << "BOUNDARYCONDITION::PERIODIC::ERROR: Invalid axis argument in call to modifyRunningIndices." << std::endl;

    for (unsigned axis_iter = 0; axis_iter < 3; axis_iter++)
    {
        for (unsigned stencil_iter = 0; stencil_iter < 2*stencil_size + 1; stencil_iter++)
        {
            if (axis_iter == axis)
                running_indices[axis_iter][stencil_iter] += 1ULL*((modulus + loc + stencil_iter - stencil_size)%modulus)*sub_array_size;
            else
                running_indices[axis_iter][stencil_iter] += 1ULL*loc*sub_array_size;
        }
    }

    return running_indices;
}


std::vector<std::vector<const float *>> Periodic::generateScalarPointers(const std::vector<std::vector<long long unsigned>> &running_indices) const
{
    unsigned array_size = running_indices[0].size();
    std::vector<std::vector<const float*>> stencil_pointers(3, std::vector<const float*>(array_size, nullptr));

    for (unsigned axis_iter = 0; axis_iter < 3; axis_iter++)
    {
        for (unsigned stencil_iter = 0; stencil_iter < array_size; stencil_iter++)
        {
            stencil_pointers[axis_iter][stencil_iter] = &this->scalarFields[1ULL*this->numScalarComponents*running_indices[axis_iter][stencil_iter]];
        }
    }

    return stencil_pointers;
}

// std::vector<std::vector<const float *>> Periodic::generateScalarPointers(const std::vector<std::vector<long long unsigned>> &running_indices,
//                                                                           const unsigned &loc, const unsigned axis) const
// {
//     unsigned stencil_size = (running_indices[0].size()-1)/2;
//     long long unsigned sub_array_size;
//     unsigned modulus;
//     if (axis == 0)
//     {
//         sub_array_size = 1ULL*this->ny*this->nz;
//         modulus = this->nx;
//     }
//     else if (axis == 1)
//     {
//         sub_array_size = 1ULL*this->nz;
//         modulus = this->ny;
//     }
//     else if (axis == 2)
//     {
//         sub_array_size = 1ULL;
//         modulus = this->nz;
//     }
//     else
//         std::cout << "BOUNDARYCONDITIONS::PERIODIC::ERROR: Invalid axis argument in call to generateScalarPointers." << std::endl;

//     std::vector<std::vector<const float*>> stencil_pointers(3, std::vector<const float*>(2*stencil_size + 1, nullptr));

//     for (unsigned axis_iter = 0; axis_iter < 3; axis_iter++)
//     {
//         for (unsigned stencil_iter = 0; stencil_iter < 2*stencil_size + 1; stencil_iter++)
//         {
//             if (axis_iter == axis)
//             {
//                 stencil_pointers[axis_iter][stencil_iter] = &this->scalarFields[
//                     1ULL*this->numScalarComponents*( running_indices[axis_iter][stencil_iter] + 1ULL*((modulus + loc + stencil_iter - stencil_size)%modulus)*sub_array_size )
//                 ];
//             }
//             else
//             {
//                 stencil_pointers[axis_iter][stencil_iter] = &this->scalarFields[
//                     1ULL*this->numScalarComponents*running_indices[axis_iter][stencil_iter] + 1ULL*loc*sub_array_size
//                 ];
//             }
//         }
//     }

//     return stencil_pointers;
// }


std::vector<std::vector<const float *>> Periodic::generateVectorPointers(const std::vector<std::vector<long long unsigned>> &running_indices,
                                                                         const unsigned &x_loc, const unsigned &y_loc, const unsigned &z_loc) const
{
    // STENCIL SIZES ARE NOT REALLY WORKING WITH THIS PROPERLY.
    // IN PRINCIPLE I CAN HAVE A DIFFERENT WILSON LOOP STENCIL TO THE DERIVATIVE STENCIL
    // THIS FUNCTION ASSUMES THAT THE WILSON LOOP STENCIL IS THE MOST BASIC TYPE

    unsigned array_size = running_indices[0].size();
    std::vector<std::vector<const float*>> stencil_pointers(3, std::vector<const float*>(array_size + 2, nullptr));

    for (unsigned axis_iter = 0; axis_iter < 3; axis_iter++)
    {
        for (unsigned stencil_iter = 0; stencil_iter < array_size; stencil_iter++)
        {
            stencil_pointers[axis_iter][stencil_iter] = &this->vectorFields[1ULL*this->numVectorComponents*running_indices[axis_iter][stencil_iter]];
        }

        // Get extra positions that are required for the wilson loops, but need to take periodicity into account.
        int sub_array_sizes[2] = {static_cast<int>(this->ny*this->nz), 1};
        int moduli[2] = {static_cast<int>(this->nx), static_cast<int>(this->nz)};
        int locs[2] = {static_cast<int>(x_loc), static_cast<int>(z_loc)};
        if (axis_iter == 0)
        {
            sub_array_sizes[0] = static_cast<int>(this->nz);
            moduli[0] = static_cast<int>(this->ny);
            locs[0] = static_cast<int>(y_loc);
        }
        else if (axis_iter == 2)
        {
            sub_array_sizes[1] = static_cast<int>(this->nz);
            moduli[1] = static_cast<int>(this->ny);
            locs[1] = static_cast<int>(y_loc);
        }

        // Subtract off what has already been adding during the running sum and deal with the periodicity when decreasing by 1.
        stencil_pointers[axis_iter][3] = &this->vectorFields[1ULL*this->numVectorComponents*(
            running_indices[axis_iter][2] + ((moduli[0] + locs[0] - 1)%moduli[0] - locs[0])*sub_array_sizes[0]
        )];
        stencil_pointers[axis_iter][4] = &this->vectorFields[1ULL*this->numVectorComponents*(
            running_indices[axis_iter][2] + ((moduli[1] + locs[1] - 1)%moduli[1] - locs[1])*sub_array_sizes[1]
        )];

    }    

    return stencil_pointers;
}


float *Periodic::getScalarFieldPointer(const unsigned &t_step, const unsigned &x_loc, const unsigned &y_loc, const unsigned &z_loc)
{
    long long unsigned array_location = ( ( (t_step*this->nx + x_loc)*this->ny + y_loc )*this->nz + z_loc)*this->numScalarComponents;
    
    return &this->scalarFields[array_location];
}


float *Periodic::getVectorFieldPointer(const unsigned &t_step, const unsigned &x_loc, const unsigned &y_loc, const unsigned &z_loc)
{
    long long unsigned array_location = ( ( (t_step*this->nx + x_loc)*this->ny + y_loc )*this->nz + z_loc)*this->numVectorComponents;
    
    return &this->vectorFields[array_location];
}

////////////////////////////////////////////////  Constructors/Destructors  //////////////////////////////////////////////////////

Periodic::Periodic(std::vector<float> &scalar_fields, std::vector<float> &vector_fields, std::vector<Analyser *> &analysers,
                   const unsigned &num_scalar_components, const unsigned &num_vector_components,
                   Model &model, const unsigned &nx, const unsigned &ny, const unsigned &nz, const double &dt,
                   const std::vector<int> bound_vector)
    : BoundaryCondition(scalar_fields, vector_fields, analysers, num_scalar_components, num_vector_components, model, nx, ny, nz, dt, bound_vector)
{
    this->configure(std::string(SOURCE_DIR) + "/Config/Periodic.cfg", false);
}

Periodic::~Periodic()
{
    this->configure(std::string(SOURCE_DIR) + "/Config/Periodic.cfg");
}

///////////////////////////////////////////////////  Public Functions  ///////////////////////////////////////////////////////////

void Periodic::evolve(const unsigned &t_now, const unsigned &stencil_size)
{
    const unsigned t_past = !t_now;

    // Save computation by calculating the indices at lowest loop levels possible.
    std::vector< std::vector<long long unsigned> > t_running_indices(
        3, std::vector<long long unsigned>(
            2*stencil_size + 1, t_now*this->nx*this->ny*this->nz
        )
    ); 

    // Loop over all grid site that this instance is responsible for.
    for (unsigned x_iter = this->loopLimits[0][0]; x_iter < this->loopLimits[0][1]; x_iter++)
    {
        // Add x contributions to the running indices
        std::vector<std::vector<long long unsigned>> x_running_indices = this->updateRunningIndices(t_running_indices, x_iter, 0);

        for (unsigned y_iter = this->loopLimits[1][0]; y_iter < this->loopLimits[1][1]; y_iter++)
        {
            // Add y contributions to the running indices
            std::vector<std::vector<long long unsigned>> y_running_indices = this->updateRunningIndices(x_running_indices, y_iter, 1);

            for (unsigned z_iter = this->loopLimits[2][0]; z_iter < this->loopLimits[2][1]; z_iter++)
            {
                // Get the pointers to this grid point at both timesteps
                // Perhaps this can be done faster using the running indices method too?
                float* local_scalar_pointers[2] = {this->getScalarFieldPointer(t_past, x_iter, y_iter, z_iter),
                                                   this->getScalarFieldPointer(t_now, x_iter, y_iter, z_iter)};

                float* local_vector_pointers[2] = {this->getVectorFieldPointer(t_past, x_iter, y_iter, z_iter),
                                                   this->getVectorFieldPointer(t_now, x_iter, y_iter, z_iter)};

                std::vector<std::vector<long long unsigned>> z_running_indices = this->updateRunningIndices(y_running_indices, z_iter, 2);

                // Generate array of pointers to locations in the grid used by the gradient stencil.
                // Also adds the final z contribution to the indices internally.
                std::vector<std::vector<const float*>> scalar_pointers = this->generateScalarPointers(z_running_indices);

                // Generate array of pointers to the vector fields at locations in the grid used by both the gradient and wilson loop stencils.
                // Also adds the final z contribution and the extras needed for the wilson loops.
                std::vector<std::vector<const float*>> vector_pointers = this->generateVectorPointers(z_running_indices, x_iter, y_iter, z_iter);

                // Process different contributions to the equations of motion.
                this->model.calcPotentialContributions(local_scalar_pointers[1]);
                this->model.calcGradientContributions(scalar_pointers, vector_pointers);
                this->model.calcYangMillsContributions(vector_pointers, local_vector_pointers);     

                // Run all analyser functions that need to happen at every location in the grid, before the fields are evolved.
                // These analysers have access to the past and present timesteps.
                for (auto analyser : this->analysers)
                    analyser->preEvolveLocationAnalysis(1ULL*((x_iter*this->ny + y_iter)*this->nz + z_iter), 
                                                        local_scalar_pointers, scalar_pointers,
                                                        local_vector_pointers, vector_pointers);

                // Calculate the fields at the next timestep.
                this->model.evolve(local_scalar_pointers, local_vector_pointers, this->dt,
                                   this->numScalarComponents, this->numVectorComponents);
            
                // Run all analyser functions that need to happen at every location in the grid, after the fields are evolved.
                // These analysers have access to the present and future timesteps.
                for (auto analyser : this->analysers)
                    analyser->postEvolveLocationAnalysis(t_now, 1ULL*((x_iter*this->ny + y_iter)*this->nz + z_iter),
                                                         local_scalar_pointers, scalar_pointers,
                                                         local_vector_pointers, vector_pointers);

            }
        }
    }
}
