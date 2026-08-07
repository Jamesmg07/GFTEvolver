#include "Lattice.hpp"
#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////////////////
//                                   Private                                              //
////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////  Initialisers  ///////////////////////////////////////////

void Lattice::initVariables()
{
    this->nx = 0;
    this->ny = 0;
    this->nz = 0;
    this->nt = 0;

    this->dx = 0;
    this->dy = 0;
    this->dz = 0;
    this->dt = 0;

    this->numScalarFieldComponents = 0;
    this->numVectorFieldComponents = 0;

    this->scalarInitialConditions = nullptr;
    this->vectorInitialConditions = nullptr;

    this->stencilSize = 0;

    this->progressReport = false;
}

void Lattice::configure(const std::string path, const bool debug)
{
    std::ifstream ifs(path);

    // Store choices for debug
    std::vector<int> initial_condition_types(2, 0);
    std::vector<int> face_boundary_condition_types(6, 0);
    std::vector<bool> analysis_choices(3, false);

    // Read in the parameter values from the configuration file
    std::string description; // Descriptions in the config file (to be dumped)
    if (ifs.is_open())
    {
        // Skip the initial configuration file description
        for(unsigned iter = 0; iter < 10; iter++) std::getline(ifs, description);

        //////////  Grid sizes and spacings  //////////
        std::getline(ifs, description, ':');
        ifs >> this->nx >> this->ny >> this->nz;

        std::getline(ifs, description, ':');
        ifs >> this->dx >> this->dy >> this->dz;

        std::getline(ifs, description, ':');
        ifs >> this->nt;

        std::getline(ifs, description, ':');
        ifs >> this->dt;


        /////////////  Field components  ///////////////
        std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> this->numScalarFieldComponents;

        std::getline(ifs, description, ':');
        ifs >> this->numVectorFieldComponents;
        this->numVectorFieldComponents = 3*this->numVectorFieldComponents; // For the three spatial dimensions.


        // Optionally overide nt by setting it based on the light-crossing time
        // in the x direction. Rounds down.
        bool light_crossing_override;
        for (unsigned iter = 0; iter < 2; iter++) std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> light_crossing_override;

        if (light_crossing_override)
        {
            float num_light_crossings;
            std::getline(ifs, description, ':');
            ifs >> num_light_crossings;

            this->nt = static_cast<int>( num_light_crossings*static_cast<float>(nx)*dx/(2*dt) );
        }


        /////////  Choose initial conditions types  //////////////
        for (unsigned iter = 0; iter < 3; iter++) std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> initial_condition_types[0] >> initial_condition_types[1];
        this->initInitialCondition(initial_condition_types);


        /////////  Choose analyses and assign analysers.  ////////////
        bool analysis_choice;

        // Output fields?
        std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> analysis_choice;
        analysis_choices[0] = analysis_choice;
        if (analysis_choice)
            this->analysers.push_back( new OutputFields(this->scalarFields, this->numScalarFieldComponents,
                                                        this->vectorFields, this->numVectorFieldComponents)
                                     );

        // Output energy?
        std::getline(ifs, description, ':');
        ifs >> analysis_choice;
        analysis_choices[1] = analysis_choice;
        if (analysis_choice)
            this->analysers.push_back( new Energy(this->model, this->dx, this->dy, this->dz, 1ULL*this->nx*this->ny*this->nz) );

        // Output gauge condition violation?
        std::getline(ifs, description, ':');
        ifs >> analysis_choice;
        analysis_choices[2] = analysis_choice;
        if (analysis_choice)
            this->analysers.push_back( new GaugeCondition(this->model, this->dx, this->dy, this->dz, 
                                                          1ULL*this->nx*this->ny*this->nz, this->numVectorFieldComponents) );


        // Choose boundary condition types. Only the six faces are assigned by
        // the user; edge and corner behaviour is inferred automatically.
        for (unsigned iter = 0; iter < 7; iter++) std::getline(ifs, description);
        std::getline(ifs, description, ':');
        for (unsigned iter = 0; iter < 6; iter++) ifs >> face_boundary_condition_types[iter];

        // Instantiate boundary conditions
        this->initBoundaryCondition(face_boundary_condition_types);


        // Choose whether to output progress throughout the simulation.
        std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> this->progressReport;

        // And how often to print out reports.
        std::getline(ifs, description, ':');
        ifs >> this->reportFrequency;

    }

    ifs.close();

    if (debug)
    {
        std::cout << "LATTICE::\n"
                  << "nx: " << this->nx << ", ny: " << this->ny << ", nz: " << this->nz << ", nt: " << this->nt << "\n"
                  << "dx: " << this->dx << ", dy: " << this->dy << ", dz: " << this->dz << ", dt: " << this->dt << "\n"
                  << "#Scalar field components: " << numScalarFieldComponents
                  << ", #Vector field components: " << numVectorFieldComponents << "\n"
                  << "Scalar initial condition type: " << initial_condition_types[0]
                  << ", Vector initial condition type: " << initial_condition_types[1] << "\n"
                  << "Face boundary condition types:";

        for (unsigned iter = 0; iter < 6; iter++)
            std::cout << " " << face_boundary_condition_types[iter];

        std::cout << "\nOutput fields?: " << analysis_choices[0]
                  << "\nOutput energy?: " << analysis_choices[1]
                  << "\nOutput gauge condition violation?: " << analysis_choices[2]
                  << std::endl;

        std::cout << "Report progress?: " << this->progressReport << ", every " << this->reportFrequency << " timesteps.\n" << std::endl;
    }
}

void Lattice::initInitialCondition(const std::vector<int> &initial_condition_types)
{
    switch (initial_condition_types[0])
    {
    case ZERO_IC:

        std::cout << "LATTICE::WARNING: Scalar fields are being initialised to zero (default behaviour).\n"
                  << "Please ensure that this is the intended initial condition.\n" << std::endl;
        // Nothing more needs to be done as the vectors are initialised this way when memory is reserved through resize method.
        break;

    case RANDOM_UNIFORM_IC:

        this->scalarInitialConditions = new RandomUniform();
        break;

    default:

        std::cout << "LATTICE::WARNING: Initial condition type (" << initial_condition_types[0] <<  ") for the scalar fields is invalid.\n"
                  << "Reverting to default behaviour (initialisation to zero)\n" << std::endl; 

    }

    switch (initial_condition_types[1])
    {
    case ZERO_IC:

        std::cout << "LATTICE::WARNING: Gauge fields are being initialised to zero (default behaviour).\n "
                  << "Please ensure that this is the intended initial condition.\n" << std::endl;
        // Nothing more needs to be done as the vectors are initialised this way when memory is reserved through resize method.
        break;

    case RANDOM_UNIFORM_IC:

        this->vectorInitialConditions = new RandomUniform();
        break;

    default:

        std::cout << "LATTICE::WARNING: Initial condition type (" << initial_condition_types[1] << ") for the gauge fields is invalid.\n"
                  << "Reverting to default behaviour (initialisation to zero)\n" << std::endl; 

    }
}

void Lattice::initBoundaryCondition(const std::vector<int> &face_boundary_condition_types)
{

    // The user-facing interface contains exactly the six faces in the order
    // (-x), (+x), (-y), (+y), (-z), (+z).
    if (face_boundary_condition_types.size() != 6){
        throw std::runtime_error(
            "LATTICE:: Exactly six face boundary conditions must be assigned."
        );
    }

    // Automatic inference currently supports the two boundary conditions that
    // are fully implemented in V1. Neumann will be enabled separately once a
    // gauge-covariant implementation exists.
    for (unsigned face_iter = 0; face_iter < 6; face_iter++)
    {
        if (face_boundary_condition_types[face_iter] != FIXED &&
            face_boundary_condition_types[face_iter] != PERIODIC)
        {
            throw std::runtime_error(
                "LATTICE:: Face boundary condition (" +
                std::to_string(face_iter) +
                ") has unsupported type " +
                std::to_string(face_boundary_condition_types[face_iter]) + "."
            );
        }
    }

    // Periodicity is a property of a complete coordinate direction, so a
    // periodic face must always be paired with its opposite face.
    for (unsigned axis = 0; axis < 3; axis++)
    {
        const bool lower_is_periodic =
            face_boundary_condition_types[2*axis] == PERIODIC;

        const bool upper_is_periodic =
            face_boundary_condition_types[2*axis + 1] == PERIODIC;

        if (lower_is_periodic != upper_is_periodic)
            throw std::runtime_error(
                "LATTICE:: Periodic boundary conditions must be assigned "
                "to both faces of an axis."
            );
    }

    for (unsigned which_boundary = 0; which_boundary < 26; which_boundary++)
    {
        
        // which_boundary encodes the position and dimensionality of each boundary.
        // It is convenient to convert this into a vector that points perpendicularly outside the grid.
        std::vector<int> bound_vector(3, 0);
        int alternator = -1 + 2*(which_boundary%2);
        unsigned skip_counter = which_boundary/2;
        int skip_alternator = 1 - 2*(skip_counter%2);

        if (which_boundary < 6)                                    // Face boundaries
        {
            bound_vector[skip_counter] = alternator;
        }
        else if (which_boundary < 18)                              // Edge boundaries
        {
            unsigned skip10_counter = 1 + which_boundary/10;
            unsigned skip14_counter = which_boundary/14;
            bound_vector[skip10_counter] = alternator;
            bound_vector[skip14_counter] = skip_alternator;
        }
        else                                                       // Corner boundaries
        {
            int skip22_alternator = -1 + 2*((which_boundary)/22);
            bound_vector[0] = skip22_alternator;
            bound_vector[1] = skip_alternator;
            bound_vector[2] = alternator;
        }


        // Work out which of the six user-assigned faces meet at this region.
        // With the current Fixed/Periodic pair, a fixed face freezes the whole
        // region. Otherwise every incident face is periodic and the region is
        // periodic.
        int boundary_condition_type = PERIODIC;

        for (unsigned axis = 0; axis < 3; axis++)
        {
            if (bound_vector[axis] == 0)
                continue;

            const unsigned face_index =
                2*axis + (bound_vector[axis] > 0 ? 1 : 0);

            if (face_boundary_condition_types[face_index] == FIXED)
            {
                boundary_condition_type = FIXED;
                break;
            }
        }




        // Instantiate a boundary condition child class based upon the chosen type and push it back into the boundary conditions vector.
        switch (boundary_condition_type)
        {
        case UNASSIGNED_BOUNDARY_CONDITION:
            
            std::cout << "LATTICE::ERROR: Boundary condition (" << bound_vector[0] << " " << bound_vector[1]
                      << " " << bound_vector[2] << ") has not been assigned." << std::endl;
            break;
        
        case FIXED:

            this->boundaryConditions.push_back( new Fixed(this->scalarFields, this->vectorFields, this->analysers,
                                                          this->numScalarFieldComponents, this->numVectorFieldComponents,
                                                          this->model, this->nx, this->ny, this->nz, this->dt,
                                                          bound_vector)
                                              );
            break;

        case NEUMANN:

            std::cout << "LATTICE::ERROR: Boundary condition (" << bound_vector[0] << " " << bound_vector[1]
                      << " " << bound_vector[2] << ") has been chosen to be neumann, but this feature has not been added yet." << std::endl;
            break;

        case PERIODIC:

            this->boundaryConditions.push_back( new Periodic(this->scalarFields, this->vectorFields, this->analysers,
                                                             this->numScalarFieldComponents, this->numVectorFieldComponents,
                                                             this->model, this->nx, this->ny, this->nz, this->dt,
                                                             bound_vector)
                                              );
            break;

        default:
            
            std::cout << "LATTICE::ERROR: Boundary condition (" << bound_vector[0] << " " << bound_vector[1]
                      << " " << bound_vector[2] << ") has been given an invalid type: " << boundary_condition_type << "." << std::endl;
        }

    }
}

void Lattice::initFields()
{
    // (long long int) 2 both converts the whole expression to long long int 
    // and allows two timesteps to be stored.
    long long unsigned scalar_array_size = 2ULL*this->nx*this->ny*this->nz*this->numScalarFieldComponents;
    this->scalarFields.resize(scalar_array_size, 0.f);

    // Same as above (numVectorFieldComponents contains 3*the amount requested to account for spatial components of each)
    long long unsigned vector_array_size = 2ULL*this->nx*this->ny*this->nz*this->numVectorFieldComponents;
    this->vectorFields.resize(vector_array_size, 0.f);

    // A zero gauge field is stored as the SU(2) identity in quaternion form.
    if (this->model.isUsingQuaternionRepresentation())
    {
        const unsigned direction_width
            = this->numVectorFieldComponents/3U;

        for (long long unsigned loc_iter = 0;
            loc_iter < 2ULL*this->nx*this->ny*this->nz;
            loc_iter++)
        {
            for (unsigned dir_iter = 0; dir_iter < 3U; dir_iter++)
            {
                const long long unsigned direction_index
                    = loc_iter*this->numVectorFieldComponents
                    + 1ULL*dir_iter*direction_width;

                // Electroweak SU(2) identity: c0 = 1.
                this->vectorFields[direction_index + 1ULL] = 1.f;

                // The 9-component layout contains a second SU(2) factor.
                if (direction_width == 9U)
                {
                    this->vectorFields[direction_index + 5ULL] = 1.f;
                }
            }
        }
    }

    // Generate the initial conditions
    if(scalarInitialConditions)
        scalarInitialConditions->setInitialFields(scalarFields);
    if(vectorInitialConditions)
        vectorInitialConditions->setInitialFields(vectorFields);

    if (this->progressReport)
    {
        std::cout << "Initial conditions generated.\n" << std::endl;
    }
}

////////////////////////////////////////////////  Private Functions  /////////////////////////////////////////////////////

unsigned Lattice::getDefaultStencilSize() const
{
    return this->model.getDefaultStencilSize();
    // May need modification to also check other stencils later on, e.g from wilson loops.
}

std::vector< std::vector<unsigned> > Lattice::determineResponsibilities(const unsigned &stencil_size) const
{
    std::vector< std::vector<unsigned> > loop_limits(3, std::vector<unsigned>(2, stencil_size));

    // Deal with edge cases.
    if (this->nx < stencil_size)
        loop_limits[0][1] = 0;
    else
        loop_limits[0][1] = this->nx - stencil_size;

    if (this->ny < stencil_size)
        loop_limits[1][1] = 0;
    else
        loop_limits[1][1] = this->ny - stencil_size;

    if (this->nz < stencil_size)
        loop_limits[2][1] = 0;
    else
        loop_limits[2][1] = this->nz - stencil_size;

    return loop_limits;
}


void Lattice::postEvolveAnalysis(
    const unsigned &t_now,
    const std::vector<std::vector<unsigned>> &loop_limits)
{
    const unsigned t_past = !t_now;

    // Analyse all dynamic boundary sites.
    for (auto bound : this->boundaryConditions)
    {
        bound->postEvolveAnalysis(
            t_now,
            this->stencilSize);
    }

    // Calculate the initial indices for the present timestep.
    std::vector<std::vector<long long unsigned>> t_running_indices(
        3,
        std::vector<long long unsigned>(
            2*this->stencilSize + 1,
            t_now*this->nx*this->ny*this->nz
        )
    );

    // Analyse all fully interior grid sites.
    for (unsigned x_iter = loop_limits[0][0];
         x_iter < loop_limits[0][1];
         x_iter++)
    {
        std::vector<std::vector<long long unsigned>>
            x_running_indices
                = this->updateRunningIndices(
                    t_running_indices,
                    x_iter,
                    0);

        for (unsigned y_iter = loop_limits[1][0];
             y_iter < loop_limits[1][1];
             y_iter++)
        {
            std::vector<std::vector<long long unsigned>>
                y_running_indices
                    = this->updateRunningIndices(
                        x_running_indices,
                        y_iter,
                        1);

            for (unsigned z_iter = loop_limits[2][0];
                 z_iter < loop_limits[2][1];
                 z_iter++)
            {
                float* local_scalar_pointers[2]
                    = {
                        this->getScalarFieldPointer(
                            t_past, x_iter, y_iter, z_iter),
                        this->getScalarFieldPointer(
                            t_now, x_iter, y_iter, z_iter)
                    };

                float* local_vector_pointers[2]
                    = {
                        this->getVectorFieldPointer(
                            t_past, x_iter, y_iter, z_iter),
                        this->getVectorFieldPointer(
                            t_now, x_iter, y_iter, z_iter)
                    };

                std::vector<std::vector<long long unsigned>>
                    z_running_indices
                        = this->updateRunningIndices(
                            y_running_indices,
                            z_iter,
                            2);

                std::vector<std::vector<const float*>>
                    scalar_pointers
                        = this->generateScalarPointers(
                            z_running_indices);

                std::vector<std::vector<const float*>>
                    vector_pointers
                        = this->generateVectorPointers(
                            z_running_indices);

                const long long unsigned index
                    = 1ULL*((x_iter*this->ny + y_iter)*this->nz
                            + z_iter);

                for (auto analyser : this->analysers)
                {
                    analyser->postEvolveLocationAnalysis(
                        t_now,
                        index,
                        local_scalar_pointers,
                        scalar_pointers,
                        local_vector_pointers,
                        vector_pointers);
                }
            }
        }
    }
}


std::vector<std::vector<long long unsigned>> Lattice::updateRunningIndices(std::vector<std::vector<long long unsigned>> running_indices,
                                                                           const unsigned &loc, const unsigned axis) const
{
    unsigned stencil_size = static_cast<unsigned>((running_indices[0].size()-1)/2ULL);
    long long unsigned sub_array_size;
    if (axis == 0)
        sub_array_size = 1ULL*this->ny*this->nz;
    else if (axis == 1)
        sub_array_size = 1ULL*this->nz;
    else if (axis == 2)
        sub_array_size = 1;
    else
        std::cout << "LATTICE::ERROR: Invalid axis argument in call to modifyRunningIndices." << std::endl;

    for (unsigned axis_iter = 0; axis_iter < 3; axis_iter++)
    {
        for (unsigned stencil_iter = 0; stencil_iter < 2*stencil_size + 1; stencil_iter++)
        {
            if (axis_iter == axis)
                running_indices[axis_iter][stencil_iter] += 1ULL*(loc + stencil_iter - stencil_size)*sub_array_size;
            else
                running_indices[axis_iter][stencil_iter] += 1ULL*loc*sub_array_size;
        }
    }

    return running_indices;
}

std::vector<std::vector<const float *>> Lattice::generateScalarPointers(const std::vector<std::vector<long long unsigned>> &running_indices) const
{
    unsigned array_size = static_cast<unsigned>(running_indices[0].size());
    std::vector<std::vector<const float*>> stencil_pointers(3, std::vector<const float*>(array_size, nullptr));

    for (unsigned axis_iter = 0; axis_iter < 3; axis_iter++)
    {
        for (unsigned stencil_iter = 0; stencil_iter < array_size; stencil_iter++)
        {
            stencil_pointers[axis_iter][stencil_iter] = &this->scalarFields[1ULL*this->numScalarFieldComponents*running_indices[axis_iter][stencil_iter]];
        }
    }

    return stencil_pointers;
}

// std::vector<std::vector<const float *>> Lattice::generateScalarPointers(const std::vector<std::vector<long long unsigned>> &running_indices,
//                                                                    const unsigned &loc, const unsigned axis) const
// {
//     unsigned stencil_size = (running_indices[0].size()-1)/2;
//     long long unsigned sub_array_size;
//     if (axis == 0)
//         sub_array_size = 1ULL*this->ny*this->nz;
//     else if (axis == 1)
//         sub_array_size = 1ULL*this->nz;
//     else if (axis == 2)
//         sub_array_size = 1ULL;
//     else
//         std::cout << "LATTICE::ERROR: Invalid axis argument in call to generateScalarPointers." << std::endl;

//     std::vector<std::vector<const float*>> stencil_pointers(3, std::vector<const float*>(2*stencil_size + 1, nullptr));

//     for (unsigned axis_iter = 0; axis_iter < 3; axis_iter++)
//     {
//         for (unsigned stencil_iter = 0; stencil_iter < 2*stencil_size + 1; stencil_iter++)
//         {
//             if (axis_iter == axis)
//             {
//                 stencil_pointers[axis_iter][stencil_iter] = &this->scalarFields[
//                     1ULL*this->numScalarFieldComponents*( running_indices[axis_iter][stencil_iter] + 1ULL*(loc + stencil_iter - stencil_size)*sub_array_size )
//                 ];
//             }
//             else
//             {
//                 stencil_pointers[axis_iter][stencil_iter] = &this->scalarFields[
//                     1ULL*this->numScalarFieldComponents*running_indices[axis_iter][stencil_iter] + loc*sub_array_size
//                 ];
//             }
//         }
//     }

//     return stencil_pointers;
// }

std::vector<std::vector<const float *>> Lattice::generateVectorPointers(const std::vector<std::vector<long long unsigned>> &running_indices) const
{

    // STENCIL SIZES ARE NOT REALLY WORKING WITH THIS PROPERLY.
    // IN PRINCIPLE I CAN HAVE A DIFFERENT WILSON LOOP STENCIL TO THE DERIVATIVE STENCIL
    // THIS FUNCTION ASSUMES THAT THE WILSON LOOP STENCIL IS THE MOST BASIC TYPE

    unsigned array_size = static_cast<unsigned>(running_indices[0].size());
    std::vector<std::vector<const float*>> stencil_pointers(3, std::vector<const float*>(array_size + 2, nullptr));

    for (unsigned axis_iter = 0; axis_iter < 3; axis_iter++)
    {
        for (unsigned stencil_iter = 0; stencil_iter < array_size; stencil_iter++)
        {
            stencil_pointers[axis_iter][stencil_iter] = &this->vectorFields[1ULL*this->numVectorFieldComponents*running_indices[axis_iter][stencil_iter]];
        }

        // Get extra positions that are required for the wilson loops.
        long long unsigned sub_array_sizes[2] = {1ULL*this->ny*this->nz, 1ULL};
        if (axis_iter == 0)
            sub_array_sizes[0] = 1ULL*this->nz;
        else if (axis_iter == 2)
            sub_array_sizes[1] = 1ULL*this->nz;

        stencil_pointers[axis_iter][3] = &this->vectorFields[1ULL*this->numVectorFieldComponents*(running_indices[axis_iter][2] - sub_array_sizes[0])];
        stencil_pointers[axis_iter][4] = &this->vectorFields[1ULL*this->numVectorFieldComponents*(running_indices[axis_iter][2] - sub_array_sizes[1])];
    }

    return stencil_pointers;
}


float *Lattice::getScalarFieldPointer(const unsigned &t_step, const unsigned &x_loc, const unsigned &y_loc, const unsigned &z_loc)
{
    long long unsigned array_location = ( ( (t_step*this->nx + x_loc)*this->ny + y_loc )*this->nz + z_loc)*this->numScalarFieldComponents;

    return &this->scalarFields[array_location];
}


float *Lattice::getVectorFieldPointer(const unsigned &t_step, const unsigned &x_loc, const unsigned &y_loc, const unsigned &z_loc)
{
    long long unsigned array_location = ( ( (t_step*this->nx + x_loc)*this->ny + y_loc )*this->nz + z_loc)*this->numVectorFieldComponents;

    return &this->vectorFields[array_location];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                    Public                                                            //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////  Constructors/Destructors  //////////////////////////////////////////////////

Lattice::Lattice()
{
    // Initialise all variables to placeholder/default values
    this->initVariables();

    // Configure the lattice.
    this->configure(std::string(SOURCE_DIR) + "/Config/Lattice.cfg", true);
    std::cout << "Lattice configure done" << std::endl;
    
    // Configure the model (potential, gradients and wilson loops).
    this->model.configure(
        std::string(SOURCE_DIR) + "/Config/Model.cfg", 
        this->numScalarFieldComponents, this->numVectorFieldComponents,
        this->nx, this->ny, this->nz,
        this->dx, this->dy, this->dz, 
        this->dt, true
    );

    // Initialise the fields that live on the lattice.
    this->initFields();
    std::cout << "Lattice initFields done" << std::endl;

    // Determine the largest default stencil size
    this->stencilSize = this->getDefaultStencilSize();

    // PERHAPS REPLACE WITH LATTICE VERSION OF DETERMINERESPONSIBILITIES AND THIS LOOP CAN BE WITHIN THAT FUNCTION...

    // Determine the responsibilities of the boundary conditions.
    // Also check for boundaries with no responsibilities (common in 1/2D sims for example) and remove them.
    for (auto it = this->boundaryConditions.begin(); it != this->boundaryConditions.end(); )
    {
        bool is_empty = (*it)->determineResponsibilities(this->stencilSize);
        if (is_empty)
        {
            delete *it;
            it = this->boundaryConditions.erase(it);
        }
        else
        {
            ++it;
        }
    }

}

Lattice::~Lattice()
{
    // Avoid memory leaks.

    delete this->scalarInitialConditions;
    delete this->vectorInitialConditions;
    this->scalarInitialConditions = nullptr;
    this->vectorInitialConditions = nullptr;

    for (auto bound : this->boundaryConditions)
    {
        delete bound;
    }
    this->boundaryConditions.clear();

    for (auto analyser : this->analysers)
    {
        delete analyser;
    }
    this->analysers.clear();
}

/////////////////////////////////////////////////  Public Functions  //////////////////////////////////////////////////////

void Lattice::initialAnalysis() const
{
    for (auto analyser : this->analysers)
        if (analyser)
            analyser->initialAnalysis();
}

void Lattice::evolve()
{
    // Get start time of evolution stage if reporting progress.
    std::chrono::steady_clock::time_point start_time;
    if (this->progressReport) 
        start_time = std::chrono::steady_clock::now();


    // Determine how much of the grid is interior and can be evolved in the default manner
    // and how many the boundary condition classes will be responsible for evolving.
    std::vector< std::vector<unsigned> > loop_limits = this->determineResponsibilities(this->stencilSize);

    if (this->progressReport)   
        std::cout << "Beginning evolution." << std::endl;

    for (unsigned time_iter = 0; time_iter < this->nt; time_iter++)
    {
        // Track which parts of the array are the current timestep and which are the previous.
        unsigned t_now = (time_iter+1)%2;
        unsigned t_past = !t_now;

        // Update model internal evolution parameters that depend upon the timestep
        model.update(time_iter);

        //MPI: send here without blocking, move boundary evolution to later and recv just before.

         // Evolve grid points that are near/on boundaries
        for (auto bound : this->boundaryConditions)
            bound->evolve(t_now, this->stencilSize);

        // Save computation by calculating the indices at lowest loop levels possible.
        std::vector< std::vector<long long unsigned> > t_running_indices(
            3, std::vector<long long unsigned>(
                2*this->stencilSize + 1, t_now*this->nx*this->ny*this->nz
            )
        ); 
        
        // Loop over all fully interior grid sites.
        for (unsigned x_iter = loop_limits[0][0]; x_iter < loop_limits[0][1]; x_iter++)
        {
            // Add x contributions to the running indices
            std::vector<std::vector<long long unsigned>> x_running_indices = this->updateRunningIndices(t_running_indices, x_iter, 0);

            for (unsigned y_iter = loop_limits[1][0]; y_iter < loop_limits[1][1]; y_iter++)
            {
                // Add y contributions to the running indices
                std::vector<std::vector<long long unsigned>> y_running_indices = this->updateRunningIndices(x_running_indices, y_iter, 1);

                for (unsigned z_iter = loop_limits[2][0]; z_iter < loop_limits[2][1]; z_iter++)
                {
                    // Get the pointers to this grid point at both timesteps
                    float* local_scalar_pointers[2] = {this->getScalarFieldPointer(t_past, x_iter, y_iter, z_iter),
                                                       this->getScalarFieldPointer(t_now, x_iter, y_iter, z_iter)};

                    float* local_vector_pointers[2] = {this->getVectorFieldPointer(t_past, x_iter, y_iter, z_iter),
                                                       this->getVectorFieldPointer(t_now, x_iter, y_iter, z_iter)};

                    // Add z contributions to the running indices
                    std::vector<std::vector<long long unsigned>> z_running_indices = this->updateRunningIndices(y_running_indices, z_iter, 2);

                    // Generate array of pointers to the scalar fields at locations in the grid used by the gradient stencil.
                    // Also adds the final z contribution to the indices internally.
                    std::vector<std::vector<const float*>> scalar_pointers = this->generateScalarPointers(z_running_indices);

                    // Generate array of pointers to the vector fields at locations in the grid used by both the gradient and wilson loop stencils.
                    // Also adds the final z contribution and the extras needed for the wilson loops.
                    std::vector<std::vector<const float*>> vector_pointers = this->generateVectorPointers(z_running_indices);


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
                                       this->numScalarFieldComponents, this->numVectorFieldComponents);
                    
                }
            }
        }

        // Run post-evolution location analyses only after the new
        // timestep is complete across the full dynamic grid.
        this->postEvolveAnalysis(t_now, loop_limits);

        // Run all continous analyser functions that don't need to happen at every location in the dynamic grid.
        for (auto analyser : this->analysers)
            analyser->timestepAnalysis(time_iter);

        // Replace with a function that has some more nice features.
        if (this->progressReport
            && (time_iter + 1)%this->reportFrequency == 0)
        {
            std::cout
                << "Timestep "
                << std::to_string(time_iter + 1)
                << " completed.\r"
                << std::flush;
        }

    }

    if (this->progressReport)
    {
        std::cout << "Timestep " << std::to_string(this->nt) << " completed.\n";

        std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_time = end_time - start_time;
        std::cout << "Evolution finished in " << elapsed_time.count() << " seconds."  << std::endl;
    }

}

void Lattice::finalAnalysis() const
{
    for (auto analyser : this->analysers)
        analyser->finalAnalysis();
}
