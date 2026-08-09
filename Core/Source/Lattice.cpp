#include "Lattice.hpp"
#include <array>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <algorithm>

#ifdef GFT_ENABLE_MPI
#include <mpi.h>
#endif

struct Lattice::HaloExchange
{
#ifdef GFT_ENABLE_MPI
    std::array<MPI_Request, 8> requests;
#endif
    int requestCount = 0;
};

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

    this->performInitialGradientFlow = false;
    this->gradientFlowStepSize = 0.0;
    this->gradientFlowTimesteps = 0U;
    this->stopGradientFlowOnResiduals = false;
    this->scalarEquationResidualTolerance = 0.0;
    this->gaugeEquationResidualTolerance = 0.0;
    this->performDynamicalEvolution = true;
    this->activeEvolutionStep = 0.0;

    this->numScalarFieldComponents = 0;
    this->numVectorFieldComponents = 0;

    this->scalarInitialConditions = nullptr;
    this->vectorInitialConditions = nullptr;

    this->stencilSize = 0;

    this->progressReport = false;

    this->rank = 0;
    this->numRanks = 1;
    this->localNx = 0;
    this->storageNx = 0;
    this->globalXStart = 0;
    this->haloDepth = 0;
    this->ownedXBegin = 0;
    this->ownedXEnd = 0;
    this->leftRank = -1;
    this->rightRank = -1;

    this->initialConditionTypes.assign(2, 0);
    this->faceBoundaryConditionTypes.assign(6, 0);
    this->analysisChoices.assign(3, false);
}

void Lattice::configure(const std::string path, const bool debug)
{
    std::ifstream ifs(path);

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

        // Optional gradient-flow phase
        std::getline(ifs, description, ':');
        ifs >> this->performInitialGradientFlow;

        std::getline(ifs, description, ':');
        ifs >> this->gradientFlowStepSize;

        std::getline(ifs, description, ':');
        ifs >> this->gradientFlowTimesteps;

        std::getline(ifs, description, ':');
        ifs >> this->stopGradientFlowOnResiduals;

        std::getline(ifs, description, ':');
        ifs >> this->scalarEquationResidualTolerance;

        std::getline(ifs, description, ':');
        ifs >> this->gaugeEquationResidualTolerance;

        // Optional physical dynamical phase
        std::getline(ifs, description, ':');
        ifs >> this->performDynamicalEvolution;

        this->activeEvolutionStep = this->dt;

        if (!std::isfinite(this->dt) || this->dt <= 0.0)
        {
            throw std::runtime_error(
                "LATTICE:: The dynamical timestep must be finite and positive."
            );
        }

        if (this->performInitialGradientFlow
            && this->gradientFlowTimesteps > 0U
            && (!std::isfinite(this->gradientFlowStepSize)
                || this->gradientFlowStepSize <= 0.0))
        {
            throw std::runtime_error(
                "LATTICE:: The gradient-flow step size must be finite and positive."
            );
        }

        if (this->performInitialGradientFlow
            && this->stopGradientFlowOnResiduals
            && (!std::isfinite(this->scalarEquationResidualTolerance)
                || this->scalarEquationResidualTolerance <= 0.0
                || !std::isfinite(this->gaugeEquationResidualTolerance)
                || this->gaugeEquationResidualTolerance <= 0.0))
        {
            throw std::runtime_error(
                "LATTICE:: Gradient-flow residual tolerances "
                "must be finite and positive."
            );
        }



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

        if (light_crossing_override && this->performDynamicalEvolution)
        {
            float num_light_crossings;
            std::getline(ifs, description, ':');
            ifs >> num_light_crossings;

            this->nt = static_cast<int>( num_light_crossings*static_cast<float>(nx)*dx/(2*dt) );
        }


        /////////  Choose initial conditions types  //////////////
        for (unsigned iter = 0; iter < 3; iter++) std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> this->initialConditionTypes[0] >> this->initialConditionTypes[1];


        /////////  Choose analyses.  ////////////
        bool analysis_choice;

        // Output fields?
        std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> analysis_choice;
        this->analysisChoices[0] = analysis_choice;

        // Output energy?
        std::getline(ifs, description, ':');
        ifs >> analysis_choice;
        this->analysisChoices[1] = analysis_choice;

        // Output gauge condition violation?
        std::getline(ifs, description, ':');
        ifs >> analysis_choice;
        this->analysisChoices[2] = analysis_choice;

        // Choose boundary condition types. Only the six faces are assigned by
        // the user; edge and corner behaviour is inferred automatically.
        for (unsigned iter = 0; iter < 7; iter++) std::getline(ifs, description);
        std::getline(ifs, description, ':');
        for (unsigned iter = 0; iter < 6; iter++)
            ifs >> this->faceBoundaryConditionTypes[iter];


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
                  << "Scalar initial condition type: " << initialConditionTypes[0]
                  << ", Vector initial condition type: " << initialConditionTypes[1] << "\n"
                  << "Face boundary condition types:";

        for (unsigned iter = 0; iter < 6; iter++){
            std::cout << " " << faceBoundaryConditionTypes[iter];
        }
        std::cout << "\n";

        std::cout
                << "Initial gradient flow?: "
                << this->performInitialGradientFlow
                << ", step size: "
                << this->gradientFlowStepSize
                << ", steps: "
                << this->gradientFlowTimesteps
                << ", residual stopping?: "
                << this->stopGradientFlowOnResiduals
                << ", scalar tolerance: "
                << this->scalarEquationResidualTolerance
                << ", gauge tolerance: "
                << this->gaugeEquationResidualTolerance
                << "\nPerform dynamical evolution?: "
                << this->performDynamicalEvolution
                << "\n";
        
        std::cout << "\nOutput fields?: " << analysisChoices[0]
                  << "\nOutput energy?: " << analysisChoices[1]
                  << "\nOutput gauge condition violation?: " << analysisChoices[2]
                  << "\n"
                  << std::endl;

        std::cout << "Report progress?: " << this->progressReport << ", every " << this->reportFrequency << " timesteps.\n" << std::endl;
    }
}

void Lattice::initAnalysers(const std::vector<bool> &analysis_choices)
{
    if (analysis_choices.size() != 3)
        throw std::runtime_error(
            "LATTICE:: Exactly three analyser choices must be assigned."
        );

    if (analysis_choices[0])
        this->analysers.push_back(
            new OutputFields(
                this->scalarFields,
                this->numScalarFieldComponents,
                this->vectorFields,
                this->numVectorFieldComponents,
                1ULL*this->storageNx*this->ny*this->nz,
                1ULL*this->ownedXBegin*this->ny*this->nz,
                1ULL*this->ownedXEnd*this->ny*this->nz,
                this->rank,
                this->numRanks)
        );

    if (analysis_choices[1])
        this->analysers.push_back(
            new Energy(
            this->model,
            this->dx,
            this->dy,
            this->dz,
            1ULL*this->storageNx*this->ny*this->nz,
            1ULL*this->ownedXBegin*this->ny*this->nz,
            1ULL*this->ownedXEnd*this->ny*this->nz,
            this->rank,
            this->numRanks)
        );

    if (analysis_choices[2])
        this->analysers.push_back(
            new GaugeCondition(
            this->model,
            this->dx,
            this->dy,
            this->dz,
            1ULL*this->storageNx*this->ny*this->nz,
            this->numVectorFieldComponents,
            1ULL*this->ownedXBegin*this->ny*this->nz,
            1ULL*this->ownedXEnd*this->ny*this->nz,
            this->rank,
            this->numRanks)
        );
}

void Lattice::initSlabGeometry()
{
    const unsigned process_count = static_cast<unsigned>(this->numRanks);
    const unsigned process_rank = static_cast<unsigned>(this->rank);

    // Divide nx as evenly as possible. The first remainder ranks receive
    // one additional physical x-plane.
    const unsigned base_planes = this->nx/process_count;
    const unsigned remainder = this->nx%process_count;

    this->localNx = base_planes + (process_rank < remainder ? 1U : 0U);
    this->globalXStart = process_rank*base_planes
                       + (process_rank < remainder ? process_rank : remainder);

    // Preserve the original serial storage layout exactly. For more than one
    // rank, reserve enough x-width for the complete operator footprint.
    this->haloDepth = (this->numRanks > 1) ? this->stencilSize : 0U;

    if (this->localNx < 2U*this->haloDepth)
    {
        throw std::runtime_error(
            "LATTICE:: MPI x-slab is too narrow for the required stencil."
        );
    }

    this->storageNx = this->localNx + 2U*this->haloDepth;
    this->ownedXBegin = this->haloDepth;
    this->ownedXEnd = this->ownedXBegin + this->localNx;

    // Direct process neighbours always exist between adjacent slabs. At the
    // two global x-ends they exist only when x is globally periodic.
    this->leftRank = -1;
    this->rightRank = -1;

    if (this->numRanks > 1)
    {
        const bool periodic_x =
            this->faceBoundaryConditionTypes[0] == PERIODIC &&
            this->faceBoundaryConditionTypes[1] == PERIODIC;

        if (this->rank > 0)
            this->leftRank = this->rank - 1;
        else if (periodic_x)
            this->leftRank = this->numRanks - 1;

        if (this->rank < this->numRanks - 1)
            this->rightRank = this->rank + 1;
        else if (periodic_x)
            this->rightRank = 0;
    }

    // Local consistency checks 
    if (this->globalXStart + this->localNx > this->nx)
        throw std::runtime_error("LATTICE:: Invalid MPI x-slab decomposition.");

    if (this->localNx > 0 &&
        (this->localToGlobalX(this->ownedXBegin) != this->globalXStart ||
         this->localToGlobalX(this->ownedXEnd - 1) !=
             this->globalXStart + this->localNx - 1))
    {
        throw std::runtime_error("LATTICE:: Invalid local-to-global x mapping.");
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
    // gauge-covariant implementation exists. Backend Interior MPI boundary should never be specified.
    bool has_fixed_boundary = false;
    for (unsigned face_iter = 0; face_iter < 6; face_iter++)
    {
        if (face_boundary_condition_types[face_iter] == INTERIOR)
        {
            throw std::runtime_error(
                "LATTICE:: INTERIOR is an implementation-only boundary type "
                "and cannot be assigned in Lattice.cfg."
            );
        }

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

        if (face_boundary_condition_types[face_iter] == FIXED)
            has_fixed_boundary = true;
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

    if (has_fixed_boundary && this->rank == 0)
    {
        std::cout
            << "BOUNDARYCONDITIONS::FIXED::Warning: "
            << "This will currently only work when the stencil size is 1.\n"
            << "Larger stencils are not implemented yet.\n"
            << std::endl;
    }

    // Convert the six global physical face assignments into the six faces seen
    // by this rank. A process neighbour replaces the corresponding local x-face
    // with an internal interface. 
    std::vector<int> local_face_boundary_condition_types = face_boundary_condition_types;

    if (this->numRanks > 1)
    {
        if (this->leftRank >= 0)
            local_face_boundary_condition_types[0] = INTERIOR;

        if (this->rightRank >= 0)
            local_face_boundary_condition_types[1] = INTERIOR;
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


        // Work out which rank-local faces meet at this region. Physical fixed
        // support takes priority over an MPI interface, which in turn takes
        // priority over periodic treatment:
        //
        //     FIXED > INTERIOR > PERIODIC.
        int boundary_condition_type = PERIODIC;

        for (unsigned axis = 0; axis < 3; axis++)
        {
            if (bound_vector[axis] == 0)
                continue;

            const unsigned face_index =
                2*axis + (bound_vector[axis] > 0 ? 1 : 0);

            const int face_type =
                local_face_boundary_condition_types[face_index];

            if (face_type == FIXED)
            {
                boundary_condition_type = FIXED;
                break;
            }

            if (face_type == INTERIOR)
                boundary_condition_type = INTERIOR;
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
                                                          this->model, this->storageNx, this->ny, this->nz, this->activeEvolutionStep,
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
                                                             this->model, this->storageNx, this->ny, this->nz, this->activeEvolutionStep,
                                                             bound_vector)
                                              );
            break;

        case INTERIOR:

            this->boundaryConditions.push_back(
                new Interior(
                    this->scalarFields,
                    this->vectorFields,
                    this->analysers,
                    this->numScalarFieldComponents,
                    this->numVectorFieldComponents,
                    this->model,
                    this->storageNx,
                    this->ny,
                    this->nz,
                    this->activeEvolutionStep,
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
    long long unsigned scalar_array_size = 2ULL*this->storageNx*this->ny*this->nz*this->numScalarFieldComponents;
    this->scalarFields.resize(scalar_array_size, 0.f);

    // Same as above (numVectorFieldComponents contains 3*the amount requested to account for spatial components of each)
    long long unsigned vector_array_size = 2ULL*this->storageNx*this->ny*this->nz*this->numVectorFieldComponents;
    this->vectorFields.resize(vector_array_size, 0.f);

    // A zero gauge field is stored as the SU(2) identity in quaternion form.
    if (this->model.isUsingQuaternionRepresentation())
    {
        const unsigned direction_width
            = this->numVectorFieldComponents/3U;

        for (long long unsigned loc_iter = 0;
            loc_iter < 2ULL*this->storageNx*this->ny*this->nz;
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

    // Initial conditions are evaluated only on physical sites owned by this
    // rank. Halo cells retain their representation-safe defaults until the
    // initial halo exchange.
    const InitialConditionGeometry ic_geometry{
        this->nx,
        this->ny,
        this->nz,
        this->localNx,
        this->storageNx,
        this->globalXStart,
        this->ownedXBegin,
        this->ownedXEnd,
        this->dx,
        this->dy,
        this->dz,
        this->dt
    };

    if (scalarInitialConditions)
        scalarInitialConditions->setInitialFields(
            scalarFields,
            ic_geometry,
            this->numScalarFieldComponents);

    if (vectorInitialConditions)
        vectorInitialConditions->setInitialFields(
            vectorFields,
            ic_geometry,
            this->numVectorFieldComponents);


    if (this->progressReport)
    {
        std::cout << "Initial conditions generated.\n" << std::endl;
    }
}

void Lattice::synchroniseTimeBuffers(const unsigned source_time_index)
{
    if (source_time_index > 1U)
    {
        throw std::runtime_error(
            "LATTICE:: Invalid source time-buffer index."
        );
    }

    const std::size_t storage_volume
        = static_cast<std::size_t>(this->storageNx)
        *this->ny*this->nz;

    const unsigned destination_time_index
        = 1U - source_time_index;

    const auto copy_buffer =
        [source_time_index, destination_time_index](
            std::vector<float> &fields,
            const std::size_t buffer_size)
    {
        if (buffer_size == 0U)
            return;

        std::copy_n(
            fields.begin()
                + source_time_index*buffer_size,
            buffer_size,
            fields.begin()
                + destination_time_index*buffer_size
        );
    };

    copy_buffer(
        this->scalarFields,
        storage_volume*this->numScalarFieldComponents
    );

    copy_buffer(
        this->vectorFields,
        storage_volume*this->numVectorFieldComponents
    );
}

////////////////////////////////////////////////  Private Functions  /////////////////////////////////////////////////////

void Lattice::beginHaloExchange(
    const unsigned &t_step,
    HaloExchange &exchange)
{
    if (t_step > 1U)
        throw std::runtime_error(
            "LATTICE:: Invalid time-buffer index for halo exchange."
        );

    if (exchange.requestCount != 0)
        throw std::runtime_error(
            "LATTICE:: Halo exchange is already in progress."
        );

    if (this->numRanks == 1 || this->haloDepth == 0U)
        return;

#ifndef GFT_ENABLE_MPI

    throw std::runtime_error(
        "LATTICE:: Multi-rank halo exchange requires an MPI build."
    );

#else

    const int left_peer
        = (this->leftRank >= 0)
        ? this->leftRank
        : MPI_PROC_NULL;

    const int right_peer
        = (this->rightRank >= 0)
        ? this->rightRank
        : MPI_PROC_NULL;

    auto post_field =
        [this, t_step, left_peer, right_peer, &exchange](
            std::vector<float> &field,
            const unsigned num_components,
            const int left_tag,
            const int right_tag)
    {
        if (num_components == 0U)
            return;

        const unsigned long long plane_width
            = 1ULL*this->ny*this->nz*num_components;

        const unsigned long long exchange_size
            = 1ULL*this->haloDepth*plane_width;

        if (exchange_size >
            static_cast<unsigned long long>(
                std::numeric_limits<int>::max()))
        {
            throw std::runtime_error(
                "LATTICE:: MPI halo message exceeds "
                "the supported MPI count range."
            );
        }

        const int exchange_count
            = static_cast<int>(exchange_size);

        const unsigned long long time_offset
            = 1ULL*t_step*this->storageNx*plane_width;

        const unsigned long long left_halo_offset
            = time_offset
            + 1ULL*(this->ownedXBegin - this->haloDepth)
              *plane_width;

        const unsigned long long left_owned_offset
            = time_offset
            + 1ULL*this->ownedXBegin*plane_width;

        const unsigned long long right_owned_offset
            = time_offset
            + 1ULL*(this->ownedXEnd - this->haloDepth)
              *plane_width;

        const unsigned long long right_halo_offset
            = time_offset
            + 1ULL*this->ownedXEnd*plane_width;


        // Post the two receives first.
        int mpi_error = MPI_Irecv(
            field.data() + right_halo_offset,
            exchange_count,
            MPI_FLOAT,
            right_peer,
            left_tag,
            MPI_COMM_WORLD,
            &exchange.requests[exchange.requestCount]);

        if (mpi_error != MPI_SUCCESS)
            throw std::runtime_error(
                "LATTICE:: Failed to post MPI halo receive."
            );

        exchange.requestCount++;


        mpi_error = MPI_Irecv(
            field.data() + left_halo_offset,
            exchange_count,
            MPI_FLOAT,
            left_peer,
            right_tag,
            MPI_COMM_WORLD,
            &exchange.requests[exchange.requestCount]);

        if (mpi_error != MPI_SUCCESS)
            throw std::runtime_error(
                "LATTICE:: Failed to post MPI halo receive."
            );

        exchange.requestCount++;


        // Now post the two sends.
        mpi_error = MPI_Isend(
            field.data() + left_owned_offset,
            exchange_count,
            MPI_FLOAT,
            left_peer,
            left_tag,
            MPI_COMM_WORLD,
            &exchange.requests[exchange.requestCount]);

        if (mpi_error != MPI_SUCCESS)
            throw std::runtime_error(
                "LATTICE:: Failed to post MPI halo send."
            );

        exchange.requestCount++;


        mpi_error = MPI_Isend(
            field.data() + right_owned_offset,
            exchange_count,
            MPI_FLOAT,
            right_peer,
            right_tag,
            MPI_COMM_WORLD,
            &exchange.requests[exchange.requestCount]);

        if (mpi_error != MPI_SUCCESS)
            throw std::runtime_error(
                "LATTICE:: Failed to post MPI halo send."
            );

        exchange.requestCount++;
    };


    post_field(
        this->scalarFields,
        this->numScalarFieldComponents,
        100,
        101);

    post_field(
        this->vectorFields,
        this->numVectorFieldComponents,
        102,
        103);

#endif
}

void Lattice::finishHaloExchange(HaloExchange &exchange)
{
    if (exchange.requestCount == 0)
        return;

#ifndef GFT_ENABLE_MPI

    throw std::runtime_error(
        "LATTICE:: Multi-rank halo exchange requires an MPI build."
    );

#else

    const int mpi_error = MPI_Waitall(
        exchange.requestCount,
        exchange.requests.data(),
        MPI_STATUSES_IGNORE);

    if (mpi_error != MPI_SUCCESS)
        throw std::runtime_error(
            "LATTICE:: MPI halo exchange failed while waiting for completion."
        );

    exchange.requestCount = 0;

#endif
}

void Lattice::exchangeHalos(const unsigned &t_step)
{
    HaloExchange exchange;

    this->beginHaloExchange(t_step, exchange);
    this->finishHaloExchange(exchange);
}

unsigned Lattice::getDefaultStencilSize() const
{
    return this->model.getDefaultStencilSize();

}

unsigned Lattice::localToGlobalX(const unsigned &local_x) const
{
    // Deliberatley exclude halo positions
    if (local_x < this->ownedXBegin || local_x >= this->ownedXEnd)
        throw std::runtime_error(
            "LATTICE:: Local-to-global x mapping requested for a non site."
        );

    return this->globalXStart + (local_x - this->ownedXBegin);
}

std::vector< std::vector<unsigned> > Lattice::determineResponsibilities(const unsigned &stencil_size) const
{
    std::vector< std::vector<unsigned> > loop_limits(3, std::vector<unsigned>(2, stencil_size));

    // The x loop operates only on physical sites owned by this rank.
    // Halo planes are support data and must never acquire evolution
    // responsibility.
    loop_limits[0][0] = this->ownedXBegin + stencil_size;

    if (this->localNx < stencil_size)
        loop_limits[0][1] = this->ownedXBegin;
    else
        loop_limits[0][1] = this->ownedXEnd - stencil_size;

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
            t_now*this->storageNx*this->ny*this->nz
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


float *Lattice::getScalarFieldPointer(
    const unsigned &t_step,
    const unsigned &x_loc,
    const unsigned &y_loc,
    const unsigned &z_loc)
{
    long long unsigned array_location =
        (((t_step*this->storageNx + x_loc)*this->ny + y_loc)
         *this->nz + z_loc)*this->numScalarFieldComponents;

    return &this->scalarFields[array_location];
}


float *Lattice::getVectorFieldPointer(
    const unsigned &t_step,
    const unsigned &x_loc,
    const unsigned &y_loc,
    const unsigned &z_loc)
{
    long long unsigned array_location =
        (((t_step*this->storageNx + x_loc)*this->ny + y_loc)
         *this->nz + z_loc)*this->numVectorFieldComponents;

    return &this->vectorFields[array_location];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                    Public                                                            //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////  Constructors/Destructors  //////////////////////////////////////////////////

Lattice::Lattice(const int rank, const int num_ranks)
{
    this->initVariables();

    if (num_ranks < 1 || rank < 0 || rank >= num_ranks)
        throw std::runtime_error("LATTICE:: Invalid MPI rank information.");

    this->rank = rank;
    this->numRanks = num_ranks;

    this->configure(std::string(SOURCE_DIR) + "/Config/Lattice.cfg", true);

    this->model.configure(
        std::string(SOURCE_DIR) + "/Config/Model.cfg",
        this->numScalarFieldComponents, this->numVectorFieldComponents,
        this->nx, this->ny, this->nz,
        this->dx, this->dy, this->dz,
        this->dt, true
    );

    // Determine the complete operator footprint before any later
    // geometry-dependent objects are instantiated.
    this->stencilSize = this->getDefaultStencilSize();

    // Determine this rank's physical ownership and the storage geometry
    this->initSlabGeometry();

    // Instantiate the choices stored while reading Lattice.cfg only after
    // the model and current lattice geometry are known.
    this->initInitialCondition(this->initialConditionTypes);
    this->initAnalysers(this->analysisChoices);
    this->initBoundaryCondition(this->faceBoundaryConditionTypes);

    this->initFields();
    if (this->performInitialGradientFlow
    && this->gradientFlowTimesteps > 0U)
    {
        // Buffer 1 is the initial current configuration.
        // Gradient flow starts with no pre-existing flow velocity.
        this->synchroniseTimeBuffers(1U);
    }
    // Both stored initial time levels must have valid neighbour support before
    // any stencil-based operation is allowed to inspect them.
    this->exchangeHalos(0U);
    this->exchangeHalos(1U);

    // Determine the responsibilities of the boundary conditions.
    // Also check for boundaries with no responsibilities (common in 1/2D sims for example) and remove them.
    for (auto it = this->boundaryConditions.begin(); it != this->boundaryConditions.end(); )
    {
        bool is_empty = (*it)->determineResponsibilities(this->stencilSize, this->ownedXBegin, this->ownedXEnd);
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
    const unsigned maximum_gradient_flow_steps
        = this->performInitialGradientFlow
        ? this->gradientFlowTimesteps
        : 0U;

    unsigned num_gradient_flow_steps
        = maximum_gradient_flow_steps;

    const unsigned num_dynamical_steps
        = this->performDynamicalEvolution
        ? this->nt
        : 0U;

    if (num_gradient_flow_steps > std::numeric_limits<unsigned>::max() - num_dynamical_steps)
    {
        throw std::runtime_error(
            "LATTICE:: Total number of evolution steps overflows unsigned."
        );
    }

    unsigned total_evolution_steps = num_gradient_flow_steps + num_dynamical_steps;
    if (total_evolution_steps == 0U)
    {
        if (this->rank == 0)
        {
            std::cout
                << "LATTICE::WARNING: Neither gradient flow nor "
                << "dynamical evolution has been requested.\n"
                << std::endl;
        }

        return;
    }



    // Get start time of evolution stage if reporting progress.
    std::chrono::steady_clock::time_point start_time;
    if (this->progressReport)
    {
        #ifdef GFT_ENABLE_MPI
        if (this->numRanks > 1 &&
            MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
        {
            throw std::runtime_error(
                "LATTICE:: MPI barrier failed before evolution timing."
            );
        }
        #endif

        if (this->rank == 0)
            start_time = std::chrono::steady_clock::now();
    }


    // Determine how much of the grid is interior and can be evolved in the default manner
    // and how many the boundary condition classes will be responsible for evolving.
    std::vector< std::vector<unsigned> > loop_limits = this->determineResponsibilities(this->stencilSize);

    if (this->progressReport && this->rank == 0)   
        std::cout << "Beginning evolution." << std::endl;

    for (unsigned time_iter = 0; time_iter < total_evolution_steps; time_iter++)
    {
        // Gradient flow?
        const bool using_gradient_flow = time_iter < num_gradient_flow_steps;
        const unsigned phase_time_iter = using_gradient_flow ? time_iter : time_iter - num_gradient_flow_steps;
        this->activeEvolutionStep = using_gradient_flow ? this->gradientFlowStepSize : this->dt;
        this->model.setEvolutionMode( using_gradient_flow, this->activeEvolutionStep);


        // Track which parts of the array are the current timestep and which are the previous.
        unsigned t_now = (time_iter+1)%2;
        unsigned t_past = !t_now;

        // Update model internal evolution parameters that depend upon the timestep
        model.update(phase_time_iter);

        if (using_gradient_flow)
            this->model.resetGradientFlowResiduals();

        // Evolve grid points that are near/on boundaries. Once this is
        // complete, the new outgoing x-interface strips in t_past are ready.
        for (auto bound : this->boundaryConditions)
            bound->evolve(t_now, this->stencilSize);

        // Start communicating those new interface strips now. The central
        // volume below is independent of the incoming t_past halos, so its
        // evolution can overlap the MPI transfer.
        HaloExchange halo_exchange;
        this->beginHaloExchange(t_past, halo_exchange);

        // Save computation by calculating the indices at lowest loop levels possible.
        std::vector< std::vector<long long unsigned> > t_running_indices(
            3, std::vector<long long unsigned>(
                2*this->stencilSize + 1, t_now*this->storageNx*this->ny*this->nz
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
                    this->model.evolve(local_scalar_pointers, local_vector_pointers, this->activeEvolutionStep,
                                       this->numScalarFieldComponents, this->numVectorFieldComponents);
                    
                }
            }
        }

        // The full owned t_past state is now evolved. The new neighbour
        // halos must be complete before post-evolution analysis or the next
        // timestep is allowed to read them.
        this->finishHaloExchange(halo_exchange);

        //Optional Grad flow residual check
        const unsigned completed_phase_steps = phase_time_iter + 1U;

        double global_scalar_residual = 0.0;
        double global_gauge_residual = 0.0;
        bool residual_tolerances_satisfied = false;

        if (using_gradient_flow
            && this->stopGradientFlowOnResiduals)
        {
            std::array<double, 2> local_residuals{
                this->model.getMaxScalarEquationResidual(),
                this->model.getMaxGaugeEquationResidual()
            };

            std::array<double, 2> global_residuals
                = local_residuals;

            #ifdef GFT_ENABLE_MPI
            if (this->numRanks > 1)
            {
                const int reduce_error = MPI_Allreduce(
                    local_residuals.data(),
                    global_residuals.data(),
                    static_cast<int>(global_residuals.size()),
                    MPI_DOUBLE,
                    MPI_MAX,
                    MPI_COMM_WORLD
                );

                if (reduce_error != MPI_SUCCESS)
                {
                    throw std::runtime_error(
                        "LATTICE:: MPI gradient-flow "
                        "residual reduction failed."
                    );
                }
            }
            #endif

            global_scalar_residual = global_residuals[0];
            global_gauge_residual = global_residuals[1];

            residual_tolerances_satisfied
                = global_scalar_residual
                    <= this->scalarEquationResidualTolerance
                && global_gauge_residual
                    <= this->gaugeEquationResidualTolerance;

            if (residual_tolerances_satisfied)
            {
                num_gradient_flow_steps
                    = completed_phase_steps;

                total_evolution_steps
                    = num_gradient_flow_steps
                    + num_dynamical_steps;
            }
        }

        // Run post-evolution location analyses only after the new
        // timestep is complete across the full dynamic grid.
        this->postEvolveAnalysis(t_now, loop_limits);

        // Synchronises timesteps at the end of a gradient flow period.
        if (using_gradient_flow
            && completed_phase_steps
                == num_gradient_flow_steps)
        {
            // t_past contains the newest relaxed configuration.
            this->synchroniseTimeBuffers(t_past);
        }

        // Run all continous analyser functions that don't need to happen at every location in the dynamic grid.
        for (auto analyser : this->analysers)
            analyser->timestepAnalysis(time_iter);

        // Report progress using the counter for the current phase.
        const unsigned phase_total_steps
            = using_gradient_flow
            ? num_gradient_flow_steps
            : num_dynamical_steps;

        if (this->progressReport
            && this->rank == 0
            && ((completed_phase_steps % this->reportFrequency == 0U)
                || completed_phase_steps == phase_total_steps))
        {
            std::cout
                << (using_gradient_flow
                    ? "GradFlow timestep: "
                    : "Dynamical timestep: ")
                << completed_phase_steps
                << " completed."
                << std::endl;

            if (using_gradient_flow
                && this->stopGradientFlowOnResiduals)
            {
                std::cout
                    << "Scalar EOM residual: "
                    << global_scalar_residual
                    << ", gauge EOM residual: "
                    << global_gauge_residual
                    << "."
                    << std::endl;
            }

            if (using_gradient_flow
                && this->stopGradientFlowOnResiduals
                && this->rank == 0)
            {
                if (residual_tolerances_satisfied)
                {
                    std::cout
                        << "Gradient flow converged after "
                        << completed_phase_steps
                        << " timesteps."
                        << std::endl;
                }
                else if (completed_phase_steps
                    == maximum_gradient_flow_steps)
                {
                    std::cout
                        << "LATTICE::WARNING: Gradient flow reached "
                        << "its maximum number of timesteps before "
                        << "both residual tolerances were satisfied."
                        << std::endl;
                }
            }
        }

    }

    if (this->progressReport)
    {
        #ifdef GFT_ENABLE_MPI
        if (this->numRanks > 1 &&
            MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
        {
            throw std::runtime_error(
                "LATTICE:: MPI barrier failed after evolution."
            );
        }
        #endif

        if (this->rank == 0)
        {
            std::chrono::steady_clock::time_point end_time
                = std::chrono::steady_clock::now();

            std::chrono::duration<double> elapsed_time
                = end_time - start_time;

            std::cout
                << "Evolution finished in "
                << elapsed_time.count()
                << " seconds."
                << std::endl;
        }
    }

}

void Lattice::finalAnalysis() const
{
    for (auto analyser : this->analysers)
        analyser->finalAnalysis();
}
