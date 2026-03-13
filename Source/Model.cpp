#include "Model.hpp"

////////////////////////////////////  Initialisers  ////////////////////////////////////////////////

void Model::initVariables()
{
    this->potential = nullptr;
    this->gradient = nullptr;
    this->wilsonLoop = nullptr;
    this->dampingFactor = 0.f;
    this->ntDamped = 0;
    this->currentDamping = 0.f;
}

void Model::initPotential(const int &potential_type, const unsigned &num_scalar_components)
{
    switch (potential_type)
    {
    case UNASSIGNED_POTENTIAL:

        std::cout << "MODEL::ERROR: A potential has not been assigned to the model.\n" << std::endl;
        break;

    case SO_N_POTENTIAL:

        this->potential = new SO_N(num_scalar_components);
        break;

    case TWO_HDM_POTENTIAL:

        std::cout << "MODEL::ERROR: The potential type has been chosen to be the 2HDM, but this feature has not been added yet.\n" << std::endl;
        break;

    default:

        std::cout << "MODEL::ERROR: The chosen potential type (" << potential_type << ") is invalid.\n" << std::endl;
        break;

    }
}

void Model::initGradient(
    const int &gradient_type,
    const unsigned &num_scalar_components, const unsigned &num_vector_components,
    const unsigned &nx, const unsigned &ny, const unsigned &nz,
    const double &dt, const double &dx, const double &dy, const double &dz
)
{
    switch (gradient_type)
    {
    case UNASSIGNED_GRADIENT:

        std::cout << "MODEL::ERROR: A gradient has not been assigned to the model.\n" << std::endl;
        break;

    case GLOBAL_GRADIENT:

        this->gradient = new Gradients::Global(num_scalar_components, num_vector_components,
                                    nx, ny, nz,
                                    dt, dx, dy, dz);
        break;

    case SM_GRADIENT:

        this->gradient = new Gradients::StandardModel(num_scalar_components, num_vector_components,
                                                   nz, ny, nz,
                                                   dt, dx, dy, dz);
        break;

    default:

        std::cout << "MODEL::ERROR: The chosen gradient type (" << gradient_type << ") is invalid.\n" << std::endl;
        break;

    }
}

void Model::initWilsonLoop(const int &wilson_loop_type, const double &dt, const double &dx, const double &dy, const double &dz)
{
    switch (wilson_loop_type)
    {
    case UNASSIGNED_WILSON_LOOP:

        // Print out a warning, but this is okay in some cases, e.g the global model.

        std::cout << "MODEL::WARNING: A wilson loop has not been assigned to the model.\n" << std::endl;
        break;

    case SM_WILSON_LOOP:

        this->wilsonLoop = new WilsonLoops::StandardModel(dt, dx, dy, dz);
        break;

    default:

        std::cout << "MODEL::ERROR: The chosen wilson loop type (" << wilson_loop_type << ") is invalid.\n" << std::endl;
        break;

    }
}

//////////////////////////////  Constructors/Destructors  //////////////////////////////////////////

Model::Model()
{
    this->initVariables();
}

Model::~Model()
{
    delete this->potential;
    delete this->gradient;
    delete this->wilsonLoop;

    this->potential = nullptr;
    this->gradient = nullptr;
    this->wilsonLoop = nullptr;
}

/////////////////////////////////////  Functions  //////////////////////////////////////////////////

void Model::configure(
    const std::string path, 
    const unsigned &num_scalar_components, const unsigned &num_vector_components,
    const unsigned &nx, const unsigned &ny, const unsigned &nz,
    const double &dx, const double &dy, const double &dz,const double &dt,
    const bool debug)
{

    // Resize and initialise the equation of motion vectors to the appropriate size
    this->potentialContributions.resize(num_scalar_components, 0.f);
    this->derivativeContributions.resize(num_scalar_components, 0.f);
    this->currentContributions.resize(num_vector_components, 0.f);
    this->wilsonLoopMagneticContributions.resize(num_vector_components, 0.f);
    this->wilsonLoopElectricContributions.resize(num_vector_components, 0.f);


    std::ifstream ifs(path);

    int potential_type = UNASSIGNED_POTENTIAL;
    int gradient_type = UNASSIGNED_GRADIENT;
    int wilson_loop_type = UNASSIGNED_WILSON_LOOP;


    // Read in the parameter values from the configuration file
    std::string description; // Descriptions in the config file (to be dumped)
    if (ifs.is_open())
    {
        // Skip the initial configuration file description
        for(int iter = 0; iter < 11; iter++) std::getline(ifs, description);

        // Load chosen potential type and set-up the chosen potential.
        std::getline(ifs, description, ':');
        ifs >> potential_type;
        this->initPotential(potential_type, num_scalar_components);

        // Load chosen gradient type and set-up the chosen gradient.
        std::getline(ifs, description,':');
        ifs >> gradient_type;
        this->initGradient(
            gradient_type,
            num_scalar_components, num_vector_components,
            nx, ny, nz,
            dt, dx, dy, dz
        );

        // Load chosen wilson loop type and set-up the chosen wilson loop.
        std::getline(ifs, description, ':');
        ifs >> wilson_loop_type;
        this->initWilsonLoop(wilson_loop_type, dt, dx, dy, dz);

        
        // Load parameters of the damping term
        std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> this->dampingFactor;

        std::getline(ifs, description, ':');
        ifs >> this->ntDamped;


        // Optionally overide ntDamped by setting it based on the physical time at which damping should end.
        bool nt_damped_override;
        std::getline(ifs, description, ':');
        ifs >> nt_damped_override;

        float t_end_damping;
        std::getline(ifs, description, ':');
        ifs >> t_end_damping;

        if (nt_damped_override)
        {
            this->ntDamped = static_cast<int>( t_end_damping/dt );
        }
    } 

    ifs.close();

    if (debug)
    {
        std::cout << "MODEL::\n"
                  << "Potential Type: " << potential_type << ", Gradient Type: " << gradient_type
                  << ", Wilson Loop Type: " << wilson_loop_type << "\n"
                  << "Damping Factor: " << this->dampingFactor << ", Number of damped timesteps: " << this->ntDamped
                  << "\n" << std::endl;
    }
}

void Model::energyPreparation() const
{
    if (wilsonLoop)
        this->wilsonLoop->energyPreparation();
}

void Model::update(unsigned time_iter)
{
    if (time_iter < this->ntDamped)
    {
        this->currentDamping = this->dampingFactor;
        this->evolveGauge = false;
    }
    else
    {
        this->currentDamping = 0.f;
        this->evolveGauge = true;
    }
}

unsigned Model::getDefaultStencilSize() const
{
    unsigned stencil_size = this->gradient->getDefaultStencilSize();
    // Same thing for wilson loops when implemented and determine which ones are the largest.
    return stencil_size;
}

void Model::calcPotentialContributions(const float* const local_scalar_fields)
{
    this->potentialContributions = this->potential->calcPotentialDerivatives(local_scalar_fields);
}

float Model::calcPotentialEnergy(const float *const local_scalar_field) const
{
    return this->potential->calcPotentialEnergy(local_scalar_field);
}

void Model::calcGradientContributions(const std::vector<std::vector<const float*>> &scalar_pointers, const std::vector<std::vector<const float*>> &vector_pointers)
{
    this->derivativeContributions = gradient->calcDerivatives(scalar_pointers, vector_pointers);
    this->currentContributions = gradient->calcCurrents(scalar_pointers, vector_pointers);
}

float Model::calcGradientEnergy(const std::vector<std::vector<const float *>> &scalar_pointers, const std::vector<std::vector<const float*>> &vector_pointers) const
{
    return this->gradient->calcGradientEnergy(scalar_pointers, vector_pointers);
}

float Model::calcKineticEnergy(const float* const local_scalar_fields[2]) const
{
    return this->gradient->calcKineticEnergy(local_scalar_fields);
}

void Model::calcYangMillsContributions(const std::vector<std::vector<const float *>> &vector_pointers, const float *const local_vector_fields[2])
{
    if (this->wilsonLoop)
    {
        this->wilsonLoopMagneticContributions = wilsonLoop->calcMagneticContributions(vector_pointers);
        this->wilsonLoopElectricContributions = wilsonLoop->calcElectricContributions(local_vector_fields);
    }
}

float Model::calcMagneticEnergy(const std::vector<std::vector<const float *>> &vector_pointers) const
{
    if (this->wilsonLoop)
        return this->wilsonLoop->calcMagneticEnergy(vector_pointers);
    else
        return 0.f;
}

float Model::calcElectricEnergy(const float *const local_vector_fields[2]) const
{
    if (this->wilsonLoop)
        return this->wilsonLoop->calcElectricEnergy(local_vector_fields);
    else
        return 0.f;
}

void Model::evolve(float *const local_scalar_fields[2], float *const local_vector_fields[2],
                   const double &dt, const unsigned &num_scalar_components, const unsigned &num_vector_components)
{
    // Evolve the scalar fields
    for (unsigned comp_iter = 0; comp_iter < num_scalar_components; comp_iter++)
    {
        double equation_right_hand_size = dt*dt*(this->derivativeContributions[comp_iter] - this->potentialContributions[comp_iter])
                                       - 0.5*this->currentDamping*dt*(static_cast<double>(local_scalar_fields[1][comp_iter]) - static_cast<double>(local_scalar_fields[0][comp_iter]));

        local_scalar_fields[0][comp_iter] =  static_cast<float>(
            2*static_cast<double>(local_scalar_fields[1][comp_iter]) - static_cast<double>(local_scalar_fields[0][comp_iter]) + equation_right_hand_size
        );
    }

    // Evolve the vector fields
    // This depends upon the wilson loop set-up, so evolution will be done within that class.
    if (this->evolveGauge)
    {
        // Approach using 3 dofs:
        //unsigned num_eq_components = num_vector_components;

        // Approach using 4 dofs:
        unsigned num_eq_components = num_vector_components - 3;
        std::vector<double> vector_equation_RHS (num_eq_components, 0.f);
        for (unsigned comp_iter = 0; comp_iter < num_eq_components; comp_iter++)
        {
            vector_equation_RHS[comp_iter] = this->wilsonLoopElectricContributions[comp_iter]
                                        + dt*dt*(this->wilsonLoopMagneticContributions[comp_iter] 
                                        + this->wilsonLoop->getSqrCouplings(comp_iter)*this->currentContributions[comp_iter]);

            //std::cout << this->wilsonLoop->getSqrCouplings(comp_iter) << " " << this->currentContributions[comp_iter] << std::endl;

            // std::cout << comp_iter << ": " << this->wilsonLoopTemporalContributions[comp_iter] << " " << this->wilsonLoopSpatialContributions[comp_iter] 
            //           << " " << this->currentContributions[comp_iter] << std::endl;
        }

        if (this->wilsonLoop)
            this->wilsonLoop->evolve(local_vector_fields, vector_equation_RHS);
    }

}
