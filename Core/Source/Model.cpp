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
    this->numVectorEqs = 0;
}

void Model::initPotential(const int &potential_type, const unsigned &num_scalar_components)
{
    switch (potential_type)
    {
    case NULL_POTENTIAL:

        this->potential = new NullPotential(num_scalar_components);
        std::cout << "MODEL::WARNING: Using null potential, all contributions will be zero.\n" << std::endl;
        break;

    case SO_N_POTENTIAL:

        this->potential = new SO_N(num_scalar_components);
        break;

    case DOUBLE_SO_N_POTENTIAL:

        this->potential = new Double_SO_N(num_scalar_components);
        break;

    default:

        throw std::runtime_error("MODEL:: The chosen potential type (" + std::to_string(potential_type) + ") is invalid.");
        break;

    }
}

void Model::initGradient(
    const int &gradient_type,
    const unsigned &num_scalar_components, const unsigned &num_vector_components,
    const unsigned &nx, const unsigned &ny, const unsigned &nz,
    const double &dt, const double &dx, const double &dy, const double &dz,
    const bool using_generator_representation
)
{
    switch (gradient_type)
    {
    case NULL_GRADIENT:

        this->gradient = new NullGradient(num_scalar_components, num_vector_components);
        std::cout << "MODEL::WARNING: Using null gradient, all contributions will be zero.\n" << std::endl;
        break;

    case GLOBAL_GRADIENT:

        this->gradient = new Gradients::Global(num_scalar_components, num_vector_components,
                                               nx, ny, nz,
                                               dt, dx, dy, dz);
        break;

    case SM_GRADIENT:

        this->gradient = new Gradients::StandardModel(num_scalar_components, num_vector_components,
                                                      nx, ny, nz,
                                                      dt, dx, dy, dz,
                                                      using_generator_representation);
        break;

    case TWOHDM_GRADIENT:

        this->gradient = new Gradients::TwoHDM(num_scalar_components, num_vector_components,
                                               nx, ny, nz,
                                               dt, dx, dy, dz,
                                               using_generator_representation);
        break;

    default:

        throw std::runtime_error("MODEL::ERROR: The chosen gradient type (" + std::to_string(gradient_type) + ") is invalid.");
        break;

    }
}

void Model::initWilsonLoop(const int &wilson_loop_type, const unsigned &num_vector_components,
                           const double &dt, const double &dx, const double &dy, const double &dz)
{
    switch (wilson_loop_type)
    {
    case NULL_WILSON_LOOP:

        this->wilsonLoop = new NullWilsonLoop(num_vector_components);
        std::cout << "MODEL::WARNING: Using null wilson loop, all contributions will be zero.\n" << std::endl;
        break;

    case SM_WILSON_LOOP:

        this->wilsonLoop = new WilsonLoops::StandardModel(num_vector_components, dt, dx, dy, dz);
        break;

    case TWOHDM_WILSON_LOOP:

        this->wilsonLoop = new WilsonLoops::TwoHDM(num_vector_components, dt, dx, dy, dz);
        break;

    default:

        throw std::runtime_error("MODEL: The chosen wilson loop type (" + std::to_string(wilson_loop_type) + ") is invalid.");
        break;

    }

    // Establish the number of vector equations.
    // In principle this can be different from the number of vector components if there is a convenient representation for them that has redundancies.
    this->numVectorEqs = this->wilsonLoop->getNumberOfEvolutionEquations();

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

    int potential_type = NULL_POTENTIAL;
    int gradient_type = NULL_GRADIENT;
    int wilson_loop_type = NULL_WILSON_LOOP;


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

        // Load chosen gradient type but delay the set-up until after the wilson loops have been set-up.
        std::getline(ifs, description);
        std::getline(ifs, description,':');
        ifs >> gradient_type;

        // Load chosen wilson loop type and set-up the chosen wilson loop.
        std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> wilson_loop_type;
        this->initWilsonLoop(wilson_loop_type, num_vector_components, dt, dx, dy, dz);

        // Set-up the gradient.
        this->initGradient(
            gradient_type,
            num_scalar_components, num_vector_components,
            nx, ny, nz,
            dt, dx, dy, dz,
            this->wilsonLoop->isUsingGeneratorRepresentation()
        );
        
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

void Model::energyPreparation(const bool store_energy) const
{
    this->wilsonLoop->energyPreparation(store_energy);
}

void Model::update(unsigned time_iter)
{
    if (time_iter < this->ntDamped)
    {
        this->currentDamping = this->dampingFactor;
        //this->evolveGauge = false;
        this->evolveGauge = true;
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

unsigned Model::getNumberOfConstraintEquations() const
{
    unsigned gradient_num_constraints = this->gradient->getNumberOfConstraintEquations();
    unsigned loop_num_constraints = this->wilsonLoop->getNumberOfConstraintEquations();

    if (gradient_num_constraints != loop_num_constraints)
        throw std::runtime_error("MODEL: The number of constraint equations expected by the gradients (" + std::to_string(gradient_num_constraints)
                                + ") does not match the number expected by the wilson loops (" + std::to_string(loop_num_constraints) + ").");

    return loop_num_constraints;
}

void Model::calcPotentialContributions(const float *const local_scalar_fields)
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
    this->wilsonLoopMagneticContributions = wilsonLoop->calcMagneticContributions(vector_pointers);
    this->wilsonLoopElectricContributions = wilsonLoop->calcElectricContributions(local_vector_fields);
}

float Model::calcMagneticEnergy(const std::vector<std::vector<const float *>> &vector_pointers) const
{
    return this->wilsonLoop->calcMagneticEnergy(vector_pointers);
}

float Model::calcElectricEnergy(const float *const local_vector_fields[2]) const
{
    return this->wilsonLoop->calcElectricEnergy(local_vector_fields);
}

std::vector<float> Model::calcConstraintViolation(const unsigned &num_equations, const long long int &t_future_index, 
                                                  const float *const local_scalar_fields[2], 
                                                  const std::vector<std::vector<const float *>> &vector_pointers) const
{
    std::vector<float> violation(num_equations, 0.f);
    std::vector<float> gradient_contribution = this->gradient->calcConstraintContributions(local_scalar_fields);
    std::vector<float> electric_contribution = this->wilsonLoop->calcConstraintContributions(t_future_index, vector_pointers);


    for (unsigned iter = 0; iter < num_equations; iter++)
        violation[iter] = electric_contribution[iter] - this->wilsonLoop->getSqrCouplings(iter)*gradient_contribution[iter];


    return violation;
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
        std::vector<double> vector_equation_RHS(this->numVectorEqs, 0.0);
        for (unsigned comp_iter = 0; comp_iter < this->numVectorEqs; comp_iter++)
        {
            vector_equation_RHS[comp_iter] = (1.0 - this->currentDamping*dt)*this->wilsonLoopElectricContributions[comp_iter]
                                        + dt*dt*(this->wilsonLoopMagneticContributions[comp_iter] 
                                        + this->wilsonLoop->getSqrCouplings(comp_iter)*this->currentContributions[comp_iter]);

            //std::cout << comp_iter << ": " << this->wilsonLoop->getSqrCouplings(comp_iter) << ", " << this->currentContributions[comp_iter] << ", " << this->wilsonLoopElectricContributions[comp_iter] << ", " << this->wilsonLoopMagneticContributions[comp_iter] << std::endl;
        }

        this->wilsonLoop->evolve(local_vector_fields, vector_equation_RHS);
    }

}
