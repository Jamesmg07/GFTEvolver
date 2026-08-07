#include "SO_N.hpp"

//////////////////////////////////////////////  Initialisers  ///////////////////////////////////////////////////////

void SO_N::initVariables()
{
    this->lambda = 0;
    this->etaSqr = 0;
    this->fieldSqrMagnitude = 0;
}

void SO_N::configure(const std::string path, const bool debug)
{
    std::ifstream ifs(path);

    // Read in the parameter values from the configuration file
    std::string description; // Descriptions in the config file (to be dumped)
    if (ifs.is_open())
    {
        // Skip the initial configuration file description
        for(int iter = 0; iter < 15; iter++) std::getline(ifs, description);

        std::getline(ifs, description, ':');
        ifs >> this->lambda;

        std::getline(ifs, description, ':');
        ifs >> this->etaSqr;

    } 

    ifs.close();

    if (debug)
    {
        std::cout << "POTENTIALS::SO_N::\n"
                  << "#Field Components: " << this->numComponents << "\n"
                  << "lambda: " << this->lambda << ", eta^2: " << this->etaSqr 
                  << "\n" << std::endl;
    }
}

/////////////////////////////////////////  Constructors/Destructors  ///////////////////////////////////////////////

SO_N::SO_N(const unsigned &num_components)
    : numComponents(num_components)
{
    this->initVariables();
    this->configure(std::string(SOURCE_DIR) + "/Config/SO_N.cfg", true);
}

SO_N::~SO_N()
{
}

//////////////////////////////////////////////  Public functions  //////////////////////////////////////////////////

std::vector<double> SO_N::calcPotentialDerivatives(const float *field)
{
    std::vector<double> potential_contributions(this->numComponents, 0.0);

    // Calculate the |field|^2 and save it for possible later use in calculating the potential energy.
    this->fieldSqrMagnitude = 0.0;
    for (int iter = 0; iter < this->numComponents; iter++)
    {
        this->fieldSqrMagnitude += std::pow(static_cast<double>(field[iter]), 2);
    }

    // Calculate potential contributions to the equations of motion.
    for (int iter = 0; iter < this->numComponents; iter++)
    {
        potential_contributions[iter] = this->lambda*( this->fieldSqrMagnitude - this->etaSqr )*static_cast<double>(field[iter]);
    }

    return potential_contributions;
}

// Calculate the potential energy directly from the supplied field.
float SO_N::calcPotentialEnergy(const float *field) const
{
    double field_sqr_magnitude = 0.0;
    for (unsigned iter = 0; iter < this->numComponents; iter++)
    {
        field_sqr_magnitude += std::pow(static_cast<double>(field[iter]), 2);
    }

    return 0.25f*this->lambda*powf(field_sqr_magnitude - this->etaSqr, 2);
}
