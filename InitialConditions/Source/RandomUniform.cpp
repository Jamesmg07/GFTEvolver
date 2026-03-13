#include "RandomUniform.hpp"

////////////////////////////////////////////  Initialisers  //////////////////////////////////////////////////////

void RandomUniform::initVariables()
{
    this->seed = 0;
    this->minRange = 0;
    this->maxRange = 0;
    this->setTimestepsEqual = false;
}

void RandomUniform::configure(const std::string path, const bool debug)
{
    std::ifstream ifs(path);

    // Read in the parameter values from the configuration file
    std::string description; // Descriptions in the config file (to be dumped)
    if (ifs.is_open())
    {
        // Skip the initial configuration file description
        for(int iter = 0; iter < 12; iter++) std::getline(ifs, description);

        std::getline(ifs, description, ':');
        ifs >> this->seed;

        std::getline(ifs, description, ':');
        ifs >> this->minRange >> this->maxRange;

        std::getline(ifs, description, ':');
        ifs >> this->setTimestepsEqual;

    } 

    ifs.close();

    if (debug)
    {
        std::cout << "INITIALCONDITIONS::RANDOMUNIFORM::\n"
                  << "Seed: " << this->seed << ", Minimum: " << this->minRange << ", Maximum: " << this->maxRange << "\n"
                  << "Set timesteps to be equal?: " << this->setTimestepsEqual 
                  << "\n" << std::endl;
    }
}

////////////////////////////////////////  Constructors/Destructors  //////////////////////////////////////////////

RandomUniform::RandomUniform()
{
    this->initVariables();
    this->configure(std::string(SOURCE_DIR) + "/Config/RandomUniform.cfg", true);
}

RandomUniform::~RandomUniform()
{
}

////////////////////////////////////////////  Public Functions  //////////////////////////////////////////////////

void RandomUniform::setInitialFields(std::vector<float> &field) const
{
    std::mt19937 generator(this->seed);
    std::uniform_real_distribution<float> distribution(this->minRange, this->maxRange);

    size_t iterMax = field.size();
    if (this->setTimestepsEqual) 
        iterMax = iterMax/2;

    for (size_t iter = 0; iter < iterMax; iter++)
    {
        field[iter] = distribution(generator);

        if (this->setTimestepsEqual)
            field[iter + iterMax] = field[iter];
    }

}
