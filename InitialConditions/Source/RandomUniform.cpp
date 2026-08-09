#include "RandomUniform.hpp"
#include <cstddef>

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

void RandomUniform::setInitialFields(std::vector<float> &field, const InitialConditionGeometry &geometry, const unsigned num_components) const
{
    if (num_components == 0U)
        return;

    const std::size_t plane_width
        = static_cast<std::size_t>(geometry.globalNy)
        * geometry.globalNz
        * num_components;

    const std::size_t local_buffer_size
        = static_cast<std::size_t>(geometry.storageNx)
        * plane_width;

    const std::size_t global_buffer_size
        = static_cast<std::size_t>(geometry.globalNx)
        * plane_width;

    const std::size_t owned_count
        = static_cast<std::size_t>(geometry.localNx)
        * plane_width;

    const std::size_t local_owned_offset
        = static_cast<std::size_t>(geometry.ownedXBegin)
        * plane_width;

    const std::size_t global_owned_offset
        = static_cast<std::size_t>(geometry.globalXStart)
        * plane_width;

    if (field.size() != 2ULL*local_buffer_size)
        throw std::runtime_error(
            "INITIALCONDITIONS::RANDOMUNIFORM:: "
            "Field size does not match IC geometry."
        );

    std::mt19937 generator(this->seed);
    std::uniform_real_distribution<float> distribution(
        this->minRange,
        this->maxRange);

    // Advance through the actual distribution.
    for (std::size_t iter = 0;
         iter < global_owned_offset;
         iter++)
    {
        (void)distribution(generator);
    }

    // First stored time buffer: generate this rank's owned block only.
    for (std::size_t iter = 0;
         iter < owned_count;
         iter++)
    {
        field[local_owned_offset + iter]
            = distribution(generator);
    }

    const std::size_t second_local_owned_offset
        = local_buffer_size + local_owned_offset;

    if (this->setTimestepsEqual)
    {
        for (std::size_t iter = 0;
             iter < owned_count;
             iter++)
        {
            field[second_local_owned_offset + iter]
                = field[local_owned_offset + iter];
        }

        return;
    }

    // Continue through the original global sequence until this rank's
    // owned block in the second time buffer is reached.
    const std::size_t first_owned_end
        = global_owned_offset + owned_count;

    const std::size_t second_global_owned_offset
        = global_buffer_size + global_owned_offset;

    for (std::size_t iter = first_owned_end;
         iter < second_global_owned_offset;
         iter++)
    {
        (void)distribution(generator);
    }

    for (std::size_t iter = 0;
         iter < owned_count;
         iter++)
    {
        field[second_local_owned_offset + iter]
            = distribution(generator);
    }
}
