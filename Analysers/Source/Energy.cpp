#include "Energy.hpp"

////////////////////////////////////////  Initialisers  /////////////////////////////////////////////////////

void Energy::configure(const std::string path, const bool debug)
{
    std::ifstream ifs(path);

    // Read in the parameter values from the configuration file
    std::string description; // Descriptions in the config file (to be dumped)
    if (ifs.is_open())
    {
        // Skip the initial configuration file description
        for(int iter = 0; iter < 13; iter++) std::getline(ifs, description);

        // Load in the output path for global quantities output
        std::getline(ifs, description, ':');
        ifs >> this->globalQuantitiesPath;

        // Load in global output options
        this->anyGlobalOptions = false;
        for (unsigned iter = 0; iter < 6; iter++)
        {
            std::getline(ifs, description, ':');
            ifs >> this->globalOptions[iter];
            if (this->globalOptions[iter])
                this->anyGlobalOptions = true;
        }
        std::getline(ifs, description, ':');
        ifs >> this->globalFrequency;

        // Load in the output path for local quantities output
        std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> this->localQuantitiesPath;

        // Load in local output options
        this->anyLocalOptions = false;
        for (unsigned iter = 0; iter < 6; iter++)
        {
            std::getline(ifs, description, ':');
            ifs >> this->localOptions[iter];
            if (this->localOptions[iter])
                this->anyLocalOptions = true;
        }
        std::getline(ifs, description, ':');
        ifs >> this->localFrequency;
    } 

    ifs.close();

    // Start a fresh file
    if (this->anyGlobalOptions)
    {
        std::ofstream ofs(std::string(DATA_DIR) + "/" + this->globalQuantitiesPath);
        ofs.close();
    }

    if (debug)
    {
        std::cout << "ANALYSERS::ENERGY::\n"
                  << "Global quantities data path: Data/" << this->globalQuantitiesPath << "\n"
                  << "Global options:: Energy: " << this->globalOptions[0] << ", Potential: " << this->globalOptions[1]
                  << ", Gradient: " << this->globalOptions[2] << ", Kinetic: " << this->globalOptions[3] 
                  << ", Magnetic: " << this->globalOptions[4] << ", Electric: " << this->globalOptions[5] << ", every " << this->globalFrequency << " timesteps\n"
                  << "Local quantities data path: Data/" << this->localQuantitiesPath << "\n"
                  << "Local options:: Energy: " << this->localOptions[0] << ", Potential: " << this->localOptions[1]
                  << ", Gradient: " << this->localOptions[2] << ", Kinetic: " << this->localOptions[3] 
                  << ", Magnetic: << " << this->localOptions[4] << ", Electric: " << this->localOptions[5] << ", every " << this->localFrequency << " timesteps\n"
                  << std::endl;
    }
}

void Energy::initVariables(const long long unsigned grid_size)
{
    this->energy = 0.f;
    this->potential = 0.f;
    this->gradient = 0.f;
    this->kinetic = 0.f;
    this->magnetic = 0.f;
    this->electric = 0.f;

    for (unsigned iter = 0; iter < 6; iter++)
    {
        if (this->localOptions[iter])
            this->density_pointers[iter]->resize(grid_size, 0.f);
    }

    this->counter = 0;
    this->globalOutput = this->anyGlobalOptions;
    this->localOutput = this->anyLocalOptions;
}

///////////////////////////////////  Constructors/Destructors  //////////////////////////////////////////////

Energy::Energy(const Model &model, const double &dx, const double &dy, const double &dz, const long long unsigned grid_size)
    : model(model), dx(dx), dy(dy), dz(dz),
      energy_pointers{&this->energy, &this->potential, &this->gradient, &this->kinetic, &this->magnetic, &this->electric},
      density_pointers{&this->energy_density, &this->potential_density, &this->gradient_density, &this->kinetic_density,
                       &this->magnetic_density, &this->electric_density}
{
    this->configure(std::string(SOURCE_DIR) + "/Config/Energy.cfg", true);
    this->initVariables(grid_size);
}

Energy::~Energy()
{
}

///////////////////////////////////////  Public Functions  //////////////////////////////////////////////////

void Energy::initialAnalysis()
{
    this->model.energyPreparation(this->globalOutput || this->localOutput);
}

void Energy::locationAnalysis(const long long unsigned index,
                              const float* const local_scalar_pointers[2], const std::vector<std::vector<const float*>> &scalar_pointers,
                              const float* const local_vector_pointers[2], const std::vector<std::vector<const float*>> &vector_pointers)
{
    if (this->localOutput || this->globalOutput)
    {



        float potential_density = this->model.calcPotentialEnergy(local_scalar_pointers[1]);
        float gradient_density = this->model.calcGradientEnergy(scalar_pointers, vector_pointers);
        float kinetic_density = this->model.calcKineticEnergy(local_scalar_pointers);
        float magnetic_density = this->model.calcMagneticEnergy(vector_pointers);
        float electric_density = this->model.calcElectricEnergy(local_vector_pointers);

        float energy_density = potential_density + gradient_density + kinetic_density + magnetic_density + electric_density;

        // So I can loop over output choices and/or contributions to the integrated quantities.
        const float* const density_pointers[6] = {&energy_density, &potential_density, &gradient_density, &kinetic_density,
                                                  &magnetic_density, &electric_density};

        if (this->localOutput)
        {
            for (unsigned iter = 0; iter < 6; iter++)
            {
                if (this->localOptions[iter])
                    (*this->density_pointers[iter])[index] = *density_pointers[iter];
            }
        }

        // Have any global output options been chosen?
        if (this->globalOutput)
        {
            // So I can loop over the integrated quantities and the densities
            for (unsigned iter = 0; iter < 6; iter++)
            {
                if (this->globalOptions[iter])
                {
                    // Only bother doing the multiplication by volume at the timestepAnalysis stage do avoid unneccessary computation.
                    *this->energy_pointers[iter] += *density_pointers[iter];
                }
            }
            
        }

    }
}

void Energy::timestepAnalysis(const unsigned &time_step)
{
    if (this->globalOutput)
    {
        std::ofstream ofs(std::string(DATA_DIR) + "/" + this->globalQuantitiesPath, std::ios::app);

        if (ofs.is_open())
        {
            for (unsigned iter = 0; iter < 6; iter++)
            {
                if (this->globalOptions[iter])
                {
                    ofs << *this->energy_pointers[iter]*this->dx*this->dy*this->dz << " ";
                    *this->energy_pointers[iter] = 0.f;
                }
            }
            ofs << "\n";
        }

        ofs.close();
    }

    if (this->localOutput)
    {
        std::ofstream ofs(std::string(DATA_DIR) + "/" + this->localQuantitiesPath + "_" + std::to_string(counter) + ".dat");

        if (ofs.is_open())
        {
            for (unsigned iter = 0; iter < 6; iter++)
            {
                if (this->localOptions[iter])
                {
                    for (float &data : *this->density_pointers[iter])
                    {
                        ofs << data << "\n";
                    }
                }
            }
        }

        ofs.close();
    }

    // Advance the counter no matter what
    this->counter++;

    // Check if conditions for global and local output are satisfied
    this->globalOutput = this->anyGlobalOptions && this->counter%this->globalFrequency == 0;
    this->localOutput = this->anyLocalOptions && this->counter%this->localFrequency == 0;

    this->model.energyPreparation(this->globalOutput || this->localOutput);
}

void Energy::finalAnalysis()
{
}
