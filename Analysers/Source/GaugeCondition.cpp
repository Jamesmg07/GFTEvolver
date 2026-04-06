#include "GaugeCondition.hpp"

////////////////////////////////////////  Initialisers  /////////////////////////////////////////////////////

void GaugeCondition::configure(const std::string path, const bool debug)
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
        for (unsigned iter = 0; iter < this->numGlobalOptions; iter++)
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
        for (unsigned iter = 0; iter < this->numLocalOptions; iter++)
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
        std::cout << "ANALYSERS::GAUGECONDITION::\n"
                  << "Global quantities data path: Data/" << this->globalQuantitiesPath << "\n"
                  << "Global options:: Integrated: " << this->globalOptions[0] << ", Max: " << this->globalOptions[1]
                  << ", every " << this->globalFrequency << " timesteps\n"
                  << "Local quantities data path: Data/" << this->localQuantitiesPath << "\n"
                  << "Local options:: Local Violation: " << this->localOptions[0] << ", every " << this->localFrequency << " timesteps\n"
                  << std::endl;
    }
}

void GaugeCondition::initVariables(const long long unsigned &grid_size)
{
    this->integratedAbsViolation.resize(this->numEquations, 0.f);
    this->maxAbsViolation.resize(this->numEquations, 0.f);

    for (unsigned iter = 0; iter < this->numLocalOptions; iter++)
    {
        if (this->localOptions[iter])
            this->localPointers[iter]->resize(grid_size, std::vector<float>(this->numEquations, 0.f));
    }

    this->counter = 0;
    this->globalOutput = this->anyGlobalOptions;
    this->localOutput = this->anyLocalOptions;
}

///////////////////////////////////  Constructors/Destructors  //////////////////////////////////////////////

GaugeCondition::GaugeCondition(const Model &model, const double &dx, const double &dy, const double &dz, 
                               const long long unsigned grid_size, const unsigned &num_vector_components)
    : model(model), dx(dx), dy(dy), dz(dz), gridSize(grid_size), numVectorComponents(num_vector_components),
      globalPointers{&this->integratedAbsViolation, &this->maxAbsViolation},
      localPointers{&this->localViolation}
{
    this->configure(std::string(SOURCE_DIR) + "/Config/GaugeCondition.cfg", true);
}

GaugeCondition::~GaugeCondition()
{
}

///////////////////////////////////////  Public Functions  //////////////////////////////////////////////////

void GaugeCondition::initialAnalysis()
{
    // Have to do these a bit later so that model has been set-up first.
    this->numEquations = this->model.getNumberOfConstraintEquations();
    this->initVariables(this->gridSize);
}

void GaugeCondition::preEvolveLocationAnalysis(const long long unsigned index,
                                       const float* const local_scalar_pointers[2], const std::vector<std::vector<const float*>> &scalar_pointers,
                                       const float* const local_vector_pointers[2], const std::vector<std::vector<const float*>> &vector_pointers)
{
}

void GaugeCondition::postEvolveLocationAnalysis(const unsigned &t_now, const long long unsigned index, const float *const local_scalar_pointers[2], const std::vector<std::vector<const float *>> &scalar_pointers, const float *const local_vector_pointers[2], const std::vector<std::vector<const float *>> &vector_pointers)
{
    if (this->localOutput || this->globalOutput)
    {

        long long int t_future_index = this->gridSize*this->numVectorComponents;
        if (t_now == 1)
            t_future_index = -t_future_index; // Need to subtract this index rather than add.

        std::vector<float> local_violation = this->model.calcConstraintViolation(this->numEquations, t_future_index, 
                                                                                 local_scalar_pointers, vector_pointers);

        // So I can loop over output choices and/or contributions to the integrated quantities.
        const std::vector<float>* const local_pointers[this->numLocalOptions] = {&local_violation};

        if (this->localOutput)
        {
            for (unsigned iter = 0; iter < this->numLocalOptions; iter++)
            {
                if (this->localOptions[iter])
                    (*this->localPointers[iter])[index] = *local_pointers[iter];
            }
        }

        // Have any global output options been chosen?
        if (this->globalOutput)
        {
             for (unsigned iter = 0; iter < this->numEquations; iter++)
             {
                // Pre-calculate the value since it might be used multiple times
                float value = std::abs((*local_pointers[0])[iter]);

                // First option is effectively an integration of the absolute value of each violation.
                // Multiplcation by the volume factor will be done later to avoid unneccessary computation.
                if (this->globalOptions[0])
                    this->integratedAbsViolation[iter] += value;

                // Second option is about finding the largest absolute violation in the grid.
                if (this->globalOptions[1] && value > this->maxAbsViolation[iter])
                    this->maxAbsViolation[iter] = value;

             }
        
        }

    }
}

void GaugeCondition::timestepAnalysis(const unsigned &time_step)
{
    if (this->globalOutput)
    {
        std::ofstream ofs(std::string(DATA_DIR) + "/" + this->globalQuantitiesPath, std::ios::app);

        // Multiply the integrated quantity by the volume factor now
        for (unsigned iter = 0; iter < this->numEquations; iter++)
            this->integratedAbsViolation[iter] *= this->dx*this->dy*this->dz;

        if (ofs.is_open())
        {
            for (unsigned option_iter = 0; option_iter < this->numGlobalOptions; option_iter++)
            {
                if (this->globalOptions[option_iter])
                {
                    for (unsigned eq_iter = 0; eq_iter < this->numEquations; eq_iter++)
                    {
                        ofs << (*this->globalPointers[option_iter])[eq_iter]<< " ";
                        (*this->globalPointers[option_iter])[eq_iter] = 0.f;
                    }
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
            for (unsigned iter = 0; iter < this->numLocalOptions; iter++)
            {
                if (this->localOptions[iter])
                {
                    for (std::vector<float> &vector : *this->localPointers[iter])
                    {
                        for (float &data : vector)
                        {
                            ofs << data << "\n";
                        }
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
}

void GaugeCondition::finalAnalysis()
{
}
