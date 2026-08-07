#include "OutputFields.hpp"

////////////////////////////////////////////////////  Initialisers  /////////////////////////////////////////////////////////////

void OutputFields::configure(const std::string path, const bool debug)
{
    std::ifstream ifs(path);

    // Read in the parameter values from the configuration file
    std::string description; // Descriptions in the config file (to be dumped)
    if (ifs.is_open())
    {
        // Skip the initial configuration file description
        for(int iter = 0; iter < 13; iter++) std::getline(ifs, description);

        std::getline(ifs, description, ':');
        ifs >> this->outputBothTimesteps;

        std::getline(ifs, description, ':');
        ifs >> this->outputInitial;

        std::getline(ifs, description, ':');
        ifs >> this->outputFinal;


        std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> this->outputContinual;

        std::getline(ifs, description, ':');
        ifs >> this->outputFrequency;


        // Output paths
        std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> this->initialAnalysisPath;

        std::getline(ifs, description, ':');
        ifs >> this->continualAnalysisPath;

        std::getline(ifs, description, ':');
        ifs >> this->finalAnalysisPath;

    } 

    ifs.close();

    if (debug)
    {
        std::cout << "ANALYSERS::OUTPUTFIELDS::\n"
                  << "Output both timesteps?: " << this->outputBothTimesteps << "\n"
                  << "Output initial fields?: " << this->outputInitial << ", Output final fields?: " << this->outputFinal << "\n"
                  << "Output continually?: " << this->outputContinual << ", Every " << this->outputFrequency << " timesteps.\n" 
                  << "Initial analysis data path: Data/" << this->initialAnalysisPath << "\n"
                  << "Continual analysis data path: Data/" << this->continualAnalysisPath << "\n"
                  << "Final analysis data path: Data/" << this->finalAnalysisPath << "\n"
                  << std::endl;
    }
}


void OutputFields::outputFields(const std::string &path,
                                const unsigned &time_step) const
{
    std::ofstream ofs(std::string(DATA_DIR) + "/" + path);

    if (ofs.is_open())
    {
        long long unsigned loop_max;
        long long unsigned start_scalar_index;
        long long unsigned start_vector_index;

        if (this->outputBothTimesteps)
        {
            loop_max
                = this->scalarFields.size()/this->numScalarComponents;

            start_scalar_index = 0;
            start_vector_index = 0;
        }
        else
        {
            const unsigned current_time_index
                = (time_step + 1)%2;

            loop_max
                = this->scalarFields.size()
                  /(2ULL*this->numScalarComponents);

            start_scalar_index
                = current_time_index
                  *loop_max
                  *this->numScalarComponents;

            start_vector_index
                = current_time_index
                  *loop_max
                  *this->numVectorComponents;
        }

        for (long long unsigned iter = 0;
             iter < loop_max;
             iter++)
        {
            // Output the scalar fields at this location.
            for (unsigned compIter = 0;
                 compIter < this->numScalarComponents;
                 compIter++)
            {
                ofs << this->scalarFields[
                    start_scalar_index
                    + iter*this->numScalarComponents
                    + compIter] << " ";
            }

            // Output the vector fields at this location.
            for (unsigned compIter = 0;
                 compIter < this->numVectorComponents;
                 compIter++)
            {
                ofs << this->vectorFields[
                    start_vector_index
                    + iter*this->numVectorComponents
                    + compIter] << " ";
            }

            ofs << std::endl;
        }
    }

    ofs.close();
}





///////////////////////////////////////////////  Constructors/Destructors  //////////////////////////////////////////////////////

OutputFields::OutputFields(const std::vector<float> &scalar_fields, const unsigned &num_scalar_components,
                           const std::vector<float> &vector_fields, const unsigned &num_vector_components)
    : Analyser(), 
      scalarFields(scalar_fields), numScalarComponents(num_scalar_components),
      vectorFields(vector_fields), numVectorComponents(num_vector_components)
{
    this->completedTimesteps = 0;
    this->configure(std::string(SOURCE_DIR) + "/Config/OutputFields.cfg", true);
}

OutputFields::~OutputFields()
{
}

//////////////////////////////////////////////////  Public Functions  ///////////////////////////////////////////////////////////

void OutputFields::initialAnalysis()
{
    if (this->outputInitial)
    {
        this->outputFields(this->initialAnalysisPath, 0);
    }
}

void OutputFields::preEvolveLocationAnalysis(const long long unsigned index,
                                             const float* const local_scalar_pointers[2], const std::vector<std::vector<const float*>> &scalar_pointers,
                                             const float* const local_vector_pointers[2], const std::vector<std::vector<const float*>> &vector_pointers)
{
}

void OutputFields::postEvolveLocationAnalysis(const unsigned &t_now, const long long unsigned index, 
                                              const float *const local_scalar_pointers[2], const std::vector<std::vector<const float *>> &scalar_pointers, 
                                              const float *const local_vector_pointers[2], const std::vector<std::vector<const float *>> &vector_pointers)
{
}

void OutputFields::timestepAnalysis(const unsigned &time_step)
{
    this->completedTimesteps = time_step + 1;

    if (this->outputContinual
        && this->completedTimesteps%this->outputFrequency == 0)
    {
        this->outputFields(
            this->continualAnalysisPath
                + "_"
                + std::to_string(this->completedTimesteps)
                + ".dat",
            this->completedTimesteps);
    }
}

void OutputFields::finalAnalysis()
{
    if (this->outputFinal)
    {
        this->outputFields(
            this->finalAnalysisPath,
            this->completedTimesteps);
    }
}

 