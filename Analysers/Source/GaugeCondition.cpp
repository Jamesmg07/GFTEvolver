#include "GaugeCondition.hpp"

#include <algorithm>
#include <filesystem>
#include <stdexcept>

#ifdef GFT_ENABLE_MPI
#include <mpi.h>
#endif

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
    if (this->anyGlobalOptions && this->rank == 0)
    {
        std::ofstream ofs(std::string(DATA_DIR) + "/" + this->globalQuantitiesPath);
        ofs.close();
    }

    if (debug && this->rank == 0)
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

std::string GaugeCondition::rankLocalPath(
    const std::string &path,
    const int file_rank) const
{
    if (this->numRanks == 1)
        return path;

    const std::string suffix
        = "_rank" + std::to_string(file_rank);

    const std::size_t separator_position
        = path.find_last_of("/\\");

    const std::size_t extension_position
        = path.find_last_of('.');

    if (extension_position == std::string::npos ||
        (separator_position != std::string::npos &&
         extension_position < separator_position))
    {
        return path + suffix;
    }

    return path.substr(0, extension_position)
         + suffix
         + path.substr(extension_position);
}

void GaugeCondition::mergeRankLocalOutput(
    const std::string &path,
    const std::vector<unsigned long long> &owned_site_counts) const
{
    namespace fs = std::filesystem;

    if (this->numRanks == 1)
        return;

    if (owned_site_counts.size()
        != static_cast<std::size_t>(this->numRanks))
    {
        throw std::runtime_error(
            "ANALYSERS::GAUGECONDITION:: "
            "Invalid MPI site-count data."
        );
    }

    const fs::path output_path
        = fs::path(DATA_DIR) / path;

    const fs::path temporary_path
        = fs::path(output_path.string() + ".merge_tmp");

    std::vector<fs::path> rank_paths;
    rank_paths.reserve(
        static_cast<std::size_t>(this->numRanks));

    for (int file_rank = 0;
         file_rank < this->numRanks;
         file_rank++)
    {
        const fs::path rank_path
            = fs::path(DATA_DIR)
            / this->rankLocalPath(path, file_rank);

        if (!fs::exists(rank_path))
        {
            throw std::runtime_error(
                "ANALYSERS::GAUGECONDITION:: "
                "Missing rank-local output file: "
                + rank_path.string()
            );
        }

        rank_paths.push_back(rank_path);
    }

    std::ofstream merged(
        temporary_path,
        std::ios::binary | std::ios::trunc);

    if (!merged.is_open())
    {
        throw std::runtime_error(
            "ANALYSERS::GAUGECONDITION:: "
            "Could not create merged local output file."
        );
    }

    // Gauge has only one local output block, so rank order
    // is already global x order.
    for (int file_rank = 0;
         file_rank < this->numRanks;
         file_rank++)
    {
        std::ifstream rank_file(
            rank_paths[
                static_cast<std::size_t>(file_rank)],
            std::ios::binary);

        if (!rank_file.is_open())
        {
            throw std::runtime_error(
                "ANALYSERS::GAUGECONDITION:: "
                "Could not read rank-local output file."
            );
        }

        const unsigned long long line_count
            = owned_site_counts[
                  static_cast<std::size_t>(file_rank)]
              * static_cast<unsigned long long>(
                    this->numEquations);

        std::string line;

        for (unsigned long long line_iter = 0;
             line_iter < line_count;
             line_iter++)
        {
            if (!std::getline(rank_file, line))
            {
                throw std::runtime_error(
                    "ANALYSERS::GAUGECONDITION:: "
                    "Rank-local output file is shorter "
                    "than expected."
                );
            }

            merged.write(
                line.data(),
                static_cast<std::streamsize>(
                    line.size()));

            merged.put('\n');
        }

        if (!merged)
        {
            throw std::runtime_error(
                "ANALYSERS::GAUGECONDITION:: "
                "Failed while merging local output."
            );
        }
    }

    merged.close();

    if (!merged)
    {
        throw std::runtime_error(
            "ANALYSERS::GAUGECONDITION:: "
            "Failed to complete merged local output."
        );
    }

    if (fs::exists(output_path))
        fs::remove(output_path);

    fs::rename(temporary_path, output_path);

    for (const fs::path &rank_path : rank_paths)
    {
        if (!fs::remove(rank_path))
        {
            throw std::runtime_error(
                "ANALYSERS::GAUGECONDITION:: "
                "Could not remove merged rank-local file: "
                + rank_path.string()
            );
        }
    }
}

///////////////////////////////////  Constructors/Destructors  //////////////////////////////////////////////

GaugeCondition::GaugeCondition(
    const Model &model,
    const double &dx,
    const double &dy,
    const double &dz,
    const long long unsigned grid_size,
    const unsigned &num_vector_components,
    const long long unsigned owned_site_begin,
    const long long unsigned owned_site_end,
    const int rank,
    const int num_ranks)
    : model(model),
      globalPointers{
          &this->integratedAbsViolation,
          &this->maxAbsViolation},
      localPointers{&this->localViolation},
      dx(dx), dy(dy), dz(dz),
        gridSize(grid_size),
        numVectorComponents(num_vector_components),
        rank(rank), numRanks(num_ranks),
        ownedSiteBegin(owned_site_begin),
        ownedSiteEnd(owned_site_end)
{
    this->configure(
        std::string(SOURCE_DIR)
            + "/Config/GaugeCondition.cfg",
        true);
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

                // First option integrates the absolute violation over the
                // evolved-site diagnostic domain.
                // Multiplcation by the volume factor will be done later to avoid unneccessary computation.
                if (this->globalOptions[0])
                    this->integratedAbsViolation[iter] += value;

                // Second option finds the largest absolute violation in the
                // evolved-site diagnostic domain.
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
        std::vector<float> global_integrated
            = this->integratedAbsViolation;

        std::vector<float> global_maximum
            = this->maxAbsViolation;

        if (this->numRanks > 1)
        {
    #ifndef GFT_ENABLE_MPI

            throw std::runtime_error(
                "ANALYSERS::GAUGECONDITION:: Multi-rank reduction "
                "requires an MPI build."
            );

    #else

            if (this->globalOptions[0])
            {
                const int reduce_error = MPI_Reduce(
                    this->integratedAbsViolation.data(),
                    global_integrated.data(),
                    static_cast<int>(this->numEquations),
                    MPI_FLOAT,
                    MPI_SUM,
                    0,
                    MPI_COMM_WORLD);

                if (reduce_error != MPI_SUCCESS)
                {
                    throw std::runtime_error(
                        "ANALYSERS::GAUGECONDITION:: "
                        "MPI integrated-violation reduction failed."
                    );
                }
            }

            if (this->globalOptions[1])
            {
                const int reduce_error = MPI_Reduce(
                    this->maxAbsViolation.data(),
                    global_maximum.data(),
                    static_cast<int>(this->numEquations),
                    MPI_FLOAT,
                    MPI_MAX,
                    0,
                    MPI_COMM_WORLD);

                if (reduce_error != MPI_SUCCESS)
                {
                    throw std::runtime_error(
                        "ANALYSERS::GAUGECONDITION:: "
                        "MPI maximum-violation reduction failed."
                    );
                }
            }

    #endif
        }

        if (this->rank == 0)
        {
            // Preserve the original serial rounding behaviour.
            if (this->globalOptions[0])
            {
                for (unsigned iter = 0;
                    iter < this->numEquations;
                    iter++)
                {
                    global_integrated[iter]
                        *= this->dx*this->dy*this->dz;
                }
            }

            std::ofstream ofs(
                std::string(DATA_DIR)
                    + "/"
                    + this->globalQuantitiesPath,
                std::ios::app);

            if (ofs.is_open())
            {
                if (this->globalOptions[0])
                {
                    for (unsigned eq_iter = 0;
                        eq_iter < this->numEquations;
                        eq_iter++)
                    {
                        ofs << global_integrated[eq_iter]
                            << " ";
                    }
                }

                if (this->globalOptions[1])
                {
                    for (unsigned eq_iter = 0;
                        eq_iter < this->numEquations;
                        eq_iter++)
                    {
                        ofs << global_maximum[eq_iter]
                            << " ";
                    }
                }

                ofs << "\n";
            }

            ofs.close();
        }

        std::fill(
            this->integratedAbsViolation.begin(),
            this->integratedAbsViolation.end(),
            0.f);

        std::fill(
            this->maxAbsViolation.begin(),
            this->maxAbsViolation.end(),
            0.f);
    }

    if (this->localOutput)
    {
        const std::string path
            = this->localQuantitiesPath
            + "_"
            + std::to_string(counter)
            + ".dat";

        std::ofstream ofs(
            std::string(DATA_DIR)
            + "/"
            + this->rankLocalPath(path, this->rank));

        if (ofs.is_open())
        {
            for (unsigned iter = 0;
                iter < this->numLocalOptions;
                iter++)
            {
                if (this->localOptions[iter])
                {
                    for (unsigned long long site_iter
                            = this->ownedSiteBegin;
                        site_iter < this->ownedSiteEnd;
                        site_iter++)
                    {
                        for (float data :
                            (*this->localPointers[iter])[site_iter])
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
    if (this->numRanks == 1 ||
        !this->anyLocalOptions ||
        this->counter == 0)
    {
        return;
    }

#ifndef GFT_ENABLE_MPI

    throw std::runtime_error(
        "ANALYSERS::GAUGECONDITION:: "
        "Multi-rank local output merging requires an MPI build."
    );

#else

    const unsigned long long local_owned_site_count
        = this->ownedSiteEnd - this->ownedSiteBegin;

    std::vector<unsigned long long> owned_site_counts(
        this->rank == 0
            ? static_cast<std::size_t>(this->numRanks)
            : 0U);

    const int gather_error = MPI_Gather(
        &local_owned_site_count,
        1,
        MPI_UNSIGNED_LONG_LONG,
        this->rank == 0
            ? owned_site_counts.data()
            : nullptr,
        1,
        MPI_UNSIGNED_LONG_LONG,
        0,
        MPI_COMM_WORLD);

    if (gather_error != MPI_SUCCESS)
    {
        throw std::runtime_error(
            "ANALYSERS::GAUGECONDITION:: "
            "Failed to gather MPI output metadata."
        );
    }

    if (this->rank == 0)
    {
        for (unsigned output_counter = 0;
             output_counter < this->counter;
             output_counter += this->localFrequency)
        {
            const std::string path
                = this->localQuantitiesPath
                + "_"
                + std::to_string(output_counter)
                + ".dat";

            this->mergeRankLocalOutput(
                path,
                owned_site_counts);
        }
    }

    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
    {
        throw std::runtime_error(
            "ANALYSERS::GAUGECONDITION:: "
            "MPI barrier failed after local-output merge."
        );
    }

#endif
}
