#include "Energy.hpp"

#include <filesystem>
#include <stdexcept>

#ifdef GFT_ENABLE_MPI
#include <mpi.h>
#endif

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

    for (unsigned iter = 0; iter < this->numLocalOptions; iter++)
    {
        if (this->localOptions[iter])
            this->densityPointers[iter]->resize(grid_size, 0.f);
    }

    this->counter = 0;
    this->globalOutput = this->anyGlobalOptions;
    this->localOutput = this->anyLocalOptions;
}

std::string Energy::rankLocalPath(
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

void Energy::mergeRankLocalOutput(
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
            "ANALYSERS::ENERGY:: Invalid MPI site-count data."
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
                "ANALYSERS::ENERGY:: "
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
            "ANALYSERS::ENERGY:: "
            "Could not create merged local output file."
        );
    }

    unsigned output_block = 0U;

    for (unsigned option_iter = 0;
         option_iter < this->numLocalOptions;
         option_iter++)
    {
        if (!this->localOptions[option_iter])
            continue;

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
                    "ANALYSERS::ENERGY:: "
                    "Could not read rank-local output file."
                );
            }

            const unsigned long long site_count
                = owned_site_counts[
                    static_cast<std::size_t>(file_rank)];

            std::string line;

            // Move to the requested quantity block
            // in this rank's file.
            for (unsigned long long line_iter = 0;
                 line_iter
                     < static_cast<unsigned long long>(
                           output_block)*site_count;
                 line_iter++)
            {
                if (!std::getline(rank_file, line))
                {
                    throw std::runtime_error(
                        "ANALYSERS::ENERGY:: "
                        "Rank-local output file is shorter "
                        "than expected."
                    );
                }
            }

            // Copy this rank's contribution to this quantity.
            for (unsigned long long line_iter = 0;
                 line_iter < site_count;
                 line_iter++)
            {
                if (!std::getline(rank_file, line))
                {
                    throw std::runtime_error(
                        "ANALYSERS::ENERGY:: "
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
                    "ANALYSERS::ENERGY:: "
                    "Failed while merging local output."
                );
            }
        }

        output_block++;
    }

    merged.close();

    if (!merged)
    {
        throw std::runtime_error(
            "ANALYSERS::ENERGY:: "
            "Failed to complete merged local output."
        );
    }

    // Do not remove the rank pieces until the merged
    // file has been successfully completed.
    if (fs::exists(output_path))
        fs::remove(output_path);

    fs::rename(temporary_path, output_path);

    for (const fs::path &rank_path : rank_paths)
    {
        if (!fs::remove(rank_path))
        {
            throw std::runtime_error(
                "ANALYSERS::ENERGY:: "
                "Could not remove merged rank-local file: "
                + rank_path.string()
            );
        }
    }
}

///////////////////////////////////  Constructors/Destructors  //////////////////////////////////////////////

Energy::Energy(const Model &model,
               const double &dx,
               const double &dy,
               const double &dz,
               const long long unsigned grid_size,
               const long long unsigned owned_site_begin,
               const long long unsigned owned_site_end,
               const int rank,
               const int num_ranks)
    : model(model),
      energyPointers{&this->energy, &this->potential, &this->gradient,
                     &this->kinetic, &this->magnetic, &this->electric},
      densityPointers{&this->energyDensity, &this->potentialDensity,
                      &this->gradientDensity, &this->kineticDensity,
                      &this->magneticDensity, &this->electricDensity},
      dx(dx), dy(dy), dz(dz),
        rank(rank), numRanks(num_ranks),
        ownedSiteBegin(owned_site_begin),
        ownedSiteEnd(owned_site_end)
{
    this->configure(
        std::string(SOURCE_DIR) + "/Config/Energy.cfg",
        true);

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

void Energy::preEvolveLocationAnalysis(const long long unsigned index,
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
        const float* const density_pointers[this->numLocalOptions] = {&energy_density, &potential_density, &gradient_density, &kinetic_density,
                                                  &magnetic_density, &electric_density};

        if (this->localOutput)
        {
            for (unsigned iter = 0; iter < this->numLocalOptions; iter++)
            {
                if (this->localOptions[iter])
                    (*this->densityPointers[iter])[index] = *density_pointers[iter];
            }
        }

        // Have any global output options been chosen?
        if (this->globalOutput)
        {
            // So I can loop over the integrated quantities and the densities
            for (unsigned iter = 0; iter < this->numGlobalOptions; iter++)
            {
                if (this->globalOptions[iter])
                {
                    // Only bother doing the multiplication by volume at the timestepAnalysis stage do avoid unneccessary computation.
                    *this->energyPointers[iter] += *density_pointers[iter];
                }
            }
            
        }

    }
}

void Energy::postEvolveLocationAnalysis(const unsigned &t_now, const long long unsigned index, const float *const local_scalar_pointers[2], const std::vector<std::vector<const float *>> &scalar_pointers, const float *const local_vector_pointers[2], const std::vector<std::vector<const float *>> &vector_pointers)
{
}

void Energy::timestepAnalysis(const unsigned &time_step)
{
    if (this->globalOutput)
    {
        std::vector<float> local_quantities(
            this->numGlobalOptions, 0.f);

        std::vector<float> global_quantities(
            this->numGlobalOptions, 0.f);

        for (unsigned iter = 0;
            iter < this->numGlobalOptions;
            iter++)
        {
            local_quantities[iter]
                = *this->energyPointers[iter];
        }

        if (this->numRanks == 1)
        {
            global_quantities = local_quantities;
        }
        else
        {
    #ifndef GFT_ENABLE_MPI

            throw std::runtime_error(
                "ANALYSERS::ENERGY:: Multi-rank reduction "
                "requires an MPI build."
            );

    #else

            const int reduce_error = MPI_Reduce(
                local_quantities.data(),
                global_quantities.data(),
                static_cast<int>(this->numGlobalOptions),
                MPI_FLOAT,
                MPI_SUM,
                0,
                MPI_COMM_WORLD);

            if (reduce_error != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    "ANALYSERS::ENERGY:: "
                    "MPI energy reduction failed."
                );
            }

    #endif
        }

        if (this->rank == 0)
        {
            std::ofstream ofs(
                std::string(DATA_DIR)
                    + "/"
                    + this->globalQuantitiesPath,
                std::ios::app);

            if (ofs.is_open())
            {
                for (unsigned iter = 0;
                    iter < this->numGlobalOptions;
                    iter++)
                {
                    if (this->globalOptions[iter])
                    {
                        ofs << global_quantities[iter]
                            *this->dx*this->dy*this->dz
                            << " ";
                    }
                }

                ofs << "\n";
            }

            ofs.close();
        }

        // All ranks must begin the next accumulation from zero.
        for (unsigned iter = 0;
            iter < this->numGlobalOptions;
            iter++)
        {
            *this->energyPointers[iter] = 0.f;
        }
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
                        ofs
                            << (*this->densityPointers[iter])[site_iter]
                            << "\n";
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
    if (this->numRanks == 1 ||
        !this->anyLocalOptions ||
        this->counter == 0)
    {
        return;
    }

#ifndef GFT_ENABLE_MPI

    throw std::runtime_error(
        "ANALYSERS::ENERGY:: "
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
            "ANALYSERS::ENERGY:: "
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
            "ANALYSERS::ENERGY:: "
            "MPI barrier failed after local-output merge."
        );
    }

#endif
}
