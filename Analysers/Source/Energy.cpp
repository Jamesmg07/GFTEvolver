#include "Energy.hpp"

#include <filesystem>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <limits>

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

        // Load in monopole maximum-position output path
        std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> this->monoSeparationEnabled;
        
        std::getline(ifs, description, ':');
        ifs >> this->monoSeparationPath;

        // Load in monopole maximum-position output frequency
        std::getline(ifs, description, ':');
        ifs >> this->monoSeparationFrequency;

        // Load in minimum physical monopole separation
        std::getline(ifs, description, ':');
        ifs >> this->monoMinimumSeparation;

    } 

    ifs.close();



    if (this->rank == 0)
    {
        std::ofstream ofs(
            std::string(DATA_DIR)
            + "/"
            + this->monoSeparationPath);

        ofs.close();
    }

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
                  << ", Magnetic: " << this->localOptions[4] << ", Electric: " << this->localOptions[5] << ", every " << this->localFrequency << " timesteps\n"
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

    // Allocate energy density if either:
    // 1. local energy-density output is enabled, or
    // 2. monopole tracking is enabled.
    //
    // Monopole tracking needs the total energy density even when
    // local energy-density output is disabled.
    if (this->localOptions[0] ||
        this->monoSeparationEnabled)
    {
        this->energyDensity.resize(grid_size, 0.f);
    }

    // The remaining local density arrays are only needed when
    // their corresponding local output is enabled.
    for (unsigned iter = 1;
        iter < this->numLocalOptions;
        iter++)
    {
        if (this->localOptions[iter])
        {
            this->densityPointers[iter]->resize(
                grid_size,
                0.f);
        }
    }

    this->counter = 0;
  

    this->globalOutput = this->anyGlobalOptions;
    this->localOutput = this->anyLocalOptions;

    this->monoSeparationOutput =
    this->monoSeparationEnabled &&
    this->counter % this->monoSeparationFrequency == 0;

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
               const double &dx, const double &dy, const double &dz,
               const long long unsigned grid_size,
               const long long unsigned owned_site_begin,
               const long long unsigned owned_site_end,
               const int rank, const int num_ranks,
               const unsigned nx, const unsigned ny, const unsigned nz,
               const unsigned global_x_start)
    : model(model), energyPointers{&this->energy, &this->potential, &this->gradient,
                     &this->kinetic, &this->magnetic, &this->electric},
      densityPointers{&this->energyDensity, &this->potentialDensity,
                      &this->gradientDensity, &this->kineticDensity,
                      &this->magneticDensity, &this->electricDensity},
      dx(dx), dy(dy), dz(dz),
      rank(rank), numRanks(num_ranks),
      ownedSiteBegin(owned_site_begin), ownedSiteEnd(owned_site_end),
      nx(nx), ny(ny), nz(nz), globalXStart(global_x_start)
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
    this->model.energyPreparation(
        this->globalOutput ||
        this->localOutput ||
        this->monoSeparationOutput);
}




void Energy::preEvolveLocationAnalysis(const long long unsigned index,
                                        const float* const local_scalar_pointers[2],
                                        const std::vector<std::vector<const float*>> &scalar_pointers,
                                        const float* const local_vector_pointers[2],
                                        const std::vector<std::vector<const float*>> &vector_pointers)
{

    
    if (this->localOutput ||
        this->globalOutput ||
        this->monoSeparationOutput)
    {
    

        float potential_density = this->model.calcPotentialEnergy(local_scalar_pointers[1]);
        float gradient_density = this->model.calcGradientEnergy(scalar_pointers, vector_pointers);
        float kinetic_density = this->model.calcKineticEnergy(local_scalar_pointers);
        float magnetic_density = this->model.calcMagneticEnergy(vector_pointers);
        float electric_density = this->model.calcElectricEnergy(local_vector_pointers);

        float energy_density = potential_density + gradient_density + kinetic_density + magnetic_density + electric_density;

        // Always store total energy density when maxima analysis
        // is enabled, even if local output is disabled.
        if (this->localOptions[0] ||
    this->monoSeparationOutput)
        {
            this->energyDensity[index] = energy_density;
        }

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
    ////////////////////////////////////////////////////////////////
    // Find the two strongest local maxima of the 3D energy density
    ////////////////////////////////////////////////////////////////

    if (this->monoSeparationOutput)
    {
        struct Candidate { float value; unsigned long long x, y, z; };

        const unsigned long long planeSize = 1ULL*this->ny*this->nz;
        const unsigned long long ownedXBeginLocal = this->ownedSiteBegin / planeSize;
        const unsigned long long ownedXEndLocal   = this->ownedSiteEnd   / planeSize;

        auto flatIndex = [this](unsigned long long xl, unsigned long long y, unsigned long long z)
        { return (xl*this->ny + y)*this->nz + z; };

        std::vector<Candidate> localCandidates;

        for (unsigned long long xl = ownedXBeginLocal; xl < ownedXEndLocal; xl++)
        {
            const unsigned long long gx = this->globalXStart + (xl - ownedXBeginLocal);
            if (gx == 0 || gx == this->nx - 1) continue;

            for (unsigned long long y = 1; y < this->ny - 1; y++)
            for (unsigned long long z = 1; z < this->nz - 1; z++)
            {
                const unsigned long long idx = flatIndex(xl, y, z);
                const float value = this->energyDensity[idx];
                bool isMax = true;

                for (int ddx = -1; ddx <= 1 && isMax; ddx++)
                for (int ddy = -1; ddy <= 1 && isMax; ddy++)
                for (int ddz = -1; ddz <= 1; ddz++)
                {
                    if (!ddx && !ddy && !ddz) continue;
                    if (this->energyDensity[flatIndex(xl+ddx, y+ddy, z+ddz)] > value)
                    { isMax = false; break; }
                }
                if (isMax) localCandidates.push_back({value, gx, y, z});
            }
        }

        // --- Assemble every rank's candidates on rank 0 ---
        std::vector<Candidate> all;

        if (this->numRanks == 1)
        {
            all = localCandidates;
        }
        else
        {
#ifndef GFT_ENABLE_MPI

            throw std::runtime_error(
                "ANALYSERS::ENERGY:: Multi-rank monopole tracking "
                "requires an MPI build."
            );

#else

            int localCount = static_cast<int>(localCandidates.size());
            std::vector<int> counts(this->rank == 0 ? this->numRanks : 0);

            if (MPI_Gather(&localCount, 1, MPI_INT,
                    this->rank == 0 ? counts.data() : nullptr,
                    1, MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    "ANALYSERS::ENERGY:: MPI candidate-count gather failed."
                );
            }

            std::vector<int> displs;
            if (this->rank == 0)
            {
                displs.resize(this->numRanks);
                int running = 0;
                for (int r = 0; r < this->numRanks; r++) { displs[r] = running; running += counts[r]; }
                all.resize(running);
            }

            MPI_Datatype candType;
            MPI_Type_contiguous(sizeof(Candidate), MPI_BYTE, &candType);
            MPI_Type_commit(&candType);

            const int gatherv_error = MPI_Gatherv(
                localCandidates.data(), localCount, candType,
                this->rank == 0 ? all.data() : nullptr,
                this->rank == 0 ? counts.data() : nullptr,
                this->rank == 0 ? displs.data() : nullptr,
                candType, 0, MPI_COMM_WORLD);

            MPI_Type_free(&candType);

            if (gatherv_error != MPI_SUCCESS)
            {
                throw std::runtime_error(
                    "ANALYSERS::ENERGY:: MPI candidate gather failed."
                );
            }

#endif
        }

        if (this->rank == 0)
        {
            auto worse = [](const Candidate &a, const Candidate &b)
            {
                if (a.value != b.value) return a.value < b.value;
                if (a.x != b.x) return a.x > b.x;
                if (a.y != b.y) return a.y > b.y;
                return a.z > b.z;
            };

            Candidate first{-std::numeric_limits<float>::infinity(), 0, 0, 0};
            for (auto &c : all) if (worse(first, c)) first = c;
            const bool firstFound = first.value > -std::numeric_limits<float>::infinity();

            Candidate second{-std::numeric_limits<float>::infinity(), 0, 0, 0};
            if (firstFound)
            {
                const double minSepSq = this->monoMinimumSeparation * this->monoMinimumSeparation;
                for (auto &c : all)
                {
                    if (c.x == first.x && c.y == first.y && c.z == first.z) continue;
                    const double ddx = (double(c.x) - double(first.x)) * this->dx;
                    const double ddy = (double(c.y) - double(first.y)) * this->dy;
                    const double ddz = (double(c.z) - double(first.z)) * this->dz;
                    if (ddx*ddx + ddy*ddy + ddz*ddz < minSepSq) continue;
                    if (worse(second, c)) second = c;
                }
            }
            const bool secondFound = second.value > -std::numeric_limits<float>::infinity();

            std::ofstream ofs(std::string(DATA_DIR) + "/" + this->monoSeparationPath, std::ios::app);
            if (ofs.is_open())
            {
                if (firstFound)
                {
                    ofs << time_step << " "
                        << first.x*this->dx << " " << first.y*this->dy << " " << first.z*this->dz << " "
                        << first.value << " ";
                    if (secondFound)
                        ofs << second.x*this->dx << " " << second.y*this->dy << " " << second.z*this->dz << " " << second.value;
                    else
                        ofs << "nan nan nan nan";
                }
                else
                {
                    ofs << time_step << " nan nan nan nan nan nan nan nan";
                }
                ofs << "\n";
            }
        }
    }

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
                ofs << time_step << " ";

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

    this->counter++;

    // Check if conditions for global and local output are satisfied.
    this->globalOutput =
        this->anyGlobalOptions &&
        this->counter % this->globalFrequency == 0;

    this->localOutput =
        this->anyLocalOptions &&
        this->counter % this->localFrequency == 0;

    // Maximum-position analysis is independent of global/local output.
    this->monoSeparationOutput =
        this->monoSeparationEnabled &&
        this->counter % this->monoSeparationFrequency == 0;

    this->model.energyPreparation(
        this->globalOutput ||
        this->localOutput ||
        this->monoSeparationOutput);
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
