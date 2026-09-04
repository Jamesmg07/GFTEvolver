#include "OutputFields.hpp"
#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <stdexcept>

#ifdef GFT_ENABLE_MPI
#include <mpi.h>
#endif

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

        std::getline(ifs, description, ':');
        ifs >> this->outputRMagnitudeEnabled;

        std::getline(ifs, description, ':');
        ifs >> this->RMagnitudeFrequency;

        std::getline(ifs, description, ':');
        ifs >> this->RMagnitudePath;

    } 

    ifs.close();

    if (debug)
    {
        std::cout << "ANALYSERS::OUTPUTFIELDS::\n"
          << "Output both timesteps?: " << this->outputBothTimesteps << "\n"
          << "Output initial fields?: " << this->outputInitial
          << ", Output final fields?: " << this->outputFinal << "\n"
          << "Output continually?: " << this->outputContinual
          << ", Every " << this->outputFrequency << " timesteps.\n"
          << "Initial analysis data path: Data/" << this->initialAnalysisPath << "\n"
          << "Continual analysis data path: Data/" << this->continualAnalysisPath << "\n"
          << "Final analysis data path: Data/" << this->finalAnalysisPath << "\n"
          << "Output R magnitude?: " << this->outputRMagnitudeEnabled
          << ", Every " << this->RMagnitudeFrequency << " timesteps.\n"
          << "R magnitude data path: Data/" << this->RMagnitudePath << "\n"
          << std::endl;
    }
}


std::string OutputFields::rankLocalPath(const std::string &path, const int file_rank) const
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


void OutputFields::mergeRankOutput(
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
            "ANALYSERS::OUTPUTFIELDS:: Invalid MPI site-count data."
        );
    }

    const fs::path output_path
        = fs::path(DATA_DIR)/path;

    const fs::path temporary_path
        = fs::path(output_path.string() + ".merge_tmp");

    std::vector<fs::path> rank_paths;
    rank_paths.reserve(static_cast<std::size_t>(this->numRanks));

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
                "ANALYSERS::OUTPUTFIELDS:: Missing rank-local field file: "
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
            "ANALYSERS::OUTPUTFIELDS:: Could not create merged field file: "
            + temporary_path.string()
        );
    }

    if (!this->outputBothTimesteps)
    {
        for (const fs::path &rank_path : rank_paths)
{
    std::ifstream rank_file(rank_path);

    if (!rank_file.is_open())
    {
        throw std::runtime_error(
            "ANALYSERS::OUTPUTFIELDS:: Could not read "
            "rank-local R magnitude file: "
            + rank_path.string()
        );
    }

    std::string line;

    while (std::getline(rank_file, line))
    {
        merged << line << '\n';

        if (!merged)
        {
            throw std::runtime_error(
                "ANALYSERS::OUTPUTFIELDS:: Failed writing "
                "merged R magnitude output."
            );
        }
    }

    if (rank_file.bad())
    {
        throw std::runtime_error(
            "ANALYSERS::OUTPUTFIELDS:: Failed reading "
            "rank-local R magnitude file: "
            + rank_path.string()
        );
    }
}
    }
    else
    {
        // Each rank file is [t0 slab][t1 slab].
        // Reconstruct [global t0][global t1].
        for (unsigned time_index = 0U;
             time_index < 2U;
             time_index++)
        {
            for (int file_rank = 0;
                 file_rank < this->numRanks;
                 file_rank++)
            {
                std::ifstream rank_file(
                    rank_paths[static_cast<std::size_t>(file_rank)],
                    std::ios::binary);

                if (!rank_file.is_open())
                {
                    throw std::runtime_error(
                        "ANALYSERS::OUTPUTFIELDS:: Could not read rank-local field file."
                    );
                }

                const unsigned long long site_count
                    = owned_site_counts[
                        static_cast<std::size_t>(file_rank)];

                std::string line;

                // For t1, pass over this rank's t0 block first.
                for (unsigned long long line_iter = 0;
                     line_iter < time_index*site_count;
                     line_iter++)
                {
                    if (!std::getline(rank_file, line))
                    {
                        throw std::runtime_error(
                            "ANALYSERS::OUTPUTFIELDS:: Rank-local field file is shorter than expected."
                        );
                    }
                }

                for (unsigned long long line_iter = 0;
                     line_iter < site_count;
                     line_iter++)
                {
                    if (!std::getline(rank_file, line))
                    {
                        throw std::runtime_error(
                            "ANALYSERS::OUTPUTFIELDS:: Rank-local field file is shorter than expected."
                        );
                    }

                    merged.write(
                        line.data(),
                        static_cast<std::streamsize>(line.size()));

                    merged.put('\n');
                }

                if (!merged)
                {
                    throw std::runtime_error(
                        "ANALYSERS::OUTPUTFIELDS:: Failed while merging field output."
                    );
                }
            }
        }
    }

    merged.close();

    if (!merged)
    {
        throw std::runtime_error(
            "ANALYSERS::OUTPUTFIELDS:: Failed to complete merged field file."
        );
    }

    // Only replace/remove files after a complete merged file exists.
    if (fs::exists(output_path))
        fs::remove(output_path);

    fs::rename(temporary_path, output_path);

    for (const fs::path &rank_path : rank_paths)
    {
        if (!fs::remove(rank_path))
        {
            throw std::runtime_error(
                "ANALYSERS::OUTPUTFIELDS:: Could not remove merged rank-local file: "
                + rank_path.string()
            );
        }
    }
}

void OutputFields::mergeRMagnitudeOutput(
    const std::string &path) const
{
    namespace fs = std::filesystem;

    if (this->numRanks == 1)
        return;

    const fs::path output_path =
        fs::path(DATA_DIR) / path;

    const fs::path temporary_path =
        fs::path(output_path.string() + ".merge_tmp");

    std::vector<fs::path> rank_paths;
    rank_paths.reserve(
        static_cast<std::size_t>(this->numRanks));

    // Find every rank-local R^2 file.
    for (int file_rank = 0;
         file_rank < this->numRanks;
         file_rank++)
    {
        const fs::path rank_path =
            fs::path(DATA_DIR)
            / this->rankLocalPath(path, file_rank);

        if (!fs::exists(rank_path))
        {
            throw std::runtime_error(
                "ANALYSERS::OUTPUTFIELDS:: Missing rank-local "
                "R magnitude file: "
                + rank_path.string()
            );
        }

        rank_paths.push_back(rank_path);
    }

    for (const fs::path &rank_path : rank_paths)
    {
        const auto file_size = fs::file_size(rank_path);
        std::cerr << "ANALYSERS::OUTPUTFIELDS:: rank-local file "
                  << rank_path.string() << " size = " << file_size << " bytes\n";

        if (file_size == 0)
        {
            std::cerr << "ANALYSERS::OUTPUTFIELDS:: WARNING - empty rank-local file: "
                      << rank_path.string() << "\n";
        }
    }

    std::ofstream merged(
        temporary_path,
        std::ios::binary | std::ios::trunc);

    if (!merged.is_open())
    {
        throw std::runtime_error(
            "ANALYSERS::OUTPUTFIELDS:: Could not create merged "
            "R magnitude file: "
            + temporary_path.string()
        );
    }

    // Each R^2 file contains exactly one value per owned site.
    // Therefore simply concatenate the rank files in rank order.
    for (const fs::path &rank_path : rank_paths)
    {
        std::ifstream rank_file(
            rank_path,
            std::ios::binary);

        if (!rank_file.is_open())
        {
            throw std::runtime_error(
                "ANALYSERS::OUTPUTFIELDS:: Could not read "
                "rank-local R magnitude file: "
                + rank_path.string()
            );
        }

               merged << rank_file.rdbuf();

        if (rank_file.bad())
        {
            throw std::runtime_error(
                "ANALYSERS::OUTPUTFIELDS:: Failed while merging "
                "R magnitude output."
            );
        }

        if (!merged && !merged.bad())
        {
            // A zero-length rdbuf insertion (empty rank file) legitimately
            // sets failbit even though nothing went wrong; clear it.
            merged.clear(merged.rdstate() & ~std::ios::failbit);
        }
    }

    if (merged.bad())
    {
        throw std::runtime_error(
            "ANALYSERS::OUTPUTFIELDS:: Failed while merging "
            "R magnitude output."
        );
    }

    merged.close();

    if (!merged)
    {
        throw std::runtime_error(
            "ANALYSERS::OUTPUTFIELDS:: Failed to complete merged "
            "R magnitude file."
        );
    }

    // Only replace the final file once the merge is complete.
    if (fs::exists(output_path))
        fs::remove(output_path);

    fs::rename(temporary_path, output_path);

    // Remove rank-local files after successful merge.
    for (const fs::path &rank_path : rank_paths)
    {
        if (!fs::remove(rank_path))
        {
            throw std::runtime_error(
                "ANALYSERS::OUTPUTFIELDS:: Could not remove merged "
                "R magnitude rank-local file: "
                + rank_path.string()
            );
        }
    }
}


void OutputFields::outputFields(const std::string &path,
                                const unsigned &time_step) const
{
    std::ofstream ofs(
        std::string(DATA_DIR) + "/" + this->rankLocalPath(path, this->rank));

    if (ofs.is_open())
    {
        const unsigned first_time_index
            = this->outputBothTimesteps
            ? 0U
            : (time_step + 1)%2;

        const unsigned time_count
            = this->outputBothTimesteps ? 2U : 1U;

        for (unsigned time_iter = 0;
             time_iter < time_count;
             time_iter++)
        {
            const unsigned time_index
                = this->outputBothTimesteps
                ? time_iter
                : first_time_index;

            const unsigned long long scalar_time_offset
                = 1ULL*time_index*this->storageVolume
                  *this->numScalarComponents;

            const unsigned long long vector_time_offset
                = 1ULL*time_index*this->storageVolume
                  *this->numVectorComponents;

            for (unsigned long long site_iter = this->ownedSiteBegin;
                 site_iter < this->ownedSiteEnd;
                 site_iter++)
            {
                const unsigned long long scalar_index
                    = scalar_time_offset
                    + site_iter*this->numScalarComponents;

                const unsigned long long vector_index
                    = vector_time_offset
                    + site_iter*this->numVectorComponents;

                for (unsigned comp_iter = 0;
                     comp_iter < this->numScalarComponents;
                     comp_iter++)
                {
                    ofs << this->scalarFields[
                        scalar_index + comp_iter] << " ";
                }

                for (unsigned comp_iter = 0;
                     comp_iter < this->numVectorComponents;
                     comp_iter++)
                {
                    ofs << this->vectorFields[
                        vector_index + comp_iter] << " ";
                }

                ofs << std::endl;
            }
        }
    }

    ofs.close();
}


void OutputFields::outputRMagnitude(
    const std::string &path,
    const unsigned &time_step) const
{
    std::ofstream ofs(
        std::string(DATA_DIR) + "/" +
        this->rankLocalPath(path, this->rank));

    if (!ofs.is_open())
    {
        throw std::runtime_error(
            "ANALYSERS::OUTPUTFIELDS:: Could not open R magnitude output file: "
            + std::string(DATA_DIR) + "/" +
            this->rankLocalPath(path, this->rank)
        );
    }

    // Use the same timestep convention as outputFields().
    const unsigned time_index =
        (time_step + 1) % 2;

    const unsigned long long scalar_time_offset =
        1ULL * time_index
        * this->storageVolume
        * this->numScalarComponents;

    for (unsigned long long site_iter = this->ownedSiteBegin;
         site_iter < this->ownedSiteEnd;
         site_iter++)
    {
        const unsigned long long scalar_index =
            scalar_time_offset
            + site_iter * this->numScalarComponents;

        // The eight real scalar fields:
        //
        // phi1 = (psi1 + i psi2, psi3 + i psi4)
        // phi2 = (psi5 + i psi6, psi7 + i psi8)

        const double psi1 =
            this->scalarFields[scalar_index + 0];

        const double psi2 =
            this->scalarFields[scalar_index + 1];

        const double psi3 =
            this->scalarFields[scalar_index + 2];

        const double psi4 =
            this->scalarFields[scalar_index + 3];

        const double psi5 =
            this->scalarFields[scalar_index + 4];

        const double psi6 =
            this->scalarFields[scalar_index + 5];

        const double psi7 =
            this->scalarFields[scalar_index + 6];

        const double psi8 =
            this->scalarFields[scalar_index + 7];


        ////////////////////////////////////////////////////////////////
        // |phi1|^2 and |phi2|^2
        ////////////////////////////////////////////////////////////////

        const double phi1_norm_sq =
              psi1*psi1 + psi2*psi2
            + psi3*psi3 + psi4*psi4;

        const double phi2_norm_sq =
              psi5*psi5 + psi6*psi6
            + psi7*psi7 + psi8*psi8;


        ////////////////////////////////////////////////////////////////
        // R^mu
        ////////////////////////////////////////////////////////////////


        const double R1 =
            2.0 * (
                  psi1*psi5 + psi2*psi6
                + psi3*psi7 + psi4*psi8
                );

        const double R2 =
            2.0 * (
                  psi1*psi6 - psi2*psi5
                + psi3*psi8 - psi4*psi7
                );

        const double R3 =
            phi1_norm_sq - phi2_norm_sq;


        ////////////////////////////////////////////////////////////////
        // R^2
        ////////////////////////////////////////////////////////////////

        const double R_squared = R1*R1
            + R2*R2
            + R3*R3;


        ofs << R_squared << "\n";

        if (!ofs)
        {
            throw std::runtime_error(
                "ANALYSERS::OUTPUTFIELDS:: Failed writing R magnitude value "
                "for site " + std::to_string(site_iter));
        }
    }

    ofs.close();

    if (ofs.fail())
    {
        throw std::runtime_error(
            "ANALYSERS::OUTPUTFIELDS:: Failed to close R magnitude output file."
        );
    }
}

///////////////////////////////////////////////  Constructors/Destructors  //////////////////////////////////////////////////////

OutputFields::OutputFields(
    const std::vector<float> &scalar_fields,
    const unsigned &num_scalar_components,
    const std::vector<float> &vector_fields,
    const unsigned &num_vector_components,
    const long long unsigned storage_volume,
    const long long unsigned owned_site_begin,
    const long long unsigned owned_site_end,
    const int rank,
    const int num_ranks)
    : Analyser(),
      scalarFields(scalar_fields),
      vectorFields(vector_fields),
      numScalarComponents(num_scalar_components),
      numVectorComponents(num_vector_components),
      storageVolume(storage_volume),
      ownedSiteBegin(owned_site_begin),
      ownedSiteEnd(owned_site_end),
      rank(rank),
      numRanks(num_ranks)
{
    this->completedTimesteps = 0;


    this->configure(
        std::string(SOURCE_DIR) + "/Config/OutputFields.cfg",
        true);
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

    if (this->outputRMagnitudeEnabled)
    {
        const std::string r_magnitude_path =
            this->RMagnitudePath.substr(
                0,
                this->RMagnitudePath.find_last_of('.'))
            + "_0.dat";

        // Every rank calculates and writes its own R^2 values,
        // exactly as timestepAnalysis() does at later timesteps.
        this->outputRMagnitude(
            r_magnitude_path,
            0);

#ifdef GFT_ENABLE_MPI

        // Make sure every rank has finished writing before rank 0
        // attempts to read the rank-local files.
        if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
        {
            throw std::runtime_error(
                "ANALYSERS::OUTPUTFIELDS:: MPI barrier failed "
                "before initial R magnitude merge."
            );
        }

        // Rank 0 reconstructs the global file.
        if (this->rank == 0)
        {
            this->mergeRMagnitudeOutput(
                r_magnitude_path);
        }

        // Do not let the other ranks continue until the merge is
        // completely finished.
        if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
        {
            throw std::runtime_error(
                "ANALYSERS::OUTPUTFIELDS:: MPI barrier failed "
                "after initial R magnitude merge."
            );
        }

#endif
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
        && this->completedTimesteps % this->outputFrequency == 0)
    {
        this->outputFields(
            this->continualAnalysisPath
                + "_"
                + std::to_string(this->completedTimesteps)
                + ".dat",
            this->completedTimesteps);
    }

    if (this->outputRMagnitudeEnabled
    && this->completedTimesteps % this->RMagnitudeFrequency == 0)
{
    const std::string r_magnitude_path =
        this->RMagnitudePath.substr(
            0,
            this->RMagnitudePath.find_last_of('.'))
        + "_"
        + std::to_string(this->completedTimesteps)
        + ".dat";

    // Every rank calculates and writes its own R^2 values.
    this->outputRMagnitude(
        r_magnitude_path,
        this->completedTimesteps);

#ifdef GFT_ENABLE_MPI

    // Make sure every rank has finished writing before rank 0
    // attempts to read the rank-local files.
    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
    {
        throw std::runtime_error(
            "ANALYSERS::OUTPUTFIELDS:: MPI barrier failed "
            "before R magnitude merge."
        );
    }

    // Rank 0 reconstructs the global file.
    if (this->rank == 0)
    {
        this->mergeRMagnitudeOutput(
            r_magnitude_path);
    }

    // Do not let the other ranks continue until the merge is
    // completely finished.
    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
    {
        throw std::runtime_error(
            "ANALYSERS::OUTPUTFIELDS:: MPI barrier failed "
            "after R magnitude merge."
        );
    }

#endif
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

    if (this->numRanks == 1)
        return;

    std::vector<std::string> output_paths;

    const auto add_output_path = [&output_paths](const std::string &path)
    {
        if (std::find(output_paths.begin(), output_paths.end(), path)
            == output_paths.end())
        {
            output_paths.push_back(path);
        }
    };

    if (this->outputInitial)
        add_output_path(this->initialAnalysisPath);

    if (this->outputContinual)
    {
        for (unsigned timestep = this->outputFrequency;
             timestep <= this->completedTimesteps;
             timestep += this->outputFrequency)
        {
            add_output_path(
                this->continualAnalysisPath
                + "_"
                + std::to_string(timestep)
                + ".dat");
        }
    }

    if (this->outputFinal)
        add_output_path(this->finalAnalysisPath);

    if (output_paths.empty())
        return;

#ifndef GFT_ENABLE_MPI

    throw std::runtime_error(
        "ANALYSERS::OUTPUTFIELDS:: Multi-rank output merging requires "
        "MPI build."
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
        this->rank == 0 ? owned_site_counts.data() : nullptr,
        1,
        MPI_UNSIGNED_LONG_LONG,
        0,
        MPI_COMM_WORLD);

    if (gather_error != MPI_SUCCESS)
    {
        throw std::runtime_error(
            "ANALYSERS::OUTPUTFIELDS:: Failed to gather MPI output metadata."
        );
    }

    if (this->rank == 0)
    {
        for (const std::string &path : output_paths)
            this->mergeRankOutput(path, owned_site_counts);
    }

    // Hold the other ranks here while rank 0 completes the disk merge.
    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
    {
        throw std::runtime_error(
            "ANALYSERS::OUTPUTFIELDS:: MPI barrier failed after field merge."
        );
    }

#endif
}

 