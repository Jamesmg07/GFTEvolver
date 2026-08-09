#include "UserDefinedProfile.hpp"

#include <cstddef>
#include <fstream>
#include <iostream>
#include <stdexcept>


////////////////////////////////////////////  Initialisers  //////////////////////////////////////////////////////

void UserDefinedProfile::configure(
    const std::string path,
    const bool debug)
{
    std::ifstream ifs(path);

    if (!ifs.is_open())
    {
        throw std::runtime_error(
            "INITIALCONDITIONS::USERDEFINEDPROFILE:: "
            "Could not open configuration file: " + path
        );
    }

    std::string description;

    // Skip the initial configuration-file description.
    for (unsigned iter = 0; iter < 14U; iter++)
    {
        std::getline(ifs, description);
    }

    std::getline(ifs, description, ':');
    ifs >> this->scalarProfileDataFilename;

    std::getline(ifs, description, ':');
    ifs >> this->gaugeProfileDataFilename;

    if (!ifs)
    {
        throw std::runtime_error(
            "INITIALCONDITIONS::USERDEFINEDPROFILE:: "
            "Could not read the profile configuration."
        );
    }

    ifs.close();

    const std::string &selected_filename
        = this->fieldType == FieldType::Scalar
        ? this->scalarProfileDataFilename
        : this->gaugeProfileDataFilename;

    /*
     * USER SECTION: OPTIONAL PROFILE-DATA LOADING
     *
     * If the profile uses a small tabulated data file, load it here once and
     * store the resulting data in class members declared in the header.
     *
     * Do not open a file inside setScalarProfile() or setGaugeProfile(), since
     * those functions are called once for every owned lattice site.
     *
     * For a path relative to InitialConditions, one possible convention is:
     *
     * const std::string data_path = std::string(SOURCE_DIR) + "/" + selected_filename;
     *
     * Each MPI rank may independently read a small shared profile table.
     */

    if (debug)
    {
        const char *field_name
            = this->fieldType == FieldType::Scalar
            ? "scalar"
            : "gauge";

        std::cout
            << "INITIALCONDITIONS::USERDEFINEDPROFILE::\n"
            << "Field type: " << field_name
            << ", profile-data filename: "
            << selected_filename;

        if (this->fieldType == FieldType::Gauge)
        {
            std::cout
                << ", quaternion representation?: "
                << this->usingQuaternionRepresentation;
        }

        std::cout << "\n" << std::endl;
    }
}


////////////////////////////////////////  Constructors/Destructors  //////////////////////////////////////////////

UserDefinedProfile::UserDefinedProfile(
    const FieldType field_type,
    const bool using_quaternion_representation)
    : fieldType(field_type),
      usingQuaternionRepresentation(using_quaternion_representation),
      scalarProfileDataFilename("NONE"),
      gaugeProfileDataFilename("NONE")
{
    this->configure(
        std::string(SOURCE_DIR)
            + "/Config/UserDefinedProfile.cfg",
        true
    );

    /*
     * This warning is deliberately retained in the student template.
     * Remove it once the corresponding profile has been implemented.
     */
    std::cout
        << "INITIALCONDITIONS::USERDEFINEDPROFILE::WARNING: "
        << "The supplied user-defined profile is a blank template.\n"
        << "Representation-safe default field values will be retained "
        << "until its profile functions are implemented.\n"
        << std::endl;
}

UserDefinedProfile::~UserDefinedProfile()
{
}



////////////////////////////////////////////  Private Functions  /////////////////////////////////////////////////

void UserDefinedProfile::setScalarProfile(
    float *previous_fields,
    float *current_fields,
    const unsigned global_x,
    const unsigned global_y,
    const unsigned global_z,
    const InitialConditionGeometry &geometry,
    const unsigned num_components) const
{
    /*
     * USER SECTION: SCALAR INITIAL CONDITIONS
     *
     * 1. Construct physical coordinates using the GLOBAL lattice indices.
     *    For a profile centred on the middle of the lattice, for example:
     *
     *     const double x = (static_cast<double>(global_x) - 0.5*static_cast<double>(geometry.globalNx - 1U))*geometry.dx;
     *
     *     const double y = (static_cast<double>(global_y) - 0.5*static_cast<double>(geometry.globalNy - 1U))*geometry.dy;
     *
     *     const double z = (static_cast<double>(global_z) - 0.5*static_cast<double>(geometry.globalNz - 1U)) *geometry.dz;
     *
     * 2. Evaluate an analytic profile, or interpolate data loaded once in
     *    configure().
     *
     * 3. Assign current_fields[component] for every required scalar
     *    component. The meaning of each component is fixed by the chosen model.
     *
     * 4. For zero initial time derivative, copy the current values:
     *
     *      for (unsigned comp = 0; comp < num_components; comp++)
     *      {
     *            previous_fields[comp] = current_fields[comp];
     *      }
     *
     *      To encode a non-zero initial derivative v using the backward time
     *      level, a common finite-difference choice is
     *
     *      previous_fields[comp] = current_fields[comp] - geometry.dt*v;
     *
     *      Check that convention against the intended discretised equations.
     */

    // Remove these casts when the corresponding arguments are used.
    (void)previous_fields;
    (void)current_fields;
    (void)global_x;
    (void)global_y;
    (void)global_z;
    (void)geometry;
    (void)num_components;
}

void UserDefinedProfile::setGaugeProfile(
    float *previous_fields,
    float *current_fields,
    const unsigned global_x,
    const unsigned global_y,
    const unsigned global_z,
    const InitialConditionGeometry &geometry,
    const unsigned num_components) const
{
    /*
     * USER SECTION: GAUGE-LINK INITIAL CONDITIONS
     *
     * The component array contains the x-, y- and z-directed links in three
     * consecutive blocks. Their common width is
     *
     *      const unsigned direction_width = num_components/3U;
     *
     * so
     *
     *      current_fields[dir*direction_width + component]
     *
     * selects direction dir = 0, 1, 2.
     *
     * The component meaning inside each block is determined by the chosen
     * Wilson-loop model and gauge representation.
     *
     * Generator representation:
     *     zero components represent trivial links.
     *
     * Quaternion representation:
     *     every SU(2) factor must be assigned as (c0,c1,c2,c3), satisfying
     *
     *     c0*c0 + c1*c1 + c2*c2 + c3*c3 = 1.
     *
     *     The arrays entering this blank function already contain identity
     *     links in the quaternion slots. Do not replace those four components
     *     by all zeros.
     *
     * After assigning current_fields, copy the complete representation into
     * previous_fields for zero initial electric field.
     *
     * Non-zero electric initial data are group- and representation-dependent
     * and should not be introduced by a naive component-wise difference.
     */

    // Remove these casts when the corresponding arguments are used.
    (void)previous_fields;
    (void)current_fields;
    (void)global_x;
    (void)global_y;
    (void)global_z;
    (void)geometry;
    (void)num_components;
}


////////////////////////////////////////////  Public Functions  //////////////////////////////////////////////////

void UserDefinedProfile::setInitialFields(
    std::vector<float> &field,
    const InitialConditionGeometry &geometry,
    const unsigned num_components) const
{
    if (num_components == 0U)
    {
        return;
    }

    if (geometry.ownedXBegin > geometry.ownedXEnd
        || geometry.ownedXEnd > geometry.storageNx
        || geometry.ownedXEnd - geometry.ownedXBegin
            != geometry.localNx
        || geometry.globalXStart + geometry.localNx
            > geometry.globalNx)
    {
        throw std::runtime_error(
            "INITIALCONDITIONS::USERDEFINEDPROFILE:: "
            "Invalid IC geometry."
        );
    }

    if (this->fieldType == FieldType::Gauge
        && num_components%3U != 0U)
    {
        throw std::runtime_error(
            "INITIALCONDITIONS::USERDEFINEDPROFILE:: "
            "Gauge components must divide into three spatial directions."
        );
    }

    const std::size_t local_site_count = static_cast<std::size_t>(geometry.storageNx)*geometry.globalNy*geometry.globalNz;

    const std::size_t buffer_size = local_site_count*num_components;

    if (field.size() != 2ULL*buffer_size)
    {
        throw std::runtime_error(
            "INITIALCONDITIONS::USERDEFINEDPROFILE:: "
            "Field size does not match IC geometry."
        );
    }

    for (unsigned local_x = geometry.ownedXBegin; local_x < geometry.ownedXEnd; local_x++)
    {
        const unsigned global_x = geometry.localToGlobalX(local_x);

        for (unsigned global_y = 0U; global_y < geometry.globalNy; global_y++)
        {
            for (unsigned global_z = 0U; global_z < geometry.globalNz; global_z++)
            {
                const std::size_t local_site_index
                    = (static_cast<std::size_t>(local_x)*geometry.globalNy + global_y)*geometry.globalNz + global_z;

                const std::size_t field_index = local_site_index*num_components;

                // Buffer 0: previous physical timestep.
                float *previous_fields = field.data() + field_index;

                // Buffer 1: current initial configuration.
                float *current_fields = field.data() + buffer_size + field_index;

                if (this->fieldType == FieldType::Scalar)
                {
                    this->setScalarProfile(
                        previous_fields,
                        current_fields,
                        global_x,
                        global_y,
                        global_z,
                        geometry,
                        num_components
                    );
                }
                else
                {
                    this->setGaugeProfile(
                        previous_fields,
                        current_fields,
                        global_x,
                        global_y,
                        global_z,
                        geometry,
                        num_components
                    );
                }
            }
        }
    }
}