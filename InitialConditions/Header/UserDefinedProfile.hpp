#pragma once

#include "InitialCondition.hpp"

#include <string>


/*
 * Template for position-dependent initial conditions.
 * The generic lattice indexing and MPI-slab mapping are implemented here;
 * users only need to edit the two profile functions in the source file.
 */
class UserDefinedProfile:
    public InitialCondition
{
public:

    enum class FieldType
    {
        Scalar,
        Gauge
    };

private:

    ////////////////////////////////////////////////  Variables  /////////////////////////////////////////////////////

    FieldType fieldType;
    bool usingQuaternionRepresentation;

    // Optional filenames read from UserDefinedProfile.cfg. Add further
    // profile parameters or loaded data containers here as required.
    std::string scalarProfileDataFilename;
    std::string gaugeProfileDataFilename;

    ///////////////////////////////////////////////  Initialisers  ///////////////////////////////////////////////////

    /*
     * Loads the user-defined profile settings from a config file.
     * Small tabulated profiles should be loaded once here, not inside the
     * site-by-site profile functions below.
     */
    void configure(
        const std::string path,
        const bool debug = false
    );

    /////////////////////////////////////////////  Profile Functions  ////////////////////////////////////////////////

    /*
     * Student-editable scalar-profile function. The two pointers address the
     * complete set of scalar components at one global physical lattice site.
     */
    void setScalarProfile(
        float *previous_fields,
        float *current_fields,
        const unsigned global_x,
        const unsigned global_y,
        const unsigned global_z,
        const InitialConditionGeometry &geometry,
        const unsigned num_components
    ) const;

    /*
     * Student-editable gauge-profile function. The two pointers address all
     * stored gauge-link components at one global physical lattice site.
     */
    void setGaugeProfile(
        float *previous_fields,
        float *current_fields,
        const unsigned global_x,
        const unsigned global_y,
        const unsigned global_z,
        const InitialConditionGeometry &geometry,
        const unsigned num_components
    ) const;

public:

    //////////////////////////////////////////  Constructors/Destructors  ////////////////////////////////////////////

    UserDefinedProfile(
        const FieldType field_type,
        const bool using_quaternion_representation = false
    );

    ~UserDefinedProfile() override;

    /////////////////////////////////////////////  Public Functions  /////////////////////////////////////////////////

    /*
     * Loops over this MPI rank's owned physical sites and dispatches to the
     * appropriate student-editable profile function.
     */
    void setInitialFields(
        std::vector<float> &field,
        const InitialConditionGeometry &geometry,
        const unsigned num_components
    ) const override;
};