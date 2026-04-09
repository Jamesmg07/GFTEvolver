#pragma once

#include "Gradient.hpp"

#include <cmath>
#include <iomanip>

// Handles all details relating to the covariant derivative, as used in the 2HDM.

namespace Gradients {

    class TwoHDM:
        public Gradient
    {
    private:

        /////////////////////////////////////////////  Variables  ///////////////////////////////////////////////

        const unsigned &numScalarComponents, &numVectorComponents;
        const double &dt, &dx, &dy, &dz;
        
        const bool usingGeneratorRepresentation;
        unsigned gaugeNum; // 4 for SM, 5 for SM x U(1) and 7 for SM x SU(2)

        double divisionByZeroTolerance;

        ////////////////////////////////////////////  Initialisers  /////////////////////////////////////////////

        /*
        * Loads in the parameters relating to the 2HDM covariant derivatives.
        * 
        * @param        string path                Path to the config file.
        * @param        bool debug                 Outputs loaded parameters if true.
        */
        void configure(const std::string path, const bool debug = false);

        /////////////////////////////////////////  Private Functions  ////////////////////////////////////////////

        /*
        * Takes in the scalar field and transforms it under the action of a given link variable.
        * The first component of the vector_field is the hypercharge phase, y_i in e^{iy_i}.
        * The next three are the isospin, w_i^a in e^{iw_i^a\sigma^a}.
        * Next can either be one or three based on whether the Higgs family has a gauged U(1) or SU(2).
        * There are three sets of these, corresponding to the three spatial directions.
        * 
        * @param        vector<double> &scalar_fields                Scalar fields at a given position.
        * @param        float* vector_pointer                        Pointer to the vector field at a given position.
        * @param        unsigned dir_index                           Selects which spatial component of the vector field to use
        * @param        bool conjugate                               Conjugates the transformation matrix if true.
        * 
        * @return       vector<double>                               The scalar fields after the transformation.
        */
        std::vector<double> transform(const std::vector<double> &scalar_fields, const float* const vector_pointer, 
                                      const unsigned dir_index, const bool conjugate) const;

        /*
         * Returns the SU(2) matrix pointed to, in a different representation.
         * The input representation is assumed to be in the form e^{iw_i^a\sigma^a}.
         * The output representation is (c^0, c^a), with U = c^0\sigma^0 + ic^a\sigma^a.
         * 
         * @param        float* vector_pointer                Pointer to the vector fields at a specific grid position.
         * @param        unsigned dir_index                   Selects which spatial component of the vector field to use.
         * @param        unsigned offset                      Specifies where the SU(2) information resides in the array
         * @param        bool conjugate                       Conjugates the matrix if true.
         * 
         * @return       vector<double>                        The coefficents of the SU(2) representation described above.
         */
        std::vector<double> getSU2Representation(const float* const vector_pointer, const unsigned dir_index, const unsigned offset, 
                                                 const bool conjugate) const;

        /*
        * Simply returns a vector containing all components of the scalar field, that are pointed to by the pointer.
        *
        * @param        float* scalar_pointer                Pointer to the scalar field at a given position.
        * 
        * @return       vector<double>                       The scalar field at this position.
        */
        std::vector<double> getScalarField(const float* const scalar_pointer) const;

    public:

        //////////////////////////////////////  Constructors/Destructors  ////////////////////////////////////////

        TwoHDM(
            const unsigned &num_scalar_components, const unsigned &num_vector_components,
            const unsigned &nx, const unsigned &ny, const unsigned &nz,
            const double &dt, const double &dx, const double &dy, const double &dz,
            const bool using_generator_representation
        );
        virtual ~TwoHDM();


        /////////////////////////////////////////  Public Functions  ////////////////////////////////////////////

        /*
        * Get the size of the default stencil (assumed the same in each direction).
        * It is assumed that the stencil is symmetric, so ,for example, 1 means that the stencil looks at neighbours
        * that are one lattice site away on both sides.
        * 
        * @return    unsigned                The stencil size.
        */
        unsigned getDefaultStencilSize() const;

        /*
        * Get the number of constraint equations (due to the fixing of the temporal gauge) that should be satisfied throughout.
        * 
        * @return        unsigned                                             The number of constraint equations.
        */
        unsigned getNumberOfConstraintEquations() const;

        /*
        * Calculates the gradient energy.
        * 
        * @param        vector<vector<float*>> &scalar_pointers                 Array of pointers to the scalar fields at the required grid positions.
        * @param        vector<vector<float*>> &vector_pointers                 Array of pointers to the vector fields at the required grid positions.
        * 
        * @return       float                                                   The gradient energy (density) at this position.
        */
        float calcGradientEnergy(const std::vector<std::vector<const float*>> &scalar_pointers, 
                                const std::vector<std::vector<const float*>> &vector_pointers) const;

        /*
        * Calculate the kinetic energy.
        * 
        * @param        float* local_scalar_fields[2]                           Pointers to scalar fields at this location, for both timesteps.
        * 
        * @return       float                                                   The kinetic energy (density) at this position.
        */
        float calcKineticEnergy(const float* const local_scalar_fields[2]) const;

        /*
        * Calculates the contributions from the gradients to the constraint equations associated with the temporal gauge fixing.
        * 
        * @param        float* local_scalar_fields[2]                           Pointers to scalar fields at this location, for both timesteps.
        * 
        * @return       vector<float>                                           Contribution to the constraint equations.
        */
        std::vector<float> calcConstraintContributions(const float* const local_scalar_fields[2]) const;

        /*
        * Calculate the contribution from the gradient energy term to the scalar equations of motion.
        * 
        * @param        vector<vector<float*>> &scalar_pointers                 Array of pointers to the scalar fields at the required grid positions.
        * @param        vector<vector<float*>> &vector_pointers                 Array of pointers to the vector fields at the required grid positions.
        * 
        * @return       vector<double>                                          Contribution to the scalar equations of motion.
        */
        std::vector<double> calcDerivatives(const std::vector<std::vector<const float*>> &scalar_pointers,
                                            const std::vector<std::vector<const float*>> &vector_pointers) const;

        /*
        * Calculate the contribution from the gradient energy term to the vector (gauge) equations of motion.
        * 
        * @param        vector<vector<float*>> &scalar_pointers                 Array of pointers to the scalar fields at the required grid positions.
        * @param        vector<vector<float*>> &vector_pointers                 Array of pointers to the vector fields at the required grid positions.
        * 
        * @return       vector<double>                                          Contribution to the vector equations of motion.
        */
        std::vector<double> calcCurrents(const std::vector<std::vector<const float*>> &scalar_pointers,
                                         const std::vector<std::vector<const float*>> &vector_pointers) const;

    };

}