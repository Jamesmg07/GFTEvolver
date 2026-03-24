#pragma once

#include "Gradient.hpp"

#include <cmath>

namespace Gradients{

    class Global:
        public Gradient
    {
    private:

        ////////////////////////////////////////////////  Variables  //////////////////////////////////////////////////////////

        const unsigned &numScalarComponents, &numVectorComponents;
        const double &dt, &dx, &dy, &dz;

        ///////////////////////////////////////////////  Initialisers  ////////////////////////////////////////////////////////

        /*
        * Load in any parameters relating to the global derivatives from the associated config file.
        * Currently there are none for this case, but an empty function and config file are still provided
        * in case they are needed in the future.
        * 
        * @param        string path                Path to the config file.
        * @param        bool debug                 Outputs loaded parameters if true.
        */
        void configure(const std::string path, const bool debug = 0);

    public:

        //////////////////////////////////////////  Constructors/Destructors  /////////////////////////////////////////////////

        Global(
            const unsigned &num_scalar_components, const unsigned &num_vector_components,
            const unsigned &nx, const unsigned &ny, const unsigned &nz,
            const double &dt, const double &dx, const double &dy, const double &dz
        );
        virtual ~Global();

        /////////////////////////////////////////////  Public Functions  //////////////////////////////////////////////////////

        /*
        * Get the size of the default stencil (assumed the same in each direction).
        * It is assumed that the stencil is symmetric, so ,for example, 1 means that the stencil looks at neighbours
        * that are one lattice site away on both sides.
        * 
        * @return    unsigned                The stencil size.
        */
        unsigned getDefaultStencilSize() const;

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
        * Calculates the kinetic energy.
        * 
        * @param        float* local_scalar_fields[2]                           Pointers to scalar fields at this location, for both timesteps.
        * 
        * @return       float                                                   The kinetic energy (density) at this position.
        */
        float calcKineticEnergy(const float* const local_scalar_fields[2]) const;

        /*
        * This function is not required because there are no constraint equations associated with the gauge fixing in the global case.
        * It still takes in parameters to match the pure virtual function that may be called.
        * 
        * @param        float* local_scalar_fields[2]                           Pointers to scalar fields at this location, for both timesteps.
        * 
        * @return       vector<float>                                           Contribution to the constraint equations.
        */
        std::vector<float> calcConstraintContributions(const float* const local_scalar_fields[2]) const;

        /*
        * Calculates the contribution to the scalar equations of motion from the gradient energy term.
        * For global derivatives, this is the only contribution.
        * 
        * @param        vector<vector<float*>> &scalar_pointers                 Array of pointers to the scalar fields at the required grid positions.
        * @param        vector<vector<float*>> &vector_pointers                 Array of pointers to the vector fields at the required grid positions.
        * 
        * @return       vector<double>                                          Contribution to the scalar equations of motion.
        */
        std::vector<double> calcDerivatives(const std::vector<std::vector<const float*>> &scalar_pointers,
                                            const std::vector<std::vector<const float*>> &vector_pointers) const;

        /*
        * There is no gauge field equations of motion, so this function does nothing because it is not required.
        * It still takes in parameters to match the pure virtual function that may be called.
        * 
        * @param        vector<vector<float*>> &scalar_pointers                 Array of pointers to the scalar fields at the required grid positions.
        * @param        vector<vector<float*>> &vector_pointers                 Array of pointers to the vector fields at the required grid positions.
        * 
        * @return       vector<double>                                          Empty vector.
        */
        std::vector<double> calcCurrents(const std::vector<std::vector<const float*>> &scalar_pointers,
                                         const std::vector<std::vector<const float*>> &vector_pointers) const;

    };

}