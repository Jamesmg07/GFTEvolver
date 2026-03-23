#pragma once

#include "WilsonLoop.hpp"

#include <cmath>
#include <iomanip>

// Handles all details related to the field strength-tensor, as used in the Standard Model.

namespace WilsonLoops{

    class StandardModel:
        public WilsonLoop
    {
    private:

        ////////////////////////////////////////////////  Variables  ///////////////////////////////////////////////////

        double divisionByZeroTolerance;
        const double &dt, &dx, &dy, &dz;
        double inverse_sqr_spacings[3], inverse_sqr_dt;

        double g_sqr, gp_sqr, inverse_g_sqr, inverse_gp_sqr;

        bool storeEnergy;
        mutable float storedMagneticEnergy, storedElectricEnergy; 

        ///////////////////////////////////////////////  Initialisers  /////////////////////////////////////////////////

        /*
        * Loads in the from the associated configuration file.
        * 
        * @param        string path             The file path to the config file.
        * @param        bool debug              Outputs the parameters read in from the config file if true.
        */
        void configure(const std::string path, const bool debug = false);

        /*
         * Initialises variables
         */
        void initVariables();


        //////////////////////////////////////////////  Private Functions  /////////////////////////////////////////////

        /*
         * Returns the SU(2) matrix pointed to, in a different representation.
         * The input representation is assumed to be in the form e^{iw_i^a\sigma^a}.
         * The output representation is (c^0, c^a), with U = c^0\sigma^0 + ic^a\sigma^a.
         * 
         * @param        float* vector_pointer                Pointer to the vector fields at a specific grid position.
         * @param        unsigned dir_index                   Selects which spatial component of the vector field to use.
         * @param        bool conjugate                       Conjugates the matrix if true.
         * 
         * @return       vector<double>                        The coefficents of the SU(2) representation described above.
         */
        std::vector<double> getSU2Representation(const float* const vector_pointer, const unsigned dir_index, const bool conjugate) const;

        /*
         * Does the opposite of above. Takes in a representation of SU(2) like (c^0, c^a) and returns w_i^a.
         * The SU(2) matrix is U = c^0\sigma^0 +_ ic^a\sigma^a = e^{iw_i^a\sigma^a}.
         * 
         * @param        vector<double> U_representation              Array of c^\mu components (4).
         * 
         * @return       vector<float>                               Array of w_i^a components (3).
         */
        std::vector<float> invertSU2Representation(const std::vector<double> U_representation) const;

        /*
         * Returns the matrix product of two SU(2) matrices. The inputs are assumed to have size 4, with the components
         * corresponding to an SU(2) matrix through: U = c^0\sigma^0 + ic^a\sigma^a.
         * In order for this to be an SU(2) matrix, they need to satisfy c^0c^0 + c^ac^a = 1, but this is not checked.
         * 
         * @param        vector<float> U1                A representation of the first SU(2) matrix, as described above.
         * @param        vector<float> U2                Same as above, but for the second SU(2) matrix.
         * @param        bool calc_trace                 If false, only bother calculating the traceless part of the result.
         * 
         * @return       vector<float>                   Returns a representation of the resulting SU(2).
         */
        std::vector<double> SU2Product(const std::vector<double> U1, const std::vector<double> U2, const bool calc_trace) const;

    public:

        //////////////////////////////////////////  Constructors/Destructors  /////////////////////////////////////////

        StandardModel(const double &dt, const double &dx, const double &dy, const double &dz);
        virtual ~StandardModel();

        //////////////////////////////////////////////  Public Functions  /////////////////////////////////////////////

        /*
        * Returns the relevant squares of the gauge couplings, for each component.
        * 
        * @param        unsigned comp_iter                                     Integer that decides which component of the gauge field.
        * 
        * @return       float                                                  Square of the relevant gauge coupling.
        */
        float getSqrCouplings(const unsigned comp_iter) const;

        /*
        * Set up for the energy calculation. Only runs if energy analyser is being used.
        *
        * @param        bool store_energy                                      Boolean that decides if energy calculations are performed during evolution.
        */
        void energyPreparation(const bool store_energy);

        /*
        * Calculates the energy density coming from the spatial part of both Yang-Mills terms.
        * 
        * @param        vector<vector<float*>> &vector_pointers                Array of pointers to the vector fields at the required grid positions.
        * 
        * @return       float                                                  Energy density at this position from the Yang-Mills term.
        */
        float calcMagneticEnergy(const std::vector<std::vector<const float*>> &vector_pointers) const;

        /*
        * Calculates the energy density coming the temporal part of both Yang-Mills terms.
        * 
        * @param       float* local_vector_fields[2]                           Pointers to the vector fields at this location, for both timesteps.
        * 
        * @return      float                                                   Energy density at this position from the (temporal part of the) Yang-Mills term.
        */
        float calcElectricEnergy(const float* const local_vector_fields[2]) const;

        /*
        * Calculates the contribution to the equations of motion coming from the spatial part of both Yang-Mills terms.
        * The contribution from the time components will be treated separately, in the evolve function.
        * 
        * @param        vector<vector<float*>> &vector_pointers                Array of pointers to the vector fields at the required grid positions.
        * 
        * @return       vector<float>                                          Contribution to the vector equations of motion.
        */
        std::vector<double> calcMagneticContributions(const std::vector<std::vector<const float*>> &vector_pointers) const;

        /*
        * Calculates the contribution to the equations of motion coming from half of the temporal part of the Yang-Mills term,
        * the half that depends upon the current timestep and the previous one.
        * The other half will be used to evolve the fields and dealt with separately.
        * 
        * @param        float* local_vector_fields[2]                          Pointers to the vector fields at this location, for both timesteps.
        * 
        * @return       vector<float>                                          Contribution to the vector equations of motion.
        */
        std::vector<double> calcElectricContributions(const float* const local_vector_fields[2]) const;


        /*
        * Determines the group element at the next timestep, given that the imaginary part of the trace of \sigma^a U_i(t+dt)U_i^\dagger(t) = RHS_i^a.
        * It is assumed that, for sufficiently small dt, U_i(t+dt)U_i^\dagger(t) is close to the identity which resolves the remaining ambiguity.
        * Finally, multiply on the right by U_i(t) to get the group element at the next timestep.
        * 
        * @param        float* local_vector_fields[2]                          Pointers to the vector fields at this location, for both timesteps.
        * @param        vector<double> equation_RHS                             Array containing the right-hand side of the equation, for all components.
        */
        void evolve(float* const local_vector_fields[2], std::vector<double> equation_RHS);

    };

}