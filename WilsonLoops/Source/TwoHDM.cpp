#include "TwoHDMLoops.hpp"

/////////////////////////////////////////////  Initialisers  ///////////////////////////////////////////////

void WilsonLoops::TwoHDM::configure(const std::string path, const bool debug)
{

    // Load in from file
    std::ifstream ifs(path);

    // Read in the parameter values from the configuration file
    std::string description; // Descriptions in the config file (to be dumped)
    double g, gp, gpp;
    if (ifs.is_open())
    {
        // Skip the initial configuration file description
        for (unsigned iter = 0; iter < 12; iter++)
            std::getline(ifs, description);

        // Gauge couplings
        std::getline(ifs, description, ':');
        ifs >> g;
        this->g_sqr = g*g;

        std::getline(ifs, description, ':');
        ifs >> gp;
        this->gp_sqr = gp*gp;

        std::getline(ifs, description, ':');
        ifs >> gpp;
        this->gpp_sqr = gpp*gpp;

        std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> this->usingGeneratorRepresentation;

        std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> this->divisionByZeroTolerance;

    } 

    ifs.close();

    if (debug)
    {
        std::cout << "WILSONLOOPS::TWOHDM::\n"
                  << "Isospin gauge coupling: " << g << ", Hypercharge gauge coupling: " << gp << ", Higgs family gauge coupling: " << gpp << "\n"
                  << "Is the generator representation being used?: " << this->usingGeneratorRepresentation << "\n"
                  << "Division by zero tolerance: " << this->divisionByZeroTolerance
                  << "\n" << std::endl;
    }
}


void WilsonLoops::TwoHDM::initVariables()
{
    // Defining convenient array of spacings that can be indexed.

    this->inverse_sqr_spacings[0] = 1.f/(this->dx*this->dx);
    this->inverse_sqr_spacings[1] = 1.f/(this->dy*this->dy);
    this->inverse_sqr_spacings[2] = 1.f/(this->dz*this->dz);

    this->inverse_sqr_dt = 1.f/(this->dt*this->dt);

    // Pre-calculate inverses of the gauge couplings (only affects the energy)
    // They will be set to zero if the squared gauge couplings are smaller than the division by zero tolerance.
    if (this->g_sqr > this->divisionByZeroTolerance)
        this->inverse_g_sqr = 1.0/this->g_sqr;
    else
        this->inverse_g_sqr = 0.0;

    if (this->gp_sqr > this->divisionByZeroTolerance)
        this->inverse_gp_sqr = 1.0/this->gp_sqr;
    else
        this->inverse_gp_sqr = 0.0;

    if (this->gpp_sqr > this->divisionByZeroTolerance)
        this->inverse_gpp_sqr = 1.0/this->gpp_sqr;
    else
        this->inverse_gpp_sqr = 0.0;

    // Set boolean output for the energy calculations to false.
    this->storeEnergy = false;

}

/////////////////////////////////////////////  Private Functions  //////////////////////////////////////////

WilsonLoops::TwoHDM::SU2Representation
WilsonLoops::TwoHDM::getSU2Representation(
    const float *const vector_pointer,
    const long long int dir_index,
    const unsigned offset,
    const bool conjugate) const
{
    SU2Representation c_representation{};

    int conj_fac = 1;
    if (conjugate)
        conj_fac = -1;

    if (this->usingGeneratorRepresentation)    // Approach when 3 dofs are stored
    {
        double w_mag = 0.0;
        for (unsigned comp_iter = offset; comp_iter < 3 + offset; comp_iter++)
            w_mag += std::pow(static_cast<double>(vector_pointer[dir_index + comp_iter]), 2);

        w_mag = std::sqrt(w_mag);
        c_representation[0] = std::cos(w_mag);

        if (w_mag > this->divisionByZeroTolerance)
        {
            for (unsigned comp_iter = 1; comp_iter < 4; comp_iter++)
                c_representation[comp_iter] = conj_fac*std::sin(w_mag)*static_cast<double>(vector_pointer[dir_index + comp_iter + offset - 1])/w_mag;
        }
        else
        {
            c_representation[1] = std::sin(w_mag); // Arbitrarily choose this direction.
        }
    }
    else                                       // Approach when 4 dofs are stored
    {
        c_representation[0] = static_cast<double>(vector_pointer[dir_index + offset]);
        for (int comp_iter = 1; comp_iter < 4; comp_iter++) 
            c_representation[comp_iter] = static_cast<double>(conj_fac*vector_pointer[dir_index + comp_iter + offset]);
    }

    return c_representation;
}


std::vector<float> WilsonLoops::TwoHDM::invertSU2Representation(const WilsonLoops::TwoHDM::SU2Representation &c_representation) const
{
    std::vector<float> w_representation(3, 0.f);

    //float phase = std::acos(c_representation[0]); // Returns results in the range 0 to pi, so that sin(phase)>=0.

    double c_mag = 0;
    for (unsigned comp_iter = 1; comp_iter < 4; comp_iter++)
        c_mag += std::pow(c_representation[comp_iter], 2);

    c_mag = std::sqrt(c_mag);

    double phase;
    if (c_mag > 1)
        //std::cout << "Panic in invert (c_mag = " << c_mag << " )" << " " << c_representation[1] << " " << c_representation[2] << " " << c_representation[3] << std::endl;

        // Should just be ever so slightly > 1 from rounding errors

        phase = M_PI_2;

    else

        phase = std::asin(c_mag);    // Will return result in range 0 to pi/2 -> sin(phase)>=0 but need to check what sign I want cos(phase) to have.

        if (c_representation[0] < 0)
            phase = M_PI - phase;

    // double phase = std::asin(c_mag);
    // if (c_representation[0] < 0)
    //     phase = M_PI - phase;

    if (c_mag > this->divisionByZeroTolerance)
    {
        for (unsigned comp_iter = 0; comp_iter < 3; comp_iter++)
            w_representation[comp_iter] = static_cast<float>(phase*c_representation[comp_iter + 1]/c_mag);
    }
    else
    {
        w_representation[0] = static_cast<float>(phase); // Arbitrarily choose this direction.
    }

    return w_representation;
}


WilsonLoops::TwoHDM::SU2Representation
WilsonLoops::TwoHDM::SU2Product(
    const WilsonLoops::TwoHDM::SU2Representation &U1,
    const WilsonLoops::TwoHDM::SU2Representation &U2,
    const bool calc_trace) const
{
    SU2Representation U_product{};

    // Remaining three are given by c^0_1c^a_2 + c^a_1c^0_2 - epsilon^{abc}c^b_1c^c_2
    U_product[1] = U1[0]*U2[1] + U1[1]*U2[0] - U1[2]*U2[3] + U1[3]*U2[2];
    U_product[2] = U1[0]*U2[2] + U1[2]*U2[0] - U1[3]*U2[1] + U1[1]*U2[3];
    U_product[3] = U1[0]*U2[3] + U1[3]*U2[0] - U1[1]*U2[2] + U1[2]*U2[1];

    if (calc_trace)
    {
        //First component is given by c^0_1c^0_2 - c^a_1c^a_2
        double dot_product = 0.0;
        for (unsigned comp_iter = 1; comp_iter < 4; comp_iter++)
            dot_product += U1[comp_iter]*U2[comp_iter];

        U_product[0] = U1[0]*U2[0] - dot_product;

        //SU2 must satisfy U[0]^2 + U[a]^2 = 1. Renormalise this to avoid slowly drifting away from SU(2) matrices.
        double norm = std::sqrt( std::pow(U_product[0], 2) + std::pow(U_product[1], 2) + std::pow(U_product[2], 2) + std::pow(U_product[3], 2) );

        for (unsigned comp_iter = 0; comp_iter < 4; comp_iter++)
            U_product[comp_iter] = U_product[comp_iter]/norm;

    }

    // if (1.0 - U_product[0] < 0.0)
    //     std::cout << std::setprecision(12) << "Panic in SU2Product: " << U_product[0] << std::endl;

    return U_product;

}

/////////////////////////////////////////  Constructors/Destructors  ///////////////////////////////////////

WilsonLoops::TwoHDM::TwoHDM(const unsigned &num_vector_components, const double &dt, const double &dx, const double &dy, const double &dz)
    : numVectorComponents(num_vector_components), dt(dt), dx(dx), dy(dy), dz(dz)
{
    this->configure(std::string(SOURCE_DIR)+"/Config/TwoHDM.cfg", true);
    this->initVariables();

    if (this->usingGeneratorRepresentation)
    {
        if (this->numVectorComponents == 12)
            this->gaugeNum = 4U; // Just SM
        else if (this->numVectorComponents == 15)
            this->gaugeNum = 5U; // SM x U(1)
        else if (this->numVectorComponents == 21)
            this->gaugeNum = 7U; // SM x SU(2)
        else
            throw std::runtime_error("WILSONLOOPS::TWOHDM:: When using the generator representation, either 12 (just SM), 15 (SM x U(1)) "
                                     "or 21 (SM x SU(2)) vector components are required but " + std::to_string(this->numVectorComponents) +
                                     " have been assigned.\n"
                                     "Note that 3 spatial dimensions are assumed so requested number of components is multiplied by 3 internally.");
    }
    else
    {
        if (this->numVectorComponents == 15)
            this->gaugeNum = 4U; // Just SM
        else if (this->numVectorComponents == 18)
            this->gaugeNum = 5U; // SM x U(1)
        else if (this->numVectorComponents == 27) 
            this->gaugeNum = 7U; // SM x SU(2)
        else
            throw std::runtime_error("WILSONLOOPS::TWOHDM:: When using the quaternion representation, either 15 (just SM), 18 (SM x U(1)) "
                                     "or 27 (SM x SU(2)) vector components are required but " + std::to_string(this->numVectorComponents) +
                                     " have been assigned.\n"
                                     "Note that 3 spatial dimensions are assumed so requested number of components is multiplied by 3 internally.");
    }

    std::cout << "WILSONLOOPS::TWOHDM:: " << this->gaugeNum << std::endl;
}


WilsonLoops::TwoHDM::~TwoHDM()
{
}

/////////////////////////////////////////////  Public Functions  //////////////////////////////////////////

unsigned WilsonLoops::TwoHDM::getDefaultStencilSize() const
{
    return 1;
}

bool WilsonLoops::TwoHDM::isUsingGeneratorRepresentation() const
{
    return this->usingGeneratorRepresentation;
}

float WilsonLoops::TwoHDM::getSqrCouplings(const unsigned comp_iter) const
{
    unsigned iter_modulo = comp_iter%this->gaugeNum;
    if (iter_modulo == 0)
        return this->gp_sqr;
    else if (iter_modulo/4 == 0)
        return this->g_sqr;
    else
        return this->gpp_sqr;
}

// Can combine the two following functions into one and move the multiplication by #dims outside.
unsigned WilsonLoops::TwoHDM::getNumberOfEvolutionEquations() const
{
    return this->gaugeNum*3U;
}

unsigned WilsonLoops::TwoHDM::getNumberOfConstraintEquations() const
{
    return this->gaugeNum;
}

void WilsonLoops::TwoHDM::energyPreparation(const bool store_energy)
{
    this->storeEnergy = store_energy;
    this->storedMagneticEnergy = 0.f;
    this->storedElectricEnergy = 0.f;
}


float WilsonLoops::TwoHDM::calcMagneticEnergy(const std::vector<std::vector<const float*>> &vector_pointers) const
{
    // Calculate real part of spatial plaquette action for the hypercharge and trace for the isospin.
    // This will have already been done within the calcMagneticContributions function and stored.
    // This function now just needs to be return the value and reset the stored energy to zero for the next sum.

    float magnetic_energy_density = this->storedMagneticEnergy;
    this->storedMagneticEnergy = 0.f;

    return magnetic_energy_density;
}


float WilsonLoops::TwoHDM::calcElectricEnergy(const float *const local_vector_fields[2]) const
{
    // Same as above but for the temporal loops. Will have already done some calculations within calcElectricContributions.
    // The only difference compared to above is that the multiplication by 1/dt^2 is common here, so has been moved outside into this function.
    float electric_energy_density = this->inverse_sqr_dt*this->storedElectricEnergy;
    this->storedElectricEnergy = 0.f;

    return electric_energy_density;
}

std::vector<float> WilsonLoops::TwoHDM::calcConstraintContributions(const long long int &t_future_index, 
                                                                           const std::vector<std::vector<const float *>> &vector_pointers) const
{
    // WARNING:: THIS CURRENTLY ASSUMES THAT THE STENCIL BEING USED IS ALWAYS 3-POINT. WILL NEED TO MAKE SOME ALTERATIONS TO GET IT TO WORK FOR
    // OTHER STENCILS!

    std::vector<float> contribution(this->gaugeNum, 0.f);

    for (unsigned dir_iter = 0; dir_iter < 3; dir_iter++)
    {

        int dir_index = dir_iter*this->numVectorComponents/3;
        long long int future_index = t_future_index + dir_index;

        //////////////////  Hypercharge U(1) calculations:  ////////////////////////

        double local_loop_angle = static_cast<double>(vector_pointers[0][1][future_index]) 
                                - static_cast<double>(vector_pointers[0][1][dir_index]);

        double neighbour_loop_angle = static_cast<double>(vector_pointers[dir_iter][0][t_future_index + dir_index])
                                    - static_cast<double>(vector_pointers[dir_iter][0][dir_index]);

        contribution[0] += ( std::sin( local_loop_angle ) - std::sin( neighbour_loop_angle ) )*this->inverse_sqr_spacings[dir_iter];


        //////////////////  Isospin SU(2) calculations:  ///////////////////////////

        // Calculate the contribution from the first Wilson loop:
        // Representation of right-hand matrix
        SU2Representation U_product1 = this->getSU2Representation(vector_pointers[0][1], future_index, 1U, false);

        // Representation of the left-hand matrix
        SU2Representation U_multiply = this->getSU2Representation(vector_pointers[0][1], dir_index, 1U, true);

        U_product1 = this->SU2Product(U_multiply, U_product1, false);


        // Repeat for the second Wilson loop:
        SU2Representation U_product2 = this->getSU2Representation(vector_pointers[dir_iter][0], dir_index, 1U, true);
        U_multiply = this->getSU2Representation(vector_pointers[dir_iter][0], future_index, 1U, false);

        U_product2 = this->SU2Product(U_multiply, U_product2, false);

        for (unsigned eq_iter = 1; eq_iter < 4; eq_iter++)
            contribution[eq_iter] += static_cast<float>( ( U_product1[eq_iter] - U_product2[eq_iter] )*this->inverse_sqr_spacings[dir_iter] );


        ///////////////////  Higgs family calculations:  ////////////////////////////

        unsigned offset;
        if (this->gaugeNum > 4U)
        {
            if (this->usingGeneratorRepresentation)
                offset = 4U;
            else
                offset = 5U;
        }

        if (this->gaugeNum == 5U)
        {
            // Higgs family U(1)

            local_loop_angle = static_cast<double>(vector_pointers[0][1][future_index + offset]) 
                             - static_cast<double>(vector_pointers[0][1][dir_index + offset]);

            neighbour_loop_angle = static_cast<double>(vector_pointers[dir_iter][0][t_future_index + dir_index + offset])
                                 - static_cast<double>(vector_pointers[dir_iter][0][dir_index + offset]);

            contribution[4] += ( std::sin( local_loop_angle ) - std::sin( neighbour_loop_angle ) )*this->inverse_sqr_spacings[dir_iter];
        }
        else if (this->gaugeNum == 7U)
        {
            // Higgs family SU(2)

            // Calculate the contribution from the first Wilson loop:
            U_product1 = this->getSU2Representation(vector_pointers[0][1], future_index, offset, false);
            U_multiply = this->getSU2Representation(vector_pointers[0][1], dir_index, offset, true);

            U_product1 = this->SU2Product(U_multiply, U_product1, false);

            // Repeat for the second Wilson loop:
            U_product2 = this->getSU2Representation(vector_pointers[dir_iter][0], dir_index, offset, true);
            U_multiply = this->getSU2Representation(vector_pointers[dir_iter][0], future_index, offset, false);

            U_product2 = this->SU2Product(U_multiply, U_product2, false);

            for (unsigned eq_iter = 4; eq_iter < 7; eq_iter++)
                contribution[eq_iter] += static_cast<float>( ( U_product1[eq_iter - 3U] - U_product2[eq_iter - 3U] )*this->inverse_sqr_spacings[dir_iter] );
        }

    }

    return contribution;
}

std::vector<double> WilsonLoops::TwoHDM::calcMagneticContributions(const std::vector<std::vector<const float *>> &vector_pointers) const
{
    // WARNING:: THIS CURRENTLY ASSUMES THAT THE STENCIL BEING USED IS ALWAYS 3-POINT. WILL NEED TO MAKE SOME ALTERATIONS TO GET IT TO WORK FOR
    // OTHER STENCILS!

    std::vector<double> contribution(this->gaugeNum*3U, 0.0);

    for (unsigned dir1_iter = 0; dir1_iter < 3; dir1_iter++) // This loops over the components which will be evolved
    {
        unsigned eq_dir_index = dir1_iter*this->gaugeNum;
        unsigned dir1_index = dir1_iter*this->numVectorComponents/3U;

        unsigned diag_index = 2;
        for (unsigned dir2_iter = 0; dir2_iter < 3; dir2_iter++) // This loops over other components which contribute through the loops
        {
            if (dir2_iter != dir1_iter)
            {
                unsigned dir2_index = dir2_iter*this->numVectorComponents/3U;

                // Also need to work out where the grid position at +dir1 - dir2 is stored.
                diag_index += 1;


                /////////////////////////  Hypercharge U(1) calculations:  //////////////////////////////////


                double local_loop_angle =     -static_cast<double>(vector_pointers[dir1_iter][2][dir2_index]) + static_cast<double>(vector_pointers[dir2_iter][2][dir1_index])
                                             + static_cast<double>(vector_pointers[dir1_iter][1][dir2_index]) - static_cast<double>(vector_pointers[dir2_iter][1][dir1_index]);
                double neighbour_loop_angle = -static_cast<double>(vector_pointers[dir1_iter][diag_index][dir2_index]) - static_cast<double>(vector_pointers[dir2_iter][0][dir1_index])
                                             + static_cast<double>(vector_pointers[dir2_iter][0][dir2_index]) + static_cast<double>(vector_pointers[dir1_iter][1][dir1_index]);

                contribution[eq_dir_index] += ( std::sin( local_loop_angle ) - std::sin( neighbour_loop_angle ) )*this->inverse_sqr_spacings[dir2_iter];


                ////////////////////////  Isospin SU(2) calculations:  ///////////////////////////////////////


                // Calculate the contribution from the first Wilson loop:
                // Representation of right-hand matrix
                SU2Representation U_product1 = this->getSU2Representation(vector_pointers[dir2_iter][1], dir1_index, 1U, true);

                // Representation of the left-hand matrix
                SU2Representation U_multiply = this->getSU2Representation(vector_pointers[dir1_iter][1], dir2_index, 1U, false);

                // Multiply them together
                U_product1 = this->SU2Product(U_multiply, U_product1, true);

                // Repeat for the rest of the loop.
                U_multiply = this->getSU2Representation(vector_pointers[dir2_iter][2], dir1_index, 1U, false);
                U_product1 = this->SU2Product(U_multiply, U_product1, true);

                U_multiply = this->getSU2Representation(vector_pointers[dir1_iter][2], dir2_index, 1U, true);

                // Keep the trace part for the last product only if energy is being calculated.
                U_product1 = this->SU2Product(U_multiply, U_product1, this->storeEnergy);



                // Now repeat for the second Wilson loop:
                SU2Representation U_product2 = this->getSU2Representation(vector_pointers[dir1_iter][diag_index], dir2_index, 1U, true);
                U_multiply = this->getSU2Representation(vector_pointers[dir2_iter][0], dir1_index, 1U, true);
                U_product2 = this->SU2Product(U_multiply, U_product2, true);

                U_multiply = this->getSU2Representation(vector_pointers[dir2_iter][0], dir2_index, 1U, false);
                U_product2 = this->SU2Product(U_multiply, U_product2, true);

                U_multiply = this->getSU2Representation(vector_pointers[dir1_iter][1], dir1_index, 1U, false);

                // Never need the trace part for this final product (it will be calculated again at other points in the grid)
                U_product2 = this->SU2Product(U_multiply, U_product2, false);

                // Add results to the contributions to the isospin equations
                for (unsigned comp_iter = 1; comp_iter < 4; comp_iter++)
                    contribution[eq_dir_index + comp_iter] += ( U_product1[comp_iter] - U_product2[comp_iter] )*this->inverse_sqr_spacings[dir2_iter];


                // If energy is being calculated, do some additional calculations now so that wilson loops don't need to be recalculated.
                // This only includes the SM gauge fields, Higgs family are done later to keep it optional and for reuse of variables.
                if (this->storeEnergy)
                    this->storedMagneticEnergy += 2.0*this->inverse_sqr_spacings[dir1_iter]*this->inverse_sqr_spacings[dir2_iter]*(
                        this->inverse_g_sqr*(1.0 - U_product1[0]) 
                      + this->inverse_gp_sqr*(1.0 - std::cos(local_loop_angle))
                    );

                
                ////////////////////  Higgs family calculations:  ////////////////////////////////////

                unsigned offset;
                if (this->gaugeNum > 4U)
                {
                    if (this->usingGeneratorRepresentation)
                        offset = 4U;
                    else
                        offset = 5U;
                }

                if (this->gaugeNum == 5U)
                {
                    // Higgs family U(1)

                    local_loop_angle =     -static_cast<double>(vector_pointers[dir1_iter][2][dir2_index + offset]) + static_cast<double>(vector_pointers[dir2_iter][2][dir1_index + offset])
                                          + static_cast<double>(vector_pointers[dir1_iter][1][dir2_index + offset]) - static_cast<double>(vector_pointers[dir2_iter][1][dir1_index + offset]);
                    neighbour_loop_angle = -static_cast<double>(vector_pointers[dir1_iter][diag_index][dir2_index + offset]) - static_cast<double>(vector_pointers[dir2_iter][0][dir1_index + offset])
                                          + static_cast<double>(vector_pointers[dir2_iter][0][dir2_index + offset]) + static_cast<double>(vector_pointers[dir1_iter][1][dir1_index + offset]);

                    contribution[eq_dir_index + 4] += ( std::sin( local_loop_angle ) - std::sin( neighbour_loop_angle ) )*this->inverse_sqr_spacings[dir2_iter];

                    if (this->storeEnergy)
                        this->storedMagneticEnergy += 2.0*this->inverse_sqr_spacings[dir1_iter]*this->inverse_sqr_spacings[dir2_iter]
                                                    *this->inverse_gpp_sqr*(1.0 - std::cos(local_loop_angle));

                }
                else if (this->gaugeNum == 7U)
                {
                    // Higgs family SU(2) (see isospin case above for more comments)

                    // Calculate the contribution from the first Wilson loop:
                    U_product1 = this->getSU2Representation(vector_pointers[dir2_iter][1], dir1_index, offset, true);
                    U_multiply = this->getSU2Representation(vector_pointers[dir1_iter][1], dir2_index, offset, false);
                    U_product1 = this->SU2Product(U_multiply, U_product1, true);

                    U_multiply = this->getSU2Representation(vector_pointers[dir2_iter][2], dir1_index, offset, false);
                    U_product1 = this->SU2Product(U_multiply, U_product1, true);

                    U_multiply = this->getSU2Representation(vector_pointers[dir1_iter][2], dir2_index, offset, true);
                    U_product1 = this->SU2Product(U_multiply, U_product1, this->storeEnergy);


                    // Now repeat for the second Wilson loop:
                    U_product2 = this->getSU2Representation(vector_pointers[dir1_iter][diag_index], dir2_index, offset, true);
                    U_multiply = this->getSU2Representation(vector_pointers[dir2_iter][0], dir1_index, offset, true);
                    U_product2 = this->SU2Product(U_multiply, U_product2, true);

                    U_multiply = this->getSU2Representation(vector_pointers[dir2_iter][0], dir2_index, offset, false);
                    U_product2 = this->SU2Product(U_multiply, U_product2, true);

                    U_multiply = this->getSU2Representation(vector_pointers[dir1_iter][1], dir1_index, offset, false);
                    U_product2 = this->SU2Product(U_multiply, U_product2, false);

                    // Add results to the contributions to the higgs family equations
                    for (unsigned comp_iter = 4; comp_iter < 7; comp_iter++)
                        contribution[eq_dir_index + comp_iter] += ( U_product1[comp_iter - 3U] - U_product2[comp_iter - 3U] )*this->inverse_sqr_spacings[dir2_iter];

                    if (this->storeEnergy)
                        this->storedMagneticEnergy += 2.0*this->inverse_sqr_spacings[dir1_iter]*this->inverse_sqr_spacings[dir2_iter]
                                                    *this->inverse_gpp_sqr*(1.0 - U_product1[0]);
                }
            }

            // Otherwise there is no contribution
        }
    }

    return contribution;
}


std::vector<double> WilsonLoops::TwoHDM::calcElectricContributions(const float *const local_vector_fields[2]) const
{
    std::vector<double> contribution(this->gaugeNum*3U, 0.0);

    for (unsigned dir_iter = 0; dir_iter < 3; dir_iter++)
    {
        unsigned eq_dir_index = dir_iter*this->gaugeNum;
        unsigned dir_index = dir_iter*this->numVectorComponents/3U;


        /////////////////////////  Hypercharge U(1) calculations:  ////////////////////////////////////


        double loop_angle = static_cast<double>(local_vector_fields[1][dir_index]) - static_cast<double>(local_vector_fields[0][dir_index]);
        contribution[eq_dir_index] = std::sin( loop_angle );


        /////////////////////////  Isospin SU(2) calculations:  ///////////////////////////////////////


        SU2Representation U_product = this->getSU2Representation(local_vector_fields[0], dir_index, 1U, true);
        SU2Representation U_multiply = this->getSU2Representation(local_vector_fields[1], dir_index, 1U, false);

        U_product = this->SU2Product(U_multiply, U_product, this->storeEnergy);

        for (unsigned comp_iter = 1; comp_iter < 4; comp_iter++)
            contribution[eq_dir_index + comp_iter] = U_product[comp_iter];


        // If energy is being calculated, do some calculations now to avoid re-evaluations of wilson loops.
        // This only includes the SM gauge fields, Higgs family will be done later to keep it optional and for reuse of variables.
        if (this->storeEnergy)
            this->storedElectricEnergy += 4.0*this->inverse_sqr_spacings[dir_iter]*(
                this->inverse_g_sqr*(1.0 - U_product[0]) 
              + this->inverse_gp_sqr*(1.0 - std::cos(loop_angle))
            );

        
        //////////////////////////  Higgs family calculations:  ///////////////////////////////////////


        unsigned offset;
        if (this->gaugeNum > 4U)
        {
            if (this->usingGeneratorRepresentation)
                offset = 4U;
            else
                offset = 5U;
        }

        if (this->gaugeNum == 5U)
        {
            // Higgs family U(1)

            loop_angle = static_cast<double>(local_vector_fields[1][dir_index + offset]) - static_cast<double>(local_vector_fields[0][dir_index + offset]);
            contribution[eq_dir_index + 4U] = std::sin( loop_angle );

            if (this->storeEnergy)
                this->storedElectricEnergy += 4.0*this->inverse_sqr_spacings[dir_iter]*this->inverse_gpp_sqr*(1.0 - std::cos(loop_angle));
        }
        else if (this->gaugeNum == 7U)
        {
            // Higgs family SU(2)

            U_product = this->getSU2Representation(local_vector_fields[0], dir_index, offset, true);
            U_multiply = this->getSU2Representation(local_vector_fields[1], dir_index, offset, false);

            U_product = this->SU2Product(U_multiply, U_product, this->storeEnergy);

            for (unsigned comp_iter = 4; comp_iter < 7; comp_iter++)
                contribution[eq_dir_index + comp_iter] = U_product[comp_iter - 3U];

            if (this->storeEnergy)
                this->storedElectricEnergy += 4.0*this->inverse_sqr_spacings[dir_iter]*this->inverse_gpp_sqr*(1.0 - U_product[0]);
        }
    }

    return contribution;
}


void WilsonLoops::TwoHDM::evolve(float *const local_vector_fields[2], std::vector<double> equation_RHS)
{
    // I will refer to u_i = U_i(t+dt)U_i^\dagger(t) throughout (no sum over i), where the group U belongs to is implicit in the sections.
    for (unsigned dir_iter = 0; dir_iter < 3; dir_iter++)
    {
        unsigned eq_dir_index = dir_iter*this->gaugeNum;
        unsigned dir_index = dir_iter*this->numVectorComponents/3U;


        //////////////////////////  Hypercharge U(1) calculations:  /////////////////////////////////////


        // Equation gives me the imaginary part of u_i. For small dt, u_i is close to 1, so take sin^-1(rhs) in the range -pi/2 to pi/2 for the phase of u_i.
        // Next timestep is U_i(t+dt) = u_iU_i(t), so I just add the phase of u_i to that of U_i(t).

        local_vector_fields[0][dir_index] = static_cast<float>(std::asin(equation_RHS[eq_dir_index]) + static_cast<double>(local_vector_fields[1][dir_index]));

        // What if RHS > 1?????
        if (std::pow(equation_RHS[eq_dir_index],2) > 1)
            std::cout << "Panic in evolve (rhs = " << equation_RHS[eq_dir_index] << " )" << std::endl;


        //////////////////////////  Isospin SU(2) calculations:  ////////////////////////////////////////


        // Equation gives me the traceless part of u_i, c_i^a, where u_i = c_i^0\sigma^0 + ic_i^a\sigma^a. 
        // For small dt, u_i is close to \sigma^0, so take c_i^0 = +\sqrt{1-c_i^ac_i^a}.
        // Then I can do the SU(2) multiplication u_iU_i(t) to get U_i(t+dt) in c_i^\mu form above.
        // The last step is to convert that into the three degrees of freedom form that is stored.

        double c_mag_sqr = 0;
        for (unsigned comp_iter = 1; comp_iter < 4; comp_iter++)
            c_mag_sqr += std::pow(equation_RHS[eq_dir_index + comp_iter], 2);

        // WHAT IF C_MAG_SQR > 1????????
        if (c_mag_sqr > 1)
            std::cout << "Panic in evolve (c_mag_sqr = " << c_mag_sqr << " )" << std::endl;

        SU2Representation u_multiply = {std::sqrt(1.0 - c_mag_sqr), equation_RHS[eq_dir_index + 1], equation_RHS[eq_dir_index + 2], equation_RHS[eq_dir_index + 3]};
        SU2Representation U_product = this->getSU2Representation(local_vector_fields[1], dir_index, 1U, false);

        // Get the SU(2) matrix uU(t) in c_i^\mu form.
        U_product = this->SU2Product(u_multiply, U_product, true);

        if (this->usingGeneratorRepresentation)    // Approach when 3 dofs are stored:
        {
            // Convert c_i^\mu form into w_i^a form.
            std::vector<float> U_next = this->invertSU2Representation(U_product);
            
            for (unsigned comp_iter = 0; comp_iter < 3; comp_iter++)
                local_vector_fields[0][dir_index + comp_iter + 1] = U_next[comp_iter];
        }
        else                                       // Approach when 4 dofs are stored:
        {
            for (unsigned comp_iter = 0; comp_iter < 4; comp_iter++)
                local_vector_fields[0][dir_index + comp_iter + 1] = U_product[comp_iter];
        }


        /////////////////////////  Higgs family calculations:  /////////////////////////////////////////


        unsigned offset;
        if (this->gaugeNum > 4U)
        {
            if (this->usingGeneratorRepresentation)
                offset = 4U;
            else
                offset = 5U;
        }

        if (this->gaugeNum == 5U)
        {
            // Higgs family U(1)

            local_vector_fields[0][dir_index + offset] = static_cast<float>(std::asin(equation_RHS[eq_dir_index + 4U]) + static_cast<double>(local_vector_fields[1][dir_index + offset]));

            // What if RHS > 1?????
            if (std::pow(equation_RHS[eq_dir_index + 4U],2) > 1)
                std::cout << "Panic in evolve (HF rhs = " << equation_RHS[eq_dir_index + 4U] << " )" << std::endl;
        }
        else if (this->gaugeNum == 7U)
        {
            // Higgs family SU(2)

            c_mag_sqr = 0;
            for (unsigned comp_iter = 4; comp_iter < 7; comp_iter++)
                c_mag_sqr += std::pow(equation_RHS[eq_dir_index + comp_iter], 2);

            // WHAT IF C_MAG_SQR > 1????????
            if (c_mag_sqr > 1)
                std::cout << "Panic in evolve (HF c_mag_sqr = " << c_mag_sqr << " )" << std::endl;

            u_multiply = {std::sqrt(1.0 - c_mag_sqr), equation_RHS[eq_dir_index + 4], equation_RHS[eq_dir_index + 5], equation_RHS[eq_dir_index + 6]};
            U_product = this->getSU2Representation(local_vector_fields[1], dir_index, offset, false);

            U_product = this->SU2Product(u_multiply, U_product, true);

            if (this->usingGeneratorRepresentation)    // Approach when 3 dofs are stored:
            {
                std::vector<float> U_next = this->invertSU2Representation(U_product);
                
                for (unsigned comp_iter = 0; comp_iter < 3; comp_iter++)
                    local_vector_fields[0][dir_index + comp_iter + offset] = U_next[comp_iter];
            }
            else                                       // Approach when 4 dofs are stored:
            {
                for (unsigned comp_iter = 0; comp_iter < 4; comp_iter++)
                    local_vector_fields[0][dir_index + comp_iter + offset] = U_product[comp_iter];
            }
        }
    }
}