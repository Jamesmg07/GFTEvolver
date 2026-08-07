#include "Double_SO_N.hpp"

//////////////////////////////////////////////  Initialisers  ///////////////////////////////////////////////////////

void Double_SO_N::initVariables()
{
    this->applyPotentialNormalisation = false;
    this->potentialEnergyOffset = 0.0;

    this->m1 = 0;
    this->m2 = 0;
    this->m12r = 0;
    this->m12i = 0;

    this->l1 = 0;
    this->l2 = 0;
    this->l3 = 0;
    this->l4p5 = 0;
    this->l4m5 = 0;
    this->l5i = 0;
    this->l6r = 0;
    this->l6i = 0;
    this->l7r = 0;
    this->l7i = 0;

    this->breakToSON = false;
    this->breakToUHalfN = false;

    this->fieldSqrMagnitude[0] = 0;
    this->fieldSqrMagnitude[1] = 0;
    this->fieldsDot = 0;
    this->fieldsAntisym = 0;
}

void Double_SO_N::configure(const std::string path, const bool debug)
{
    std::ifstream ifs(path);

    // Read in the parameter values from the configuration file
    std::string description; // Descriptions in the config file (to be dumped)
    if (ifs.is_open())
    {
        // Skip the initial configuration file description
        for(int iter = 0; iter < 40; iter++) std::getline(ifs, description);

        // Load number of components in f1.
        std::getline(ifs, description, ':');
        ifs >> this->numComponents[0];

        // Load the optional user-supplied potential-energy shift.
        std::getline(ifs, description, ':');
        ifs >> this->applyPotentialNormalisation;

        std::getline(ifs, description, ':');
        ifs >> this->potentialEnergyOffset;

        if (this->numComponents[0] > this->totalNumComponents)
            throw std::runtime_error("POTENTIAL::DOUBLE_SO_N:: The number of field components in f1 has been set to " + std::to_string(this->numComponents[0])
                                    + " and the total number has been set to " + std::to_string(this->totalNumComponents)
                                    + ".\n The total number must be greater or equal to the number in f1."); 
        
        this->numComponents[1] = this->totalNumComponents - this->numComponents[0];

        //////////////////////// Parameters that respect the full SO(N) x SO(M) symmetry //////////////////////////////
        std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> this->m1 >> this->m2;

        std::getline(ifs, description, ':');
        ifs >> this->l1 >> this->l2;

        std::getline(ifs, description, ':');
        ifs >> this->l3;

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Turn on parameters that explicitly break the symmetry down to a single SO(N)?
        for (unsigned iter = 0; iter < 4; iter++) std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> this->breakToSON;

        if (this->breakToSON)
        {
            // Check that the number of field components are equal.
            if (this->numComponents[0] != this->numComponents[1])
                throw std::runtime_error("POTENTIAL::DOUBLE_SO_N:: Attempting to use the potential parameters that explicitly break the symmetry down"
                "to a single SO(N) requires that the number of components in f1 and f2 are the same.\n"
                "The total number of scalar field components assigned (" + std::to_string(this->totalNumComponents) + ") should be twice the "
                "the number of components assigned to f1 (" + std::to_string(this->numComponents[0]) + ").");

            std::getline(ifs, description);
            std::getline(ifs, description, ':');
            ifs >> this->m12r;

            std::getline(ifs, description, ':');
            ifs >> this->l4p5;

            std::getline(ifs, description, ':');
            ifs >> this->l6r >> this->l7r;

            // Turn on parameters that explicitly break the symmetry down to a single U(N/2)?
            for (unsigned iter = 0; iter < 5; iter++) std::getline(ifs, description);
            std::getline(ifs, description, ':');
            ifs >> this->breakToUHalfN;

            if (this->breakToUHalfN)
            {
                // Check that the number of field components in f1 (=f2) is even.
                if (this->numComponents[0]%2 != 0)
                    throw std::runtime_error("POTENTIAL::DOUBLE_SO_N:: Attempting to use the potential parameters that explicitly break the symmetry down"
                    "to a single U(N/2) requires that the number of components in f1 (equal to f2) is even.\n"
                    "The number of components assigned to f1 is " + std::to_string(this->numComponents[0]) );

                std::getline(ifs, description);
                std::getline(ifs, description, ':');
                ifs >> this->m12i;

                std::getline(ifs, description, ':');
                ifs >> this->l4m5 >> this->l5i;

                std::getline(ifs, description, ':');
                ifs >> this->l6i >> this->l7i;
            }
        }

    } 

    ifs.close();

    if (debug)
    {
        std::cout << "POTENTIALS::DOUBLE_SO_N::\n"
                  << "#Components in f1: " << this->numComponents[0] << ", #Components in f2:" << this->numComponents[1] << "\n"
                  << "Apply optional potential-energy normalisation: " << this->applyPotentialNormalisation << ", potential-energy offset added when enabled: " << this->potentialEnergyOffset << "\n"
                  << "m1: " << this->m1 << ", m2: " << this->m2 << ", l1: " << this->l1 << ", l2: " << this->l2 << ", l3: " << this->l3
                  << "\n m12r: " << this->m12r << ", l4p5: " << this->l4p5 << ", l6r: " << this->l6r << ", l7r: " << this->l7r
                  << "\n m12i: " << this->m12i << ", l4m5: " << this->l4m5 << ", l5i: " << this->l5i << ", l6i: " << this->l6i << ", l7i: " << this->l7i
                  << "\n" << std::endl;
    }
}

/////////////////////////////////////////  Constructors/Destructors  ///////////////////////////////////////////////

Double_SO_N::Double_SO_N(const unsigned &total_num_components)
    : totalNumComponents(total_num_components)
{
    this->initVariables();
    this->configure(std::string(SOURCE_DIR) + "/Config/Double_SO_N.cfg", true);
}

Double_SO_N::~Double_SO_N()
{
}

//////////////////////////////////////////////  Public functions  //////////////////////////////////////////////////

std::vector<double> Double_SO_N::calcPotentialDerivatives(const float *field)
{
    // Calculate the field magnitudes, the dot product of the two fields and the antisymmetric product (as defined in config) of the two fields.
    for (unsigned field_iter = 0; field_iter < 2; field_iter++)
    {
        this->fieldSqrMagnitude[field_iter] = 0.0;
        for (unsigned comp_iter = 0; comp_iter < this->numComponents[field_iter]; comp_iter++)
        {
            this->fieldSqrMagnitude[field_iter] += std::pow(static_cast<double>(field[field_iter*this->numComponents[0] + comp_iter]), 2);
        }

        // Below assumes the fields have the same number of components
        if (this->breakToSON)
        {
            this->fieldsDot = 0.0;
            for (unsigned comp_iter = 0; comp_iter < this->numComponents[0]; comp_iter++)
            {
                this->fieldsDot += static_cast<double>(field[comp_iter])*static_cast<double>(field[this->numComponents[0] + comp_iter]);
            }

            // Below additionally assumes the number of components is even
            if (this->breakToUHalfN)
            {
                this->fieldsAntisym = 0.0;
                for (unsigned comp_iter = 0; comp_iter < this->numComponents[0]/2; comp_iter++)
                {
                    this->fieldsAntisym += static_cast<double>(field[2*comp_iter])*static_cast<double>(field[this->numComponents[0] + 2*comp_iter + 1])
                                         - static_cast<double>(field[2*comp_iter + 1])*static_cast<double>(field[this->numComponents[0] + 2*comp_iter]);
                }
            }
        }
    }



    // Now that all field quantities have been updated, can calculate the contribution from the terms that do not explicitly break the symmetry.
    std::vector<double> potential_contributions(this->totalNumComponents, 0.0);

    // Start with f1:
    for (unsigned comp_iter = 0; comp_iter < this->numComponents[0]; comp_iter++)
    {
        // Full symmetry contribution
        potential_contributions[comp_iter] += ( -2.0*this->m1 + 4.0*this->l1*this->fieldSqrMagnitude[0] + 2.0*this->l3*this->fieldSqrMagnitude[1] 
                                              )*static_cast<double>(field[comp_iter]);

        // Single SO(N) contribution
        if (this->breakToSON)
        {
            unsigned f2index = this->numComponents[0] + comp_iter;

            potential_contributions[comp_iter] += ( -this->m12r + 2.0*this->l4p5*this->fieldsDot + this->l6r*this->fieldSqrMagnitude[0] 
                                                + this->l7r*this->fieldSqrMagnitude[1] )*static_cast<double>(field[f2index])
                                                + 2.0*this->l6r*this->fieldsDot*static_cast<double>(field[comp_iter]);

            // U(N/2) contribution
            if (this->breakToUHalfN)
            {
                int alternator = 1 - 2*(comp_iter%2); // +1 or -1
                potential_contributions[comp_iter] += ( this->m12i + 2.0*this->l4m5*this->fieldsAntisym - 2.0*this->l5i*this->fieldsDot 
                                                      - this->l6i*this->fieldSqrMagnitude[0] - this->l7i*this->fieldSqrMagnitude[1]
                                                      )*alternator*static_cast<double>(field[f2index + alternator])
                                                    - 2.0*this->l5i*this->fieldsAntisym*static_cast<double>(field[f2index])
                                                    - 2.0*this->l6i*this->fieldsAntisym*static_cast<double>(field[comp_iter]);
            }
        }
    }

    // Now repeat for f2:
    for (unsigned comp_iter = 0; comp_iter < this->numComponents[1]; comp_iter++)
    {
        unsigned f2index = this->numComponents[0] + comp_iter;

        // Full symmetry contribution
        potential_contributions[f2index] += ( -2.0*this->m2 + 4.0*this->l2*this->fieldSqrMagnitude[1] + 2.0*this->l3*this->fieldSqrMagnitude[0] )
                                            * static_cast<double>(field[f2index]);

        // Single SO(N) contribution
        if (this->breakToSON)
        {
            potential_contributions[f2index] += ( -this->m12r + 2.0*this->l4p5*this->fieldsDot + this->l6r*this->fieldSqrMagnitude[0]
                                                + this->l7r*this->fieldSqrMagnitude[1] )*static_cast<double>(field[comp_iter])
                                                + 2.0*this->l7r*this->fieldsDot*static_cast<double>(field[f2index]);

            // U(N/2) contribution
            if (this->breakToUHalfN)
            {
                int alternator = 1 - 2*(comp_iter%2);
                potential_contributions[f2index] += -( this->m12i + 2.0*this->l4m5*this->fieldsAntisym - 2.0*this->l5i*this->fieldsDot
                                                    - this->l6i*this->fieldSqrMagnitude[0] - this->l7i*this->fieldSqrMagnitude[1] 
                                                    )*alternator*static_cast<double>(field[comp_iter + alternator])
                                                    - 2.0*this->l5i*this->fieldsAntisym*static_cast<double>(field[comp_iter])
                                                    - 2.0*this->l7i*this->fieldsAntisym*static_cast<double>(field[f2index]);
            }
        }
    }

    return potential_contributions;
}

// V = V_full + V_single + V_U, where

// V_full = -m1*|f1|^2 - m2*|f2|^2 + l1*|f1|^4 + l2*|f2|^4 + l3*|f1|^2|f2|^2,
// which is the part of the potential that is symmetric under SO(N) x SO(M) symmetry. f1 is an N component field and f2 is an M component field.

// In order for V_single to be used, f1 and f2 need to have the same number of components, as there will now just be a single SO(N) symmetry that
// acts the same way on both fields. Potential parameters below will only be used if a boolean is set to true.

// V_single = -m12r*f1^T*f2 + l4p5*(f1^T*f2)^2 + l6r*|f1|^2(f1^T*f2) + l7r*|f2|^2(f1^T*f2),
// which are all of the terms that can be formed from dot products of the two field vectors (v^T*u should be understood to be the dot product of v and u).

// The next part of the potential will further explicitly break the symmetry down to U(N/2), as it includes terms like f1^T*J*f2,
// where J = I_{N/2} \otimes \epsilon and \epsilon is the 2x2 antisymmetric tensor. This seems like a slightly odd construction,
// but in the language of complex (N/2)-tuples, it is simply Im(\Phi_1^\dagger\Phi_2), whereas f1^T*f2 is Re(\Phi_1^\dagger\Phi_2).

// V_U = +m12i*f1^T*J*f2 + l4m5*(f1^T*Jf2)^2 - 2*l5i*(f1^T*f2)*(f1^T*J*f2) - l6i*|f1|^2*(f1^T*J*f2) - l7i*|f2|^2*(f1^T*J*f2).

// Calculate the potential energy directly from the supplied field.
float Double_SO_N::calcPotentialEnergy(const float *field) const
{
    // Calculate the squared magnitude of each field.
    double field_sqr_magnitude[2] = {0.0, 0.0};

    for (unsigned field_iter = 0; field_iter < 2; field_iter++)
    {
        for (unsigned comp_iter = 0;
             comp_iter < this->numComponents[field_iter];
             comp_iter++)
        {
            field_sqr_magnitude[field_iter]
                += std::pow(
                    static_cast<double>(
                        field[field_iter*this->numComponents[0] + comp_iter]),
                    2);
        }
    }

    double fields_dot = 0.0;
    double fields_antisym = 0.0;

    // These quantities are only needed when the potential breaks
    // SO(N) x SO(N) to a single SO(N).
    if (this->breakToSON)
    {
        for (unsigned comp_iter = 0;
             comp_iter < this->numComponents[0];
             comp_iter++)
        {
            fields_dot
                += static_cast<double>(field[comp_iter])
                 * static_cast<double>(
                       field[this->numComponents[0] + comp_iter]);
        }

        // This quantity is only needed when the symmetry is further
        // broken to U(N/2).
        if (this->breakToUHalfN)
        {
            for (unsigned comp_iter = 0;
                 comp_iter < this->numComponents[0]/2;
                 comp_iter++)
            {
                fields_antisym
                    += static_cast<double>(field[2*comp_iter])
                     * static_cast<double>(
                           field[this->numComponents[0]
                               + 2*comp_iter + 1])
                     - static_cast<double>(field[2*comp_iter + 1])
                     * static_cast<double>(
                           field[this->numComponents[0]
                               + 2*comp_iter]);
            }
        }
    }

    // Full SO(N) x SO(M) contribution.
    float potential
        = -this->m1*field_sqr_magnitude[0]
          -this->m2*field_sqr_magnitude[1]
          +this->l1*std::pow(field_sqr_magnitude[0], 2)
          +this->l2*std::pow(field_sqr_magnitude[1], 2)
          +this->l3*field_sqr_magnitude[0]*field_sqr_magnitude[1];

    // Single SO(N) contribution.
    if (this->breakToSON)
    {
        potential
            += -this->m12r*fields_dot
               +this->l4p5*std::pow(fields_dot, 2)
               +this->l6r*field_sqr_magnitude[0]*fields_dot
               +this->l7r*field_sqr_magnitude[1]*fields_dot;

        // U(N/2) contribution.
        if (this->breakToUHalfN)
        {
            potential
                += this->m12i*fields_antisym
                   +this->l4m5*std::pow(fields_antisym, 2)
                   -2.0*this->l5i*fields_dot*fields_antisym
                   -this->l6i*field_sqr_magnitude[0]*fields_antisym
                   -this->l7i*field_sqr_magnitude[1]*fields_antisym;
        }
    }

    // A constant shift changes the energy zero but not the equations of motion.
    if (this->applyPotentialNormalisation)
    {
        potential += static_cast<float>(
            this->potentialEnergyOffset);
    }

    return potential;
}
