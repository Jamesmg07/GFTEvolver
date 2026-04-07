#include "Double_SO_N.hpp"

//////////////////////////////////////////////  Initialisers  ///////////////////////////////////////////////////////

void Double_SO_N::initVariables()
{
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
        bool break_to_SO_N = false;
        for (unsigned iter = 0; iter < 4; iter++) std::getline(ifs, description);
        std::getline(ifs, description, ':');
        ifs >> break_to_SO_N;

        if (break_to_SO_N)
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
            bool break_to_U_half_N = false;
            for (unsigned iter = 0; iter < 5; iter++) std::getline(ifs, description);
            std::getline(ifs, description, ':');
            ifs >> break_to_U_half_N;

            if (break_to_U_half_N)
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

float Double_SO_N::calcPotentialEnergy(const float *field) const
{
    //return 0.25f*this->lambda*powf(this->fieldSqrMagnitude - this->etaSqr, 2);
    return 0.f;
}

std::vector<double> Double_SO_N::calcPotentialDerivatives(const float *field)
{
    std::vector<double> potential_contributions(this->totalNumComponents, 0.f);

    // // Calculate the |field|^2 and save it for possible later use in calculating the potential energy.
    // this->fieldSqrMagnitude = 0.f;
    // for (int iter = 0; iter < this->numComponents; iter++)
    // {
    //     this->fieldSqrMagnitude += std::pow(static_cast<double>(field[iter]), 2);
    // }

    // // Calculate potential contributions to the equations of motion.
    // for (int iter = 0; iter < this->numComponents; iter++)
    // {
    //     potential_contributions[iter] = this->lambda*( this->fieldSqrMagnitude - this->etaSqr )*static_cast<double>(field[iter]);
    // }

    return potential_contributions;
}
