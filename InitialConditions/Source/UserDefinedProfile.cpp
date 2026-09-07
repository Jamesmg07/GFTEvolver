#include "UserDefinedProfile.hpp"

#include <cstddef>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <complex>
#include <vector>

using namespace std;



namespace
{
    double interpProfile(
        const std::vector<double>& profile,
        double rPos,
        double farValue)
    {
        const int nProfile =
            static_cast<int>(profile.size());

        const int rC =
            static_cast<int>(std::round(rPos));

        const int rM = rC - 1;
        const int rP = rC + 1;

        if (rP >= nProfile)
        {
            return farValue;
        }

        if (rC == 0)
        {
            return
                (-(rC - rPos) * profile[rP])
                + ((rP - rPos) * profile[rC]);
        }

        return
            (((rM - rPos) * (rC - rPos) * profile[rP]) / 2.0)
            - (((rM - rPos) * (rP - rPos) * profile[rC]))
            + (((rC - rPos) * (rP - rPos) * profile[rM]) / 2.0);
    }
}

namespace
{
static inline void build_u_matrices(
    double x_1,
    double y_1,
    double z_1,
    double r_1,
    double x_2,
    double y_2,
    double z_2,
    double r_2,
    complex<double> u_1[2][2],
    complex<double> u_2[2][2]
) {
    double EPS = 1e-8;

    if (fabs(z_1 - r_1) < EPS) {
        u_1[0][0] = complex<double>(1.0, 0.0);
        u_1[0][1] = complex<double>(0.0, 0.0);
        u_1[1][0] = complex<double>(0.0, 0.0);
        u_1[1][1] = complex<double>(1.0, 0.0);

        u_2[0][0] = complex<double>(0.0, 0.0);
        u_2[0][1] = complex<double>(-1.0, 0.0);
        u_2[1][0] = complex<double>(1.0, 0.0);
        u_2[1][1] = complex<double>(0.0, 0.0);
        return;
    }

    if (fabs(z_2 + r_2) < EPS) {
        u_1[0][0] = complex<double>(0.0, 0.0);
        u_1[0][1] = complex<double>(-1.0, 0.0);
        u_1[1][0] = complex<double>(1.0, 0.0);
        u_1[1][1] = complex<double>(0.0, 0.0);

        u_2[0][0] = complex<double>(-1.0, 0.0);
        u_2[0][1] = complex<double>(0.0, 0.0);
        u_2[1][0] = complex<double>(0.0, 0.0);
        u_2[1][1] = complex<double>(-1.0, 0.0);
        return;
    }
   
    if (fabs(z_1 - r_1) >= EPS && fabs(z_2 + r_2) >= EPS && fabs(z_2 - r_2) < EPS) {
        u_1[0][0] = complex<double>(0.0, 0.0);
        u_1[0][1] = complex<double>(-1.0, 0.0);
        u_1[1][0] = complex<double>(1.0, 0.0);
        u_1[1][1] = complex<double>(0.0, 0.0);

        u_2[0][0] = complex<double>(0.0, 0.0);
        u_2[0][1] = complex<double>(-1.0, 0.0);
        u_2[1][0] = complex<double>(1.0, 0.0);
        u_2[1][1] = complex<double>(0.0, 0.0);
        return;
    }

    const double cos_1 = pow(0.5 * (1 + (z_1 / r_1)), 0.5);
    u_1[0][0] = complex<double>(cos_1, 0.0);
    u_1[0][1] = complex<double>(-x_1 / (2 * r_1 * cos_1), y_1 / (2 * r_1 * cos_1));
    u_1[1][0] = complex<double>(x_1 / (2 * r_1 * cos_1), y_1 / (2 * r_1 * cos_1));
    u_1[1][1] = complex<double>(cos_1, 0.0);

    const double sin_2 = pow(0.5 * (1 - (z_2 / r_2)), 0.5);
    u_2[0][0] = complex<double>(-sin_2, 0.0);
    u_2[0][1] = complex<double>(-x_2 / (2 * r_2 * sin_2), y_2 / (2 * r_2 * sin_2));
    u_2[1][0] = complex<double>(x_2 / (2 * r_2 * sin_2), y_2 / (2 * r_2 * sin_2));
    u_2[1][1] = complex<double>(-sin_2, 0.0);
    }
}

namespace
{
    void buildScalarField(
    complex<double> u_1[2][2],
    complex<double> u_2[2][2],
    const double phi_both[4],
    complex<double> (&phi)[4],
    double gamma1,
    double gamma2
) {
    const double pi = 4.0 * atan(1.0);
    const double gamma_param_1 = gamma1 * pi;
    const double gamma_param_2 = gamma2 * pi;  
    const double monopolePrefactor = pow(2, -1.5);
        
    const double sqrt2_inv = 1.0 / sqrt(2.0);
    complex<double> M[2][2] = {
        {sqrt2_inv, sqrt2_inv},
        {-sqrt2_inv, sqrt2_inv}
    };

    complex<double> T1[2][2] = {
        {exp(complex<double>(0, 0.5 * gamma_param_1)), complex<double>(0.0, 0.0)},
        {complex<double>(0.0, 0.0), exp(complex<double>(0, -0.5 * gamma_param_1))}
    };
    complex<double> T2[2][2] = {
        {exp(complex<double>(0, 0.5 * gamma_param_2)), complex<double>(0.0, 0.0)},
        {complex<double>(0.0, 0.0), exp(complex<double>(0, -0.5 * gamma_param_2))}
    };

    complex<double> C1[2][2], C2[2][2];
    for (int row = 0; row < 2; ++row) {
        for (int col = 0; col < 2; ++col) {
            C1[row][col] = complex<double>(0.0, 0.0);
            C2[row][col] = complex<double>(0.0, 0.0);
            for (int index = 0; index < 2; ++index) {
                C1[row][col] += T1[row][index] * u_2[index][col];
                C2[row][col] += T2[row][index] * u_2[index][col];
            }
        }
    }

    complex<double> A1[2][2], A2[2][2];
    for (int row = 0; row < 2; ++row) {
        for (int col = 0; col < 2; ++col) {
            A1[row][col] = complex<double>(0.0, 0.0);
            A2[row][col] = complex<double>(0.0, 0.0);
            for (int index = 0; index < 2; ++index) {
                A1[row][col] += u_1[row][index] * C1[index][col];
                A2[row][col] += u_1[row][index] * C2[index][col];
            }
        }
    }

    complex<double> B1[2][2], B2[2][2];
    for (int row = 0; row < 2; ++row) {
        for (int col = 0; col < 2; ++col) {
            B1[row][col] = complex<double>(0.0, 0.0);
            B2[row][col] = complex<double>(0.0, 0.0);
            for (int index = 0; index < 2; ++index) {
                B1[row][col] += A1[row][index] * M[index][col];
                B2[row][col] += A2[row][index] * M[index][col];
            }
        }
    }

    complex<double> TP[4][4];
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 2; ++c) {
            for (int s = 0; s < 2; ++s) {
                for (int t = 0; t < 2; ++t) {
                    TP[2 * r + s][2 * c + t] = B1[r][c] * B2[s][t];
                }
            }
        }
    }

    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 4; ++c) {
            TP[r][c] *= monopolePrefactor;
        }
    }



    for (int r = 0; r < 4; ++r) {
        phi[r] = complex<double>(0.0, 0.0);
        for (int c = 0; c < 4; ++c) {
            phi[r] += TP[r][c] * phi_both[c];
        }
    }


    }
}

namespace
{
    void buildGaugeFields(
    double x, double y, double z,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double gamma1, double gamma2,
    double V_profile, double W_profile,
    double g, double gpp,
    double V_inf_amp, double W_inf_amp,
    double V_ai[3][3], double W_ai[3][3])
    {
        const double pi = 4.0 * atan(1.0);
        const double gamma_param_1 = gamma1 * pi;
        const double gamma_param_2 = gamma2 * pi;

        const double R_gamma1[3][3] = {
            {cos(gamma_param_1),  sin(gamma_param_1), 0.0},
            {-sin(gamma_param_1), cos(gamma_param_1), 0.0},
            {0.0, 0.0, 1.0}
        };

        const double R_gamma2[3][3] = {
            {cos(gamma_param_2),  sin(gamma_param_2), 0.0},
            {-sin(gamma_param_2), cos(gamma_param_2), 0.0},
            {0.0, 0.0, 1.0}
        };
    
        const double EPS_LOCAL = 1.0e-12;
        
        double r_sq = x * x + y * y + z * z;
        if (r_sq < EPS_LOCAL * EPS_LOCAL) {
            r_sq = EPS_LOCAL * EPS_LOCAL;
        }
        double r1_sq = x1 * x1 + y1 * y1 + z1 * z1;
        if (r1_sq < EPS_LOCAL * EPS_LOCAL) {
            r1_sq = EPS_LOCAL * EPS_LOCAL;
        }
        double r2_sq = x2 * x2 + y2 * y2 + z2 * z2;
        if (r2_sq < EPS_LOCAL * EPS_LOCAL) {
            r2_sq = EPS_LOCAL * EPS_LOCAL;
        }
    
    
        double rho_1 = hypot(x1, y1);
        double rho_2 = hypot(x2, y2);

        double phi1   = atan2(y1, x1);
    
        double phi2   = atan2(y2, x2);
    
        double sintheta1 = rho_1 / sqrt(r1_sq);
        double costheta1 = z1 / sqrt(r1_sq);
        double sintheta2 = rho_2 / sqrt(r2_sq);
        double costheta2 = z2 / sqrt(r2_sq);

        double zhat[3] = {0.0, 0.0, 1.0};
   
        double halfsintheta1 = sqrt((sqrt(r1_sq) - z1) / (2 * sqrt(r1_sq)));
        double halfsintheta2 = sqrt((sqrt(r2_sq) - z2) / (2 * sqrt(r2_sq)));
    
        double halfcostheta1 = sqrt((sqrt(r1_sq) + z1) / (2 * sqrt(r1_sq)));
        double halfcostheta2 = sqrt((sqrt(r2_sq) + z2) / (2 * sqrt(r2_sq)));
    
        double P[3][3] = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 0.0}
        };
    
   
    
    
    
        double R_m[3][3] = {
            {halfcostheta1 * halfcostheta1 - cos(2 * phi1) *halfsintheta1 *halfsintheta1, -1 * halfsintheta1 * halfsintheta1 * sin(2 * phi1), cos(phi1) * sintheta1},
            {-1 * halfsintheta1 * halfsintheta1 * sin(2 * phi1), halfcostheta1 * halfcostheta1 + cos(2 * phi1) * halfsintheta1 * halfsintheta1, sintheta1 * sin(phi1)},
            {-1 * cos(phi1) * sintheta1, -1 * sintheta1 * sin(phi1), costheta1}
        };
        
        double R_mbar[3][3] = {
            {1.0 / 2.0 *(1 - costheta2 - 2 * halfcostheta2 * halfcostheta2 *  cos(2 * phi2)), -1 * halfcostheta2 * halfcostheta2 * sin(2 * phi2), -1 * cos(phi2) * sintheta2},
            {-1 * halfcostheta2 * halfcostheta2 * sin(2 * phi2), halfcostheta2 * halfcostheta2 * cos(2 * phi2) + halfsintheta2 * halfsintheta2, -1 * sintheta2 * sin(phi2)},
            {cos(phi2) * sintheta2, sintheta2 * sin(phi2), -1 * costheta2}
        };
        // MISSING BOTTOM ROW
    
        double gLr_ai_m[3][3] = {
            {-2 * (cos(phi1) * halfsintheta1 * halfsintheta1 * sin(phi1)), cos(phi1) * cos(phi1) + costheta1 * sin(phi1) * sin(phi1),-1 * sintheta1 * sin(phi1)},
            {-1 * costheta1 * cos(phi1) * cos(phi1) - sin(phi1) * sin(phi1), halfsintheta1 * halfsintheta1 * sin(2 * phi1), cos(phi1) * sintheta1},
            {0,0,0}
        };
        double gLr_ai_mbar[3][3] = {
            {halfcostheta2 * halfcostheta2 * sin(2 * phi2), -1 * cos(phi2) * cos(phi2) + costheta2 * sin(phi2) * sin(phi2), -1 * sintheta2 * sin(phi2)},
            {-1 * costheta2 * cos(phi2) * cos(phi2) + sin(phi2) * sin(phi2), -2 * halfcostheta2 * halfcostheta2 * cos(phi2) * sin(phi2), cos(phi2) * sintheta2},
            {0,0,0}
        };
    
        double L_ai_m[3][3] = {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        };
        double L_ai_mbar[3][3] = {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        };
    
        for (int a_idx = 0; a_idx < 3; ++a_idx) {
            for (int i_idx = 0; i_idx < 3; ++i_idx) {
                L_ai_m[a_idx][i_idx] = gLr_ai_m[a_idx][i_idx] / (g * sqrt(r1_sq));
                L_ai_mbar[a_idx][i_idx] = gLr_ai_mbar[a_idx][i_idx] / (g * sqrt(r2_sq));
            }
        }
    
        double theta1 = atan2(rho_1, z1);
        double theta2 = atan2(rho_2, z2);
    
        double sin_term = 0.0;
        if (abs(theta1) < 1e-7 && abs(theta2) < 1e-7) {
            sin_term = 1.0;
        }else if (abs(theta1 - pi) < 1e-7 && abs(theta2) < 1e-7) {
            sin_term = 1.0;
        }else if (abs(theta1 - pi) < 1e-7 && abs(theta2 - pi) < 1e-7) {
            sin_term = 1.0;
        }else {
            sin_term = sintheta2 / sintheta1;
        }
    
        double gr_div_term_m[3][3] = {
            {-1 * cos(phi1) * (1 - costheta1) * sin(phi1), cos(phi1) * cos(phi1) * (1 - costheta1), 0},
            {-1 * sin(phi1) * (1 - costheta1) * sin(phi1), sin(phi1) * cos(phi1) * (1 - costheta1), 0},
            {0, 0, 0}
        };
    
    
        double div_term_m[3][3] = {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        };
    
        for (int a_idx = 0; a_idx < 3; ++a_idx) {
            for (int i_idx = 0; i_idx < 3; ++i_idx) {
                div_term_m[a_idx][i_idx] = gr_div_term_m[a_idx][i_idx] * sin_term / (g * sqrt(r1_sq));
            }
        }
    
        double Q_H[3][3] = {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        };
        double S_H[3][3] = {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        };
    
        double Q_L[3][3] = {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        };
        double S_L[3][3] = {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        };
    

        auto matMul3 = [](const double A[3][3], const double B[3][3], double out[3][3]) {
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c) {
                out[r][c] = 0.0;
                for (int k_idx = 0; k_idx < 3; ++k_idx)
                    out[r][c] += A[r][k_idx] * B[k_idx][c];
            }
        };

        double PR1[3][3], PR1M[3][3], X1[3][3]; // X1 = P * R_gamma1 * R_mbar * P
        matMul3(P, R_gamma1, PR1);
        matMul3(PR1, R_mbar, PR1M);
        matMul3(PR1M, P, X1);

        double PR2[3][3], PR2M[3][3], X2[3][3]; // X2 = P * R_gamma2 * R_mbar * P
        matMul3(P, R_gamma2, PR2);
        matMul3(PR2, R_mbar, PR2M);
        matMul3(PR2M, P, X2);

        for (int a_idx = 0; a_idx < 3; ++a_idx) {
            for (int i_idx = 0; i_idx < 3; ++i_idx) {
                S_H[a_idx][i_idx] = 0.0;
                S_L[a_idx][i_idx] = 0.0;
                for (int e_idx = 0; e_idx < 3; ++e_idx) {
                    S_H[a_idx][i_idx] += (g / gpp) * L_ai_m[e_idx][i_idx] * X1[e_idx][a_idx];
                    S_L[a_idx][i_idx] += L_ai_m[e_idx][i_idx] * X2[e_idx][a_idx];
                }
            }
        }


        for (int a_idx = 0; a_idx < 3; ++a_idx) {
            for (int i_idx = 0; i_idx < 3; ++i_idx) {
                for (int b_idx = 0; b_idx < 3; ++b_idx) {
                    Q_H[a_idx][i_idx] += (g / gpp) * P[a_idx][b_idx] * L_ai_mbar[b_idx][i_idx];
                    Q_L[a_idx][i_idx] += P[a_idx][b_idx] * L_ai_mbar[b_idx][i_idx];
                    
    
                }
            }
        }
    
        double L_ai[3][3] = {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        };
        double H_ai[3][3] = {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        };
    
        for (int a_idx = 0; a_idx < 3; ++a_idx) {
            for (int i_idx = 0; i_idx < 3; ++i_idx) {
                L_ai[a_idx][i_idx] = L_ai_mbar[a_idx][i_idx];
                H_ai[a_idx][i_idx] = (g / gpp) * L_ai_mbar[a_idx][i_idx];
                for (int b_idx = 0; b_idx < 3; ++b_idx) {
                    for (int c_idx = 0; c_idx < 3; ++c_idx) {
                        L_ai[a_idx][i_idx] += L_ai_m[c_idx][i_idx] * R_gamma2[c_idx][b_idx] * R_mbar[b_idx][a_idx];
                        H_ai[a_idx][i_idx] += (g / gpp) * L_ai_m[c_idx][i_idx] * R_gamma1[c_idx][b_idx] * R_mbar[b_idx][a_idx];
                    }
                }
            }
        }
    
        double A_i[3] = {0.0, 0.0, 0.0};
    
        for (int i_idx = 0; i_idx < 3; ++i_idx) {
            for (int a_idx = 0; a_idx < 3; ++a_idx) {
                A_i[i_idx] += (-1.0/2.0) * (L_ai[a_idx][i_idx] * zhat[a_idx] - (gpp / g) * H_ai[a_idx][i_idx] * zhat[a_idx]);
            }
        }
    
    
        double R_1[3][3] = {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        };
    
        double R_2[3][3] = {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        };

        double Rm1[3][3], Rm2[3][3];
        matMul3(R_m, R_gamma1, Rm1);
        matMul3(Rm1, R_mbar, R_1);
        matMul3(R_m, R_gamma2, Rm2);
        matMul3(Rm2, R_mbar, R_2);
            
        // for (int a_idx = 0; a_idx < 3; ++a_idx) {
        //     for (int b_idx = 0; b_idx < 3; ++b_idx) {
        //         for (int c_idx = 0; c_idx < 3; ++c_idx) {
        //             for (int d_idx = 0; d_idx < 3; ++d_idx) {
        //                 R_1[a_idx][b_idx] += R_m[a_idx][d_idx] * R_gamma1[d_idx][c_idx] * R_mbar[c_idx][b_idx];
        //                 R_2[a_idx][b_idx] += R_m[a_idx][d_idx] * R_gamma2[d_idx][c_idx] * R_mbar[c_idx][b_idx];
        //             }
        //         }
        //     }
        // }
    
    
        double V_local[3][3] = {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        };
    
        double W_local[3][3] = {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        };
    
        for (int a_idx = 0; a_idx < 3; ++a_idx) {
            for (int i_idx = 0; i_idx < 3; ++i_idx) {
                for (int b_idx = 0; b_idx < 3; ++b_idx) {
                    V_local[a_idx][i_idx] += R_1[a_idx][b_idx] * (-1 * (g / gpp) * A_i[i_idx] * zhat[b_idx] - Q_H[b_idx][i_idx] - S_H[b_idx][i_idx] - (g / gpp) * div_term_m[b_idx][i_idx]);
                    W_local[a_idx][i_idx] += R_2[a_idx][b_idx] * ( A_i[i_idx] * zhat[b_idx] - Q_L[b_idx][i_idx] - S_L[b_idx][i_idx] - div_term_m[b_idx][i_idx]);
                }
            }
        }
    
    
        for (int a_idx = 0; a_idx < 3; ++a_idx) {
            for (int i_idx = 0; i_idx < 3; ++i_idx) {
                V_ai[a_idx][i_idx] = V_inf_amp * V_profile * V_local[a_idx][i_idx];
                W_ai[a_idx][i_idx] = W_inf_amp * W_profile * W_local[a_idx][i_idx];
            }
        }
    }
}



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

    // Read profile filenames.
    std::getline(ifs, description, ':');
    ifs >> this->scalarProfileDataFilename;

    std::getline(ifs, description, ':');
    ifs >> this->gaugeProfileDataFilename;

    std::getline(ifs, description, ':');
    ifs >> this->sorProfileDataFilename;

    // Read monopole positions.
    std::getline(ifs, description, ':');
    ifs >> this->monopole1X;

    std::getline(ifs, description, ':');
    ifs >> this->monopole1Y;

    std::getline(ifs, description, ':');
    ifs >> this->monopole1Z;

    std::getline(ifs, description, ':');
    ifs >> this->monopole2X;

    std::getline(ifs, description, ':');
    ifs >> this->monopole2Y;

    std::getline(ifs, description, ':');
    ifs >> this->monopole2Z;

    // Read asymptotic field amplitudes.
    std::getline(ifs, description, ':');
    ifs >> this->WInfAmp;

    std::getline(ifs, description, ':');
    ifs >> this->VInfAmp;

    // Read gamma parameters.
    std::getline(ifs, description, ':');
    ifs >> this->gamma1;

    std::getline(ifs, description, ':');
    ifs >> this->gamma2;

    // Read couplings.
    std::getline(ifs, description, ':');
    ifs >> this->g;

    std::getline(ifs, description, ':');
    ifs >> this->gpp;

    // Read monopole grid spacing.
    std::getline(ifs, description, ':');
    ifs >> this->monopoleGridSpacing;

    if (!ifs)
    {
        throw std::runtime_error(
            "INITIALCONDITIONS::USERDEFINEDPROFILE:: "
            "Could not read the profile configuration."
        );
    }

    ifs.close();

    /*
     * Select the appropriate profile-data file according to the field type.
     */if (this->sorProfileDataFilename != "NONE")
    {
        const std::string data_path
            = std::string(SOURCE_DIR) + "/" + this->sorProfileDataFilename;

        std::ifstream profileFile(data_path);

        if (!profileFile.is_open())
        {
            throw std::runtime_error(
                "INITIALCONDITIONS::USERDEFINEDPROFILE:: "
                "Could not open profile file: " + data_path
            );
        }

        double kValue;
        double kPrimeValue;
        double hWValue;
        double hVValue;

        while (profileFile >> kValue
                        >> kPrimeValue
                        >> hWValue
                        >> hVValue)
        {
            this->k.push_back(kValue);
            this->kPrime.push_back(kPrimeValue);
            this->hW.push_back(hWValue);
            this->hV.push_back(hVValue);
        }

        if (this->k.empty()
            || this->kPrime.size() != this->k.size()
            || this->hW.size() != this->k.size()
            || this->hV.size() != this->k.size())
        {
            throw std::runtime_error(
                "INITIALCONDITIONS::USERDEFINEDPROFILE:: "
                "Invalid or empty profile file: " + data_path
            );
        }
    }

    /*
     * SOR profile is separate from the normal scalar/gauge profile files.
     * Load/use sorProfileDataFilename here when you implement the SOR fields.
     */
}   /*
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
     */
    const double x = (static_cast<double>(global_x) - 0.5*static_cast<double>(geometry.globalNx - 1U))*geometry.dx;
    
    const double y = (static_cast<double>(global_y) - 0.5*static_cast<double>(geometry.globalNy - 1U))*geometry.dy;
    
    const double z = (static_cast<double>(global_z) - 0.5*static_cast<double>(geometry.globalNz - 1U)) *geometry.dz;

    const double x1 = x - this->monopole1X;
    const double y1 = y - this->monopole1Y;
    const double z1 = z - this->monopole1Z;
    
    const double x2 = x - this->monopole2X;
    const double y2 = y - this->monopole2Y;
    const double z2 = z - this->monopole2Z;

    const double r_1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    const double r_2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);

    const double rPos1 = r_1 / this->monopoleGridSpacing;
    const double rPos2 = r_2 / this->monopoleGridSpacing;

    const double k_1 = interpProfile(this->k, rPos1, 1.0);
    const double k_1_p = interpProfile(this->kPrime, rPos1, 0.0);


    const double k_2 = interpProfile(this->k, rPos2, 1.0);
    const double k_2_p = interpProfile(this->kPrime, rPos2, 0.0);

    const double g_1_p = (k_1 + k_1_p);
    const double g_1 = (k_1 - k_1_p);
    const double g_2_p = (k_2 + k_2_p);
    const double g_2 = (k_2 - k_2_p);

    complex<double> phi[4];

    const double phi_both[4] = {(-g_1_p * g_2_p), (g_1 * g_2), (-g_1 * g_2), (g_1_p * g_2_p)};

    complex<double> u_1[2][2];
    complex<double> u_2[2][2];
    build_u_matrices(x1, y1, z1, r_1, x2, y2, z2, r_2, u_1, u_2);
    
    buildScalarField(u_1, u_2, phi_both, phi, this->gamma1, this->gamma2);


        
    for (unsigned comp = 0; comp < 4U; ++comp)
    {
        current_fields[2*comp]     = static_cast<float>(phi[comp].real());
        current_fields[2*comp + 1] = static_cast<float>(phi[comp].imag());
        previous_fields[2*comp]     = current_fields[2*comp];
        previous_fields[2*comp + 1] = current_fields[2*comp + 1];
    }
     /*
     * 2. Evaluate an analytic profile, or interpolate data loaded once in
     *    configure().
     *
     * 3. Assign current_fields[component] for every required scalar
     *    component. The meaning of each component is fixed by the chosen model.
     *
     * 4. For zero initial time derivative, copy the current values:
     */

     
     /*      To encode a non-zero initial derivative v using the backward time
     *      level, a common finite-difference choice is
     *
     *      previous_fields[comp] = current_fields[comp] - geometry.dt*v;
     *
     *      Check that convention against the intended discretised equations.
     */

    // Remove these casts when the corresponding arguments are used.
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
    const double x =
        (static_cast<double>(global_x)
        - 0.5 * static_cast<double>(geometry.globalNx - 1U))
        * geometry.dx;

    const double y =
        (static_cast<double>(global_y)
        - 0.5 * static_cast<double>(geometry.globalNy - 1U))
        * geometry.dy;

    const double z =
        (static_cast<double>(global_z)
        - 0.5 * static_cast<double>(geometry.globalNz - 1U))
        * geometry.dz;


    const double x1 = x - this->monopole1X;
    const double y1 = y - this->monopole1Y;
    const double z1 = z - this->monopole1Z;
    
    const double x2 = x - this->monopole2X;
    const double y2 = y - this->monopole2Y;
    const double z2 = z - this->monopole2Z;

    
    const double r_1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);
    const double r_2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);



    const double rPos1 = r_1 / this->monopoleGridSpacing;
    const double rPos2 = r_2 / this->monopoleGridSpacing;

    const double W1 = interpProfile(this->hW, rPos1, 1.0);
    const double W2 = interpProfile(this->hW, rPos2, 1.0);

    const double V1 = interpProfile(this->hV, rPos1, 1.0);
    const double V2 = interpProfile(this->hV, rPos2, 1.0);


    const double WProfile = W1 * W2;
    const double VProfile = V1 * V2;


    // 4. Construct continuum gauge fields
    double V_ai[3][3];
    double W_ai[3][3];


    buildGaugeFields(
    x, y, z,
    x1, y1, z1,
    x2, y2, z2,
    this->gamma1, this->gamma2,
    VProfile, WProfile,
    this->g, this->gpp,
    this->VInfAmp, this->WInfAmp,
    V_ai, W_ai
    );



    double final_w_ai[3][3];
    double final_v_ai[3][3];
    double d_var[3] = {geometry.dx, geometry.dy, geometry.dz};
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            final_w_ai[i][j] = 1.0/2.0 * g * d_var[j] * W_ai[i][j];
            final_v_ai[i][j] = 1.0/2.0 * gpp * d_var[j] * V_ai[i][j];
        }
    }

    double final_y_ai[3]= {0.0, 0.0, 0.0};

    const unsigned direction_width = num_components / 3U;

    for (unsigned dir = 0; dir < 3; ++dir)
    {
        // Hypercharge U(1)
        current_fields[dir * direction_width + 0] = final_y_ai[dir];

        // SM SU(2)
        for (unsigned a = 0; a < 3; ++a)
            current_fields[dir * direction_width + 1 + a] = -final_w_ai[a][dir];

        // Higgs-family SU(2)
        for (unsigned a = 0; a < 3; ++a)
            current_fields[dir * direction_width + 4 + a] = -final_v_ai[a][dir];
    }



    //     for (int i = 0; i < 3; ++i)
    // {
    //     for (int a = 0; a < 3; ++a)
    //     {
    //         current_fields[i * 3 + a]
    //             = static_cast<float>(final_w_ai[a][i]);
    //     }
    // }
    // 5. Convert V_ai/W_ai into whatever
    //    GFTEvolver expects in current_fields[]


    // 6. Zero initial electric field
    for (unsigned comp = 0; comp < num_components; ++comp)
    {
        previous_fields[comp] = current_fields[comp];
    }
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
