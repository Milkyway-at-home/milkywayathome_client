#include "bfe.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "basis.h"
#include "milkyway/milkyway_math.h"

// --- Forward Declarations ---
static inline int get_coeff_index(int n, int l, int m, int nmax, int lmax);
static void bfe_compute_radial_basis(int n, int l, double s, double* phi_nl_out, double* dphi_nl_ds_out);
static void bfe_sum_potential_and_derivatives(const BFEModel* model, double s,
                                              const double* Plm, const double* dPlm_dtheta,
                                              const double* cos_m_phi, const double* sin_m_phi,
                                              double* Phi, double* dPhi_ds, double* dPhi_dtheta, double* dPhi_dphi);
static void bfe_convert_force_to_cartesian(double r, double cos_theta, double sin_theta,
                                           double cos_phi, double sin_phi,
                                           double dPhi_dr, double dPhi_dtheta, double dPhi_dphi,
                                           double* force_out);

// --- File I/O and Lifecycle Functions (Unchanged) ---
// bfe_create_from_file, bfe_destroy, bfe_evolve_coeffs are identical to the previous version.

BFEModel* bfe_create_from_file(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        perror("Error opening BFE coefficient file");
        return NULL;
    }
    BFEModel* model = (BFEModel*)malloc(sizeof(BFEModel));
    if (!model) {
        fclose(file);
        return NULL;
    }
    if (fscanf(file, "# nmax lmax scale_radius\n%d %d %lf\n", &model->nmax, &model->lmax, &model->scale_radius) != 3) {
        fprintf(stderr, "Invalid BFE coefficient file header format.\n");
        fclose(file);
        free(model);
        return NULL;
    }
    int num_coeffs = (model->nmax + 1) * (model->lmax + 1) * (model->lmax + 1);
    model->S_coeffs = (double*)calloc(num_coeffs, sizeof(double));
    model->T_coeffs = (double*)calloc(num_coeffs, sizeof(double));
    if (!model->S_coeffs || !model->T_coeffs) {
        fprintf(stderr, "Failed to allocate memory for BFE coefficients\n");
        fclose(file);
        free(model->S_coeffs);
        free(model);
        return NULL;
    }
    int n, l, m;
    double s_val, t_val;
    char line_buffer[256];
    fgets(line_buffer, sizeof(line_buffer), file);
    while (fgets(line_buffer, sizeof(line_buffer), file) != NULL) {
        if (line_buffer[0] == '#' || strlen(line_buffer) < 5) continue;
        if (sscanf(line_buffer, "%d %d %d %lf %lf", &n, &l, &m, &s_val, &t_val) == 5) {
            int index = get_coeff_index(n, l, m, model->nmax, model->lmax);
            model->S_coeffs[index] = s_val;
            model->T_coeffs[index] = t_val;
        }
    }
    fclose(file);
    return model;
}

void bfe_destroy(BFEModel* model) {
    if (!model) return;
    free(model->S_coeffs);
    free(model->T_coeffs);
    free(model);
}

void bfe_evolve_coeffs(BFEModel* model, const Particle particles[], int num_particles, double dt) {
    (void)model; (void)particles; (void)num_particles; (void)dt;
    return;
}

// --- Core Physics Functions ---

void bfe_calculate_potential_and_force(double pos[3], const BFEModel* model, double* potential_out, double* force_out) {
    *potential_out = 0.0;
    force_out[0] = 0.0; force_out[1] = 0.0; force_out[2] = 0.0;
    if (!model) return;

    double r = mw_sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
    if (r < 1e-9) return;

    double phi = mw_atan2(pos[1], pos[0]);
    double cos_theta = pos[2] / r;
    if (cos_theta > 1.0) cos_theta = 1.0;
    if (cos_theta < -1.0) cos_theta = -1.0;
    double s = r / model->scale_radius;

    int lmax = model->lmax;
    double* cos_m_phi = malloc((lmax + 1) * sizeof(double));
    double* sin_m_phi = malloc((lmax + 1) * sizeof(double));
    int legendre_size = (lmax + 1) * (lmax + 2) / 2;
    double* Plm = malloc(legendre_size * sizeof(double));
    double* dPlm_dtheta = malloc(legendre_size * sizeof(double));
    if (!cos_m_phi || !sin_m_phi || !Plm || !dPlm_dtheta) {
        free(cos_m_phi); free(sin_m_phi); free(Plm); free(dPlm_dtheta);
        return;
    }

    basis_sincos(lmax, phi, cos_m_phi, sin_m_phi);
    basis_legendre_deriv(lmax, cos_theta, Plm, dPlm_dtheta);

    double Phi = 0.0, dPhi_ds = 0.0, dPhi_dtheta = 0.0, dPhi_dphi = 0.0;
    bfe_sum_potential_and_derivatives(model, s, Plm, dPlm_dtheta, cos_m_phi, sin_m_phi,
                                      &Phi, &dPhi_ds, &dPhi_dtheta, &dPhi_dphi);
    
    free(cos_m_phi); free(sin_m_phi); free(Plm); free(dPlm_dtheta);

    const double G_CONST = 4.30091e-6;
    const double GALAXY_MASS = 1.0e11;
    const double POTENTIAL_SCALE = G_CONST * GALAXY_MASS;

    *potential_out = Phi * POTENTIAL_SCALE;

    double dPhi_dr = dPhi_ds / model->scale_radius;
    dPhi_dr     *= POTENTIAL_SCALE;
    dPhi_dtheta *= POTENTIAL_SCALE;
    dPhi_dphi   *= POTENTIAL_SCALE;

    double sin_theta = mw_sqrt(1.0 - cos_theta*cos_theta);
    bfe_convert_force_to_cartesian(r, cos_theta, sin_theta, pos[0]/r*sin_theta, pos[1]/r*sin_theta,
                                   dPhi_dr, dPhi_dtheta, dPhi_dphi, force_out);
}

// Wrapper function for backward compatibility
void bfe_calculate_force(double pos[3], const BFEModel* model, double* force_out, 
                         double* cos_m_phi, double* sin_m_phi, double* Plm, double* dPlm_dtheta) {
    double potential; // Dummy variable
    bfe_calculate_potential_and_force(pos, model, &potential, force_out);
}

// --- Static Helper Function Implementations ---

static inline int get_coeff_index(int n, int l, int m, int nmax, int lmax) {
    return n * (lmax + 1) * (lmax + 1) + l * (lmax + 1) + m;
}

static void bfe_compute_radial_basis(int n, int l, double s, double* phi_nl_out, double* dphi_nl_ds_out) {
    double alpha = 2.0 * l + 1.5;
    double xi = (s - 1.0) / (s + 1.0);
    double s_plus_1 = s + 1.0;
    double C_n_values[n + 2];
    basis_gegenbauer(n + 1, alpha, xi, C_n_values);
    double C_n = C_n_values[n];
    double K_nl;
    if (n == 0) {
        K_nl = sqrt( (2.0*l + 1.5) * (l + 1.0) / 2.0 );
    } else {
        K_nl = sqrt( (n + 2.0*l + 1.5) * (n + l + 1.0) * n / 2.0 );
    }
    K_nl *= 2.0 / sqrt(M_PI);
    double envelope = mw_pow(s, l) * mw_pow(s_plus_1, -2.0 * l - 1.0);
    *phi_nl_out = -K_nl * envelope * C_n;
    double C_n_plus_1 = C_n_values[n+1];
    double dphi_term = (n + 2.0 * l + 1.5) * C_n - (n + 1.0) * C_n_plus_1;
    *dphi_nl_ds_out = K_nl * mw_pow(s, l) * mw_pow(s_plus_1, -2.0*l-2.0) * dphi_term;
}

static void bfe_sum_potential_and_derivatives(const BFEModel* model, double s,
                                          const double* Plm, const double* dPlm_dtheta,
                                          const double* cos_m_phi, const double* sin_m_phi,
                                          double* Phi, double* dPhi_ds, double* dPhi_dtheta, double* dPhi_dphi) {
    int nmax = model->nmax;
    int lmax = model->lmax;

    for (int n = 0; n <= nmax; ++n) {
        for (int l = 0; l <= lmax; ++l) {
            double phi_nl, dphi_nl_ds;
            bfe_compute_radial_basis(n, l, s, &phi_nl, &dphi_nl_ds);

            for (int m = 0; m <= l; ++m) {
                int coeff_idx = get_coeff_index(n, l, m, nmax, lmax);
                double S_nlm = model->S_coeffs[coeff_idx];
                double T_nlm = model->T_coeffs[coeff_idx];

                if (S_nlm == 0 && T_nlm == 0) continue;

                int legendre_idx = l * (l + 1) / 2 + m;
                double ST_term = S_nlm * cos_m_phi[m] + T_nlm * sin_m_phi[m];

                *Phi         += phi_nl * Plm[legendre_idx] * ST_term;
                *dPhi_ds     += dphi_nl_ds * Plm[legendre_idx] * ST_term;
                *dPhi_dtheta += phi_nl * dPlm_dtheta[legendre_idx] * ST_term;
                *dPhi_dphi   += phi_nl * Plm[legendre_idx] * m * (-S_nlm * sin_m_phi[m] + T_nlm * cos_m_phi[m]);
            }
        }
    }
}

static void bfe_convert_force_to_cartesian(double r, double cos_theta, double sin_theta,
                                           double cos_phi, double sin_phi,
                                           double dPhi_dr, double dPhi_dtheta, double dPhi_dphi,
                                           double* force_out) {
    double r_inv = 1.0 / r;
    double sin_theta_inv = (sin_theta > 1e-9) ? 1.0 / sin_theta : 0.0;
    double F_r     = -dPhi_dr;
    double F_theta = -r_inv * dPhi_dtheta;
    double F_phi   = -r_inv * sin_theta_inv * dPhi_dphi;
    force_out[0] = sin_theta * cos_phi * F_r + cos_theta * cos_phi * F_theta - sin_phi * F_phi;
    force_out[1] = sin_theta * sin_phi * F_r + cos_theta * sin_phi * F_theta + cos_phi * F_phi;
    force_out[2] = cos_theta * F_r - sin_theta * F_theta;
}