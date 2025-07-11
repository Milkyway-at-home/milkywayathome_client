#include "basis.h"
#include <math.h>
#include <string.h>
#include "milkyway/milkyway_math.h"

void basis_sincos(int mmax, double phi, double* c, double* s) {
    c[0] = 1.0;
    s[0] = 0.0;
    if (mmax == 0) return;

    c[1] = mw_cos(phi);
    s[1] = mw_sin(phi);

    for (int m = 2; m <= mmax; m++) {
        c[m] = 2.0 * c[1] * c[m-1] - c[m-2];
        s[m] = 2.0 * c[1] * s[m-1] - s[m-2];
    }
}

void basis_legendre_deriv(int lmax, double x, double* p, double* dp) {
    // This function calculates Associated Legendre Polynomials P_l^m(x)
    // and their derivatives with respect to theta, d(P_l^m)/d(theta),
    // where x = cos(theta). It uses a packed storage format for efficiency.

    // Calculate sin(theta) from x=cos(theta), ensuring it's non-negative.
    double sin_theta = mw_sqrt(1.0 - x*x);

    // --- Part 1: Calculate P_l^m(x) using recurrence relations ---
    // Base case P_0^0 = 1
    p[0] = 1.0;
    if (lmax == 0) {
        if (dp != NULL) dp[0] = 0.0;
        return;
    }

    // Calculate diagonal terms P_m^m using recurrence
    // P_m^m = -(2m-1) * sin(theta) * P_{m-1}^{m-1}
    // Note: p_prev_mm holds P_{m-1}^{m-1}
    double p_prev_mm = 1.0;
    for (int m = 1; m <= lmax; m++) {
        int idx = m * (m + 1) / 2 + m;
        p[idx] = -(2.0 * m - 1.0) * sin_theta * p_prev_mm;
        p_prev_mm = p[idx];
    }

    // Calculate P_l^m for m < l using recurrence
    // (l-m)P_l^m = (2l-1)x * P_{l-1}^m - (l+m-1)P_{l-2}^m
    for (int l = 1; l <= lmax; l++) {
        int idx_lm1_m = (l - 1) * l / 2 + l-1; // P_{l-1}^{l-1}
        p[l*(l+1)/2 + l-1] = (2.0*l - 1.0) * x * p[idx_lm1_m];

        for (int m = l - 2; m >= 0; m--) {
             int idx_lm = l * (l + 1) / 2 + m;
             int idx_lm1_m = (l - 1) * l / 2 + m;
             int idx_lm2_m = (l - 2) * (l-1) / 2 + m;
             p[idx_lm] = ((2.0 * l - 1.0) * x * p[idx_lm1_m] - (l + m - 1.0) * p[idx_lm2_m]) / (l - m);
        }
    }

    // --- Part 2: Calculate d(P_l^m)/d(theta) ---
    if (dp == NULL) return;

    // Use the stable recurrence relation:
    // sin(theta) * d(P_l^m)/d(theta) = -l*x*P_l^m + (l+m)*P_{l-1}^m
    double sin_theta_inv = (sin_theta > 1e-9) ? 1.0 / sin_theta : 0.0;

    for (int l = 0; l <= lmax; l++) {
        for (int m = 0; m <= l; m++) {
            int idx_lm = l * (l + 1) / 2 + m;
            int idx_lm1_m = (l > 0) ? (l - 1) * l / 2 + m : 0;

            double term1 = -l * x * p[idx_lm];
            double term2 = (l > 0 && m < l) ? (l + m) * p[idx_lm1_m] : 0.0;
            
            dp[idx_lm] = (term1 + term2) * sin_theta_inv;
        }
    }
}

void basis_gegenbauer(int nmax, double alpha, double x, double* C_out) {
    if (nmax < 0) return;
    C_out[0] = 1.0;
    if (nmax == 0) return;

    C_out[1] = 2.0 * alpha * x;
    if (nmax == 1) return;

    for (int n = 2; n <= nmax; ++n) {
        double term1 = 2.0 * (n + alpha - 1.0) * x * C_out[n - 1];
        double term2 = (n + 2.0 * alpha - 2.0) * C_out[n - 2];
        C_out[n] = (term1 - term2) / n;
    }
}