#include <stdio.h>
#include "nbody_defaults.h"
#include "nbody_types.h"

int main() {
    int n = 11;  // number of dwarfs
    printf("=== DEBUG: Expanded NBody Defaults ===\n");

    printf("B start:\n");
    for (int i = 0; i < n; i++) {
        printf("%d: %f\n", i, i < 2 ? defaultNBodyCtx.b[i] : -1.0);
    }

    printf("R start:\n");
    for (int i = 0; i < n; i++) {
        printf("%d: %f\n", i, i < 2 ? defaultNBodyCtx.r[i] : -1.0);
    }

    printf("VX start:\n");
    for (int i = 0; i < n; i++) {
        printf("%d: %f\n", i, i < 2 ? defaultNBodyCtx.vx[i] : -1.0);
    }

    printf("VY start:\n");
    for (int i = 0; i < n; i++) {
        printf("%d: %f\n", i, i < 2 ? defaultNBodyCtx.vy[i] : -1.0);
    }

    printf("VZ start:\n");
    for (int i = 0; i < n; i++) {
        printf("%d: %f\n", i, i < 2 ? defaultNBodyCtx.vz[i] : -1.0);
    }

    printf("=======================================\n");
    return 0;
}