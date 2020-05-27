/* ==========================================================================
======================= n-Step Gauss-Seidel Solver ==========================
========================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>


int GS_nstep(double** f, double** phi, int Nx, int Ny, double D_x, double D_y, double epsilon, int nGS, int H_cells, int Bj, int Bi) {

    /* Initilaizations */
    int i, j;
    double nx = Nx;
    double ny = Ny;

    double f_norm;
    double integral;
    double lambda = pow(D_x, -2);
    double RHS;    
    double laplace_phi_minus_f_norm;
   
    int step = 1;
    int max_num_steps = nGS;

    double** laplace_phi = (double**)calloc(Ny, sizeof(double*));
    for (j = 0; j < Ny; j++) {
        laplace_phi[j] = (double*)calloc(Nx, sizeof(double));
    }


    /* Solving -------------------------------------------------------------------------------- */

    if (Nx != 1 && Ny != 1) {
        f_norm = 0; /* compute f_norm */
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {

                if ( ((j < Bj) || ( j > (Bj + H_cells - 1))) || ((i < Bi) || (i > (Bi + H_cells - 1)))  ) {     /* if not inside the void... */
                    f_norm = f_norm + pow(f[j][i], 2);
                }

            }
        }
       
        f_norm = sqrt(f_norm);
       


    }


    do {

        /* ================================================= CALCULATING PHI ================================================= */
        if (Nx == 1 || Ny == 1) { /* ===================== SKIP ALL CONDITION ===================== */
           
            /* Skip everything */

        }
        else if (Nx == 2 || Ny == 2) { /* ===================== THIN DOMAIN CONDITION ===================== */
           
            if (Nx == 2 && Ny == 2) { /* ======= If domain is a 2x2 square ======= */
                /* Only corner calculations are needed */

                /* 1st corner calculations */
                if ( ((0 >= Bj) && (0 <= (Bj + H_cells - 1))) && ((0 >= Bi) && (0 <= (Bi + H_cells - 1)))  ) {
                    /* Inside void - skip this position */
                }
                else if (0 == Bi - 1 && 0 >= Bj && 0 <= Bj + H_cells - 1) {
                    /* At left void boundary (to the left of the void) */
                    phi[0][0] = phi[1][0] - f[0][0] / lambda;
                }
                else if (0 == Bj - 1 && 0 >= Bi && 0 <= Bi + H_cells - 1) {
                    /* At bottom void boundary (beneath the void) */
                    phi[0][0] = phi[0][1] - f[0][0] / lambda;
                }
                else {
                    /* Not in void or at void boundary - calculate normally */
                    phi[0][0] = (phi[1][0] + phi[0][1]) / 2 - f[0][0] / (2 * lambda);

                }

                /* 2nd corner calculations */
                if ( ((0 >= Bj) && (0 <= (Bj + H_cells - 1))) && ((Nx - 1 >= Bi) && (Nx - 1 <= (Bi + H_cells - 1)))  ) {
                    /* Inside void - skip this position */
                }
                else if (Nx - 1 == Bi + H_cells && 0 >= Bj && 0 <= Bj + H_cells - 1) {
                    /* At right void boundary */
                    phi[0][Nx - 1] = phi[1][Nx - 1] - f[0][Nx - 1] / lambda;
                }
                else if (0 == Bj - 1 && Nx - 1 >= Bi && Nx - 1 <= Bi + H_cells - 1) {
                    /* At bottom void boundary */
                    phi[0][Nx - 1] = phi[0][Nx - 2] - f[0][Nx - 1] / lambda;
                }
                else {
                    /* Not in void or at void boundary - calculate normally */
                    phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) / 2 - f[0][Nx - 1] / (2 * lambda);

                }

                /* 3rd corner calculations */
                if ( ((Ny - 1 >= Bj) && (Ny - 1 <= (Bj + H_cells - 1))) && ((0 >= Bi) && (0 <= (Bi + H_cells - 1)))  ) {
                    /* Inside void - skip this position */
                }
                else if (0 == Bi - 1 && Ny - 1 >= Bj && Ny - 1 <= Bj + H_cells - 1) {
                    /* At left void boundary */
                    phi[Ny - 1][0] = phi[Ny - 2][0] - f[Ny - 1][0] / lambda;
                }
                else if (Ny - 1 == Bj + H_cells && 0 >= Bi && 0 <= Bi + H_cells - 1) {
                    /* At top void boundary */
                    phi[Ny - 1][0] = phi[Ny - 1][1] - f[Ny - 1][0] / lambda;
                }
                else {
                    /* Not in void or at void boundary - calculate normally */
                    phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) / 2 - f[Ny - 1][0] / (2 * lambda);

                }

                /* 4th corner calculations */
                if ( ((Ny - 1 >= Bj) && (Ny - 1 <= (Bj + H_cells - 1))) && ((Nx - 1 >= Bi) && (Nx - 1 <= (Bi + H_cells - 1)))  ) {
                    /* Inside void - skip this position */
                }
                else if (Nx - 1 == Bi + H_cells && Ny - 1 >= Bj && Ny - 1 <= Bj + H_cells - 1) {
                    /* At right void boundary */
                    phi[Ny - 1][Nx - 1] = phi[Ny - 2][Nx - 1] - f[Ny - 1][Nx - 1] / lambda;
                }
                else if (Ny - 1 == Bj + H_cells && Nx - 1 >= Bi && Nx - 1 <= Bi + H_cells - 1) {
                    /* At top void boundary */
                    phi[Ny - 1][Nx - 1] = phi[Ny - 1][Nx - 2] - f[Ny - 1][Nx - 1] / lambda;
                }
                else {
                    /* Not in void or at void boundary - calculate normally */
                    phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) / 2 - f[Ny - 1][Nx - 1] / (2 * lambda);

                }

            }
            else { /* ======= If domain is not square ======= */
                /* Loop through entire space */
                for (j = 0; j < Ny; j++) {
                    for (i = 0; i < Nx; i++) {

                        if (i == 0 && j == 0) {
                            /* 1st corner calculations */
                            if ( ((0 >= Bj) && (0 <= (Bj + H_cells - 1))) && ((0 >= Bi) && (0 <= (Bi + H_cells - 1)))  ) {
                                /* Inside void - skip this position */
                            }
                            else if (0 == Bi - 1 && 0 >= Bj && 0 <= Bj + H_cells - 1) {
                                /* At left void boundary (to the left of the void) */
                                phi[0][0] = phi[1][0] - f[0][0] / lambda;
                            }
                            else if (0 == Bj - 1 && 0 >= Bi && 0 <= Bi + H_cells - 1) {
                                /* At bottom void boundary (beneath the void) */
                                phi[0][0] = phi[0][1] - f[0][0] / lambda;
                            }
                            else {
                                /* Not in void or at void boundary - calculate normally */
                                phi[0][0] = (phi[1][0] + phi[0][1]) / 2 - f[0][0] / (2 * lambda);

                            }

                        }
                        else if (i == Nx - 1 && j == 0) {
                            /* 2nd corner calculations */
                            if ( ((0 >= Bj) && (0 <= (Bj + H_cells - 1))) && ((Nx - 1 >= Bi) && (Nx - 1 <= (Bi + H_cells - 1)))  ) {
                                /* Inside void - skip this position */
                            }
                            else if (Nx - 1 == Bi + H_cells && 0 >= Bj && 0 <= Bj + H_cells - 1) {
                                /* At right void boundary */
                                phi[0][Nx - 1] = phi[1][Nx - 1] - f[0][Nx - 1] / lambda;
                            }
                            else if (0 == Bj - 1 && Nx - 1 >= Bi && Nx - 1 <= Bi + H_cells - 1) {
                                /* At bottom void boundary */
                                phi[0][Nx - 1] = phi[0][Nx - 2] - f[0][Nx - 1] / lambda;
                            }
                            else {
                                /* Not in void or at void boundary - calculate normally */
                                phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) / 2 - f[0][Nx - 1] / (2 * lambda);

                            }

                        }
                        else if (i == 0 && j == Ny - 1) {
                            /* 3rd corner calculations */
                            if ( ((Ny - 1 >= Bj) && (Ny - 1 <= (Bj + H_cells - 1))) && ((0 >= Bi) && (0 <= (Bi + H_cells - 1)))  ) {
                                /* Inside void - skip this position */
                            }
                            else if (0 == Bi - 1 && Ny - 1 >= Bj && Ny - 1 <= Bj + H_cells - 1) {
                                /* At left void boundary */
                                phi[Ny - 1][0] = phi[Ny - 2][0] - f[Ny - 1][0] / lambda;
                            }
                            else if (Ny - 1 == Bj + H_cells && 0 >= Bi && 0 <= Bi + H_cells - 1) {
                                /* At top void boundary */
                                phi[Ny - 1][0] = phi[Ny - 1][1] - f[Ny - 1][0] / lambda;
                            }
                            else {
                                /* Not in void or at void boundary - calculate normally */
                                phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) / 2 - f[Ny - 1][0] / (2 * lambda);

                            }

                        }
                        else if (i == Nx - 1 && j ==  Ny - 1) {
                            /* 4th corner calculations */
                            if ( ((Ny - 1 >= Bj) && (Ny - 1 <= (Bj + H_cells - 1))) && ((Nx - 1 >= Bi) && (Nx - 1 <= (Bi + H_cells - 1))) ) {
                                /* Inside void - skip this position */
                            }
                            else if (Nx - 1 == Bi + H_cells && Ny - 1 >= Bj && Ny - 1 <= Bj + H_cells - 1) {
                                /* At right void boundary */
                                phi[Ny - 1][Nx - 1] = phi[Ny - 2][Nx - 1] - f[Ny - 1][Nx - 1] / lambda;
                            }
                            else if (Ny - 1 == Bj + H_cells && Nx - 1 >= Bi && Nx - 1 <= Bi + H_cells - 1) {
                                /* At top void boundary */
                                phi[Ny - 1][Nx - 1] = phi[Ny - 1][Nx - 2] - f[Ny - 1][Nx - 1] / lambda;
                            }
                            else {
                                /* Not in void or at void boundary - calculate normally */
                                phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) / 2 - f[Ny - 1][Nx - 1] / (2 * lambda);

                            }

                        }
                        else if (Nx == 2) { /* ======= For a tall domain ======= */
                            /* Edge calculations for Nx = 2 - there is no iterior */
                            if (i == 0) { /* Left */
                                if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                                    /* Inside void - skip this position */
                                }
                                else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells - 1) {
                                    /* At bottom void boundary */
                                    phi[j][i] = (phi[j][i + 1] + phi[j - 1][i]) / 2 - f[j][i] / (2 * lambda);
                                }
                                else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells - 1) {
                                    /* At left void boundary */
                                    phi[j][i] = (phi[j - 1][i] + phi[j + 1][i]) / 2 - f[j][i] / (2 * lambda);
                                }
                                else if (j == Bj + H_cells && i >= Bi && i <= Bi + H_cells - 1) {
                                    /* At top void boundary */
                                    phi[j][i] = (phi[j][i + 1] + phi[j + 1][i]) / 2 - f[j][i] / (2 * lambda);
                                }
                                else {
                                    /* Not in void or at void boundary - calculate normally */
                                    phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

                                }
                            }
                            else if (i == Nx - 1) { /* Right */
                                if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                                    /* Inside void - skip this position */
                                }
                                else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells - 1) {
                                    /* At bottom void boundary */
                                    phi[j][i] = (phi[j][i - 1] + phi[j - 1][i]) / 2 - f[j][i] / (2 * lambda);
                                }
                                else if (i == Bi + H_cells && j >= Bj && j <= Bj + H_cells - 1) {
                                    /* At right void boundary */
                                    phi[j][i] = (phi[j - 1][i] + phi[j + 1][i]) / 2 - f[j][i] / (2 * lambda);
                                }
                                else if (j == Bj + H_cells && i >= Bi && i <= Bi + H_cells - 1) {
                                    /* At top void boundary */
                                    phi[j][i] = (phi[j][i - 1] + phi[j + 1][i]) / 2 - f[j][i] / (2 * lambda);
                                }
                                else {
                                    /* Not in void or at void boundary - calculate normally */
                                    phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

                                }
                               
                            }

                        }
                        else if (Ny == 2) { /* ======= For a flat domain ======= */
                            /* Edge calculations for Ny = 2 - there is no iterior */
                            if (j == 0) { /* Bottom */
                                if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                                    /* Inside void - skip this position */
                                }
                                else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells - 1) {
                                    /* At bottom void boundary */
                                    phi[j][i] = (phi[j][i + 1] + phi[j][i - 1]) / 2 - f[j][i] / (2 * lambda);
                                }
                                else if (i == Bi + H_cells && j >= Bj && j <= Bj + H_cells - 1) {
                                    /* At right void boundary */
                                    phi[j][i] = (phi[j][i + 1] + phi[j + 1][i]) / 2 - f[j][i] / (2 * lambda);
                                }
                                else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells - 1) {
                                    /* At left void boundary */
                                    phi[j][i] = (phi[j][i - 1] + phi[j + 1][i]) / 2 - f[j][i] / (2 * lambda);
                                }
                                else {
                                    /* Not in void or at void boundary - calculate normally */
                                    phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

                                }
                           
                            }
                            else if (j == Ny - 1) { /* Top */
                                if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                                    /* Inside void - skip this position */
                                }
                                else if (j == Bj + H_cells && i >= Bi && i <= Bi + H_cells - 1) {
                                    /* At top void boundary */
                                    phi[j][i] = (phi[j][i + 1] + phi[j][i - 1]) / 2 - f[j][i] / (2 * lambda);
                                }
                                else if (i == Bi + H_cells && j >= Bj && j <= Bj + H_cells - 1) {
                                    /* At right void boundary */
                                    phi[j][i] = (phi[j][i + 1] + phi[j - 1][i]) / 2 - f[j][i] / (2 * lambda);
                                }
                                else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells - 1) {
                                    /* At left void boundary */
                                    phi[j][i] = (phi[j][i - 1] + phi[j - 1][i]) / 2 - f[j][i] / (2 * lambda);
                                }
                                else {
                                    /* Not in void or at void boundary - calculate normally */
                                    phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) / 3 - f[j][i] / (3 * lambda);

                                }
                           
                            }

                        }

                    }
                }
            }
        }
        else if (Nx >= 3 && Ny >= 3) { /* ===================== NORMAL DOMAIN CONDITION ===================== */

            /* Loop through entire space */
            for (j = 0; j < Ny; j++) {
                for (i = 0; i < Nx; i++) {

                    if ((i > 0 && i < Nx - 1) && (j > 0 && j < Ny - 1)) {
                        /* Interior calculations */
                        if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                            /* Inside void - skip this position */
                        }
                        else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At left void boundary */
                            phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                        }
                        else if (i == Bi + H_cells && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At right void boundary */
                            phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                        }
                        else if (j == Bj + H_cells && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At top void boundary */
                            phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                        }
                        else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At bottom void boundary */
                            phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) / 3 - f[j][i] / (3 * lambda);
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 4 - f[j][i] / (4 * lambda);
                       
                        }

                    }
                    else if (j == 0 && i != 0 && i != Nx - 1) {
                        /* Edge calcualtions on bottom */
                        if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                            /* Inside void - skip this position */
                        }
                        else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At bottom void boundary */
                            phi[j][i] = (phi[j][i + 1] + phi[j][i - 1]) / 2 - f[j][i] / (2 * lambda);
                        }
                        else if (i == Bi + H_cells && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At right void boundary */
                            phi[j][i] = (phi[j][i + 1] + phi[j + 1][i]) / 2 - f[j][i] / (2 * lambda);
                        }
                        else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At left void boundary */
                            phi[j][i] = (phi[j][i - 1] + phi[j + 1][i]) / 2 - f[j][i] / (2 * lambda);
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                            
                        }

                    }
                    else if (j == Ny - 1 && i != 0 && i != Nx - 1) {
                        /* Edge calcualtions on top */
                        if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                            /* Inside void - skip this position */
                        }
                        else if (j == Bj + H_cells && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At top void boundary */
                            phi[j][i] = (phi[j][i + 1] + phi[j][i - 1]) / 2 - f[j][i] / (2 * lambda);
                        }
                        else if (i == Bi + H_cells && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At right void boundary */
                            phi[j][i] = (phi[j][i + 1] + phi[j - 1][i]) / 2 - f[j][i] / (2 * lambda);
                        }
                        else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At left void boundary */
                            phi[j][i] = (phi[j][i - 1] + phi[j - 1][i]) / 2 - f[j][i] / (2 * lambda);
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) / 3 - f[j][i] / (3 * lambda);
                            
                        }

                    }
                    else if (i == 0 && j != 0 && j != Ny - 1) {
                        /* Edge calcualtions on left */
                        if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                            /* Inside void - skip this position */
                        }
                        else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At bottom void boundary */
                            phi[j][i] = (phi[j][i + 1] + phi[j - 1][i]) / 2 - f[j][i] / (2 * lambda);
                        }
                        else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At left void boundary */
                            phi[j][i] = (phi[j - 1][i] + phi[j + 1][i]) / 2 - f[j][i] / (2 * lambda);
                        }
                        else if (j == Bj + H_cells && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At top void boundary */
                            phi[j][i] = (phi[j][i + 1] + phi[j + 1][i]) / 2 - f[j][i] / (2 * lambda);
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                            
                           
                        }

                    }
                    else if (i == Nx - 1 && j != 0 && j != Ny - 1) {
                        /* Edge calcualtions on right */
                        if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                            /* Inside void - skip this position */
                        }
                        else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At bottom void boundary */
                            phi[j][i] = (phi[j][i - 1] + phi[j - 1][i]) / 2 - f[j][i] / (2 * lambda);
                        }
                        else if (i == Bi + H_cells && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At right void boundary */
                            phi[j][i] = (phi[j - 1][i] + phi[j + 1][i]) / 2 - f[j][i] / (2 * lambda);
                        }
                        else if (j == Bj + H_cells && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At top void boundary */
                            phi[j][i] = (phi[j][i - 1] + phi[j + 1][i]) / 2 - f[j][i] / (2 * lambda);
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                            
                        }

                    }
                    else if (i == 0 && j == 0) {
                        /* 1st corner calculation */
                        if ( ((0 >= Bj) && (0 <= (Bj + H_cells - 1))) && ((0 >= Bi) && (0 <= (Bi + H_cells - 1)))  ) {
                            /* Inside void - skip this position */
                        }
                        else if (0 == Bi - 1 && 0 >= Bj && 0 <= Bj + H_cells - 1) {
                            /* At left void boundary (to the left of the void) */
                            phi[0][0] = phi[1][0] - f[0][0] / lambda;
                        }
                        else if (0 == Bj - 1 && 0 >= Bi && 0 <= Bi + H_cells - 1) {
                            /* At bottom void boundary (beneath the void) */
                            phi[0][0] = phi[0][1] - f[0][0] / lambda;
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            phi[0][0] = (phi[1][0] + phi[0][1]) / 2 - f[0][0] / (2 * lambda);
                        
                        }

                    }
                    else if (i == Nx - 1 && j == 0) {
                        /* 2nd corner calculations */
                        if ( ((0 >= Bj) && (0 <= (Bj + H_cells - 1))) && ((Nx - 1 >= Bi) && (Nx - 1 <= (Bi + H_cells - 1)))  ) {
                            /* Inside void - skip this position */
                        }
                        else if (Nx - 1 == Bi + H_cells && 0 >= Bj && 0 <= Bj + H_cells - 1) {
                            /* At right void boundary */
                            phi[0][Nx - 1] = phi[1][Nx - 1] - f[0][Nx - 1] / lambda;
                        }
                        else if (0 == Bj - 1 && Nx - 1 >= Bi && Nx - 1 <= Bi + H_cells - 1) {
                            /* At bottom void boundary */
                            phi[0][Nx - 1] = phi[0][Nx - 2] - f[0][Nx - 1] / lambda;
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) / 2 - f[0][Nx - 1] / (2 * lambda);
                            
                        }

                    }
                    else if (i == 0 && j == Ny - 1) {
                        /* 3rd corner calculations */
                        if ( ((Ny - 1 >= Bj) && (Ny - 1 <= (Bj + H_cells - 1))) && ((0 >= Bi) && (0 <= (Bi + H_cells - 1)))  ) {
                            /* Inside void - skip this position */
                        }
                        else if (0 == Bi - 1 && Ny - 1 >= Bj && Ny - 1 <= Bj + H_cells - 1) {
                            /* At left void boundary */
                            phi[Ny - 1][0] = phi[Ny - 2][0] - f[Ny - 1][0] / lambda;
                        }
                        else if (Ny - 1 == Bj + H_cells && 0 >= Bi && 0 <= Bi + H_cells - 1) {
                            /* At top void boundary */
                            phi[Ny - 1][0] = phi[Ny - 1][1] - f[Ny - 1][0] / lambda;
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) / 2 - f[Ny - 1][0] / (2 * lambda);
                            
                        }

                    }
                    else if (i == Nx - 1 && j ==  Ny - 1) {
                        /* 4th corner calculations */
                            if ( ((Ny - 1 >= Bj) && (Ny - 1 <= (Bj + H_cells - 1))) && ((Nx - 1 >= Bi) && (Nx - 1 <= (Bi + H_cells - 1))) ) {
                                /* Inside void - skip this position */
                            }
                            else if (Nx - 1 == Bi + H_cells && Ny - 1 >= Bj && Ny - 1 <= Bj + H_cells - 1) {
                                /* At right void boundary */
                                phi[Ny - 1][Nx - 1] = phi[Ny - 2][Nx - 1] - f[Ny - 1][Nx - 1] / lambda;
                            }
                            else if (Ny - 1 == Bj + H_cells && Nx - 1 >= Bi && Nx - 1 <= Bi + H_cells - 1) {
                                /* At top void boundary */
                                phi[Ny - 1][Nx - 1] = phi[Ny - 1][Nx - 2] - f[Ny - 1][Nx - 1] / lambda;
                            }
                            else {
                                /* Not in void or at void boundary - calculate normally */
                                phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) / 2 - f[Ny - 1][Nx - 1] / (2 * lambda);
                                
                            }

                    }

                }
            }  
        }
        else {
            /* This can only happen in the event of an error in the values of Nx and Ny.  Print information to user */
            printf("ERROR IN VALUES OF Nx AND Ny AT PHI CALCULATIONS: Either at least one is not an integer,\n at least one is 0, or at least one is negative.\n");

        }



        /* ================================================= CALCULATING LAPLACE PHI ================================================= */
        if (Nx == 1 || Ny == 1) { /* ===================== SKIP ALL CONDITION ===================== */
           
            /* Skip everything */

        }
        else if (Nx == 2 || Ny == 2) { /* ===================== THIN DOMAIN CONDITION ===================== */
           
            if (Nx == 2 && Ny == 2) { /* ======= If domain is a 2x2 square ======= */
                /* Only corner calculations are needed */

                /* 1st corner calculation */
                if ( ((0 >= Bj) && (0 <= (Bj + H_cells - 1))) && ((0 >= Bi) && (0 <= (Bi + H_cells - 1)))  ) {
                    /* Inside void - skip this position */
                }
                else if (0 == Bi - 1 && 0 >= Bj && 0 <= Bj + H_cells - 1) {
                    /* At left void boundary (to the left of the void) */
                    laplace_phi[j][i] = (phi[j + 1][i]) * lambda - lambda * phi[j][i];
                }
                else if (0 == Bj - 1 && 0 >= Bi && 0 <= Bi + H_cells - 1) {
                    /* At bottom void boundary (beneath the void) */
                    laplace_phi[j][i] = (phi[j][i + 1]) * lambda - lambda * phi[j][i];
                }
                else {
                    /* Not in void or at void boundary - calculate normally */
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];

                }

                /* 2nd corner calculations */
                if ( ((0 >= Bj) && (0 <= (Bj + H_cells - 1))) && ((Nx - 1 >= Bi) && (Nx - 1 <= (Bi + H_cells - 1)))  ) {
                    /* Inside void - skip this position */
                }
                else if (Nx - 1 == Bi + H_cells && 0 >= Bj && 0 <= Bj + H_cells - 1) {
                    /* At right void boundary */
                    laplace_phi[j][i] = (phi[j + 1][i]) * lambda - lambda * phi[j][i];
                }
                else if (0 == Bj - 1 && Nx - 1 >= Bi && Nx - 1 <= Bi + H_cells - 1) {
                    /* At bottom void boundary */
                    laplace_phi[j][i] = (phi[j][i - 1]) * lambda - lambda * phi[j][i];
                }
                else {
                    /* Not in void or at void boundary - calculate normally */
                    laplace_phi[j][i] = (phi[j][i - 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];

                }

                /* 3rd corner calculations */
                if ( ((Ny - 1 >= Bj) && (Ny - 1 <= (Bj + H_cells - 1))) && ((0 >= Bi) && (0 <= (Bi + H_cells - 1)))  ) {
                    /* Inside void - skip this position */
                }
                else if (0 == Bi - 1 && Ny - 1 >= Bj && Ny - 1 <= Bj + H_cells - 1) {
                    /* At left void boundary */
                }
                else if (Ny - 1 == Bj + H_cells && 0 >= Bi && 0 <= Bi + H_cells - 1) {
                    /* At top void boundary */
                }
                else {
                    /* Not in void or at void boundary - calculate normally */

                }

                /* 4th corner calculations */
                if ( ((Ny - 1 >= Bj) && (Ny - 1 <= (Bj + H_cells - 1))) && ((Nx - 1 >= Bi) && (Nx - 1 <= (Bi + H_cells - 1)))  ) {
                    /* Inside void - skip this position */
                }
                else if (Nx - 1 == Bi + H_cells && Ny - 1 >= Bj && Ny - 1 <= Bj + H_cells - 1) {
                    /* At right void boundary */
                }
                else if (Ny - 1 == Bj + H_cells && Nx - 1 >= Bi && Nx - 1 <= Bi + H_cells - 1) {
                    /* At top void boundary */
                }
                else {
                    /* Not in void or at void boundary - calculate normally */

                }

            }
            else { /* ======= If domain is not square ======= */
                /* Loop through entire space */
                for (j = 0; j < Ny; j++) {
                    for (i = 0; i < Nx; i++) {

                        if (i == 0 && j == 0) {
                            /* 1st corner calculation */
                            if ( ((0 >= Bj) && (0 <= (Bj + H_cells - 1))) && ((0 >= Bi) && (0 <= (Bi + H_cells - 1)))  ) {
                                /* Inside void - skip this position */
                            }
                            else if (0 == Bi - 1 && 0 >= Bj && 0 <= Bj + H_cells - 1) {
                                /* At left void boundary (to the left of the void) */
                                laplace_phi[j][i] = (phi[j + 1][i]) * lambda - lambda * phi[j][i];
                            }
                            else if (0 == Bj - 1 && 0 >= Bi && 0 <= Bi + H_cells - 1) {
                                /* At bottom void boundary (beneath the void) */
                                laplace_phi[j][i] = (phi[j][i + 1]) * lambda - lambda * phi[j][i];
                            }
                            else {
                                /* Not in void or at void boundary - calculate normally */
                                laplace_phi[j][i] = (phi[j][i + 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];

                            }

                        }
                        else if (i == Nx - 1 && j == 0) {
                            /* 2nd corner calculations */
                            if ( ((0 >= Bj) && (0 <= (Bj + H_cells - 1))) && ((Nx - 1 >= Bi) && (Nx - 1 <= (Bi + H_cells - 1)))  ) {
                                /* Inside void - skip this position */
                            }
                            else if (Nx - 1 == Bi + H_cells && 0 >= Bj && 0 <= Bj + H_cells - 1) {
                                /* At right void boundary */
                                laplace_phi[j][i] = (phi[j + 1][i]) * lambda - lambda * phi[j][i];
                            }
                            else if (0 == Bj - 1 && Nx - 1 >= Bi && Nx - 1 <= Bi + H_cells - 1) {
                                /* At bottom void boundary */
                                laplace_phi[j][i] = (phi[j][i - 1]) * lambda - lambda * phi[j][i];
                            }
                            else {
                                /* Not in void or at void boundary - calculate normally */
                                laplace_phi[j][i] = (phi[j][i - 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];

                            }

                        }
                        else if (i == 0 && j == Ny - 1) {
                            /* 3rd corner calculations */
                        if ( ((Ny - 1 >= Bj) && (Ny - 1 <= (Bj + H_cells - 1))) && ((0 >= Bi) && (0 <= (Bi + H_cells - 1)))  ) {
                            /* Inside void - skip this position */
                        }
                        else if (0 == Bi - 1 && Ny - 1 >= Bj && Ny - 1 <= Bj + H_cells - 1) {
                            /* At left void boundary */
                            laplace_phi[j][i] = (phi[j - 1][i]) * lambda - lambda * phi[j][i];
                        }
                        else if (Ny - 1 == Bj + H_cells && 0 >= Bi && 0 <= Bi + H_cells - 1) {
                            /* At top void boundary */
                            laplace_phi[j][i] = (phi[j][i + 1]) * lambda - lambda * phi[j][i];
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i]) * lambda - 2 * lambda * phi[j][i];

                        }

                        }
                        else if (i == Nx - 1 && j ==  Ny - 1) {
                            /* 4th corner calculations */
                            if ( ((Ny - 1 >= Bj) && (Ny - 1 <= (Bj + H_cells - 1))) && ((Nx - 1 >= Bi) && (Nx - 1 <= (Bi + H_cells - 1))) ) {
                                /* Inside void - skip this position */
                            }
                            else if (Nx - 1 == Bi + H_cells && Ny - 1 >= Bj && Ny - 1 <= Bj + H_cells - 1) {
                                /* At right void boundary */
                                laplace_phi[j][i] = (phi[j - 1][i]) * lambda - lambda * phi[j][i];
                            }
                            else if (Ny - 1 == Bj + H_cells && Nx - 1 >= Bi && Nx - 1 <= Bi + H_cells - 1) {
                                /* At top void boundary */
                                laplace_phi[j][i] = (phi[j][i - 1]) * lambda - lambda * phi[j][i];
                            }
                            else {
                                /* Not in void or at void boundary - calculate normally */
                                laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i]) * lambda - 2 * lambda * phi[j][i];

                            }

                        }
                        else if (Nx == 2) { /* ======= For a tall domain ======= */
                            /* Edge calculations for Nx = 2 - there is no iterior */
                            if (i == 0) { /* Left */
                                if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                                    /* Inside void - skip this position */
                                }
                                else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells - 1) {
                                    /* At bottom void boundary */
                                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];
                                }
                                else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells - 1) {
                                    /* At left void boundary */
                                    laplace_phi[j][i] = (phi[j - 1][i] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];
                                }
                                else if (j == Bj + H_cells && i >= Bi && i <= Bi + H_cells - 1) {
                                    /* At top void boundary */
                                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];
                                }
                                else {
                                    /* Not in void or at void boundary - calculate normally */
                                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];
                                   
                                }
                            }
                            else if (i == Nx - 1) { /* Right */
                                if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                                    /* Inside void - skip this position */
                                }
                                else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells - 1) {
                                    /* At bottom void boundary */
                                    laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i]) * lambda - 2 * lambda * phi[j][i];
                                }
                                else if (i == Bi + H_cells && j >= Bj && j <= Bj + H_cells - 1) {
                                    /* At right void boundary */
                                    laplace_phi[j][i] = (phi[j - 1][i] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];
                                }
                                else if (j == Bj + H_cells && i >= Bi && i <= Bi + H_cells - 1) {
                                    /* At top void boundary */
                                    laplace_phi[j][i] = (phi[j][i - 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];
                                }
                                else {
                                    /* Not in void or at void boundary - calculate normally */
                                    laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];

                                }
                               
                            }

                        }
                        else if (Ny == 2) { /* ======= For a flat domain ======= */
                            /* Edge calculations for Ny = 2 - there is no iterior */
                            if (j == 0) { /* Bottom */
                                if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                                    /* Inside void - skip this position */
                                }
                                else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells - 1) {
                                    /* At bottom void boundary */
                                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1]) * lambda - 2 * lambda * phi[j][i];
                                }
                                else if (i == Bi + H_cells && j >= Bj && j <= Bj + H_cells - 1) {
                                    /* At right void boundary */
                                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];
                                }
                                else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells - 1) {
                                    /* At left void boundary */
                                    laplace_phi[j][i] = (phi[j][i - 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];
                                }
                                else {
                                    /* Not in void or at void boundary - calculate normally */
                                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];

                                }
                           
                            }
                            else if (j == Ny - 1) { /* Top */
                                if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                                    /* Inside void - skip this position */
                                }
                                else if (j == Bj + H_cells && i >= Bi && i <= Bi + H_cells - 1) {
                                    /* At top void boundary */
                                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1]) * lambda - 2 * lambda * phi[j][i];
                                }
                                else if (i == Bi + H_cells && j >= Bj && j <= Bj + H_cells - 1) {
                                    /* At right void boundary */
                                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i]) * lambda - 2 * lambda * phi[j][i];
                                }
                                else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells - 1) {
                                    /* At left void boundary */
                                    laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i]) * lambda - 2 * lambda * phi[j][i];
                                }
                                else {
                                    /* Not in void or at void boundary - calculate normally */
                                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - 3 * lambda * phi[j][i];

                                }
                           
                            }

                        }

                    }
                }
            }
        }
        else if (Nx >= 3 && Ny >= 3) { /* ===================== NORMAL DOMAIN CONDITION ===================== */

            /* Loop through entire space */
            for (j = 0; j < Ny; j++) {
                for (i = 0; i < Nx; i++) {

                    if ((i > 0 && i < Nx - 1) && (j > 0 && j < Ny - 1)) {
                        /* Interior calculations */
                        if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                            /* Inside void - skip this position */
                        }
                        else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At left void boundary */
                            laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];
                        }
                        else if (i == Bi + H_cells && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At right void boundary */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];
                        }
                        else if (j == Bj + H_cells && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At top void boundary */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];
                        }
                        else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At bottom void boundary */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - 3 * lambda * phi[j][i];
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 4 * lambda * phi[j][i];
                            
                        }

                    }
                    else if (j == 0 && i != 0 && i != Nx - 1) {
                        /* Edge calcualtions on bottom */
                        if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                            /* Inside void - skip this position */
                        }
                        else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At bottom void boundary */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1]) * lambda - 2 * lambda * phi[j][i];
                        }
                        else if (i == Bi + H_cells && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At right void boundary */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];
                        }
                        else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At left void boundary */
                            laplace_phi[j][i] = (phi[j][i - 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];
                        }

                    }
                    else if (j == Ny - 1 && i != 0 && i != Nx - 1) {
                        /* Edge calcualtions on top */
                        if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                            /* Inside void - skip this position */
                        }
                        else if (j == Bj + H_cells && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At top void boundary */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1]) * lambda - 2 * lambda * phi[j][i];
                        }
                        else if (i == Bi + H_cells && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At right void boundary */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i]) * lambda - 2 * lambda * phi[j][i];
                        }
                        else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At left void boundary */
                            laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i]) * lambda - 2 * lambda * phi[j][i];
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - 3 * lambda * phi[j][i];

                        }

                    }
                    else if (i == 0 && j != 0 && j != Ny - 1) {
                        /* Edge calcualtions on left */
                        if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                            /* Inside void - skip this position */
                        }
                        else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At bottom void boundary */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];
                        }
                        else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At left void boundary */
                            laplace_phi[j][i] = (phi[j - 1][i] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];
                        }
                        else if (j == Bj + H_cells && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At top void boundary */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];
                           
                        }

                    }
                    else if (i == Nx - 1 && j != 0 && j != Ny - 1) {
                        /* Edge calcualtions on right */
                        if ( ((j >= Bj) && (j <= (Bj + H_cells - 1))) && ((i >= Bi) && (i <= (Bi + H_cells - 1))) ) {
                            /* Inside void - skip this position */
                        }
                        else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At bottom void boundary */
                            laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i]) * lambda - 2 * lambda * phi[j][i];
                        }
                        else if (i == Bi + H_cells && j >= Bj && j <= Bj + H_cells - 1) {
                            /* At right void boundary */
                            laplace_phi[j][i] = (phi[j - 1][i] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];
                        }
                        else if (j == Bj + H_cells && i >= Bi && i <= Bi + H_cells - 1) {
                            /* At top void boundary */
                            laplace_phi[j][i] = (phi[j][i - 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];

                        }

                    }
                    else if (i == 0 && j == 0) {
                        /* 1st corner calculation */
                        if ( ((0 >= Bj) && (0 <= (Bj + H_cells - 1))) && ((0 >= Bi) && (0 <= (Bi + H_cells - 1)))  ) {
                            /* Inside void - skip this position */
                        }
                        else if (0 == Bi - 1 && 0 >= Bj && 0 <= Bj + H_cells - 1) {
                            /* At left void boundary (to the left of the void) */
                            laplace_phi[j][i] = (phi[j + 1][i]) * lambda - lambda * phi[j][i];
                        }
                        else if (0 == Bj - 1 && 0 >= Bi && 0 <= Bi + H_cells - 1) {
                            /* At bottom void boundary (beneath the void) */
                            laplace_phi[j][i] = (phi[j][i + 1]) * lambda - lambda * phi[j][i];
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];

                        }

                    }
                    else if (i == Nx - 1 && j == 0) {
                        /* 2nd corner calculations */
                        if ( ((0 >= Bj) && (0 <= (Bj + H_cells - 1))) && ((Nx - 1 >= Bi) && (Nx - 1 <= (Bi + H_cells - 1)))  ) {
                            /* Inside void - skip this position */
                        }
                        else if (Nx - 1 == Bi + H_cells && 0 >= Bj && 0 <= Bj + H_cells - 1) {
                            /* At right void boundary */
                            laplace_phi[j][i] = (phi[j + 1][i]) * lambda - lambda * phi[j][i];
                        }
                        else if (0 == Bj - 1 && Nx - 1 >= Bi && Nx - 1 <= Bi + H_cells - 1) {
                            /* At bottom void boundary */
                            laplace_phi[j][i] = (phi[j][i - 1]) * lambda - lambda * phi[j][i];
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            laplace_phi[j][i] = (phi[j][i - 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];

                        }

                    }
                    else if (i == 0 && j == Ny - 1) {
                        /* 3rd corner calculations */
                        if ( ((Ny - 1 >= Bj) && (Ny - 1 <= (Bj + H_cells - 1))) && ((0 >= Bi) && (0 <= (Bi + H_cells - 1)))  ) {
                            /* Inside void - skip this position */
                        }
                        else if (0 == Bi - 1 && Ny - 1 >= Bj && Ny - 1 <= Bj + H_cells - 1) {
                            /* At left void boundary */
                            laplace_phi[j][i] = (phi[j - 1][i]) * lambda - lambda * phi[j][i];
                        }
                        else if (Ny - 1 == Bj + H_cells && 0 >= Bi && 0 <= Bi + H_cells - 1) {
                            /* At top void boundary */
                            laplace_phi[j][i] = (phi[j][i + 1]) * lambda - lambda * phi[j][i];
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i]) * lambda - 2 * lambda * phi[j][i];

                        }

                    }
                    else if (i == Nx - 1 && j ==  Ny - 1) {
                        /* 4th corner calculations */
                        if ( ((Ny - 1 >= Bj) && (Ny - 1 <= (Bj + H_cells - 1))) && ((Nx - 1 >= Bi) && (Nx - 1 <= (Bi + H_cells - 1))) ) {
                            /* Inside void - skip this position */
                        }
                        else if (Nx - 1 == Bi + H_cells && Ny - 1 >= Bj && Ny - 1 <= Bj + H_cells - 1) {
                            /* At right void boundary */
                            laplace_phi[j][i] = (phi[j - 1][i]) * lambda - lambda * phi[j][i];
                        }
                        else if (Ny - 1 == Bj + H_cells && Nx - 1 >= Bi && Nx - 1 <= Bi + H_cells - 1) {
                            /* At top void boundary */
                            laplace_phi[j][i] = (phi[j][i - 1]) * lambda - lambda * phi[j][i];
                        }
                        else {
                            /* Not in void or at void boundary - calculate normally */
                            laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i]) * lambda - 2 * lambda * phi[j][i];

                        }

                    }

                }
            }  
        }
        else {
            /* This can only happen in the event of an error in the values of Nx and Ny.  Print information to user */
            printf("ERROR IN VALUES OF Nx AND Ny AT LAPLACE PHI CALCULATIONS:\nEither at least one is not an integer,\n at least one is 0, or at least one is negative.\n");

        }
   



        /* compute the norm */
        laplace_phi_minus_f_norm = 0;

        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {

               
                if ( ((j < Bj) || ( j > (Bj + H_cells))) || ((i < Bi) || (i > (Bi + H_cells)))  ) {     /* if not inside the void... */
                    laplace_phi_minus_f_norm += pow((laplace_phi[j][i] - f[j][i]), 2);
                }

            }
        }

       

        laplace_phi_minus_f_norm = sqrt(laplace_phi_minus_f_norm);



        if (f_norm == 0) {
            RHS = epsilon;
        }
        else {
            RHS = epsilon * f_norm;
        }


        /* break out of the loop if the max number of steps has been reached */
        if (step == max_num_steps) {
            break;
        }


        step += 1;

    } while (laplace_phi_minus_f_norm > RHS && Nx != 1 && Ny != 1);



    if (Nx != 1 && Ny != 1) {

        /* Impose condition that the integral over the domain is equal to zero */
        integral = 0;
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
               
                if ( ((j < Bj) || ( j > (Bj + H_cells))) || ((i < Bi) || (i > (Bi + H_cells)))  ) {     /* if not inside the void... */
                    integral += phi[j][i] * D_x * D_y;
                }

            }
        }

        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {

                if ( ((j < Bj) || ( j > (Bj + H_cells))) || ((i < Bi) || (i > (Bi + H_cells)))  ) {     /* if not inside the void... */
                    phi[j][i] -= integral / (Nx * Ny - pow(H_cells, 2));
                }
            
            }
        }

    }


    FILE* fid_laplace_phi = fopen("laplace_p.txt", "w");
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            fprintf(fid_laplace_phi, "%.20lf ", laplace_phi[j][i]);
        }
        fprintf(fid_laplace_phi, "\n");
    }
    fclose(fid_laplace_phi);



    /* Freeing arrays made in this function */
    for (i = 0; i < Ny; i++) {
        free(laplace_phi[i]);
    }

    free(laplace_phi);

    return 0;

}