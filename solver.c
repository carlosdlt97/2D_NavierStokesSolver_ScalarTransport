#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>


double* TriDiag_GaussElim(int size, double dx_or_dy, double dt, double Re, double* d, int right_outlet, int solving_for_du_s) {

    ////input size should be Nx
    ////should pass in d of whole row except for first point (should be size:(size - 1) )
    ////right_outlet = 1 if the right ourlet boundary is to the right of what we are looking at

    int i;

    double* b = (double*)calloc(size-1, sizeof(double));
    double* a = (double*)calloc(size-1, sizeof(double));
    double* c = (double*)calloc(size-1, sizeof(double));


    //initialize tridiagonal matrix
    for (i = 0; i < size - 1; i++) { 
        b[i] = 1 + dt / (Re * pow(dx_or_dy, 2) );
        a[i] = -dt / (2 * Re * pow(dx_or_dy,2) );
        c[i] = -dt / (2 * Re * pow(dx_or_dy,2) );
    }


    //if we have the outlet boundary in our domain
    if (right_outlet == 1)  {
        b[size - 2] = 1 + dt / (2 * Re * pow(dx_or_dy,2));
    }

    //if we are solving for du_s 
    if (solving_for_du_s == 1) {
        b[0] = 1 + 3 * dt / (2 * Re * pow(dx_or_dy, 2));
        b[size - 2] = 1 + 3 * dt / (2 * Re * pow(dx_or_dy, 2));
    }



    //forward elimination
    for (i = 1; i < size - 1; i++) {
        b[i] = b[i] - c[i-1] * a[i] / b[i - 1];
        
        d[i] = d[i] - d[i-1] * a[i] / b[i-1]; 
    }


    //back substitution
    d[size - 2] = d[size - 2] / b[size - 2];

    for (i = size - 3; i > -1; i--) {
        d[i] = ( d[i] - c[i] * d[i+1] ) / b[i];
    }



    free(a);
    free(b);
    free(c);


    return d;

}




/* ==========================================================================
======================= n-Step Gauss-Seidel Solver ==========================
========================================================================== */

int GS_nstep(double** f, double** phi, int Nx, int Ny, int D_x, int D_y, double epsilon, int nGS, int H_cells, int Bj, int Bi) {

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

                if ( ((j < Bj) || ( j > (Bj + H_cells - 1))) && ((i < Bi) || (i > (Bi + H_cells - 1)))  ) {     /* if not inside the void... */
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

                    if ((i > 0 && i < Nx - 1) && (j > 0 && j < Ny + 1)) {
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

                    if ((i > 0 && i < Nx - 1) && (j > 0 && j < Ny + 1)) {
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

                
                if ( ((j < Bj) || ( j > (Bj + H_cells))) && ((i < Bi) || (i > (Bi + H_cells)))  ) {     /* if not inside the void... */
                    laplace_phi_minus_f_norm = laplace_phi_minus_f_norm + pow((laplace_phi[j][i] - f[j][i]), 2);
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
                
                if ( ((j < Bj) || ( j > (Bj + H_cells))) && ((i < Bi) || (i > (Bi + H_cells)))  ) {     /* if not inside the void... */
                    integral += phi[j][i] * D_x * D_y;
                }

            }
        }

        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {

                if ( ((j < Bj) || ( j > (Bj + H_cells))) && ((i < Bi) || (i > (Bi + H_cells)))  ) {     /* if not inside the void... */
                    phi[j][i] -= integral / (Nx * Ny - pow(H_cells, 2));
                }
            
            }
        }

    }


    /* Freeing arrays made in this function */
    for (i = 0; i < nx; i++) {
        free(laplace_phi[i]);
    }

    free(laplace_phi);

    return 0;

}







/* ================================================================================
= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
======================================= MAIN ======================================
= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
================================================================================ */


int main() {
    /* -------------------------------------------------------------------------
    -------------------------- Initializing Variables --------------------------
    ------------------------------------------------------------------------- */
    int i = 0, j = 0, k = 0; /* Counters for all "for loops" */
    int iter = 0;

    double Re = 20;  /* Problem parameters */
    double dt = 0.0001;
    double H = 1;
    double U_inlet = 1;
    double Diff = pow(10,-4);
    double nGS = 15000; /* Number of Gauss-Seidel steps to take per relaxation in mutligrid acceleration */
    
    char convective_method[] = "quick";  // options are "upwind", "centraldiff", or "quick"

    double epsilon = pow(10, -3);
    int max_time_steps = 1;

    int Nx = 500;
    int Ny = 100;
    double nx = Nx;
    double ny = Ny;
    int H_cells = 20;
    double dx = H / H_cells;
    double dy = H / H_cells;

    int Bj = 40;
    int Bi = 200;



    /* Initializing variables for this step 1 */
    double Hx;
    double Hy;
    double u_cc; /* cc = Cell cented velocity */
    double v_cc;
    double u_cc_im1; /* im1 = "i minus 1"; cell centered velocity at 1 position back in the i direction */
    double v_cc_jm1; /* jm1 = "j minus 1"; cell centered velocity at 1 position back in the j direction */
    double u_s;  /* s = staggered velocity (i-1/2, j-1/2) from the cc values */
    double v_s;
    double u_s_ip1; /* ip1 = "i plus 1"; staggered velocity at 1 position forward in the i direction */
    double u_s_jp1; /* jp1 = "j plus 1"; staggered velocity at 1 position forward in the j direction */
    double v_s_ip1;
    double v_s_jp1;

    /* For scalar transport */
    double convec;
    double diffu;


    /* -------------------------------------------------------------------------
    -------------------------- Initializing Arrays -----------------------------
    ------------------------------------------------------------------------- */

    double** u = (double**)calloc(ny, sizeof(double*));
    double** v = (double**)calloc(ny, sizeof(double*));

    double** Hx_old = (double**)calloc(ny, sizeof(double*));
    double** Hy_old = (double**)calloc(ny, sizeof(double*));

    double** p = (double**)calloc(ny, sizeof(double*));             //cell centered
    double** grad_u_star_over_dt = (double**)calloc(ny, sizeof(double*));
    
    double** phi = (double**)calloc(ny, sizeof(double*));
    double** phi_new = (double**)calloc(ny, sizeof(double*));

    double** step1_mat_x = (double**)calloc(ny, sizeof(double*));
    double** step1_mat_y = (double**)calloc(ny, sizeof(double*));


    double** du_s = (double**)calloc(ny, sizeof(double*));
    double** dv_s = (double**)calloc(ny, sizeof(double*));
    double** du_ss = (double**)calloc(ny, sizeof(double*));
    double** dv_ss = (double**)calloc(ny, sizeof(double*));


    // check the size of these to make sure they align with the code

    double* temp_vec_x_small = (double*)calloc(99, sizeof(double));
    double* temp_vec_x_medium = (double*)calloc(19 * 20 - 1, sizeof(double));
    double* temp_vec_x_medium2 = (double*)calloc(19 * 20, sizeof(double));
    double* temp_vec_x_long = (double*)calloc(nx - 1, sizeof(double));

    double* temp_vec_y_small = (double*)calloc(40, sizeof(double));
    double* temp_vec_y_small2 = (double*)calloc(39, sizeof(double));
    double* temp_vec_y_long = (double*)calloc(ny, sizeof(double));
    double* temp_vec_y_long2 = (double*)calloc(ny - 1, sizeof(double));



    for (i = 0; i < Ny; i++) {
        u[i] = (double*)calloc(nx, sizeof(double));
        v[i] = (double*)calloc(nx, sizeof(double));

        Hx_old[i] = (double*)calloc(nx, sizeof(double));
        Hy_old[i] = (double*)calloc(nx, sizeof(double));

        p[i] = (double*)calloc(nx, sizeof(double));
        grad_u_star_over_dt[i] = (double*)calloc(nx, sizeof(double));
        
        phi[i] = (double*)calloc(nx, sizeof(double));
        phi_new[i] = (double*)calloc(nx, sizeof(double));

        du_s[i] = (double*)calloc(nx, sizeof(double));
        dv_s[i] = (double*)calloc(nx, sizeof(double));
        du_ss[i] = (double*)calloc(nx, sizeof(double));
        dv_ss[i] = (double*)calloc(nx, sizeof(double));

        step1_mat_x[i] = (double*)calloc(nx, sizeof(double));
        step1_mat_y[i] = (double*)calloc(nx, sizeof(double));
        
    }

    //these variables must be a little larger for the multigrid to work
    double** u_star = (double**)calloc(ny, sizeof(double*));
    double** v_star = (double**)calloc(ny + 1, sizeof(double*));

    for (i = 0; i < Ny; i++) {
        u_star[i] = (double*)calloc(nx + 1, sizeof(double));
    }

    for (i = 0; i < Ny + 1; i++) {
        v_star[i] = (double*)calloc(nx, sizeof(double));
    }

    //Initialize velocity at the left boundary
    for (j = 0; j < Ny; j++) {
        u[j][0] = U_inlet;
        u_star[j][0] = U_inlet;    
    }
 

    
    /* -------------------------------------------------------------------------
    ----------------------------------- Step 1 ---------------------------------
    ------------------------------------------------------------------------- */


    


    /* -------------------------------------------------------------------------
    ---------------------------------- Step 2 ----------------------------------
    ------------------------------------------------------------------------- */

    /* Calculating the cell-centered divergence of velocity_star divided by dt*/
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            
            grad_u_star_over_dt[j][i] = ((u_star[j][i + 1] - u_star[j][i]) / dx + (v_star[j + 1][i] - v_star[j][i]) / dy) / dt;


        }
    }

    /* Solving for cell-centered pressure using multigrid acceleration method */
    GS_nstep(grad_u_star_over_dt, p, Nx, Ny, dx, dy, epsilon, nGS, H_cells, Bj, Bi);




    /* -------------------------------------------------------------------------
    ---------------------------------- Step 3 ----------------------------------
    ------------------------------------------------------------------------- */

    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {

            if ((j >= Bj && j <= Bj + H_cells - 1) && (i >= Bi && i <= Bi + H_cells - 1)) {
                /* Inside the void - set velocity equal to 0 */
                u[j][i] = 0;
                v[j][i] = 0;
            }
            else if (i == 0 && j !=0 && j != Ny - 1) {
                /* On the domain's left edge */
                u[j][i] = U_inlet;
                v[j][i] = 0;
            }
            else if (i == Nx - 1 && j !=0 && j != Ny - 1) {
                /* On the domain's right edge */
                u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;
            }
            else if (j == 0) {
                /* On the domain's bottom edge */
                u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                v[j][i] = 0;
            }
            else if (j == Ny - 1) {
                /* On the domain's top edge */
                u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;
            }
            else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells) {
                /* On the void's left edge */
                u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;
            }
            else if (i == Bi + H_cells && j >= Bj && j <= Bj + H_cells - 1) {
                /* On the void's right edge */
                u[j][i] = 0;
                v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;
            }
            else if (j == Bj + H_cells && i >= Bi && i <= Bi + H_cells - 1) {
                /* On the void's top edge */
                u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                v[j][i] = 0;
            }
            else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells) {
                /* On the void's bottom edge */
                u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;
            }
            else {
                /* In normal space - calculate normally */
                u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;

            }

        }
    }








    /* -------------------------------------------------------------------------
    ----------------------------- Scalar Transport -----------------------------
    ------------------------------------------------------------------------- */


    
    

    /* -------------------------------------------------------------------------
    ---------------------------- Freeing Variables -----------------------------
    ------------------------------------------------------------------------- */
    
    for (i = 0; i < Ny; i++) {
        free(u[i]);
        free(v[i]);

        free(Hx_old[i]);
        free(Hy_old[i]);

        free(p[i]);
        free(grad_u_star_over_dt[i]);
        
        free(phi[i]);
        free(phi_new[i]);

        free(step1_mat_x[i]);
        free(step1_mat_y[i]);

        free(du_s[i]);
        free(dv_s[i]);
        free(du_ss[i]);
        free(dv_ss[i]);

    }

    free(u);
    free(v);

    free(Hx_old);
    free(Hy_old);

    free(p);
    free(grad_u_star_over_dt);
    
    free(phi);
    free(phi_new);

    free(step1_mat_x);
    free(step1_mat_y);


    free(du_s);
    free(dv_s);
    free(du_ss);
    free(dv_ss);

    free(temp_vec_x_small);
    free(temp_vec_x_medium);
    free(temp_vec_x_medium2);
    free(temp_vec_x_long);

    free(temp_vec_y_small);
    free(temp_vec_y_small2);
    free(temp_vec_y_long);
    free(temp_vec_y_long2);

    for (i = 0; i < Ny; i++) {
        free(u_star[i]);
    }

    for (i = 0; i < Ny + 1; i++) {
        free(v_star[i]);
    }

    free(u_star);
    free(v_star);

    printf("\n End of script");

    return 0;
    

}