#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>



/* ==========================================================================
======================= n-Step Gauss-Seidel Solver ==========================
========================================================================== */

int GS_nstep(double** f, double** phi, int Nx, int Ny, double epsilon, int nGS) {

    /* Initilaizations */
    int i, j;
    double nx = Nx;
    double ny = Ny;
    double D_x = 1 / nx;
    double D_y = 1 / ny;

    double  f_norm;
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

    f_norm = 0; /* compute f_norm */
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            f_norm = f_norm + pow(f[j][i], 2);
        }
    }

    f_norm = sqrt(f_norm);


    if (Nx > 1 && Ny > 1) {

        do {

            if (Nx > 2 && Ny > 2) {

                for (j = 1; j < Ny - 1; j++) {
                    for (i = 1; i < Nx - 1; i++) {
                        phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 4 - f[j][i] / (4 * lambda);

                    }
                }
                

                //update left boundary values
                for (j = 1; j < Ny - 1; j++) {
                    i = 0;
                    phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

                }

                //update right boundary values
                for (j = 1; j < Ny - 1; j++) {
                    i = Nx - 1;
                    phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

                }

                //update bottom boundary values 
                for (i = 1; i < Nx - 1; i++) {
                    j = 0;
                    phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                }

                //update top boundary values
                for (i = 1; i < Nx - 1; i++) {
                    j = Ny - 1;
                    phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) / 3 - f[j][i] / (3 * lambda);
                }

            }


            //update corner points
            phi[0][0] = (phi[1][0] + phi[0][1]) / 2 - f[0][0] / (2 * lambda);
            phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) / 2 - f[0][Nx - 1] / (2 * lambda);
            phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) / 2 - f[Ny - 1][0] / (2 * lambda);
            phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) / 2 - f[Ny - 1][Nx - 1] / (2 * lambda);

            
            if (Nx > 2 && Ny > 2) {

                //compute laplace_p matrix
                for (j = 1; j < Ny - 1; j++) {
                    for (i = 1; i < Nx - 1; i++) {
                        laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 4 * lambda * phi[j][i];
                    }
                }
                for (j = 1; j < Ny - 1; j++) {
                    i = 0;
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
                }
                for (j = 1; j < Ny - 1; j++) {
                    i = Nx - 1;
                    laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
                }
                for (i = 1; i < Nx - 1; i++) {
                    j = 0;
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
                }
                for (i = 1; i < Nx - 1; i++) {
                    j = Ny - 1;
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - phi[j][i] * (3 * lambda);
                }

            }

            laplace_phi[0][0] = (phi[1][0] + phi[0][1]) * lambda - phi[0][0] * (2 * lambda);
            laplace_phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) * lambda - phi[0][Nx - 1] * (2 * lambda);
            laplace_phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) * lambda - phi[Ny - 1][0] * (2 * lambda);
            laplace_phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) * lambda - phi[Ny - 1][Nx - 1] * (2 * lambda);



            //compute the norm
            laplace_phi_minus_f_norm = 0;

            for (j = 0; j < Ny; j++) {
                for (i = 0; i < Nx; i++) {
                    laplace_phi_minus_f_norm = laplace_phi_minus_f_norm + pow((laplace_phi[j][i] - f[j][i]), 2);
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

        } while (laplace_phi_minus_f_norm > RHS);
    }


    /* Impose condition that the integral over the domain is equal to zero */
    integral = 0;
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            integral += phi[j][i] * D_x * D_y;
        }
    }

    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            phi[j][i] = phi[j][i] - integral / (Nx * Ny);
        }
    }


    for (i = 0; i < nx; i++) {
        free(laplace_phi[i]);
    }

    free(laplace_phi);

    return 0;

}



/* ====================================================================================================================
=================================== MG_recursion ======================================================================
==================================================================================================================== */


int MG_recursion(double** f, double** phi, int Nx, int Ny, double epsilon, int nGS) {

    /* Initializing variables */
    int finished;
    int i, j;

    double nx = Nx;
    double ny = Ny;
    double D_x = 1 / nx;
    double D_y = 1 / ny;
    double lambda = pow(D_x, -2);

    double half_nx = ceil(nx / 2);
    double half_ny = ceil(ny / 2);
    double half_D_x = 1 / half_nx;
    double half_D_y = 1 / half_nx;
    double half_lambda = pow(half_D_x, -2);
    int half_Nx = half_nx;
    int half_Ny = half_ny;

    double** residual = (double**)calloc(Ny, sizeof(double*));
    double** residual2 = (double**)calloc(half_Ny, sizeof(double*));
    double** error = (double**)calloc(Ny, sizeof(double*));
    double** error2 = (double**)calloc(half_Ny, sizeof(double*));
    double** laplace_phi = (double**)calloc(Ny, sizeof(double*));
    
    for (i = 0; i < Nx; i++) {
        residual[i] = (double*)calloc(Nx, sizeof(double));
        error[i] = (double*)calloc(Nx, sizeof(double));
        laplace_phi[i] = (double*)calloc(Nx, sizeof(double));
    }
    for (i = 0; i < half_Nx; i++) {
        residual2[i] = (double*)calloc(half_Nx, sizeof(double));
        error2[i] = (double*)calloc(half_Nx, sizeof(double));
    }
    
    /* ================================================================================================*/
    
    /* ------- GETTING RESIDUAL ------- */
    /* compute laplace_p matrix */
    if (Nx > 1 && Ny > 1) {
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 4 * lambda * phi[j][i];
            }
        }
        for (j = 1; j < Ny - 1; j++) {
            i = 0;
            laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        for (j = 1; j < Ny - 1; j++) {
            i = Nx - 1;
            laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        for (i = 1; i < Nx - 1; i++) {
            j = 0;
            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        for (i = 1; i < Nx - 1; i++) {
            j = Ny - 1;
            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }

        laplace_phi[0][0] = (phi[1][0] + phi[0][1]) * lambda - phi[0][0] * (2 * lambda);
        laplace_phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) * lambda - phi[0][Nx - 1] * (2 * lambda);
        laplace_phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) * lambda - phi[Ny - 1][0] * (2 * lambda);
        laplace_phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) * lambda - phi[Ny - 1][Nx - 1] * (2 * lambda);

        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                residual[j][i] = f[j][i] - laplace_phi[j][i];

            }
        }
    }

    

    /* ------- RESTRICTION ------- */
    /* Restricting residual by taking average of a point and its surounding 4 points and placing it into residual2 (adjusted for edges & corners) */
    if (Nx > 2 && Ny > 2) { /* For residual sizes greater than 2x2 */

        /* For interior points */
        for (j = 1; j < (half_Ny - 1); j++){
            for(i = 1; i < (half_Nx - 1);  i++) {

                residual2[j][i] = (residual[j * 2][i * 2] + residual[j * 2][i * 2 + 1] + residual[j * 2 + 1][i * 2] + residual[j * 2][i * 2 - 1] + residual[j * 2 - 1][i * 2]) / 5;
            
            }
        }

        /* For the left boundary */
        for (j = 1; j < (half_Ny - 1); j++) {
            i = 0;
            residual2[j][i] = (residual[j * 2][i * 2] + residual[j * 2][i * 2 + 1] + residual[j * 2 + 1][i * 2] + residual[j * 2 - 1][i * 2]) / 4;

        }

        /* For the bottom boundary */
        for (i = 1; j < (half_Nx - 1); j++) {
            j = 0;
            residual2[j][i] = (residual[j * 2][i * 2] + residual[j * 2 + 1][i * 2] + residual[j * 2][i * 2 + 1] + residual[j * 2][i * 2 - 1]) / 4;

        }

        /* ODD */
        if (nx / half_nx != 2 && ny / half_ny != 2) { /* For odd-sized square meshes (7x7, 23x23, etc.) */
            /* For the right boundary */
            for (j = 1; j < (half_Ny - 1); j++) {
                i = half_Nx - 1;
                residual2[j][i] = (residual[j * 2][i * 2] + residual[j * 2][i * 2 + 1] + residual[j * 2 + 1][i * 2] + residual[j * 2 - 1][i * 2]) / 4;

            }

            /* For the top boundary */
            for (i = 1; j < (half_Nx - 1); j++) {
                j = half_Ny - 1;
                residual2[j][i] = (residual[j * 2][i * 2] + residual[j * 2 + 1][i * 2] + residual[j * 2][i * 2 + 1] + residual[j * 2][i * 2 - 1]) / 4;

            }

            /* For the corners */
            residual2[0][0] = (residual[0][0] + residual[0][1] + residual[1][0]) / 3; /* bottom left */
            residual2[0][half_Nx - 1] = (residual[0][(half_Nx - 1) * 2] + residual[0][((half_Nx - 1) * 2) - 1] + residual[1][(half_Nx - 1) * 2]) / 3; /* bottom right */
            residual2[half_Ny - 1][0] = (residual[(half_Ny - 1) * 2][0] + residual[((half_Ny - 1) * 2) - 1][0] + residual[(half_Ny - 1) * 2][1]) / 3; /* top left */
            residual2[half_Ny - 1][half_Nx - 1] = (residual[(half_Ny - 1) * 2][(half_Nx - 1) * 2] +  + residual[(half_Ny - 1) * 2][(half_Nx - 1) * 2 - 1] + residual[(half_Ny - 1) * 2 - 1][(half_Nx - 1) * 2]) / 3; /* top right */
            
        }
        /* EVEN */
        else if (nx / half_nx == 2 && ny / half_ny == 2){ /* For even-sized square meshes (6x6, 24x24, etc.) */
                /* For the right boundary */
            for (j = 1; j < (half_Ny - 1); j++) {
                i = half_Nx - 1;
                residual2[j][i] = (residual[j * 2][i * 2] + residual[j * 2][i * 2 + 1] + residual[j * 2 + 1][i * 2] + residual[j * 2][i * 2 - 1] + residual[j * 2 - 1][i * 2]) / 5;

            }

            /* For the top boundary */
            for (i = 1; j < (half_Nx - 1); j++) {
                j = half_Ny - 1;
                residual2[j][i] = (residual[j * 2][i * 2] + residual[j * 2][i * 2 + 1] + residual[j * 2 + 1][i * 2] + residual[j * 2][i * 2 - 1] + residual[j * 2 - 1][i * 2]) / 5;
            
            }

            /* For the corners */
            residual2[0][0] = (residual[0][0] + residual[0][1] + residual[1][0]) / 3; /* bottom left */
            residual2[0][half_Nx - 1] = (residual[0][(half_Nx - 1) * 2] + residual[0][((half_Nx - 1) * 2) + 1] + residual[0][((half_Nx - 1) * 2) - 1] + residual[1][(half_Nx - 1) * 2]) / 4; /* bottom right */
            residual2[half_Ny - 1][0] = (residual[(half_Ny - 1) * 2][0] + residual[((half_Ny - 1) * 2) + 1][0] + residual[((half_Ny - 1) * 2) - 1][0] + residual[(half_Ny - 1) * 2][1]) / 4; /* top left */
            residual2[half_Ny - 1][half_Nx - 1] = (residual[(half_Ny - 1) * 2][(half_Nx - 1) * 2] + residual[(half_Ny - 1) * 2][(half_Nx - 1) * 2 + 1] + residual[(half_Ny - 1) * 2 + 1][(half_Nx - 1) * 2] + residual[(half_Ny - 1) * 2][(half_Nx - 1) * 2 - 1] + residual[(half_Ny - 1) * 2 - 1][(half_Nx - 1) * 2]) / 5; /* top right */

        }
        else { /* Making sure the mesh is square, and terminating recursion otherwise (also frees memory) */
            printf("\n\nYOU DID NOT INPUT A SQUARE MESH.\n\n");

            for (j = 0; j < Ny; j++) {
                free(residual[j]);
                free(error[j]);
                free(laplace_phi[j]);
            }
            for (j = 0; j < half_Ny; j++) {
                free(residual2[j]);
                free(error2[j]);
            }
            free(residual);
            free(residual2);
            free(error);
            free(error2);
            free(laplace_phi);                

            return 0;
        }

    }
    else if (Nx == 2 && Ny == 2){ /* For phi & f sizes of 2x2 */
        residual2[0][0] = (residual[0][0] + residual[0][1] + residual[1][0]) / 3;
    }
    /* residual size of 1x1 is not restricted */


    /* ------- RELAXATION ON error ------- */
    GS_nstep(residual2, error2, half_Nx, half_Ny, epsilon, nGS);



    /* ------- RECURSION ------- */
    if (half_Nx > 1 && half_Ny > 1) {
        MG_recursion(residual2, error2, half_Nx, half_Ny, epsilon, nGS);
    }

    
    /* ------- PROLONGATION OF ERROR ------- */
    /* EVEN */
    if (nx / half_nx == 2 && ny / half_ny == 2){ /* For even-sized square meshes (6x6, 24x24, etc.) */
        for (j = 0; j < half_Ny; j++){
            for(i = 0; i < half_Nx; i++) {
                error[j * 2][i * 2] = error2[j][i];
                error[j * 2 + 1][i * 2] = error2[j][i];
                error[j * 2][i * 2 + 1] = error2[j][i];
                error[j * 2 + 1][i * 2 + 1] = error2[j][i];
            }
        }
    }
    /* ODD */
    else { /* For odd-sized square meshes (7x7, 23x23, etc.) */
        /* interior points */
        for (j = 0; j < half_Ny - 1; j++){
            for(i = 0; i < half_Nx - 1; i++) {
                error[j * 2][i * 2] = error2[j][i];
                error[j * 2 + 1][i * 2] = error2[j][i];
                error[j * 2][i * 2 + 1] = error2[j][i];
                error[j * 2 + 1][i * 2 + 1] = error2[j][i];
            }
        }
        
        /* edges */
        /* right edge */
        for (j = 0; j < half_Ny - 1; j++) {
            i = half_Nx - 1;
            error[j * 2][i * 2] = error2[j][i];
            error[j * 2 + 1][i * 2] = error2[j][i];
        }

        /* top edge */
        for (i = 0; i < half_Nx - 1; i++) {
            j = half_Ny - 1;
            error[j * 2][i * 2] = error2[j][i];
            error[j * 2][i * 2 + 1] = error2[j][i];
        }

        /* top-right corner */
        error[Ny - 1][Nx - 1] = error2[half_Ny - 1][half_Nx - 1];
    }

    

    /* ------- CORRECTION ------- */
    for (j = 0; j < Ny; j++){
        for(i = 0; i < Nx; i++) {

            phi[j][i] += error[j][i];

        }
    }

    /* Exit GS solving */
    finished = GS_nstep(f, phi, Nx, Ny, epsilon, nGS);


    for (j = 0; j < Ny; j++) {
        free(residual[j]);
        free(error[j]);
        free(laplace_phi[j]);
    }
    for (j = 0; j < half_Ny; j++) {
        free(residual2[j]);
        free(error2[j]);
    }
    free(residual);
    free(residual2);
    free(error);
    free(error2);
    free(laplace_phi);  

    return 0;
}



/* ====================================================================================================================
=================================== Multigrid solver ==================================================================
==================================================================================================================== */


int Multigrid_solver(double** f, double** phi, int Nx, int Ny, int D_x, int D_y, double epsilon, int nGS) {
    
    int i, j;
    double nx = Nx;
    double ny = Ny;
    double lambda = pow(D_x, -2);
    int step = 1;
    int max_num_steps = 2000;
    double laplace_phi_minus_f_norm;
    double RHS;
    double integral;
    double f_norm = 0;
    double epsilon2 = epsilon * pow(10, -2);

    
    double** laplace_phi = (double**)calloc(Ny, sizeof(double*));

    for (i = 0; i < Nx; i++) {
        laplace_phi[i] = (double*)calloc(Nx, sizeof(double));
    }

    /* ------- Computing f_norm ------- */
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            f_norm += pow(f[j][i], 2);
        }
    }

    f_norm = sqrt(f_norm);


    /* ------- RELAXATION ON phi ------- */
    GS_nstep(f, phi, Nx, Ny, epsilon, nGS);

    do {
        

        /* ------- RECURSION ------- */
        MG_recursion(f, phi, Nx, Ny, epsilon2, nGS);
        
        /* compute laplace_p matrix */
        for (j = 1; j < Ny - 1; j++) {
            for (i = 1; i < Nx - 1; i++) {
                laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 4 * lambda * phi[j][i];
            }
        }
        for (j = 1; j < Ny - 1; j++) {
            i = 0;
            laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        for (j = 1; j < Ny - 1; j++) {
            i = Nx - 1;
            laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        for (i = 1; i < Nx - 1; i++) {
            j = 0;
            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }
        for (i = 1; i < Nx - 1; i++) {
            j = Ny - 1;
            laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - phi[j][i] * (3 * lambda);
        }

        laplace_phi[0][0] = (phi[1][0] + phi[0][1]) * lambda - phi[0][0] * (2 * lambda);
        laplace_phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) * lambda - phi[0][Nx - 1] * (2 * lambda);
        laplace_phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) * lambda - phi[Ny - 1][0] * (2 * lambda);
        laplace_phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) * lambda - phi[Ny - 1][Nx - 1] * (2 * lambda);


        /* Compute the norm */
        laplace_phi_minus_f_norm = 0;

        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                laplace_phi_minus_f_norm = laplace_phi_minus_f_norm + pow((laplace_phi[j][i] - f[j][i]), 2);
            }
        }

        
        laplace_phi_minus_f_norm = sqrt(laplace_phi_minus_f_norm);


        if (f_norm == 0) {
            RHS = epsilon;
        }
        else {
            RHS = epsilon * f_norm;
        }


        if (step == max_num_steps) {
            break;
        }


        step += 1;

    } while (laplace_phi_minus_f_norm > RHS);


    /* Impose condition that the integral over the domain is equal to zero */
    integral = 0;
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            integral += phi[j][i] * D_x * D_y;
        }
    }

    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            phi[j][i] = phi[j][i] - integral / (Nx * Ny);
        }
    }

    /* Freeing memory */
    for (i = 0; i < Nx; i++) {
        free(laplace_phi[i]);
    }
    free(laplace_phi);

    
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
    double nGS = 7; /* Number of Gauss-Seidel steps to take per relaxation in mutligrid acceleration */

    double epsilon = pow(10, -3);

    int Nx = 500;
    int Ny = 100;
    double nx = Nx;
    double ny = Ny;
    double dx = H / 20;
    double dy = H / 20;



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

    


    /* -------------------------------------------------------------------------
    -------------------------- Initializing Arrays -----------------------------
    ------------------------------------------------------------------------- */

    double** u = (double**)calloc(ny, sizeof(double*));
    double** v = (double**)calloc(ny, sizeof(double*));
    double** u_star = (double**)calloc(ny, sizeof(double*));
    double** v_star = (double**)calloc(ny, sizeof(double*));

    double** Hx_old = (double**)calloc(ny, sizeof(double*));
    double** Hy_old = (double**)calloc(ny, sizeof(double*));

    double** p = (double**)calloc(ny, sizeof(double*));             //cell centered
    double** grad_u_star_over_dt = (double**)calloc(ny, sizeof(double*));

    double** step1_mat_x = (double**)calloc(ny, sizeof(double*));
    double** step1_mat_y = (double**)calloc(ny, sizeof(double*));


    double* du_s = (double*)calloc(nx*ny, sizeof(double));
    double* dv_s = (double*)calloc(nx*ny, sizeof(double));
    double* du_ss = (double*)calloc(nx*ny, sizeof(double));
    double* dv_ss = (double*)calloc(nx*ny, sizeof(double));

    double* step1_mat_x_vec = (double*)calloc(nx*ny, sizeof(double));
    double* step1_mat_y_vec = (double*)calloc(nx*ny, sizeof(double));


    for (i = 0; i < ny; i++) {
        u[i] = (double*)calloc(nx, sizeof(double));
        v[i] = (double*)calloc(nx, sizeof(double));
        u_star[i] = (double*)calloc(nx, sizeof(double));
        v_star[i] = (double*)calloc(nx, sizeof(double));

        Hx_old[i] = (double*)calloc(nx, sizeof(double));
        Hy_old[i] = (double*)calloc(nx, sizeof(double));

        p[i] = (double*)calloc(nx, sizeof(double));
        grad_u_star_over_dt[i] = (double*)calloc(nx, sizeof(double));

        step1_mat_x[i] = (double*)calloc(nx, sizeof(double));
        step1_mat_y[i] = (double*)calloc(nx, sizeof(double));
        

    }

    

    
    //sample for loop
    for(i = 0; i < Nx; i++) {
        for (j = 0; j < Ny; j++) {
            u[j][i] = 0;
            printf("(%d, %d)", j, i);
        }
        printf("\n");
    }



    /* -------------------------------------------------------------------------
    ----------------------------------- Step 1 ---------------------------------
    ------------------------------------------------------------------------- */

    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            /* Calculating step1_mat_x */
            if (i > 0 && j > 0 && i < (Nx - 1) && j < (Ny - 1)) {  /* At interior points except at the edges  */
                u_cc = (u[j][i] + u[j][i + 1]) / 2;
                u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
            }
            else if (i == 0) {      /* At the left wall  */
                u_cc = (u[j][i] + u[j][i + 1]) / 2;
                u_cc_im1 = U_inlet;

                u_s = U_inlet;
                v_s = 0;
                u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                v_s_jp1 = 0;

                Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
            }
            else if (i == Nx - 1) {     /* At the right wall  */
                u_cc = u[j][i];
                u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
            }
            else if (j == 0) {      /* At the bottom wall  */
                u_cc = (u[j][i] + u[j][i + 1]) / 2;
                u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                u_s = u[j][i];     //this could also be set to zero
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j][i]) / pow(dy, 2)));
            }
            else if (j == Ny - 1) {    /* At the top lid  */
                u_cc = (u[j][i] + u[j][i + 1]) / 2;
                u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_jp1 = 0;
                v_s_jp1 = 0;

                Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
            }






            /* Calculating step1_mat_y */
            if (i > 0 && j > 0 && i < (Nx - 1) && j < (Ny - 1)) {  /* At interior points except at the edges  */
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
            }
            else if (j == 0) {    /* At the bottom wall  */
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                v_cc_jm1 = 0;

                u_s = 0;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_ip1 = 0;
                v_s_ip1 = 0;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + (v[j + 1][i] - 2 * v[j][i] + v[j][i]) / pow(dy, 2)));
            }
            else if (j == Ny - 1) {     /* At the top wall  */
                v_cc = (v[j][i] + 0) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + (v[j][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
            }
            else if (i == 0) {      /* At left wall  */
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = U_inlet;
                v_s = 0;
                u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i]) / pow(dx, 2) + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
            }
            else if (i == Nx - 1) {     /* At the right wall  */
                v_cc = (v[j][i] + v[j + 1][i]) / 2;
                v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                u_s = (u[j][i] + u[j - 1][i]) / 2;
                v_s = (v[j][i] + v[j][i - 1]) / 2;
                u_s_ip1 = u_s;
                v_s_ip1 = v_s;

                Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
            }
        }
    }

    //deal with corner points for step1_mat_x and step1_mat_y


    //vectorize step1_mat_x and step1_mat_y


    //tridiagonal solve for du_ss and dv_ss


    //tridiagonal solve for du_s and dv_s


    //solve for u_star



    /* -------------------------------------------------------------------------
    ---------------------------------- Step 2 ----------------------------------
    ------------------------------------------------------------------------- */

    /* Calculating the cell-centered divergence of velocity_star divided by dt*/
    for (j = 0; j < (Ny - 1); j++) {
        for (i = 0; i < (Nx - 1); i++) {

            grad_u_star_over_dt[j][i] = ((u_star[j][i + 1] - u_star[j][i]) / dx + (v_star[j + 1][i] - v_star[j][i]) / dy) / dt;


        }
    }
    /* Right edge */
    for (j = 0; j < Ny; j++) {
        i = Nx - 1;
    }

    /* Solving for cell-centered pressure using multigrid acceleration method */
    Multigrid_solver(grad_u_star_over_dt, p, Nx, Ny, dx, dy, epsilon, nGS);




    /* -------------------------------------------------------------------------
    ---------------------------------- Step 3 ----------------------------------
    ------------------------------------------------------------------------- */










    /* -------------------------------------------------------------------------
    ----------------------------- Scalar Transport -----------------------------
    ------------------------------------------------------------------------- */



    /* -------------------------------------------------------------------------
    ---------------------------- Freeing Variables -----------------------------
    ------------------------------------------------------------------------- */
    for (i = 0; i < Ny; i++) {
        free(u[i]);
        free(v[i]);
        free(u_star[i]);
        free(v_star[i]);

        free(Hx_old[i]);
        free(Hy_old[i]);

        free(p[i]);
        free(grad_u_star_over_dt[i]);

        free(step1_mat_x[i]);
        free(step1_mat_y[i]);
    }

    free(u);
    free(v);
    free(u_star);
    free(v_star);

    free(Hx_old);
    free(Hy_old);

    free(p);
    free(grad_u_star_over_dt);

    free(step1_mat_x);
    free(step1_mat_y);


    free(du_s);
    free(dv_s);
    free(du_ss);
    free(dv_s);

    free(step1_mat_x_vec);
    free(step1_mat_y_vec);

}