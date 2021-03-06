#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <string.h>


//function used to print the data to a text file
int print_current_data(int step, double** u, double** v, double** phi, int NX, int NY, char method[]) {
    int i, j;


    char filename1[40] = "step_";
    itoa(step, filename1 + 4, 10);
    strcat(filename1, "_");
    strcat(filename1, method);
    strcat(filename1, "_u_data.txt");

    FILE* fpointer1 = fopen(filename1, "w");
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            fprintf(fpointer1, "%.20lf ", u[j][i]);
        }
        fprintf(fpointer1, "\n");
    }
    fclose(fpointer1);



    char filename2[25] = "step_";
    itoa(step, filename2 + 4, 10);
    strcat(filename2, "_");
    strcat(filename2, method);
    strcat(filename2, "_v_data.txt");

    FILE* fpointer2 = fopen(filename2, "w");
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            fprintf(fpointer2, "%.20lf ", v[j][i]);
        }
        fprintf(fpointer2, "\n");
    }
    fclose(fpointer2);



    char filename3[30] = "step_";
    itoa(step, filename3 + 4, 10);
    strcat(filename3, "_");
    strcat(filename3, method);
    strcat(filename3, "_phi_data.txt");

    FILE* fpointer3 = fopen(filename3, "w");
    for (j = 0; j < NY; j++) {
        for (i = 0; i < NX; i++) {
            fprintf(fpointer3, "%.20lf ", phi[j][i]);
        }
        fprintf(fpointer3, "\n");
    }
    fclose(fpointer3);



    return 1;
}


double* TriDiag_GaussElim(int size, double dx_or_dy, double dt, double Re, double* d, int right_outlet, int solving_for_du_s) {

    ////input size should be Nx
    ////should pass in d of whole row except for first point (should be size:(size - 1) )
    ////right_outlet = 1 if the right outlet boundary is to the right of what we are looking at

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

int GS_nstep(double** f, double** phi, int Nx, int Ny, double D_x, double D_y, double epsilon, int nGS, int H_cells, int Bj, int Bi) {

    /* Initilaizations */
    int i, j;
    double nx = Nx;
    double ny = Ny;

    double f_norm;
    double lambda = pow(D_x, -2);
    double RHS;    
    double laplace_phi_minus_f_norm;
    double integral = 0;
    double area = 0;
    double phi_avg;
   
    int step = 1;
    int max_num_steps = nGS;

    double** laplace_phi = (double**)calloc(Ny, sizeof(double*));
    for (j = 0; j < Ny; j++) {
        laplace_phi[j] = (double*)calloc(Nx, sizeof(double));
    }


    /* Solving -------------------------------------------------------------------------------- */

    /* compute f_norm */
    f_norm = 0; 
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {

            if ( ((j < Bj) || ( j > (Bj + H_cells - 1))) || ((i < Bi) || (i > (Bi + H_cells - 1)))  ) {     /* if not inside the void... */
                f_norm = f_norm + pow(f[j][i], 2);
            }

        }
    }
    
    f_norm = sqrt(f_norm);
      
    



    do {

        /* ================================================= CALCULATING PHI ================================================= */

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
                    phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);
                                            
                }
                else if (j == Ny - 1 && i != 0 && i != Nx - 1) {
                    /* Edge calcualtions on top */
                    phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) / 3 - f[j][i] / (3 * lambda);
                        
                }
                else if (i == 0 && j != 0 && j != Ny - 1) {
                    /* Edge calcualtions on left */
                    phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

                }
                else if (i == Nx - 1 && j != 0 && j != Ny - 1) {
                    /* Edge calcualtions on right */
                    phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) / 3 - f[j][i] / (3 * lambda);

                }
                else if (i == 0 && j == 0) {
                    /* 1st corner calculation */
                    phi[0][0] = (phi[1][0] + phi[0][1]) / 2 - f[0][0] / (2 * lambda);

                }
                else if (i == Nx - 1 && j == 0) {
                    /* 2nd corner calculations */
                    phi[0][Nx - 1] = (phi[1][Nx - 1] + phi[0][Nx - 2]) / 2 - f[0][Nx - 1] / (2 * lambda);

                }
                else if (i == 0 && j == Ny - 1) {
                    /* 3rd corner calculations */
                    phi[Ny - 1][0] = (phi[Ny - 1][1] + phi[Ny - 2][0]) / 2 - f[Ny - 1][0] / (2 * lambda);

                }
                else if (i == Nx - 1 && j ==  Ny - 1) {
                    /* 4th corner calculations */
                    phi[Ny - 1][Nx - 1] = (phi[Ny - 1][Nx - 2] + phi[Ny - 2][Nx - 1]) / 2 - f[Ny - 1][Nx - 1] / (2 * lambda);

                }

            }
        }



        /* ================================================= CALCULATING LAPLACE PHI ================================================= */

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
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];

                }
                else if (j == Ny - 1 && i != 0 && i != Nx - 1) {
                    /* Edge calcualtions on top */
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j][i - 1] + phi[j - 1][i]) * lambda - 3 * lambda * phi[j][i];

                }
                else if (i == 0 && j != 0 && j != Ny - 1) {
                    /* Edge calcualtions on left */
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];

                }
                else if (i == Nx - 1 && j != 0 && j != Ny - 1) {
                    /* Edge calcualtions on right */
                    laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i] + phi[j + 1][i]) * lambda - 3 * lambda * phi[j][i];

                }
                else if (i == 0 && j == 0) {
                    /* 1st corner calculation */
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];

                }
                else if (i == Nx - 1 && j == 0) {
                    /* 2nd corner calculations */
                    laplace_phi[j][i] = (phi[j][i - 1] + phi[j + 1][i]) * lambda - 2 * lambda * phi[j][i];

                }
                else if (i == 0 && j == Ny - 1) {
                    /* 3rd corner calculations */
                    laplace_phi[j][i] = (phi[j][i + 1] + phi[j - 1][i]) * lambda - 2 * lambda * phi[j][i];

                }
                else if (i == Nx - 1 && j ==  Ny - 1) {
                    /* 4th corner calculations */
                    laplace_phi[j][i] = (phi[j][i - 1] + phi[j - 1][i]) * lambda - 2 * lambda * phi[j][i];

                }

            }
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

    } while (laplace_phi_minus_f_norm > RHS);





    /* Impose condition that the integral over the domain is equal to zero */
    integral = 0;
    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {
            
            if ( ((j < Bj) || ( j > (Bj + H_cells))) || ((i < Bi) || (i > (Bi + H_cells)))  ) {     /* if not inside the void... */
                integral += phi[j][i] * D_x * D_y;
                area += D_x * D_y;
            }

        }
    }

    phi_avg = integral / area;

    for (j = 0; j < Ny; j++) {
        for (i = 0; i < Nx; i++) {

            if ( ((j < Bj) || ( j > (Bj + H_cells))) || ((i < Bi) || (i > (Bi + H_cells)))  ) {     /* if not inside the void... */
                phi[j][i] -= phi_avg;
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
    double dt = 0.001;
    double H = 1;
    double U_inlet = 1;
    double Diff = pow(10,-4);
    double nGS = 15000; /* Number of Gauss-Seidel steps to take per relaxation in mutligrid acceleration */

    char convective_method[] = "centraldiff";  // options are "upwind", "centraldiff", or "quick"

    double epsilon = pow(10, -3);
    int max_time_steps = 1000;

    int Nx = 500;
    int Ny = 100;
    double nx = Nx;
    double ny = Ny;
    int H_cells = 20;
    double dx = H / H_cells;
    double dy = H / H_cells;

    int Bj = 40;
    int Bi = 100;
   
    printf("\n \n");

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

    /* Initializing variable for step 3 */
    double u_out_total = 0;
    double u_offset;
    double area_out = ny * dx * dy;
    double u_out_avg;

    /* For scalar transport */
    double convec;
    double diffu;

    int print_now;


    /* -------------------------------------------------------------------------
    -------------------------- Initializing Arrays -----------------------------
    ------------------------------------------------------------------------- */

    double** u = (double**)calloc(ny, sizeof(double*));
    double** v = (double**)calloc(ny, sizeof(double*));

    double** Hx_old = (double**)calloc(ny, sizeof(double*));
    double** Hy_old = (double**)calloc(ny, sizeof(double*));

    double** p = (double**)calloc(ny, sizeof(double*));             //cell centered
    double** phi = (double**)calloc(ny, sizeof(double*));
    double** phi_new = (double**)calloc(ny, sizeof(double*));
    double** grad_u_star_over_dt = (double**)calloc(ny, sizeof(double*));

    double** step1_mat_x = (double**)calloc(ny, sizeof(double*));
    double** step1_mat_y = (double**)calloc(ny, sizeof(double*));


    double** du_s = (double**)calloc(ny, sizeof(double*));
    double** dv_s = (double**)calloc(ny, sizeof(double*));
    double** du_ss = (double**)calloc(ny, sizeof(double*));
    double** dv_ss = (double**)calloc(ny, sizeof(double*));



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

    for (iter = 0; iter < max_time_steps; iter++) {

        /* -------------------------------------------------------------------------
        ----------------------------------- Step 1 ---------------------------------
        ------------------------------------------------------------------------- */


        //--//--//--//--//-- Solve for step1_mat_x and step1_mat_y along boundaries and interior points --//--//--//--//--//

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
               
                    Hx_old[j][i] = Hx;
                }
                else if (i == 0 && j != 0 && j != Ny - 1) {      /* At the left wall  */ //////checked
                    u_cc = (u[j][i] + u[j][i + 1]) / 2;
                    u_cc_im1 = U_inlet;     ///technically don't use this value below

                    u_s = U_inlet;
                    v_s = 0;
                    u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                    v_s_jp1 = 0;

                    Hx = ((pow(u_cc, 2) - pow(U_inlet, 2)) * 2 / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                    step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * (          ( (u[j][i+1] - u[j][i]) / dx - (u_cc - U_inlet) * 2 / dx )/dx                    + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
               
                    Hx_old[j][i] = Hx;
                }
                else if (i == Nx - 1 && j != 0 && j != Ny - 1) {     /* At the right wall  */  //////checked
                    u_cc = u[j][i];
                    u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                    u_s = (u[j][i] + u[j - 1][i]) / 2;
                    v_s = (v[j][i] + v[j][i - 1]) / 2;
                    u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                    v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                    Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                    step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
                    
                    Hx_old[j][i] = Hx;
                }
                else if (j == 0 && i != 0 && i != Nx - 1) {      /* At the bottom wall  */   ///////checked
                    u_cc = (u[j][i] + u[j][i + 1]) / 2;
                    u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                    u_s = 0;    
                    v_s = 0;
                    u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
                    v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

                    Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                    step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) +    ( (u[j+1][i] - u[j][i]) / dy - (u[j][i] - u_s) *(2/dy) ) / dy   ));
               
                    Hx_old[j][i] = Hx;
                }
                else if (j == Ny - 1 && i != 0 && i != Nx - 1) {    /* At the top wall  */    ///////checked
                    u_cc = (u[j][i] + u[j][i + 1]) / 2;
                    u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

                    u_s = (u[j][i] + u[j - 1][i]) / 2;
                    v_s = (v[j][i] + v[j][i - 1]) / 2;
                    u_s_jp1 = 0;
                    v_s_jp1 = 0;

                    Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

                    step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) +         ( (u_s_jp1 - u[j][i]) *(2/dy) - (u[j][i] - u[j-1][i]) / dy) / dy        ));
               
                    Hx_old[j][i] = Hx;
                }






                /* Calculating step1_mat_y */
                if (i > 0 && j > 0 && i < (Nx - 1) && j < (Ny - 1)) {  // At interior points except at the edges  
                    v_cc = (v[j][i] + v[j + 1][i]) / 2;
                    v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                    u_s = (u[j][i] + u[j - 1][i]) / 2;
                    v_s = (v[j][i] + v[j][i - 1]) / 2;
                    u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                    v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                    Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                    step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
               
                    Hy_old[j][i] = Hy;
                }
                else if (j == 0 && i != 0 && i != Nx - 1) {    /* At the bottom wall  */    /////corrected
                    v_cc = (v[j][i] + v[j + 1][i]) / 2;
                    v_cc_jm1 = 0;

                    u_s = 0;
                    v_s = 0;
                    u_s_ip1 = 0;
                    v_s_ip1 = 0;

                    Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2))*2 / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                    step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + ( (v[j+1][i] - v[j][i])/dy - (v_cc - v[j][i])*2/dy ) /dy        ));
               
                    Hy_old[j][i] = Hy;
                }
                else if (j == Ny - 1 && i != 0 && i != Nx - 1) {     /* At the top wall  */    ///////corrected
                    v_cc = (v[j][i] + 0) / 2;
                    v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                    u_s = (u[j][i] + u[j - 1][i]) / 2;
                    v_s = (v[j][i] + v[j][i - 1]) / 2;
                    u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                    v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                    Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                    step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) +     (0 - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
               
                    Hy_old[j][i] = Hy;
                }
                else if (i == 0 && j != 0 && j != Ny - 1) {      /* At left wall  */  ///////corrected
                    v_cc = (v[j][i] + v[j + 1][i]) / 2;
                    v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                    u_s = U_inlet;
                    v_s = 0;
                    u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
                    v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

                    Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                    step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) *          (     ( (v[j][i+1] - v[j][i])/dx - (v[j][i] - v_s)*2/dx )  / dx         + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
               
                    Hy_old[j][i] = Hy;
                }
                else if (i == Nx - 1 && j != 0 && j != Ny - 1) {     /* At the right wall  */    /////corrected
                    v_cc = (v[j][i] + v[j + 1][i]) / 2;
                    v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

                    u_s = (u[j][i] + u[j - 1][i]) / 2;
                    v_s = (v[j][i] + v[j][i - 1]) / 2;
                    u_s_ip1 = u_s;
                    v_s_ip1 = v_s;

                    Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

                    step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (     (v[j][i] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2)        + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
               
                    Hy_old[j][i] = Hy;
                }
            }
        }



        //--//--//--//--//-- deal with corner points for step1_mat_x and step1_mat_y --//--//--//--//--//--//

        /////////////////top left corner/////////////////////
        i = 0;
        j = Ny - 1;

        //x data
        u_cc = (u[j][i] + u[j][i + 1]) / 2;
        u_cc_im1 = U_inlet;    

        u_s = U_inlet;
        v_s = 0;
        u_s_jp1 = 0;
        v_s_jp1 = 0;

        Hx = ((pow(u_cc, 2) - pow(U_inlet, 2)) * 2 / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

        step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * (          ( (u[j][i+1] - u[j][i]) / dx - (u_cc - U_inlet) * 2 / dx )/dx                    +  ( (u_s_jp1 - u[j][i]) *(2/dy) - (u[j][i] - u[j-1][i]) / dy) / dy    ));

        Hx_old[j][i] = Hx;
        

        //y data
        v_cc = (v[j][i] + 0) / 2;
        v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

        u_s = U_inlet;
        v_s = 0;
        u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
        v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

        Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

        step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (    ( (v[j][i+1] - v[j][i])/dx - (v[j][i] - v_s)*2/dx )  / dx +     (0 - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));

        Hy_old[j][i] = Hy;



        /////////////////top right corner/////////////////
        i = Nx - 1;
        j = Ny - 1;

        //x data
        u_cc = u[j][i];
        u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

        u_s = (u[j][i] + u[j - 1][i]) / 2;
        v_s = (v[j][i] + v[j][i - 1]) / 2;
        u_s_jp1 = 0;
        v_s_jp1 = 0;

        Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

        step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * (       (u[j][i] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2)   +     ( (u_s_jp1 - u[j][i]) *(2/dy) - (u[j][i] - u[j-1][i]) / dy) / dy      ));

        Hx_old[j][i] = Hx;
               

        //y data
        v_cc = (v[j][i] + 0) / 2;
        v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

        u_s = (u[j][i] + u[j - 1][i]) / 2;
        v_s = (v[j][i] + v[j][i - 1]) / 2;
        u_s_ip1 = u_s;
        v_s_ip1 = v_s;

        Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

        step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (      (v[j][i] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2)       +     (0 - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));

        Hy_old[j][i] = Hy;



        /////////////////bottom left corner/////////////////
        i = 0;
        j = 0;

        //x data
        u_cc = (u[j][i] + u[j][i + 1]) / 2;
        u_cc_im1 = U_inlet;    

        u_s = 0;
        v_s = 0;
        u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
        v_s_jp1 = 0;

        Hx = ((pow(u_cc, 2) - pow(U_inlet, 2)) * 2 / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

        step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * (          ( (u[j][i+1] - u[j][i]) / dx - (u_cc - U_inlet) * 2 / dx )/dx          +   ( (u[j+1][i] - u[j][i]) / dy - (u[j][i] - u_s) *(2/dy) ) / dy  ));

        Hx_old[j][i] = Hx;
               

        //y data
        v_cc = (v[j][i] + v[j + 1][i]) / 2;
        v_cc_jm1 = 0;

        u_s = 0;
        v_s = 0;
        u_s_ip1 = 0;
        v_s_ip1 = 0;

        Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2))*2 / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

        step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (    ( (v[j][i+1] - v[j][i])/dx - (v[j][i] - v_s)*2/dx )  / dx      + ( (v[j+1][i] - v[j][i])/dy - (v_cc - v[j][i])*2/dy ) /dy        ));

        Hy_old[j][i] = Hy;


        /////////////////bottom right corner/////////////////
        i = Nx - 1;
        j = 0;

        //x data
        u_cc = u[j][i];
        u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

        u_s = 0;
        v_s = 0;
        u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
        v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

        Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

        step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * (   (u[j][i] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2)    +       ( (u[j+1][i] - u[j][i]) / dy - (u[j][i] - u_s) *(2/dy) ) / dy   ));

        Hx_old[j][i] = Hx;
               

        //y data
        v_cc = (v[j][i] + v[j + 1][i]) / 2;
        v_cc_jm1 = 0;

        u_s = 0;
        v_s = 0;
        u_s_ip1 = 0;
        v_s_ip1 = 0;

        Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2))*2 / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

        step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (    (v[j][i] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) + ( (v[j+1][i] - v[j][i])/dy - (v_cc - v[j][i])*2/dy ) /dy        ));

        Hy_old[j][i] = Hy;
       

       
        //--//--//--//--//-- Deal with boundary conditions for square obstacle --//--//--//--//--//

        /////////////////////Top boundary of square/////////////////////////

        // u barycenter

        j = 60;
        for (i = 100; i < 120; i++) {
            u_cc = (u[j][i] + u[j][i + 1]) / 2;
            u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

            u_s = 0;    
            v_s = 0;
            u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
            v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

            Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

            step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) +    ( (u[j+1][i] - u[j][i]) / dy - (u[j][i] - u_s) *(2/dy) ) / dy   ));
               
            Hx_old[j][i] = Hx;
        }

        /////////////////////Right boundary of square/////////////////////////

        // v barycenter
        i = 120;
        for (j = 40; j < 60; j++) {
            v_cc = (v[j][i] + v[j + 1][i]) / 2;
            v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

            u_s = 0;
            v_s = 0;
            u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
            v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

            Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

            step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) *          (     ( (v[j][i+1] - v[j][i])/dx - (v[j][i] - v_s)*2/dx )  / dx         + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
       
            Hy_old[j][i] = Hy;
        }

        /////////////////////Bottom boundary of square/////////////////////////

        //u barycenter
        j = 39;
        for (i = 100; i < 120; i++) {
            u_cc = (u[j][i] + u[j][i + 1]) / 2;
            u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

            u_s = (u[j][i] + u[j - 1][i]) / 2;
            v_s = (v[j][i] + v[j][i - 1]) / 2;
            u_s_jp1 = 0;
            v_s_jp1 = 0;

            Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

            step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) +         ( (u_s_jp1 - u[j][i]) *(2/dy) - (u[j][i] - u[j-1][i]) / dy) / dy        ));
           
            Hx_old[j][i] = Hx;
               
        }

        //v barycenter
        j = 39;
        for (i = 100; i < 120; i++) {
            v_cc = (v[j][i] + 0) / 2;
            v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

            u_s = (u[j][i] + u[j - 1][i]) / 2;
            v_s = (v[j][i] + v[j][i - 1]) / 2;
            u_s_ip1 = (u[j][i + 1] + u[j - 1][i + 1]) / 2;
            v_s_ip1 = (v[j][i + 1] + v[j][i]) / 2;

            Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

            step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * ((v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2) +     (0 - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
       
            Hy_old[j][i] = Hy;
        }



        /////////////////////Left boundary of square/////////////////////////

        //u barycenter
        i = 99;
        for (j = 40; j < 60; j++) {
            u_cc = (u[j][i] + 0 ) / 2;
            u_cc_im1 = (u[j][i - 1] + u[j][i]) / 2;

            u_s = (u[j][i] + u[j - 1][i]) / 2;
            v_s = (v[j][i] + v[j][i - 1]) / 2;
            u_s_jp1 = (u[j + 1][i] + u[j][i]) / 2;
            v_s_jp1 = (v[j + 1][i] + v[j + 1][i - 1]) / 2;

            Hx = ((pow(u_cc, 2) - pow(u_cc_im1, 2)) / dx) + (u_s_jp1 * v_s_jp1 - u_s * v_s) / dy;

            step1_mat_x[j][i] = dt * (Hx * 3 / 2 - Hx_old[j][i] / 2 + (1 / Re) * ((u[j][i + 1] - 2 * u[j][i] + u[j][i - 1]) / pow(dx, 2) + (u[j + 1][i] - 2 * u[j][i] + u[j - 1][i]) / pow(dy, 2)));
       
            Hx_old[j][i] = Hx;
        }


        //v barycenter
        i = 99;
        for (j = 40; j < 60; j++) {
            v_cc = (v[j][i] + v[j + 1][i]) / 2;
            v_cc_jm1 = (v[j - 1][i] + v[j][i]) / 2;

            u_s = (u[j][i] + u[j - 1][i]) / 2;
            v_s = (v[j][i] + v[j][i - 1]) / 2;
            u_s_ip1 = 0;
            v_s_ip1 = 0;

            Hy = ((pow(v_cc, 2) - pow(v_cc_jm1, 2)) / dy) + (u_s_ip1 * v_s_ip1 - u_s * v_s) / dx;

            step1_mat_y[j][i] = dt * (Hy * 3 / 2 - Hy_old[j][i] / 2 + (1 / Re) * (     (v[j][i + 1] - 2 * v[j][i] + v[j][i - 1]) / pow(dx, 2)        + (v[j + 1][i] - 2 * v[j][i] + v[j - 1][i]) / pow(dy, 2)));
       
            Hy_old[j][i] = Hy;
        }


        //////////// set values inside the square equal to zero for good measure

        for (j = 40; j < 60; j++) {
            for (i = 100; i < 120; i++) {
                step1_mat_x[j][i] = 0;
                step1_mat_y[j][i] = 0;
            }
        }
       
       


        //--//--//--//--//-- Tridiagonal solve for du_ss --//--//--//--//--//

        //rows below the bottom of the square
        for (j = 0; j < 40; j++) {

            //fill in temp_vec with row of pixels except first pixel
            for(i = 0; i < Nx - 1; i++) {
                temp_vec_x_long[i] = step1_mat_x[j][i+1];
            }
   
            //gaussian elimination
            TriDiag_GaussElim(Nx, dx, dt, Re, temp_vec_x_long, 1, 0);

            //update du_ss
            for(i = 0; i < Nx - 1; i++) {
                du_ss[j][i+ 1] = temp_vec_x_long[i];
            }
           
        }


        //rows to left and right of square
        for (j = 40; j < 60; j++) {

            /////////left
           
            for(i = 0; i < 99; i++) {
                temp_vec_x_small[i] = step1_mat_x[j][i+1];
            }
           
            TriDiag_GaussElim(100, dx, dt, Re, temp_vec_x_small, 0, 0);

            for(i = 0; i < 99; i++) {
                du_ss[j][i+ 1] = temp_vec_x_small[i];
            }


            ////////right

            for(i = 0; i < 19 * 20 - 1; i++) {
                temp_vec_x_medium[i] = step1_mat_x[j][i+121];
            }
   
            TriDiag_GaussElim(19*20, dx, dt, Re, temp_vec_x_medium, 1, 0);

            for(i = 0; i < 19 * 20 - 1; i++) {
                du_ss[j][i+ 121] = temp_vec_x_medium[i];
            }
           

        }

       

        //rows above the top of the square
        for (j = 60; j < Ny; j++) {

            for(i = 0; i < Nx - 1; i++) {
                temp_vec_x_long[i] = step1_mat_x[j][i+1];
            }

            TriDiag_GaussElim(Nx, dx, dt, Re, temp_vec_x_long, 1, 0);

            for(i = 0; i < Nx - 1; i++) {
                du_ss[j][i+ 1] = temp_vec_x_long[i];
            }

        }

       
       

        //--//--//--//--//-- Tridiagonal solve for dv_ss --//--//--//--//--//
       
        //rows below the bottom of the square
        for (j = 1; j < 40; j++) {

            //fill in temp_vec with row of pixels except first pixel
            for(i = 0; i < Nx - 1; i++) {
                temp_vec_x_long[i] = step1_mat_y[j][i+1];
            }
   
            //gaussian elimination
            TriDiag_GaussElim(Nx, dx, dt, Re, temp_vec_x_long, 1, 0);

            //update dv_ss
            for(i = 0; i < Nx - 1; i++) {
                dv_ss[j][i+ 1] = temp_vec_x_long[i];
            }
           
        }


        //rows to left and right of square
        for (j = 40; j < 60; j++) {

            /////////left
           
            for(i = 0; i < 99; i++) {
                temp_vec_x_small[i] = step1_mat_y[j][i+1];
            }
           
            TriDiag_GaussElim(100, dx, dt, Re, temp_vec_x_small, 0, 0);

            for(i = 0; i < 99; i++) {
                dv_ss[j][i+ 1] = temp_vec_x_small[i];
            }


            ////////right (this one was changed to incorporate the values at the boundary)

            for(i = 0; i < 19 * 20; i++) {
                temp_vec_x_medium2[i] = step1_mat_y[j][i+120];
            }
   
            TriDiag_GaussElim(19*20 + 1, dx, dt, Re, temp_vec_x_medium2, 1, 0);

            for(i = 0; i < 19 * 20; i++) {
                dv_ss[j][i+ 120] = temp_vec_x_medium2[i];
            }
           

        }


        //rows above the top of the square
        for (j = 60; j < Ny; j++) {

            for(i = 0; i < Nx - 1; i++) {
                temp_vec_x_long[i] = step1_mat_y[j][i+1];
            }

            TriDiag_GaussElim(Nx, dx, dt, Re, temp_vec_x_long, 1, 0);

            for(i = 0; i < Nx - 1; i++) {
                dv_ss[j][i+ 1] = temp_vec_x_long[i];
            }

        }



       
        //--//--//--//--//-- Tridiagonal solve for du_s --//--//--//--//--//

        //columns to the left of the square
        for (i = 1; i < 100; i++) {

            for (j = 0; j < Ny; j++) {
                temp_vec_y_long[j] = du_ss[j][i];
            }

            TriDiag_GaussElim(Ny + 1, dy, dt, Re, temp_vec_y_long, 0, 1);

            for (j = 0; j < Ny; j++) {
                du_s[j][i] = temp_vec_y_long[j];
            }

        }


       

        //columns above and below the square

        for (i = 100; i < 120; i++) {

            //below
            for (j = 0; j < 40; j++) {
                temp_vec_y_small[j] = du_ss[j][i];
            }

            TriDiag_GaussElim(41, dy, dt, Re, temp_vec_y_small, 0, 1);

            for (j = 0; j < 40; j++) {
                du_s[j][i] = temp_vec_y_small[j];
            }

            //above
            for (j = 0; j < 40; j++) {
                temp_vec_y_small[j] = du_ss[j + 60][i];
            }

            TriDiag_GaussElim(41, dy, dt, Re, temp_vec_y_small, 0, 1);

            for (j = 0; j < 40; j++) {
                du_s[j + 60][i] = temp_vec_y_small[j];
            }

        }




        //columns to the right of the square
        for (i = 120; i < Nx; i++) {

            for (j = 0; j < Ny; j++) {
                temp_vec_y_long[j] = du_ss[j][i];
            }

            TriDiag_GaussElim(Ny + 1, dy, dt, Re, temp_vec_y_long, 0, 1);

            for (j = 0; j < Ny; j++) {
                du_s[j][i] = temp_vec_y_long[j];
            }

        }

       

        //--//--//--//--//-- Tridiagonal solve for dv_s --//--//--//--//--//

        //columns to the left of the square
        for (i = 1; i < 100; i++) {

            for (j = 0; j < Ny - 1; j++) {
                temp_vec_y_long2[j] = dv_ss[j + 1][i];
            }

            TriDiag_GaussElim(Ny, dy, dt, Re, temp_vec_y_long2, 0, 0);

            for (j = 0; j < Ny - 1; j++) {
                dv_s[j + 1][i] = temp_vec_y_long2[j];
            }

        }


        //columns above and below the square

        for (i = 100; i < 120; i++) {

            //below
            for (j = 0; j < 39; j++) {
                temp_vec_y_small2[j] = dv_ss[j + 1][i];
            }

            TriDiag_GaussElim(40, dy, dt, Re, temp_vec_y_small2, 0, 0);

            for (j = 0; j < 39; j++) {
                dv_s[j + 1][i] = temp_vec_y_small2[j];
            }



            //above
            for (j = 0; j < 39; j++) {
                temp_vec_y_small2[j] = dv_ss[j + 61][i];
            }

            TriDiag_GaussElim(40, dy, dt, Re, temp_vec_y_small2, 0, 0);

            for (j = 0; j < 39; j++) {
                dv_s[j + 61][i] = temp_vec_y_small2[j];
            }


        }


        //columns to the right of the square
        for (i = 120; i < Nx; i++) {

            for (j = 0; j < Ny - 1; j++) {
                temp_vec_y_long2[j] = dv_ss[j + 1][i];
            }

            TriDiag_GaussElim(Ny, dy, dt, Re, temp_vec_y_long2, 0, 0);

            for (j = 0; j < Ny - 1; j++) {
                dv_s[j + 1][i] = temp_vec_y_long2[j];
            }

        }

       
       
       
        //--//--//--//--//-- Update u_star and v_star --//--//--//--//--//

        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                u_star[j][i] = u[j][i] + du_s[j][i];
                v_star[j][i] = v[j][i] + dv_s[j][i];
            }
        }

        // update right boundary of u_star
        for (j = 0; j < Ny; j++) {
            u_star[j][Nx] = u_star[j][Nx - 1];
        }




        /* -------------------------------------------------------------------------
        ---------------------------------- Step 2 ----------------------------------
        ------------------------------------------------------------------------- */

        /* Calculating the cell-centered divergence of velocity_star divided by dt*/
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
               
                grad_u_star_over_dt[j][i] = ((u_star[j][i + 1] - u_star[j][i]) / dx + (v_star[j + 1][i] - v_star[j][i]) / dy) / dt;


            }
        }

        FILE* fid_gusod = fopen("gusod.txt", "w");
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                fprintf(fid_gusod, "%.20lf ", grad_u_star_over_dt[j][i]);
            }
            fprintf(fid_gusod, "\n");
        }
        fclose(fid_gusod);

        FILE* fid_p = fopen("pressure.txt", "w");
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                fprintf(fid_p, "%.20lf ", p[j][i]);
            }
            fprintf(fid_p, "\n");
        }
        fclose(fid_p);


        /* Solving for cell-centered pressure using multigrid acceleration method */
        GS_nstep(grad_u_star_over_dt, p, Nx, Ny, dx, dy, epsilon, nGS, H_cells, Bj, Bi);



        /* -------------------------------------------------------------------------
        ---------------------------------- Step 3 ----------------------------------
        ------------------------------------------------------------------------- */

        FILE* fid_u = fopen("u.txt", "w");
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                fprintf(fid_u, "%.20lf ", u[j][i]);
            }
            fprintf(fid_u, "\n");
        }
        fclose(fid_u);

        FILE* fid_v = fopen("v.txt", "w");
        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                fprintf(fid_v, "%.20lf ", v[j][i]);
            }
            fprintf(fid_v, "\n");
        }
        fclose(fid_v);


        /* ------- Updating all velocities ------- */
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
                else if (j == 0 && i !=0 && i != Nx - 1) {
                    /* On the domain's bottom edge */
                    u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                    v[j][i] = 0;
                }
                else if (j == Ny - 1 && i !=0 && i != Nx - 1) {
                    /* On the domain's top edge */
                    u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                    v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;
                }
                else if (i == 0 && j == 0) {
                    /* On the domain's bottom left corner */
                    u[j][i] = 1;
                    v[j][i] = 0;
                }
                else if (i == Nx - 1 && j == 0) {
                    /* On the domain's bottom right corner */
                    u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                    v[j][i] = 0;
                }
                else if (i == 0 && j == Ny - 1) {
                    /* On the domain's top left corner */
                    u[j][i] = 1;
                    v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;
                }
                else if (i == Nx - 1 && j == Ny - 1) {
                    /* On the domain's top right corner */
                    u[j][i] = u_star[j][i] - dt * (p[j][i] - p[j][i - 1]) / dx;
                    v[j][i] = v_star[j][i] - dt * (p[j][i] - p[j - 1][i]) / dy;
                }
                else if (i == Bi - 1 && j >= Bj && j <= Bj + H_cells - 1) {
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
                else if (j == Bj - 1 && i >= Bi && i <= Bi + H_cells - 1) {
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


        
        /* ------- Imposing Conservation of Mass ------- */
        for (j = 0; j < Ny; j++) {
            u_out_total += u[j][Nx - 1] * dx * dy;
        }

        u_out_avg = u_out_total / area_out;
        
        u_offset = U_inlet - u_out_avg;

        for (j = 0; j < Ny; j++) {
            u[j][Nx - 1] += u_offset;
        }
        printf("u_out_avg = %f\n", u_out_avg);
        printf("u_offset = %f\n", u_offset);



        /* -------------------------------------------------------------------------
        ----------------------------- Scalar Transport -----------------------------
        ------------------------------------------------------------------------- */


        if ( strcmp(convective_method, "upwind") == 0 ) {
            printf("Step %d of upwind method to deal with the scalar transport \n", iter);

            //interior points

            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {

                    diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                    if (j > Ny / 2) {
                        convec =  ( u[j][i+1] * phi[j][i] - u[j][i] * phi[j][i-1]) / dx   +   (v[j + 1][i] * phi[j][i] - v[j][i] * phi[j - 1][i]) / dy;
                    }
                    else {
                        convec =  ( u[j][i+1] * phi[j][i] - u[j][i] * phi[j][i-1]) / dx   +   (v[j + 1][i] * phi[j+1][i] - v[j][i] * phi[j][i]) / dy;
                    }


                    phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

                }
            }

           
            //top boundary points except corners
            j = Ny - 1;

            for (i = 1; i < Nx - 1; i++) {

                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   ( - phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                convec =  ( u[j][i+1] * phi[j][i] - u[j][i] * phi[j][i-1]) / dx   +   (0 * phi[j][i] - v[j][i] * phi[j - 1][i]) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

            }


            //right boundary points except corners
            i = Nx - 1;

            for (j = 1; j < Ny - 1; j++) {

                diffu = Diff * (   (phi[j][i] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                convec =  ( u[j][i] * phi[j][i] - u[j][i] * phi[j][i-1]) / dx   +   (v[j + 1][i] * phi[j][i] - v[j][i] * phi[j - 1][i]) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

            }


            //left boundary points
            i = 0;

            for (j = 0; j < Ny; j++) {
                phi_new[j][i] = exp( -1 * pow(j * dy - 2.5 * H + 0.5 * dy ,2) );
            }


            //bottom boundary points
            j = 0;

            for (i = 1; i < Nx - 1; i++) {

                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );

                convec =  ( u[j][i+1] * phi[j][i] - u[j][i] * phi[j][i-1]) / dx   +   (v[j + 1][i] * phi[j][i] - v[j][i] * 0) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

            }
   

            //top right corner

            j = Ny - 1;
            i = Nx - 1;

            diffu = Diff * (   (- phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (- phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
            convec =  ( u[j][i] * phi[j][i] - u[j][i] * phi[j][i-1]) / dx   +   (0 * phi[j][i] - v[j][i] * phi[j - 1][i]) / dy;
            phi_new[j][i] = phi[j][i] + dt * (diffu - convec);


            //bottom right corner

            j = 0;
            i = Nx - 1;

            diffu = Diff * (   ( - phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );
            convec =  ( u[j][i] * phi[j][i] - u[j][i] * phi[j][i-1]) / dx   +   (v[j + 1][i] * phi[j][i] - v[j][i] * 0) / dy;
            phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

           


            //boundary of square

            //bottom
            j = 39;
            for (i = 100; i < 120; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   ( - phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * phi[j][i] - u[j][i] * phi[j][i-1]) / dx   +   (0 * phi[j][i] - v[j][i] * phi[j - 1][i]) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

            //top
            j = 60;
            for (i = 100; i < 120; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );
                convec =  ( u[j][i+1] * phi[j][i] - u[j][i] * phi[j][i-1]) / dx   +   (v[j + 1][i] * phi[j][i] - v[j][i] * 0) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

           

            //left
            i = 99;
            for (j = 40; j < 60; j++) {
                diffu = Diff * (   ( - phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( 0 * phi[j][i] - u[j][i] * phi[j][i-1]) / dx   +   (v[j + 1][i] * phi[j][i] - v[j][i] * phi[j - 1][i]) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

           

            //right
            i = 120;
            for (j = 40; j < 60; j++) {
                diffu = Diff * (   (phi[j][i+1] - phi[j][i] ) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * phi[j][i] - u[j][i] * phi[j][i-1]) / dx   +   (v[j + 1][i] * phi[j][i] - v[j][i] * 0) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

           
            //inside square, set phi_new = 0

            for (j = 40; j < 60; j++) {
                for (i = 100; i < 120; i++) {
                    phi_new[j][i] = 0;
                }
            }


        }





        if ( strcmp(convective_method, "centraldiff") == 0 ) {
            printf("Step %d of contral diff method to deal with the scalar transport \n", iter);
            //interior points

            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {

                    diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                    convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;

                    phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

                }
            }

           
            //top boundary points except corners
            j = Ny - 1;

            for (i = 1; i < Nx - 1; i++) {

                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   ( - phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (0 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

            }


            //right boundary points except corners
            i = Nx - 1;

            for (j = 1; j < Ny - 1; j++) {

                diffu = Diff * (   (phi[j][i] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                convec =  ( u[j][i] * (phi[j][i] + phi[j][i]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

            }


            //left boundary points
            i = 0;

            for (j = 0; j < Ny; j++) {
                phi_new[j][i] = exp( -1 * pow(j * dy - 2.5 * H + 0.5 * dy ,2) );
            }


            //bottom boundary points
            j = 0;

            for (i = 1; i < Nx - 1; i++) {

                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );

                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - 0 ) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

            }
   

            //top right corner

            j = Ny - 1;
            i = Nx - 1;

            diffu = Diff * (   (- phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (- phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
            convec =  ( u[j][i] * (phi[j][i] + phi[j][i]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (  0  - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
            phi_new[j][i] = phi[j][i] + dt * (diffu - convec);


            //bottom right corner

            j = 0;
            i = Nx - 1;

            diffu = Diff * (   ( - phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );
            convec =  ( u[j][i] * (phi[j][i] + phi[j][i]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - 0 ) / dy;
            phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

           


            //boundary of square

            //bottom
            j = 39;
            for (i = 100; i < 120; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   ( - phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (0 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

            //top
            j = 60;
            for (i = 100; i < 120; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - 0 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }
           

            //left
            i = 99;
            for (j = 40; j < 60; j++) {
                diffu = Diff * (   ( - phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  (   0  - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

           

            //right
            i = 120;
            for (j = 40; j < 60; j++) {
                diffu = Diff * (   (phi[j][i+1] - phi[j][i] ) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - 0 ) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

           
            //inside square, set phi_new = 0

            for (j = 40; j < 60; j++) {
                for (i = 100; i < 120; i++) {
                    phi_new[j][i] = 0;
                }
            }


        }





        if ( strcmp(convective_method, "quick") == 0 ) {
            printf("Now commencing the quick method to deal with the scalar transport \n");

            //interior points

            for (j = 1; j < Ny - 1; j++) {
                for (i = 1; i < Nx - 1; i++) {

                    diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                    if (j > Ny / 2) {
                        convec =  ( u[j][i+1] * (- phi[j][i-1]/6 + phi[j][i] * 5 /6 + phi[j][i+1]*2/6) - u[j][i] * (- phi[j][i-2]/6 + phi[j][i-1] * 5 /6 + phi[j][i]*2/6)  ) / dx   +   (v[j + 1][i] * ( phi[j+1][i]*2/6 + phi[j][i]*5/6 - phi[j-1][i]/6 ) - v[j][i] *  (phi[j][i]*2/6 + phi[j-1][i]*5/6 - phi[j-2][i]/6)  ) / dy;
                    }
                    else {
                        convec =  ( u[j][i+1] * (- phi[j][i-1]/6 + phi[j][i] * 5 /6 + phi[j][i+1]*2/6) - u[j][i] * (- phi[j][i-2]/6 + phi[j][i-1] * 5 /6 + phi[j][i]*2/6)  ) / dx   +   (v[j + 1][i] * (-phi[j+2][i]/6 + phi[j+1][i]*5/6 + phi[j][i]*2/6) - v[j][i] * (-phi[j+1][i]/6 + phi[j][i]*5/6 + phi[j-1][i]*2/6)  ) / dy;
                    }


                    phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

                }
            }


            //////////////////just replaced with second order central differencing for code below

           
            //top boundary points except corners
            j = Ny - 1;

            for (i = 1; i < Nx - 1; i++) {

                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   ( - phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (0 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

            }


            //right boundary points except corners
            i = Nx - 1;

            for (j = 1; j < Ny - 1; j++) {

                diffu = Diff * (   (phi[j][i] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );

                convec =  ( u[j][i] * (phi[j][i] + phi[j][i]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

            }


            //left boundary points
            i = 0;

            for (j = 0; j < Ny; j++) {
                phi_new[j][i] = exp( -1 * pow(j * dy - 2.5 * H + 0.5 * dy ,2) );
            }


            //bottom boundary points
            j = 0;

            for (i = 1; i < Nx - 1; i++) {

                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );

                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - 0 ) / dy;

                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

            }
   

            //top right corner

            j = Ny - 1;
            i = Nx - 1;

            diffu = Diff * (   (- phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (- phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
            convec =  ( u[j][i] * (phi[j][i] + phi[j][i]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (  0  - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
            phi_new[j][i] = phi[j][i] + dt * (diffu - convec);


            //bottom right corner

            j = 0;
            i = Nx - 1;

            diffu = Diff * (   ( - phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );
            convec =  ( u[j][i] * (phi[j][i] + phi[j][i]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - 0 ) / dy;
            phi_new[j][i] = phi[j][i] + dt * (diffu - convec);

           


            //boundary of square

            //bottom
            j = 39;
            for (i = 100; i < 120; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   ( - phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (0 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

            //top
            j = 60;
            for (i = 100; i < 120; i++) {
                diffu = Diff * (   (phi[j][i+1] - 2*phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - phi[j][i] )/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - 0 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }
           

            //left
            i = 99;
            for (j = 40; j < 60; j++) {
                diffu = Diff * (   ( - phi[j][i] + phi[j][i-1]) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  (   0  - u[j][i] * (phi[j][i-1] + phi[j][i])/2) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

           

            //right
            i = 120;
            for (j = 40; j < 60; j++) {
                diffu = Diff * (   (phi[j][i+1] - phi[j][i] ) / (pow(dx,2))   +   (phi[j + 1][i] - 2*phi[j][i] + phi[j-1][i])/(pow(dy,2))    );
                convec =  ( u[j][i+1] * (phi[j][i] + phi[j][i+1]) / 2 - 0 ) / dx   +   (v[j + 1][i] * (phi[j][i] + phi[j + 1][i])/2 - v[j][i] * (phi[j - 1][i] + phi[j][i])/2 ) / dy;
                phi_new[j][i] = phi[j][i] + dt * (diffu - convec);
            }

           
            //inside square, set phi_new = 0

            for (j = 40; j < 60; j++) {
                for (i = 100; i < 120; i++) {
                    phi_new[j][i] = 0;
                }
            }


        }


       

        for (j = 0; j < Ny; j++) {
            for (i = 0; i < Nx; i++) {
                phi[j][i] = phi_new[j][i];
            }
        }





    }


    print_now = print_current_data(iter, u, v, phi, Nx, Ny, convective_method);
    printf("Data was printed for Point Jacobi method");


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


    free(du_s);
    free(dv_s);
    free(du_ss);
    free(dv_ss);

    free(step1_mat_x);
    free(step1_mat_y);

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

    printf("\n End of script, congrats!");
   
    return 0;

}