#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main() {

    int i, j;
    int Ni = 10;
    int Nj = 3;
    double ni = Ni;
    double nj = Nj;

    double** num = (double**)calloc(nj, sizeof(double*));

    for (j = 0; j < nj; j++) {
        num[j] = (double*)calloc(ni, sizeof(double));
    }

    for (j = 1; j < Nj - 1; j++) {
        printf("check\n");
    }

    return 0;

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