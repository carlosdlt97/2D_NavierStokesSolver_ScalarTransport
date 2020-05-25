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

    

}