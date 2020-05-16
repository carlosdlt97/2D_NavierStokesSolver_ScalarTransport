#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main() {

    int i;

    double** num = (double**)calloc(12, sizeof(double*));

    for (i = 0; i < 12; i++) {
        num[i] = (double*)calloc(20, sizeof(double));
    }

    printf("%f\n", num[5][2]);

    return 0;
}