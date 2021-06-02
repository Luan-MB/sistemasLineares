#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int main ()
{
    
    int i = 1;

    // Le sistemas lineares ate encontrar EOF
    while (!feof(stdin)) {
        SistLinear_t *sistema;
        
        sistema = lerSistLinear();
        if (!sistema)
            return -1;

        real_t *result = malloc(sistema->n * sizeof(real_t));
        real_t *res;
        double tempo;
        real_t norma;
        int iter;

        printf("***** Sistema %u --> n = %u, erro: %g\n",i,sistema->n,sistema->erro);
        if ((eliminacaoGauss(sistema,result,&tempo)) == 0) {
            printf("Eliminacao Gauss: %f ms\n",tempo);
            printf("X: ");
            prnVetor(result,sistema->n);
            norma = normaL2Residuo(sistema,result,res);
            printf("Norma L2 do residuo: %g\n",norma);
            printf("\n");
        }
        
        iter = gaussJacobi(sistema,result,&tempo);

        printf("Jacobi: %f ms --- Iteracoes -> %d\n",tempo,iter);
        printf("X: ");
        prnVetor(result,sistema->n);
        norma = normaL2Residuo(sistema,result,res);
        printf("Norma L2 do residuo: %g\n",norma);

        if (norma > 5.0f) {
            iter = refinamento(sistema,result,&tempo);
            printf("\n");
            printf("Refinamento: %f ms --- Iteracoes -> %d\n",tempo,iter);
            printf("X: ");
            prnVetor(result,sistema->n);
            norma = normaL2Residuo(sistema,result,res);
            printf("Norma L2 do residuo: %g\n",norma);
        }
        printf("\n");

        iter = gaussSeidel(sistema,result,&tempo);
        printf("Gauss-Seidel: %f ms --- Iteracoes -> %u\n",tempo, iter);
        printf("X: ");
        prnVetor(result,sistema->n);
        norma = normaL2Residuo(sistema,result,res);
        printf("Norma L2 do residuo: %g\n",norma);
        printf("\n");
        
        liberaSistLinear(sistema);
        free(result);
        scanf("\n");
        
        i++;
    }
    return 1;
}

