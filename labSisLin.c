#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "SistemasLineares.h"

int main ()
{
    
    int i = 1;

    // Le sistemas lineares ate encontrar EOF
    do {
        SistLinear_t *sistema;
        
        sistema = lerSistLinear();
        if (!sistema) {
            perror("Falha ao alocar na memoria");
            return -1;
        }

        real_t *result = malloc(sistema->n * sizeof(real_t));
        real_t *res = malloc(sistema->n * sizeof(real_t));
        double tempo;
        real_t norma;
        unsigned int iter;

        printf("***** Sistema %u --> n = %u, erro: %g\n",i,sistema->n,sistema->erro);
        eliminacaoGauss(sistema,result,&tempo);
        printf("Eliminacao Gauss: %f ms\n",tempo);
        printf("X: ");
        prnVetor(result,sistema->n);
        norma = normaL2Residuo(sistema,result,res);
        printf("Norma L2 do residuo: %g\n",norma);
        printf("\n");
        
        iter = gaussJacobi(sistema,result,&tempo);
        printf("Jacobi: %f ms --- Iteracoes -> %u\n",tempo,iter);
        printf("X: ");
        prnVetor(result,sistema->n);
        norma = normaL2Residuo(sistema,result,res);
        printf("Norma L2 do residuo: %g\n",norma);
        printf("\n");
        
        iter = gaussSeidel(sistema,result,&tempo);
        printf("Seidel: %f ms --- Iteracoes -> %u\n",tempo, iter);
        printf("X: ");
        prnVetor(result,sistema->n);
        norma = normaL2Residuo(sistema,result,res);
        printf("Norma L2 do residuo: %g\n",norma);
        
        printf("\n");
        liberaSistLinear(sistema);
        free(result);
        free(res);
        scanf("\n");
        
        i++;
    } while (!feof(stdin));
    return 1;
}

