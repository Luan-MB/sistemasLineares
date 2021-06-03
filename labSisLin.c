// Luan Machado Bernardt | GRR20190363

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
        real_t *residuo = malloc(sistema->n * sizeof(real_t));
        double tempo;
        real_t norma;
        int iter;

        // Eliminação de Gauss-Jordan
        printf("***** Sistema %u --> n = %u, erro: %g\n",i,sistema->n,sistema->erro);
        if ((eliminacaoGauss(sistema,result,&tempo)) == 0) {
            printf("Eliminacao Gauss: %f ms\n",tempo);
            printf("X: ");
            prnVetor(result,sistema->n);
            norma = normaL2Residuo(sistema,result,residuo);
            printf("Norma L2 do residuo: %g\n",norma);
            printf("\n");
        }
        
        // Metodo de Jacobi
        iter = gaussJacobi(sistema,result,&tempo);
        if (iter == -1)
            printf("Jacobi: %f ms --- Nao converge\n",tempo);
        else if (iter == -2)
            printf("Jacobi: %f ms --- Sistema impossivel\n",tempo);
        else
            printf("Jacobi: %f ms --- Iteracoes -> %u\n",tempo, iter);
        printf("X: ");
        prnVetor(result,sistema->n);
        norma = normaL2Residuo(sistema,result,residuo);
        printf("Norma L2 do residuo: %g\n",norma);
        
        if (norma > 5.0f) { // Refinamento
            iter = refinamento(sistema,result,&tempo);
            printf("\n");
            if (iter == -1)
                printf("Refinamento: %f ms --- Nao converge\n",tempo);
            else if (iter == -2)
                printf("Refinamento: %f ms --- Sistema impossivel\n",tempo);
            else
                printf("Refinamento: %f ms --- Iteracoes -> %u\n",tempo, iter);
            printf("X: ");
            prnVetor(result,sistema->n);
            norma = normaL2Residuo(sistema,result,residuo);
            printf("Norma L2 do residuo: %g\n",norma);
        }
        printf("\n");

        // Metodo de Gauss-Seidel
        iter = gaussSeidel(sistema,result,&tempo);
        if (iter == -1)
            printf("Gauss-Seidel: %f ms --- Nao converge\n",tempo);
        else if (iter == -2)
            printf("Gauss-Seidel: %f ms --- Sistema impossivel\n",tempo);
        else
            printf("Gauss-Seidel: %f ms --- Iteracoes -> %u\n",tempo, iter);
        printf("X: ");
        prnVetor(result,sistema->n);
        norma = normaL2Residuo(sistema,result,residuo);
        printf("Norma L2 do residuo: %g\n",norma);
        
        if (norma > 5.0f) { // Refinamento
            iter = refinamento(sistema,result,&tempo);
            printf("\n");
            if (iter == -1)
                printf("Refinamento: %f ms --- Nao converge\n",tempo);
            else if (iter == -2)
                printf("refinamento: %f ms --- Sistema impossivel\n",tempo);
            else
                printf("Refinamento: %f ms --- Iteracoes -> %u\n",tempo, iter);
            printf("X: ");
            prnVetor(result,sistema->n);
            norma = normaL2Residuo(sistema,result,residuo);
            printf("Norma L2 do residuo: %g\n",norma);
        }
        printf("\n");
        
        liberaSistLinear(sistema);
        free(result);
        free(residuo);
        scanf("\n");
        
        i++;
    }
    return 1;
}

