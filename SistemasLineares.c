// Luan Machado Bernardt | GRR20190363

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "utils.h"
#include "SistemasLineares.h"

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
  \param res Valor do resíduo

  \return Norma L2 do resíduo.
*/

real_t normaL2Residuo(SistLinear_t *SL, real_t *x, real_t *res)
{
    real_t sum = 0.0f;
    
    for (int i=0; i<SL->n; i++) {
        res[i] = 0.0f;
        for (int j=0; j<SL->n; j++)
            res[i] += SL->A[i][j] * x[j];
        res[i] = SL->b[i] - res[i];
        sum += powf(res[i],2.0f);
    }
    return (sqrtf(sum));
}

// Copia o sistema linear base
static SistLinear_t *copySL(SistLinear_t *base) {

    SistLinear_t *copia = alocaSistLinear(base->n);
    if (!copia)
        return NULL;

    copia->n = base->n;
    copia->erro = base->erro;

    for (int i=0; i<base->n; i++) {
        copia->b[i] = base->b[i];
        for (int j=0; j<base->n; j++)
            copia->A[i][j] = base->A[i][j];
    }
    return copia;
}

static unsigned int maxValue (SistLinear_t *SL, unsigned int i) {

    unsigned int max = i;
    
    for (int j=i+1; j<SL->n; j++) {
        if (fabs(SL->A[j][i]) > fabs(SL->A[max][i]))
          max = j;
    }
    
    return max;
}

static void trocaLinha (SistLinear_t* SL, unsigned int i, unsigned int j) {

    real_t bAux;
    real_t *aAux;

    aAux = SL->A[i];
    bAux = SL->b[i];
    SL->A[i] = SL->A[j];
    SL->b[i] = SL->b[j];
    SL->A[j] = aAux;
    SL->b[j] = bAux;
}

/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTotal tempo gasto pelo método

  \return código de erro. 0 em caso de sucesso.
*/

int eliminacaoGauss (SistLinear_t *SL, real_t *x, double *tTotal)
{
    SistLinear_t *copia = copySL(SL);
    if (!copia) {
        perror("Erro Eliminacao Gauss: falha ao copiar sistema linear");
        return -1;
    }
    
    unsigned int pivo;
    *tTotal = timestamp();
    
    // Transforma a matriz em uma triangular com pivoteamento parcial
    for (int i=0; i<SL->n; i++) {
        pivo = maxValue(copia,i);
        if (pivo != i)
            trocaLinha(copia,i,pivo);
        
        for (int j=i+1; j<SL->n; j++) {
            double m = copia->A[j][i] / copia->A[i][i];
            copia->A[j][i] = 0.0f;
            for (int k=i+1; k<SL->n; k++)
                copia->A[j][k] -= copia->A[i][k] * m;
            copia->b[j] -= copia->b[i] * m;
        }
    }
    
    // Retrossubstituicao
    for (int i=SL->n-1; i>=0; i--) {
        x[i] = copia->b[i];
        for (int j=i+1; j<SL->n; j++)
            x[i] -= copia->A[i][j] * x[j];
        x[i] /= copia->A[i][i];
    }
    *tTotal = timestamp() - *tTotal;
    liberaSistLinear(copia);
    return 0;
}

// Funcao que recebe dois vetores de mesmo tamanho n
// Retorna a maior diferenca entre elementos com mesmo indice
real_t maxDiff (real_t *a, real_t *b, unsigned int n) {

    double diff, max = 0.0f;

    for (int i=0; i<n; i++) {
        diff = fabs(a[i] - b[i]);
        if (diff > max)
            max = diff;
    }

    return max;
}

// Funcao que recebe um sistema linear
// Retorna 1 se for diagonalmente dominante e 0 se nao

int diagonalDominante(SistLinear_t *SL) {

    real_t sum;
    
    for (int i=0; i<SL->n; i++) {
        sum = 0.0f;
        for (int j=0; j<SL->n; j++) {
            if (i != j)
                sum += SL->A[i][j];
        }
        if (sum > SL->A[i][i])
            return 0;
    }
    return 1;
}

/*!
  \brief Método de Jacobi

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
*/
int gaussJacobi (SistLinear_t *SL, real_t *x, double *tTotal)
{
        
    double sum, aux, diff;
    real_t oldX[SL->n];

    // Preenche o vetor oldX de zeros
    for (int k=0; k<SL->n; k++)
        oldX[k] = 0.0f;

    int iter = 0;
    *tTotal = timestamp();
    
    do {
        for (int k=0; k<SL->n; k++) {
            aux = SL->A[k][k];
            sum = SL->b[k];
            
            for (int j=0; j<SL->n; j++)
                if (k != j)
                    sum -= SL->A[k][j] * oldX[j];
            x[k] = sum / aux;
        }
        
        diff = maxDiff(oldX,x,SL->n);

        // Passa x para oldX
        for (int k=0; k<SL->n; k++)
            oldX[k] = x[k];
        
        iter++;
    } while ((iter < MAXIT) && (diff > SL->erro));
    
    *tTotal = timestamp() - *tTotal;
    
    if ((iter == 50) && (!diagonalDominante(SL)))
        return -1; // Nao converge

    for (int i=0; i<SL->n; i++)
        if ((isinf(x[i])) || (isnan(x[i])))
            return -2; // Sistema impossivel
    return iter;
}

/*!
  \brief Método de Gauss-Seidel

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int gaussSeidel (SistLinear_t *SL, real_t *x, double *tTotal)
{

    double sum, aux, diff;
    real_t oldX[SL->n];


    // Preenche o vetor oldX de zeros
    for (int k=0; k<SL->n; k++)
        oldX[k] = 0.0f;

    int iter = 0;
    *tTotal = timestamp();

    do {
        for (int k=0; k<SL->n; k++) {
            aux = SL->A[k][k];
            sum = SL->b[k];
            
            for (int j=0; j<SL->n; j++) {
                if (k != j)
                    if (j < k)
                        sum -= SL->A[k][j] * x[j];
                    else
                        sum -= SL->A[k][j] * oldX[j];
            }

            x[k] = sum / aux;
        }
        
        diff = maxDiff(oldX,x,SL->n);

        // Passa x para oldX
        for (int k=0; k<SL->n; k++)
            oldX[k] = x[k];
        
        iter++;
    } while ((iter < MAXIT) && (diff > SL->erro));
    
    *tTotal = timestamp() - *tTotal;
    
    if ((iter == 50) && (!diagonalDominante(SL)))
        return -1; // Nao converge

    for (int i=0; i<SL->n; i++)
        if ((isinf(x[i])) || (isnan(x[i])))
            return -2; // Sistema impossivel
    return iter;
}


/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução. Ao iniciar função contém
            valor inicial para início do refinamento
  \param tTotal tempo gasto pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro:
          -1 (não converge) -2 (sem solução)
  */
int refinamento (SistLinear_t *SL, real_t *x, double *tTotal)
{
    SistLinear_t *copiaSL = copySL(SL);
    if (!copiaSL) {
        perror("Erro Refinamento: falha ao copiar sistema linear");
        return 0;
    }
    real_t res[SL->n], w[SL->n], xLinha[SL->n];
    real_t norma;
    double tempo;
    int iter=0;
    
    *tTotal = timestamp();

    while ((iter < MAXIT) && (SL->erro < (norma = normaL2Residuo(SL,x,res)))) {
        for (int i=0; i<SL->n; i++) {
            xLinha[i] = x[i];
            copiaSL->b[i] = res[i];
        }
       
        eliminacaoGauss(copiaSL,w,&tempo);
        
        for (int i=0; i<SL->n; i++)
            x[i] += w[i];

        if (SL-> erro > maxDiff(x,xLinha,SL->n))
            break;

        iter++;
    }
    
    *tTotal = timestamp() - *tTotal;
    liberaSistLinear(copiaSL);
    
    if ((iter == 50) && (!diagonalDominante(SL)))
        return -1; // Nao converge

    for (int i=0; i<SL->n; i++)
        if ((isinf(x[i])) || (isnan(x[i])))
            return -2; // Sistema impossivel
    return iter;
}

/*!
  \brief Alocaçao de memória 

  \param n tamanho do SL

  \return ponteiro para SL. NULL se houve erro de alocação
  */
SistLinear_t* alocaSistLinear (unsigned int n)
{
  
    SistLinear_t *newSL = malloc(sizeof(SistLinear_t));
    
    if (!newSL)
        return NULL;
    
    // Aloca um vetor de ponteiros para real_t com tamanho n
    newSL->A = malloc(n * sizeof(real_t *));
    if (!newSL->A)
      return NULL;
    // Aloca n vetores de real_t de tamanho n
    for (int i=0; i<n; i++) {
      newSL->A[i] = malloc(n * sizeof(real_t));
      if (!newSL->A[i])
          return NULL;
    }
    
    // Aloca um vetor de real_t com tamanho n
    newSL->b = malloc(n * sizeof(real_t));
    if (!newSL->b)
        return NULL;
    
    return newSL;
}

/*!
  \brief Liberaçao de memória 

  \param sistema linear SL
  */
void liberaSistLinear (SistLinear_t *SL)
{
    
    for (int i=0; i<SL->n; i++)
        free(SL->A[i]);
    free(SL->A);
    free(SL->b);
    free(SL);
}

/*!
  \brief Leitura de SL a partir de Entrada padrão (stdin).

  \return sistema linear SL. NULL se houve erro (leitura ou alocação)
  */
SistLinear_t *lerSistLinear ()
{
    
    unsigned int n;
    real_t erro, value;

    // Recebe n e aloca o sistema linear
    scanf("%u",&n);
    if (!n) {
        perror("Ordem da matriz invalida");
        return NULL;
    }

    SistLinear_t *sisLin = alocaSistLinear(n);
    if (!sisLin)
        return NULL;
        
    sisLin->n = n;

    // Recebe o erro
    scanf("%f",&erro);
    if (!erro) {
        perror("Erro invalido");
        return NULL;
    }
    
    sisLin->erro = erro;

    // Le os termos do SL
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            scanf("%f",&value);
            sisLin->A[i][j] = value; 
        }
    }
    
    // Le os termos independentes
    for (int k=0; k<n; k++) {
        scanf("%f",&value);
        sisLin->b[k] = value;
    }

    return sisLin;
}

// Exibe SL na saída padrão
void prnSistLinear (SistLinear_t *SL)
{

    unsigned int n = SL->n;
    
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++)
            printf("%f ",SL->A[i][j]);
        printf("\n");
    }
}

// Exibe um vetor na saída padrão
void prnVetor (real_t *v, unsigned int n)
{
    for (int i=0; i<n; i++)
      printf("%g ",v[i]);
    printf("\n");
}

