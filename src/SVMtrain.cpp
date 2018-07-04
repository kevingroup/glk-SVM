//
// Created by cyhong on 10/14/2017.
//

#include "SVMtrain.h"
#include "global.h"

SVMtrain::SVMtrain(void)
{
    niter20 = 5;
}


SVMtrain::~SVMtrain(void)
{
}

void SVMtrain::initLambdas(double *x, int n)
{
    for (int i=0;i<n; i++)
    {
        x[i] = myrandom(1000000)/1000000.1;
    }
}


// using iterative method of : A discriminative framework for detecting remote protein homologies. Jaakkola T, Diekhans M, Haussler D., 2000
void SVMtrain::train(double **kernel, int npos, int nneg, double *lambdas) // calculates lambdas using iterative method
{
    int niter20 = this->niter20; // *20 iterations
    int N = npos+nneg;
    initLambdas(lambdas, N);

    int Nr = 20*N;
    int *ri=new int[Nr];
    int i;
    for(i=0;i<Nr;i++)
    {
        ri[i] = i % N;
    }
    randomPermute(ri, Nr);

    for(int iter=0; iter<niter20; iter++)
    {
        for(int ii=0;ii<Nr;ii++)
        {
            int j= ri[ii];

            double lx = 0;
            double *lambdasi = lambdas;
            double *kernelji = kernel[j];

            for(i=0;i<npos;i++)
            {
                lx+=(*kernelji++)*(*lambdasi++);
            }
            for(i=0;i<nneg;i++)
            {
                lx-=(*kernelji++)*(*lambdasi++);
            }
            double z;
            if (j<npos)
            {
                z = 1-lx+lambdas[j];
            }
            else
            {
                z = 1+lx+lambdas[j];
            }
            if (z<0)
            {
                lambdas[j] = 0;
            }
            else
            if (z>1)
            {
                lambdas[j] = 1;
            }
            else
            {
                lambdas[j] = z;
            }
        }
    }

}

double SVMtrain::evaluateObjFunc(double **kernel, int npos, int nneg, double *lambdas) // calculates objective func value for given lambdas
{
    int N = npos+nneg;
    int i;
    double objf=0;
    for(int j=0;j<N;j++)
    {
        double lx = 0;
        double *lambdasi = lambdas;
        double *kernelji = kernel[j];

        for(i=0;i<npos;i++)
        {
            lx+=(*kernelji++)*(*lambdasi++);
        }

        for(i=0;i<nneg;i++)
        {
            lx-=(*kernelji++)*(*lambdasi++);
        }
        double z;
        if (j<npos)
        {
            z = 2-lx;
        }
        else
        {
            z = 2+lx;
        }
        objf += z*lambdas[j];
    }
    return(objf);
}