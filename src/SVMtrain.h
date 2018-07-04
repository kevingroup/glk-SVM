//
// Created by cyhong on 10/14/2017.
//

#ifndef GLK_SVMTRAIN_H
#define GLK_SVMTRAIN_H


class SVMtrain
{
public:
    SVMtrain(void);

    void train(double **kernel, int npos, int nneg, double *lambdas); // calculates labdas using iterative method

    double evaluateObjFunc(double **kernel, int npos, int nneg, double *lambdas); // returns the quadratic obj func value. This is not used by train. it is added to be used to optimize parameters used to define kernel
    ~SVMtrain(void);
    int niter20;
private:
    void initLambdas(double *x, int n);

};


#endif //GLK_SVMTRAIN_H
