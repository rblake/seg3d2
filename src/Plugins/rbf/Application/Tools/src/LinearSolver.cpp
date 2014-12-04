#include <vector>
#include "SparseMatrix.h"
#include "LinearSolver.h"
#include <cmath>
#include <cstdio>

using std::vector;

LinearSolver::LinearSolver()
{

}

LinearSolver::LinearSolver(SparseMatrix *myMat)
{
	A = myMat;
}

void LinearSolver::setMatrix(SparseMatrix *myMat)
{
	A = myMat;
}

vector<double> LinearSolver::biCGStab(vector<double> &b)
{
	vector<double> x;
	biCGStab(b,x);
	return x;
}

void LinearSolver::biCGStab(vector<double> &b, vector<double> &x)
{
    int i,j,k,l,m,n, iter=0;
    vector<double> r,rhat,v,p,t,s, diag, y, z, pret, pres;
    double rho,rhoold,alpha,omega,omegaold,beta;

    n = b.size();
    x.resize(n);
    r.resize(n);
    rhat.resize(n);
    v.resize(n);
    p.resize(n);
    t.resize(n);
    s.resize(n);
    diag.resize(n);
    y.resize(n);
    z.resize(n);
    pret.resize(n);
    pres.resize(n);

    for(i=0; i<n; i++)
    {
        x[i]=v[i]=p[i]=0;
        r[i]=rhat[i]=b[i];
        diag[i] = 1.0;///A.val[i][A.length[i]-1];
    }
    rho = rhoold = alpha = omega = omegaold=1;

    if(norm(r)<n*1e-10)
        return;
    //printf("norm(r)=%lf\n", norm(r)); fflush(stdout);
    while(iter<9000)
    {
        //std::cout << "iter: " << iter << " norm(r)= " << norm(r) << std::endl;
        //printf("Iteration %d: Residual norm = %lf\n", iter, norm(r)); fflush(stdout);
        iter++;
        rho = norm(rhat,r);
        //printf("rho: %lf\n", rho);
        if(iter==1)
        {
            for(i=0; i<n; i++)
            {
                p[i]=r[i];
                //printf("p[%d]=%lf\n", i, p[i]);
            }
        }
        else
        {
            beta= (rho/rhoold)*(alpha/omegaold);
            //printf("beta: %lf\n", beta);
            for(i=0; i<n; i++)
            {
                p[i]=r[i]+beta*(p[i] - omegaold*v[i]);
                //printf("p[%d]=%lf\n", i, p[i]);
            }
        }
        for(i=0; i<n; i++)
            y[i]=diag[i]*p[i];
        SpMV(y,v);
        alpha = rho/norm(rhat,v);
        //printf("alpha: %lf\n", alpha);
        for(i=0; i<n; i++)
        {
            s[i]=r[i]-alpha*v[i];
            //printf("s[%d]=%lf\n", i, s[i]);
        }
        for(i=0; i<n; i++)
            z[i]=diag[i]*s[i];
        SpMV(z,t);

        for(i=0; i<n; i++)
            pret[i]=diag[i]*t[i];
        for(i=0; i<n; i++)
            pres[i]=diag[i]*s[i];

        omega = norm(pret,pres)/norm(pret,pret);
        //printf("omega: %lf\n", omega);
        for(i=0; i<n; i++)
            x[i]+=alpha*p[i]+omega*s[i];
        for(i=0; i<n; i++)
            r[i]=s[i]-omega*t[i];
        if(norm(r)<n*1e-15)
            break;

        rhoold=rho; omegaold=omega;
    }
    //printf("Iteration: %d: Residual norm = %lf\n", iter, norm(r)); fflush(stdout);
    /*SpMV(x,v);
    for(i=0; i<n; i++)
    {
        printf("x[%d] = %.10f r[%d]=%.10lf\n", i, x[i],i,v[i]-b[i]);
    }*/
}

double LinearSolver::norm(vector<double> &a)
{
    int i, n=a.size();
    double ret=0;
    for (i=0; i<n; i++)
        ret += a[i]*a[i];
    return sqrt(ret);	
}

double LinearSolver::norm(vector<double> &a, vector<double> &b)
{
    int i, n=a.size();
    double ret=0;
    for (i=0; i<n; i++)
        ret += a[i]*b[i];
    return ret;
}

void LinearSolver::SpMV(vector<double> &a, vector<double> &b)
{
	A->multiply(a,b);
}
