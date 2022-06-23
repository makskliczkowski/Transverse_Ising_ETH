#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>

/* number of data points to fit */
#define N 40
#define FIT(i) gsl_vector_get(w->x, i)
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

struct data 
{
    size_t n;
    double * y;
};

int expb_f (const gsl_vector * x, void *data, gsl_vector * f)
{
    size_t n = ((struct data *)data)->n;
    double *y = ((struct data *)data)->y;

    double A_1 = gsl_vector_get (x, 0);
    double A_2 = gsl_vector_get (x, 1);

    size_t i;

    for (i = 0; i < n; i++)
    {
        /* Model Yi = A_1 / x + A_2 / x**2 */
        double t = i;
        double Yi = A_1 / (t + 0.1) +A_2 / (t*t + 0.2*t + 0.01) ;
        gsl_vector_set (f, i, Yi - y[i]);
    }

    return GSL_SUCCESS;
}

int expb_df (const gsl_vector * x, void *data, gsl_matrix * J)
{
    size_t n = ((struct data *)data)->n;

    double A_1 = gsl_vector_get (x, 0);
    double A_2 = gsl_vector_get (x, 1);

    size_t i;

    for (i = 0; i < n; i++)
    {
        /* Jacobian matrix J(i,j) = dfi / dxj, */
        /* where fi = (Yi - yi)/sigma[i],      */
        /*       Yi = A_1 / (t + 0.1) +A_2 / (t*t + 0.2*t + 0.01) */
        /* and the xj are the parameters (A_1,A_2) */
        double t = i;
        double e = 1 / (t + 0.1);
        double e1 = 1 / (t*t + 0.2*t + 0.01);
        gsl_matrix_set (J, i, 0, e); 
        gsl_matrix_set (J, i, 1, e1);
    }
    return GSL_SUCCESS;
}

void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
    gsl_vector *f = gsl_multifit_nlinear_residual(w);
    gsl_vector *x = gsl_multifit_nlinear_position(w);
    double rcond;

    /* compute reciprocal condition number of J(x) */
    gsl_multifit_nlinear_rcond(&rcond, w);

    fprintf(stderr, "iter %2zu: A_1 = % e A_2 = % e cond(J) = % e, |f(x)| = % e \n", iter, gsl_vector_get(x, 0), gsl_vector_get(x, 1), 1.0 / rcond, gsl_blas_dnrm2(f));
}

int main (void)
{
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w;
    gsl_multifit_nlinear_fdf fdf;
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();
    const size_t n = N;
    const size_t p = 2;

    gsl_vector *f;
    gsl_matrix *J;
    gsl_matrix *covar = gsl_matrix_alloc (p, p);
    double y[N], weights[N];
    struct data d = { n, y };
    double x_init[2] = { 1.0, 1.0 }; /* starting values */
    gsl_vector_view x = gsl_vector_view_array (x_init, p);
    gsl_vector_view wts = gsl_vector_view_array(weights, n);
    gsl_rng * r;
    double chisq, chisq0;
    int status, info;
    size_t i;

    const double xtol = 1e-8;
    const double gtol = 1e-8;
    const double ftol = 0.0;

    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_default);

    /* define the function to be minimized */
    fdf.f = expb_f;
    fdf.df = expb_df;   /* set to NULL for finite-difference Jacobian */
    fdf.fvv = NULL;     /* not using geodesic acceleration */
    fdf.n = n;
    fdf.p = p;
    fdf.params = &d;

    /* this is the data to be fitted */
    for (i = 0; i < n; i++)
    {
        double t = i;
        double yi = (0.1 + 3.2/(t + 0.1))/(t + 0.1);
        double si = 0.1 * yi;
        double dy = gsl_ran_gaussian(r, si);

        weights[i] = 1.0 / (si * si);
        y[i] = yi + dy;
        printf ("% e % e \n",t + 0.1, y[i]);
    };

    /* allocate workspace with default parameters */
    w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);

    /* initialize solver with starting point and weights */
    gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w);

    /* compute initial cost function */
    f = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(f, f, &chisq0);

    /* solve the system with a maximum of 20 iterations */
    status = gsl_multifit_nlinear_driver(20, xtol, gtol, ftol, callback, NULL, &info, w);

    /* compute covariance of best fit parameters */
    J = gsl_multifit_nlinear_jac(w);
    gsl_multifit_nlinear_covar (J, 0.0, covar);

    /* compute final cost */
    gsl_blas_ddot(f, f, &chisq);

    fprintf(stderr, "summary from method '%s/%s'\n", gsl_multifit_nlinear_name(w), gsl_multifit_nlinear_trs_name(w));
    fprintf(stderr, "number of iterations: %zu \n", gsl_multifit_nlinear_niter(w));
    fprintf(stderr, "function evaluations: %zu \n", fdf.nevalf);
    fprintf(stderr, "Jacobian evaluations: %zu \n", fdf.nevaldf);
    fprintf(stderr, "reason for stopping: %s \n", (info == 1) ? "small step size" : "small gradient");
    fprintf(stderr, "initial |f(x)| = % e \n", sqrt(chisq0));
    fprintf(stderr, "final   |f(x)| = % e \n", sqrt(chisq));

    { 
        double dof = n - p;
        double c = GSL_MAX_DBL(1, sqrt(chisq / dof));

        fprintf(stderr, "chisq/dof = % e \n", chisq / dof);

        fprintf (stderr, "A_1      = % f +/- % f \n", FIT(0), c*ERR(0));
        fprintf (stderr, "A_2 = % f +/- % f \n", FIT(1), c*ERR(1));
    }

    fprintf (stderr, "status = %s \n", gsl_strerror (status));

    gsl_multifit_nlinear_free (w);
    gsl_matrix_free (covar);
    gsl_rng_free (r);

    return 0;
}