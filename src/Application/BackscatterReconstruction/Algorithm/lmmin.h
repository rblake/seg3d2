/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     lmmin.c
 *
 * Contents: Public interface to the Levenberg-Marquardt core implementation.
 *
 * Author:   Joachim Wuttke 2004-8 
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 *
 * Licence:  Public domain.
 */
 
#ifndef LMMIN_H
#define LMMIN_H

#ifdef __cplusplus
extern "C" {
#endif


/** Default data type for passing y(t) data to lm_evaluate **/

typedef struct {
  double *xvec;
  double *fvec;
  double (*f) (double x, double *par);
} lm_data_1d;

typedef struct {
  double *xvec;
  double *yvec;
  double *fvec;
  double (*f) (double x, double y, double *par);
} lm_data_2d;


/** User-supplied subroutines. **/
double lm_gaussian_fit_1d(double x, double *p);
void lm_evaluation_1d(double *par, int m_dat, double *fvec, void *data, int *info);
void lm_print_1d(int n_par, double *par, int m_dat, double *fvec, void *data, int iflag, int iter, int nfev);

double lm_gaussian_fit_2d(double x, double y, double *p);
void lm_evaluation_2d(double *par, int m_dat, double *fvec, void *data, int *info);
void lm_print_2d(int n_par, double *par, int m_dat, double *fvec, void *data, int iflag, int iter, int nfev);

double lm_gaussian_fit_2dg(double x, double y, double *p);


/** Compact high-level interface. **/

/* Collection of control parameters. */
typedef struct {
  double ftol;      /* relative error desired in the sum of squares. */
  double xtol;      /* relative error between last two approximations. */
  double gtol;      /* orthogonality desired between fvec and its derivs. */
  double epsilon;   /* step used to calculate the jacobian. */
  double stepbound; /* initial bound to steps in the outer loop. */
  double fnorm;     /* norm of the residue vector fvec. */
  int maxcall;      /* maximum number of iterations. */
  int nfev;	      /* actual number of iterations. */
  int info;	      /* status of minimization. */
} lm_control_type;

/* Initialize control parameters with default values. */
void lm_initialize_control(lm_control_type * control);

/* Refined calculation of Eucledian norm, typically used in printout routine. */
double lm_enorm(int, double *);

/* The actual minimization. */
void lm_minimize(int m_dat, int n_par, double *par,
		 void (*evaluate) (double *par, int m_dat, double *fvec,
                                   void *data, int *info),
                 void (*printout) (int n_par, double *par, int m_dat,
                                   double *fvec, void *data, int iflag,
                                   int iter, int nfev),
		 void *data, lm_control_type * control);


/** Legacy low-level interface. **/

/* Alternative to lm_minimize, allowing full control, and read-out
   of auxiliary arrays. For usage, see implementation of lm_minimize. */
void lm_lmdif(int m, int n, double *x, double *fvec, double ftol,
	      double xtol, double gtol, int maxfev, double epsfcn,
	      double *diag, int mode, double factor, int *info, int *nfev,
	      double *fjac, int *ipvt, double *qtf, double *wa1,
	      double *wa2, double *wa3, double *wa4,
              void (*evaluate) (double *par, int m_dat, double *fvec,
                                void *data, int *info),
              void (*printout) (int n_par, double *par, int m_dat,
                                double *fvec, void *data, int iflag,
                                int iter, int nfev),
	      void *data);

extern const char *lm_infmsg[];
extern const char *lm_shortmsg[];

#ifdef __cplusplus
}
#endif

#endif /* LMMIN_H */
