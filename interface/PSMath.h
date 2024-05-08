#ifndef PSMATH_H_
#define PSMATH_H_

#include <Rtypes.h>

class PSMath {
 public:
  PSMath();
  static int PSfit(int iloop, int &iter, int &method, int &mode, bool &noNewtonShifts, int printlevel, int np,
    std::vector<double>& a,
    std::vector<double>& astart,
    std::vector<std::vector<double>>& alimit,
    const std::vector<double>& aprec,
    std::vector<double>& daN,
    std::vector<double>& h,
    std::vector<std::vector<double>>& aMemory,
	  double chi2,
    std::vector<double>& chi2iter,
    std::vector<double>& g,
    std::vector<double>& H,
    std::vector<double>& Hinv);

  static void PSNewtonLimitShift(int sign, int np,
    std::vector<double>& a,
    const std::vector<std::vector<double>>& alimit,
    const std::vector<double>& aprec,
		std::vector<double>& daN,
    std::vector<double>& h,
    std::vector<double>& g,
    const std::vector<double>& H);

  static double PSNewtonAnalyzer(int np,
    std::vector<double>& a,
    const std::vector<std::vector<double>>& alimit,
    const std::vector<double>& aprec,
		std::vector<double>& daN,
    std::vector<double>& h,
		const std::vector<double>& g,
    const std::vector<double>& H,
    std::vector<double>& Hinv,
    double chi2, bool noNewtonShifts, int printlevel=1);

  static void PSfitShow(int iloop, int convergence, int iter, int method, int mode, int printlevel, int graphiklevel,
		int np,
    const std::vector<double>& a,
    const std::vector<double>& astart,
    const std::vector<std::vector<double>>& alimit,
		const std::vector<double>& aprec,
    const std::vector<double>& daN,
    const std::vector<double>& h,
    double chi2,
		const std::vector<double>& g,
    const std::vector<double>& H);

  static double PSLineSearch(int & mode, double hh,
    const std::vector<double>& xlimit,
		double epsx, double epsf,
    std::vector<double>& x,
    std::vector<double>& f,
		double chi2, int printlevel);

  static void PSLineLimit(int np,
    const std::vector<double>& astart,
    const std::vector<double>& daN,
    const std::vector<std::vector<double>>& alimit,
    std::vector<double>& xlimit);

  static double  PSVnorm(const std::vector<double>& x, int n);

  static void PSVprint(const char* text, const std::vector<double>& x, int n);
  static void PSMprint(const char* text, const std::vector<double>& A, int ni, int nj);
  static void PSM2print(const char* text, const std::vector<std::vector<double>>& A, int ni);

  static double PSMinverse(const std::vector<double>& H, std::vector<double>& Hinv, int p);

  static double PSMCholtest();

  static double PSMmultiply(
    std::vector<double>& A,
    const std::vector<double>& B,
    const std::vector<double>& C,int n1,int n2);

  static double PSMmultiplyMRRT(std::vector<double>& A, int n1, int n2);

  static double PSMmultiplyMT(std::vector<double>& A, const std::vector<double>& B, int n1, int n2);

  static double PSMRTrianInvert2(std::vector<double>& R, int n);

  static double PSMRTrianInvert(const std::vector<double>& R, std::vector<double>& Rinv,  int n);

  static double PSMCholesky(const std::vector<double>& M, std::vector<double>& R,  int n);

  static double PSfitCheckLimits(int np, std::vector<double>& a, std::vector<double>& h,
			    const std::vector<std::vector<double>>& alimit, const std::vector<double>& aprecision,
			    std::vector<double>& daN, const std::vector<double>& g, double d);

  static double PSminIterate(std::vector<double>& a, std::vector<double>& daN, std::vector<double>& h, int p,
			const std::vector<double>& g, const std::vector<double>& H, const std::vector<double>& Hinv, double X0);

  static double PSfuncQuadratic(const std::vector<double>& a, const std::vector<double>& amean, double F0,
			   const std::vector<double>& g, const std::vector<double>& H, int np);



  static int PSderivative(int icall, int np, std::vector<double>& a, const std::vector<double>& h,
			  double chi2, std::vector<double>& chi2iter,
			  std::vector<double>& g, std::vector<double>& H, int printlevel);

  static int PSderivative1(int icall, std::vector<double>& a, const std::vector<double>& h,
			       double F, std::vector<double>& g, std::vector<double>& H);

  static double PSfitMinStep(int np, std::vector<double>& a, std::vector<double>& h,
			const std::vector<double>& chi2iter,
			const std::vector<double>& g, const std::vector<double>& H, std::vector<double>& Hinv, std::vector<double>& daN);

  static int PSfitconstrain0(double F, double g, double H, double Fix, std::vector<double>& aix);

};

#endif
