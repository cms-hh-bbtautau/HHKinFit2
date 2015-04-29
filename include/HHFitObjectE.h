/*
 * class for fit objects
 */

#ifndef HHFitObjectE_
#define HHFitObjectE_

#include "HHLorentzVector.h"
#include "TMatrixD.h"
#include "HHFitObject.h"

namespace HHKinFit2{
class HHFitObjectE : public HHFitObject {
 public:
  HHFitObjectE(HHLorentzVector const& initial4vector);
  
  double getE() const;
  virtual HHLorentzVector changeE(double E) const = 0;
  void changeEandSave(double E);
  virtual HHLorentzVector scaleE(double scale) const;
  void scaleEandSave(double scale);
  virtual HHLorentzVector constrainEtoMinv(double minv, HHLorentzVector const& other4vector) const =0;
  void constrainEtoMinvandSave(double minv, HHLorentzVector const& other4vector);

  double getUpperFitLimitE() const;
  double getLowerFitLimitE() const;

  void setFitLimitsE(double const lowerlimit, double const upperlimit);
  void setFitLimitsE(HHLorentzVector const& own4vectorMin, double const minv, HHLorentzVector const& other4vectorMin);
  void setUpperFitLimitE(double const upperlimit);
  void setUpperFitLimitE(double const minv, HHLorentzVector const& other4vectorMin);
  void setLowerFitLimitE(double const lowerlimit);
  void setLowerFitLimitE(HHLorentzVector const& other4vectorMin);

  void setCovMatrix(double dE);

  virtual void print() const;
  void printLimits() const;

 private:
  double m_upperLimitE;
  double m_lowerLimitE;
};
}
#endif /* HHFitObjectE_ */
