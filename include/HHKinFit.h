/*
 * class for general fit constraint
 */

#ifndef HHKinFit_
#define HHKinFit_

#include "TGraph.h"
#include <vector>
#include "HHFitObjectE.h"
#include "HHFitConstraint.h"

namespace HHKinFit2{
class HHKinFit {
 public:
  HHKinFit();

  void fit();
  double getChi2() const;
  double getL() const;
  std::vector<HHFitObjectE*> getListOfFitObjects() const;
  std::vector<HHFitConstraint*> getListOfConstraints() const;
  
  void addFitObjectE(HHFitObjectE* fitobject);
  void addConstraint(HHFitConstraint* constraint);
  
  TGraph* getChi2Function(int steps);
  TGraph* getLFunction(int steps);
  int getConvergence() const;


  void printChi2() const;


 private:
  std::vector<HHFitObjectE*> m_fitobjects;
  std::vector<HHFitConstraint*> m_constraints;

  double m_chi2;
  int m_convergence;
  int m_printlevel;
  int m_maxloops;
};
}
#endif /* HHKinFit_ */
