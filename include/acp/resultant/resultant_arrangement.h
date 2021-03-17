#include "acp/poly/root2.h"
#include "acp/resultant/resultant.h"

class ResultantArrangement {
 public:
  ResultantArrangement()
      : g_xl(Parameter::constant(-1)),
        g_xh(Parameter::constant(1)),
        g_yl(Parameter::constant(-1)),
        g_yh(Parameter::constant(1)) {}

  ResultantArrangement(Parameter xl, Parameter xh, Parameter yl, Parameter yh)
      : g_xl(xl), g_xh(xh), g_yl(yl), g_yh(yh) {}

  vector<PTR<Object<PV2>>> getIntersections(PTR<Object<Poly2>> f,
                                            PTR<Object<Poly2>> g);
  void alternates_double(vector<PV2<Parameter>>& boxes, PTR<Object<Poly2>> f,
                         PTR<Object<Poly2>> g,
                         vector<PV2<Parameter>>& verified);
  bool alternates(Parameter& yp, Parameter& x_f, Parameter& x_g,
                  PTR<Object<Poly2>> f, PTR<Object<Poly2>> g);

  Parameter g_xl;
  Parameter g_xh;
  Parameter g_yl;
  Parameter g_yh;
};
