#ifndef EEWW_ANALYSIS_GETELEMENT_V2_H
#define EEWW_ANALYSIS_GETELEMENT_V2_H

#include <limits>   //  ADDED
#include <cmath>    //  ADDED

float getElement_v2(const ROOT::VecOps::RVec<float>& vec,
                 size_t idx,
                 float default_val = std::numeric_limits<float>::quiet_NaN()) //  CHANGED
{
    if (vec.size() > idx) {
        float val = vec[idx];

        //  NEW: if sentinel value, return NaN instead
        if (val <= -998.f || !std::isfinite(val))
            return std::numeric_limits<float>::quiet_NaN();

        return val;
    }

    //  CHANGED: return NaN instead of -999
    return std::numeric_limits<float>::quiet_NaN();
}

#endif //EEWW_ANALYSIS_GETELEMENT_V2_H
