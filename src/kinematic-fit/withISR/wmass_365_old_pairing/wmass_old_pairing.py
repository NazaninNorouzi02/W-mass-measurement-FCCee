"""
FCC-ee hadronic WW mass reconstruction with optional ISR recovery.

Physics workflow
----------------
- Standard raw, 4C and 5C reconstruction follows Chapter 7.
- Jet pairing uses Eq. (7.40).
- Lost ISR is represented by one collinear photon as in Sec. 7.3.2:
  pT_gamma = 0 and E_gamma = |pz_gamma|.
- y -> pz_gamma follows Eq. (7.42); b follows Chapter 2, Eq. (2.9).
- A standard-fit probability below 3% triggers an ISR refit
  (Sec. 7.3.2.2).
- No 3% event rejection is applied in this producer.  The hybrid probability
  and explicit pass_final_*_p03 flags are saved for plotting-time rejection.
- The photon model is used for both 4C+ISR and 5C+ISR. For 4C+ISR
  the equal-mass constraint is disabled; for 5C+ISR it is retained.


Numerical workflow
------------------
The ISR fit is tried from y=0 and from a value inferred from the raw missing
longitudinal momentum. The best valid constrained solution is retained.

Key helpers
-----------
isrBetaParameter()       Eq. (2.9): ISR radiator parameter.
isrEmax()                ISR photon-energy limit.
pzGammaFromGaussianY()   Eq. (7.42): y -> pz_gamma.
AugLagFunctorISR         4C+ISR / 5C+ISR objective and constraints.
runISRFitCandidate()     One ISR fit from one starting point.
runBestISRFit()          Select the best valid ISR candidate.
reconstructWW()          Pairing and all raw/4C/5C reconstruction paths.

Legacy output indices 0..14 are preserved; 4C+ISR outputs are appended at
15..19 for downstream compatibility.
"""

import ROOT
import numpy as np
import json
import os

ROOT.gROOT.SetBatch(True)

# One producer is used for all three thesis energies.  The default is 365 GeV.
# Run with WMASS_ECM=162.5, WMASS_ECM=240 or WMASS_ECM=365.
_ENERGY_CONFIGS = {
    162.5: {
        "process": "p8_ee_WW_ecm160",
        "tag": "162p5",
        "json": "all_results_ecm162p5.json",
    },
    240.0: {
        "process": "p8_ee_WW_ecm240",
        "tag": "240",
        "json": "all_results_ecm240.json",
    },
    365.0: {
        "process": "p8_ee_WW_ecm365",
        "tag": "365",
        "json": "all_results_ecm365.json",
    },
}

_requested_ecm = float(os.environ.get("WMASS_ECM", "365"))
ECM = next(
    (value for value in _ENERGY_CONFIGS if abs(value - _requested_ecm) < 1.0e-6),
    None,
)
if ECM is None:
    raise ValueError(
        f"Unsupported WMASS_ECM={_requested_ecm}. "
        "Use 162.5, 240 or 365."
    )

_ENERGY_CONFIG = _ENERGY_CONFIGS[ECM]


# =========================================================
# 1) C++ inline code: raw, convergent 4C/5C/ISR fits and pulls
# =========================================================

ROOT.gInterpreter.Declare(r'''
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <iostream>
#include <limits>  // ISRFitResult validity/default diagnostics
#include <Math/Vector4D.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/ProbFuncMathCore.h>
#include "ROOT/RVec.hxx"

using namespace ROOT::VecOps;
using P4 = ROOT::Math::PxPyPzEVector;

namespace MyFit {
    static double V[16][16];
    static double Vinv[16][16];
    void setCovariance16(const std::vector<double>& vflat,
                         const std::vector<double>& invflat){
        for(int i=0;i<16;i++){
            for(int j=0;j<16;j++){
                V[i][j] = vflat[i*16+j];
                Vinv[i][j] = invflat[i*16+j];
            }
        }
    }
}

inline double wrapPhi(double phi){
    while(phi <= -M_PI) phi += 2.0*M_PI;
    while(phi >  M_PI)  phi -= 2.0*M_PI;
    return phi;
}

inline double clampTheta(double theta){
    return std::clamp(theta, 1e-5, M_PI-1e-5);
}

struct AugLagFunctor {
    double y0[16];   
    double m0[4];    
    float sqrts;
    int a, b, c, d;  
    bool do5C;
    mutable double lambda[5];
    mutable double mu;

    AugLagFunctor(bool do5C_=true) : do5C(do5C_) {
        for(int i=0; i<5; i++) lambda[i] = 0.0;
        mu = 10.0; 
    }

    // Computes the 4 (or 5) constraint values for a given parameter vector
    void constraints(const double *par, double f[5]) const {
        P4 j[4];
        for(int i=0; i<4; i++){
            double alpha = par[4*i+0];         
            double theta = clampTheta(par[4*i+1]); 
            double phi   = wrapPhi(par[4*i+2]);     
            double x     = par[4*i+3];         
            double ex = std::exp(x);
            double m = m0[i];
            if(m < 0.1) m = 0.1; 

            double px = alpha * m * ex * std::sin(theta) * std::cos(phi);
            double py = alpha * m * ex * std::sin(theta) * std::sin(phi);
            double pz = alpha * m * ex * std::cos(theta);
            double E  = alpha * m * std::sqrt(ex*ex + 1.0);
            j[i].SetPxPyPzE(px, py, pz, E);
        }

        f[0]=f[1]=f[2]=f[3]=f[4]=0.0;
        for(int i=0; i<4; i++){
            f[0] += j[i].Px(); f[1] += j[i].Py(); f[2] += j[i].Pz(); f[3] += j[i].E();
        }
        f[3] -= sqrts; 

        if(do5C){
            f[4] = (j[a]+j[b]).M() - (j[c]+j[d]).M(); 
        }
    }

    double operator()(const double *par) const {
        double chi2 = 0.0;
        double dpar[16];
        for(int i=0; i<16; i++) dpar[i] = par[i] - y0[i];
        
        for(int i=0; i<16; i++){
            double acc = 0.0;
            for(int j=0; j<16; j++){
                acc += MyFit::Vinv[i][j] * dpar[j];
            }
            chi2 += dpar[i] * acc;
        }

        double f[5];
        constraints(par, f);

        int nC = do5C ? 5 : 4;
        double L = chi2;
        double f2 = 0.0;
        for(int k=0; k<nC; k++){
            L += lambda[k] * f[k];
            f2 += f[k] * f[k];
        }
        L += 0.5 * mu * f2;
        return L;
    }
};

// =========================================================
// ISR UPDATE — thesis Secs. 7.3.2 and 7.3.2.1
//
// Lost ISR / beamstrahlung photons are approximated as one
// effective photon emitted along the beam axis:
//
//      pT_gamma = 0
//      E_gamma  = |pz_gamma|
//
// The photon is introduced through one Gaussian fit parameter y,
// with y ~ N(0,1), and mapped to pz_gamma using Eq. 7.42:
//
//      pz_gamma(y) = sign(y) * Emax *
//                    [ erf(|y|/sqrt(2)) ]^(1/b)
//
// b is not a free tuning parameter. It is the ISR radiator
// parameter from thesis Eq. 2.9:
//
//      b = 2 alpha / pi * ( log(s / me^2) - 1 )
//
// Emax is the maximum photon energy allowed by WW production
// kinematics, Eq. 7.41.
// =========================================================

// Thesis Chapter 2, Eq. (2.9): ISR radiator parameter b.
inline double isrBetaParameter(double sqrts){
    const double alpha_em = 1.0 / 137.035999084;
    const double me       = 0.000510998950; // electron mass [GeV]
    const double s        = sqrts * sqrts;

    return (2.0 * alpha_em / M_PI) * (std::log(s / (me * me)) - 1.0);
}

// Thesis Sec. 7.3.2, Eq. (7.41): maximum allowed ISR energy.
inline double isrEmax(double sqrts, double mW_ref){
    const double s = sqrts * sqrts;
    const double x = mW_ref * mW_ref / s;

    if (x >= 1.0) return 0.0;

    // UPDATED: use thesis Eq. 7.41 as written (1 - x, not 1 - 4x)
    return 0.5 * sqrts * (1.0 - x);
}

// Thesis Sec. 7.3.2.1, Eq. (7.42): Gaussian fit parameter -> photon pz.
inline double pzGammaFromGaussianY(double y, double Emax, double b){
    if (Emax <= 0.0 || b <= 0.0) return 0.0;

    const double sign = (y >= 0.0) ? 1.0 : -1.0;
    const double u = std::erf(std::abs(y) / std::sqrt(2.0));

    return sign * Emax * std::pow(u, 1.0 / b);
}

// ISR UPDATE — photon-augmented kinematic fit.
// Sec. 7.3.2 defines the photon object. The same object is used here with
// do5C=false for 4C+ISR and do5C=true for 5C+ISR.
// Parameters: par[0..15] = same jet parameters as AugLagFunctor,
//             par[16]    = y (Gaussian N(0,1)) controlling the photon pz.
//
// The same functor is used for:
//   do5C = false : 4C + ISR
//   do5C = true  : 5C + ISR
struct AugLagFunctorISR {
    double y0_raw[16];   // raw jet parameters — chi2 residual is measured from here
    double y0_start[17]; // minimizer starting point
    double m0[4];
    float sqrts;
    double Emax;
    double bISR;
    int a, b, c, d;
    bool do5C;
    mutable double lambda[5];
    mutable double mu;

    AugLagFunctorISR(bool do5C_=true) : do5C(do5C_) {
        for(int i=0; i<5; i++) lambda[i] = 0.0;
        mu = 10.0;
        Emax = 0.0;
        bISR = 0.0;
    }

    void constraints(const double *par, double f[5]) const {
        P4 j[4];
        for(int i=0; i<4; i++){
            double alpha = par[4*i+0];
            double theta = clampTheta(par[4*i+1]);
            double phi   = wrapPhi(par[4*i+2]);
            double x     = par[4*i+3];
            double ex = std::exp(x);
            double m = m0[i];
            if(m < 0.1) m = 0.1;

            double px = alpha * m * ex * std::sin(theta) * std::cos(phi);
            double py = alpha * m * ex * std::sin(theta) * std::sin(phi);
            double pz = alpha * m * ex * std::cos(theta);
            double E  = alpha * m * std::sqrt(ex*ex + 1.0);
            j[i].SetPxPyPzE(px, py, pz, E);
        }

        const double y_isr = par[16];
        const double pzg = pzGammaFromGaussianY(y_isr, Emax, bISR);
        const double Eg  = std::abs(pzg);

        f[0]=f[1]=f[2]=f[3]=f[4]=0.0;
        for(int i=0; i<4; i++){
            f[0] += j[i].Px();
            f[1] += j[i].Py();
            f[2] += j[i].Pz();
            f[3] += j[i].E();
        }

        // Lost ISR photon is assumed collinear with the beam:
        // pT_gamma = 0 and E_gamma = |pz_gamma|.
        f[2] += pzg;
        f[3] += Eg;
        f[3] -= sqrts;

        if(do5C){
            f[4] = (j[a]+j[b]).M() - (j[c]+j[d]).M();
        }
    }

    double operator()(const double *par) const {
        double chi2 = 0.0;
        double dpar[16];
        for(int i=0; i<16; i++) dpar[i] = par[i] - y0_raw[i];

        for(int i=0; i<16; i++){
            double acc = 0.0;
            for(int j=0; j<16; j++){
                acc += MyFit::Vinv[i][j] * dpar[j];
            }
            chi2 += dpar[i] * acc;
        }

        // y_isr is a measured Gaussian parameter with mean 0 and sigma 1.
        const double y_isr = par[16];
        chi2 += y_isr * y_isr;

        double f[5];
        constraints(par, f);

        const int nC = do5C ? 5 : 4;
        double L = chi2;
        double f2 = 0.0;
        for(int k=0; k<nC; k++){
            L += lambda[k] * f[k];
            f2 += f[k] * f[k];
        }
        L += 0.5 * mu * f2;
        return L;
    }
};

// =========================================================
// PULL UPDATE
//
// For a constrained least-squares fit, the pull denominator is
//
//   sqrt(V_meas - V_fit)_ii
// = sqrt([V A^T (A V A^T)^-1 A V]_ii),
//
// evaluated with the constraint Jacobian A at the actual fitted solution.
// The ISR fit uses the extended measured vector (16 jet parameters, y_ISR)
// with covariance diag(V, 1).  Therefore ISR treatment is applied to the
// pull calculation itself, not merely used as an event-selection cut.
// =========================================================

inline bool invertSmallMatrix(const double input[5][5],
                              int n,
                              double inverse[5][5]){
    double aug[5][10] = {{0.0}};

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++) aug[i][j] = input[i][j];
        for(int j=0; j<n; j++) aug[i][j+n] = (i == j) ? 1.0 : 0.0;
    }

    for(int col=0; col<n; col++){
        int pivot = col;
        double pivotAbs = std::abs(aug[col][col]);
        for(int row=col+1; row<n; row++){
            const double candidate = std::abs(aug[row][col]);
            if(candidate > pivotAbs){
                pivotAbs = candidate;
                pivot = row;
            }
        }

        if(!std::isfinite(pivotAbs) || pivotAbs < 1.0e-18) return false;

        if(pivot != col){
            for(int j=0; j<2*n; j++) std::swap(aug[col][j], aug[pivot][j]);
        }

        const double diagonal = aug[col][col];
        for(int j=0; j<2*n; j++) aug[col][j] /= diagonal;

        for(int row=0; row<n; row++){
            if(row == col) continue;
            const double factor = aug[row][col];
            for(int j=0; j<2*n; j++) aug[row][j] -= factor * aug[col][j];
        }
    }

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++) inverse[i][j] = aug[i][j+n];
    }
    return true;
}

inline void numericalJacobianStandard(const AugLagFunctor& fit,
                                      const double par[16],
                                      double A[5][17]){
    const int nC = fit.do5C ? 5 : 4;
    for(int k=0; k<5; k++) for(int i=0; i<17; i++) A[k][i] = 0.0;

    for(int ip=0; ip<16; ip++){
        double plus[16], minus[16];
        for(int j=0; j<16; j++) plus[j] = minus[j] = par[j];

        const double h = 1.0e-5 * std::max(1.0, std::abs(par[ip]));
        plus[ip] += h;
        minus[ip] -= h;
        if(ip % 4 == 1){ plus[ip] = clampTheta(plus[ip]); minus[ip] = clampTheta(minus[ip]); }
        if(ip % 4 == 2){ plus[ip] = wrapPhi(plus[ip]); minus[ip] = wrapPhi(minus[ip]); }

        double fPlus[5], fMinus[5];
        fit.constraints(plus, fPlus);
        fit.constraints(minus, fMinus);
        for(int k=0; k<nC; k++) A[k][ip] = (fPlus[k] - fMinus[k]) / (2.0*h);
    }
}

inline void numericalJacobianISR(const AugLagFunctorISR& fit,
                                 const double par[17],
                                 double A[5][17]){
    const int nC = fit.do5C ? 5 : 4;
    for(int k=0; k<5; k++) for(int i=0; i<17; i++) A[k][i] = 0.0;

    for(int ip=0; ip<17; ip++){
        double plus[17], minus[17];
        for(int j=0; j<17; j++) plus[j] = minus[j] = par[j];

        const double h = 1.0e-5 * std::max(1.0, std::abs(par[ip]));
        plus[ip] += h;
        minus[ip] -= h;
        if(ip < 16 && ip % 4 == 1){ plus[ip] = clampTheta(plus[ip]); minus[ip] = clampTheta(minus[ip]); }
        if(ip < 16 && ip % 4 == 2){ plus[ip] = wrapPhi(plus[ip]); minus[ip] = wrapPhi(minus[ip]); }

        double fPlus[5], fMinus[5];
        fit.constraints(plus, fPlus);
        fit.constraints(minus, fMinus);
        for(int k=0; k<nC; k++) A[k][ip] = (fPlus[k] - fMinus[k]) / (2.0*h);
    }
}

inline bool computePullsStandard(const AugLagFunctor& fit,
                                 const double fitted[16],
                                 float pulls[16]){
    for(int i=0; i<16; i++) pulls[i] = std::numeric_limits<float>::quiet_NaN();

    const int nC = fit.do5C ? 5 : 4;
    double A[5][17];
    numericalJacobianStandard(fit, fitted, A);

    double VAtranspose[17][5] = {{0.0}};
    for(int i=0; i<16; i++){
        for(int k=0; k<nC; k++){
            for(int j=0; j<16; j++) VAtranspose[i][k] += MyFit::V[i][j] * A[k][j];
        }
    }

    double AVAtranspose[5][5] = {{0.0}};
    for(int k=0; k<nC; k++){
        for(int l=0; l<nC; l++){
            for(int i=0; i<16; i++) AVAtranspose[k][l] += A[k][i] * VAtranspose[i][l];
        }
    }

    double inverse[5][5] = {{0.0}};
    if(!invertSmallMatrix(AVAtranspose, nC, inverse)) return false;

    for(int i=0; i<16; i++){
        double denominator2 = 0.0;
        for(int k=0; k<nC; k++){
            for(int l=0; l<nC; l++){
                denominator2 += VAtranspose[i][k] * inverse[k][l] * VAtranspose[i][l];
            }
        }
        if(std::isfinite(denominator2) && denominator2 > 1.0e-18){
            pulls[i] = (float)((fit.y0[i] - fitted[i]) / std::sqrt(denominator2));
        }
    }
    return true;
}

inline bool computePullsISR(const AugLagFunctorISR& fit,
                            const double fitted[17],
                            float pulls[17]){
    for(int i=0; i<17; i++) pulls[i] = std::numeric_limits<float>::quiet_NaN();

    const int nC = fit.do5C ? 5 : 4;
    double A[5][17];
    numericalJacobianISR(fit, fitted, A);

    double VAtranspose[17][5] = {{0.0}};
    for(int i=0; i<17; i++){
        for(int k=0; k<nC; k++){
            if(i < 16){
                for(int j=0; j<16; j++) VAtranspose[i][k] += MyFit::V[i][j] * A[k][j];
            } else {
                VAtranspose[i][k] = A[k][16]; // Var(y_ISR)=1
            }
        }
    }

    double AVAtranspose[5][5] = {{0.0}};
    for(int k=0; k<nC; k++){
        for(int l=0; l<nC; l++){
            for(int i=0; i<17; i++) AVAtranspose[k][l] += A[k][i] * VAtranspose[i][l];
        }
    }

    double inverse[5][5] = {{0.0}};
    if(!invertSmallMatrix(AVAtranspose, nC, inverse)) return false;

    for(int i=0; i<17; i++){
        double denominator2 = 0.0;
        for(int k=0; k<nC; k++){
            for(int l=0; l<nC; l++){
                denominator2 += VAtranspose[i][k] * inverse[k][l] * VAtranspose[i][l];
            }
        }
        const double measured = (i < 16) ? fit.y0_raw[i] : 0.0;
        if(std::isfinite(denominator2) && denominator2 > 1.0e-18){
            pulls[i] = (float)((measured - fitted[i]) / std::sqrt(denominator2));
        }
    }
    return true;
}

// ISR UPDATE — numerical inverse of thesis Eq. (7.42).
// Used only to seed the second minimization from the observed missing pz;
// it does not modify the ISR probability model.
inline double gaussianYFromPzGamma(double pzg, double Emax, double bISR){
    if(Emax <= 0.0 || bISR <= 0.0 || std::abs(pzg) <= 1.0e-12){
        return 0.0;
    }

    const double sign = (pzg >= 0.0) ? 1.0 : -1.0;
    const double maxPz = std::abs(pzGammaFromGaussianY(5.0, Emax, bISR));
    const double target = std::min(std::abs(pzg), maxPz * (1.0 - 1.0e-12));

    double lo = 0.0;
    double hi = 5.0;
    for(int iter=0; iter<80; iter++){
        const double mid = 0.5 * (lo + hi);
        const double pzMid = std::abs(pzGammaFromGaussianY(mid, Emax, bISR));
        if(pzMid < target) lo = mid;
        else               hi = mid;
    }

    return sign * 0.5 * (lo + hi);
}

// Container for one photon-augmented fit candidate and its validation data.
struct ISRFitResult {
    bool valid;
    double par[17];
    double chi2;
    double maxConstraint;
    double m1;
    double m2;
    double EPhoton;
    double yPhoton;
    int minimizerStatus;

    ISRFitResult()
        : valid(false),
          chi2(std::numeric_limits<double>::infinity()),
          maxConstraint(std::numeric_limits<double>::infinity()),
          m1(0.0),
          m2(0.0),
          EPhoton(0.0),
          yPhoton(0.0),
          minimizerStatus(-1) {
        for(int i=0; i<17; i++) par[i] = 0.0;
    }
};

// Run one 4C+ISR or 5C+ISR candidate from a specified photon start value.
inline ISRFitResult runISRFitCandidate(
        bool do5C,
        const double rawPar[16],
        const double jetStart[16],
        const double jetMass[4],
        float sqrts,
        double Emax,
        double bISR,
        int a, int b, int c, int d,
        double yStart,
        const char *tag){

    ISRFitResult result;
    AugLagFunctorISR fISR(do5C);
    fISR.sqrts = sqrts;
    fISR.a = a;
    fISR.b = b;
    fISR.c = c;
    fISR.d = d;
    fISR.Emax = Emax;
    fISR.bISR = bISR;

    for(int i=0; i<16; i++){
        fISR.y0_raw[i] = rawPar[i];
        fISR.y0_start[i] = jetStart[i];
    }
    fISR.y0_start[16] = std::clamp(yStart, -5.0, 5.0);

    for(int i=0; i<4; i++){
        fISR.m0[i] = jetMass[i];
    }

    ROOT::Math::Functor functorISR([&fISR](const double *par){ return fISR(par); }, 17);
    auto minISR = std::unique_ptr<ROOT::Math::Minimizer>(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad")
    );
    minISR->SetFunction(functorISR);
    minISR->SetStrategy(1);
    minISR->SetTolerance(1e-3);
    minISR->SetMaxFunctionCalls(10000);

    for(int i=0; i<4; i++){
        const int o = 4*i;
        minISR->SetLimitedVariable(
            o+0, Form("alpha_%s_%d", tag, i),
            fISR.y0_start[o+0], 0.02, 0.1, 2.0
        );
        minISR->SetVariable(
            o+1, Form("theta_%s_%d", tag, i),
            fISR.y0_start[o+1], 0.005
        );
        minISR->SetVariable(
            o+2, Form("phi_%s_%d", tag, i),
            fISR.y0_start[o+2], 0.005
        );
        minISR->SetLimitedVariable(
            o+3, Form("x_%s_%d", tag, i),
            fISR.y0_start[o+3], 0.01,
            fISR.y0_raw[o+3]-1.5,
            fISR.y0_raw[o+3]+1.5
        );
    }

    minISR->SetLimitedVariable(
        16, Form("y_isr_%s", tag),
        fISR.y0_start[16], 0.05, -5.0, 5.0
    );

    const int nC = do5C ? 5 : 4;
    // UPDATED: strict convergence requirements
    const int maxOuterSteps = 30;
    const double constraintTolerance = 1.0e-5;
    for(int step=0; step<maxOuterSteps; step++){
        minISR->Minimize();
        const double *parTmp = minISR->X();
        double fvals[5];
        fISR.constraints(parTmp, fvals);

        double maxConstraintStep = 0.0;
        for(int k=0; k<nC; k++){
            maxConstraintStep = std::max(maxConstraintStep, std::abs(fvals[k]));
        }
        if(!std::isfinite(maxConstraintStep)) break;
        if(step > 0 && minISR->Status() == 0 &&
           maxConstraintStep < constraintTolerance) break;

        for(int k=0; k<nC; k++){
            fISR.lambda[k] += fISR.mu * fvals[k];
        }
        fISR.mu = std::min(fISR.mu * 2.0, 1.0e8);
    }

    const double *parISR = minISR->X();
    for(int i=0; i<17; i++){
        result.par[i] = parISR[i];
    }

    double fvals[5];
    fISR.constraints(parISR, fvals);
    result.maxConstraint = 0.0;
    for(int k=0; k<nC; k++){
        result.maxConstraint = std::max(result.maxConstraint, std::abs(fvals[k]));
    }

    double dpar[16];
    for(int i=0; i<16; i++){
        dpar[i] = parISR[i] - fISR.y0_raw[i];
    }

    result.chi2 = 0.0;
    for(int i=0; i<16; i++){
        double acc = 0.0;
        for(int j=0; j<16; j++){
            acc += MyFit::Vinv[i][j] * dpar[j];
        }
        result.chi2 += dpar[i] * acc;
    }
    result.chi2 += parISR[16] * parISR[16];

    P4 jISR[4];
    for(int i=0; i<4; i++){
        const double alpha = parISR[4*i+0];
        const double theta = clampTheta(parISR[4*i+1]);
        const double phi   = wrapPhi(parISR[4*i+2]);
        const double x     = parISR[4*i+3];
        const double ex    = std::exp(x);
        double m = fISR.m0[i];
        if(m < 0.1) m = 0.1;

        jISR[i].SetPxPyPzE(
            alpha*m*ex*std::sin(theta)*std::cos(phi),
            alpha*m*ex*std::sin(theta)*std::sin(phi),
            alpha*m*ex*std::cos(theta),
            alpha*m*std::sqrt(ex*ex+1.0)
        );
    }

    result.m1 = (jISR[a] + jISR[b]).M();
    result.m2 = (jISR[c] + jISR[d]).M();
    result.yPhoton = parISR[16];
    result.EPhoton = std::abs(
        pzGammaFromGaussianY(parISR[16], fISR.Emax, fISR.bISR)
    );
    result.minimizerStatus = minISR->Status();

    // UPDATED: strict validity requirements
    result.valid =
        result.minimizerStatus == 0 &&
        std::isfinite(result.chi2) &&
        std::isfinite(result.m1) &&
        std::isfinite(result.m2) &&
        std::isfinite(result.EPhoton) &&
        result.maxConstraint < 1.0e-5;

    return result;
}

// Numerical safeguard: test y=0 and a missing-pz-derived start, then retain
// the valid solution with the smaller physical chi2.
inline ISRFitResult runBestISRFit(
        bool do5C,
        const double rawPar[16],
        const double jetStart[16],
        const double jetMass[4],
        float sqrts,
        double Emax,
        double bISR,
        int a, int b, int c, int d,
        double rawMissingPz,
        const char *tagZero,
        const char *tagMissingPz){

    // Starting only at y=0 can trap a gradient minimizer at the no-photon
    // solution because Eq. 7.42 is very flat around zero for the physical b.
    const ISRFitResult fromZero = runISRFitCandidate(
        do5C, rawPar, jetStart, jetMass,
        sqrts, Emax, bISR, a, b, c, d,
        0.0, tagZero
    );

    const double yFromMissingPz =
        gaussianYFromPzGamma(rawMissingPz, Emax, bISR);

    const ISRFitResult fromMissingPz = runISRFitCandidate(
        do5C, rawPar, jetStart, jetMass,
        sqrts, Emax, bISR, a, b, c, d,
        yFromMissingPz, tagMissingPz
    );

    if(fromZero.valid && fromMissingPz.valid){
        return (fromMissingPz.chi2 < fromZero.chi2)
            ? fromMissingPz
            : fromZero;
    }
    if(fromMissingPz.valid) return fromMissingPz;
    if(fromZero.valid)      return fromZero;

    // Keep the numerically closer failed candidate for diagnostics,
    // but the caller will not mark ISR as applied unless valid=true.
    if(fromMissingPz.maxConstraint < fromZero.maxConstraint){
        return fromMissingPz;
    }
    if(fromZero.maxConstraint < fromMissingPz.maxConstraint){
        return fromZero;
    }
    return (fromMissingPz.chi2 < fromZero.chi2)
        ? fromMissingPz
        : fromZero;
}

ROOT::VecOps::RVec<float>
reconstructWW(const ROOT::VecOps::RVec<float>& px,
              const ROOT::VecOps::RVec<float>& py,
              const ROOT::VecOps::RVec<float>& pz,
              const ROOT::VecOps::RVec<float>& E,
              float sqrts) {

    const int pairings[3][4] = { {0,1,2,3}, {0,2,1,3}, {0,3,1,2} };
    
    P4 j_raw[4];
    for(int i=0; i<4; i++) j_raw[i].SetPxPyPzE(px[i], py[i], pz[i], E[i]);

    // Jet pairing, Eq. (7.40) of the thesis:
    //   chi2_pairing = (m1 - mW)^2 + (m2 - mW)^2
    // with mW the *fixed* reference W mass (80.385 GeV/c^2), NOT the
    // smallest |m1-m2|. Minimizing |m1-m2| is only equivalent to this
    // near/above 240 GeV where both Ws are on-shell; at threshold
    // (162.5 GeV) one W is genuinely off-shell, so the two raw masses
    // are expected to differ, and minimizing |m1-m2| picks the wrong
    // permutation, dragging both reconstructed masses down (exactly
    // the symptom seen: m_small/m_large means ~65/~70 instead of ~80).
    const double mW_ref = 80.385;
    int best_ip = 0;
    double min_chi2 = 1e12;
    for(int ip=0; ip<3; ip++) {
        double m1 = (j_raw[pairings[ip][0]] + j_raw[pairings[ip][1]]).M();
        double m2 = (j_raw[pairings[ip][2]] + j_raw[pairings[ip][3]]).M();
        double chi2p = (m1-mW_ref)*(m1-mW_ref) + (m2-mW_ref)*(m2-mW_ref);
        if(chi2p < min_chi2) {
            min_chi2 = chi2p;
            best_ip = ip;
        }
    }

    int a = pairings[best_ip][0]; int b = pairings[best_ip][1];
    int c = pairings[best_ip][2]; int d = pairings[best_ip][3];

    float m1_raw = (j_raw[a] + j_raw[b]).M();
    float m2_raw = (j_raw[c] + j_raw[d]).M();


    // --- 4C Fit ---
    AugLagFunctor f4C(false);
    f4C.sqrts = sqrts; f4C.a = a; f4C.b = b; f4C.c = c; f4C.d = d;

    for(int i=0; i<4; i++){
        double p = j_raw[i].P(); double m = j_raw[i].M();
        double theta = (p > 0) ? std::acos(std::clamp(j_raw[i].Pz()/p, -1.0, 1.0)) : M_PI/2.0;
        double phi = std::atan2(j_raw[i].Py(), j_raw[i].Px());
        double x_init = (m > 0.1 && p > 0.1) ? std::log(std::clamp(p/m, 0.1, 100.0)) : 0.0;

        f4C.m0[i] = m;
        f4C.y0[4*i+0] = 1.0; f4C.y0[4*i+1] = clampTheta(theta); f4C.y0[4*i+2] = wrapPhi(phi); f4C.y0[4*i+3] = std::clamp(x_init, -4.0, 4.0);
    }

    // Implementation fix: capture by reference so multiplier updates in the
    // augmented-Lagrangian loop are seen by the minimizer.
    ROOT::Math::Functor functor4C([&f4C](const double *par){ return f4C(par); }, 16);
    auto min4C = std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    min4C->SetFunction(functor4C);
    min4C->SetStrategy(1); min4C->SetTolerance(1e-3); min4C->SetMaxFunctionCalls(10000);

    for(int i=0; i<4; i++){
        int o = 4*i;
        min4C->SetLimitedVariable(o+0, Form("alpha%d",i), 1.0, 0.02, 0.1, 2.0);
        min4C->SetVariable(o+1, Form("theta%d",i), f4C.y0[o+1], 0.005);
        min4C->SetVariable(o+2, Form("phi%d",i), f4C.y0[o+2], 0.005);
        min4C->SetLimitedVariable(o+3, Form("x%d",i), f4C.y0[o+3], 0.01, f4C.y0[o+3]-1.5, f4C.y0[o+3]+1.5);
    }

    // UPDATED: strict convergence logic with Status()==0 requirement
    for(int step=0; step<30; step++) {
        min4C->Minimize();
        const double* parTmp = min4C->X();
        double fvals[5];
        f4C.constraints(parTmp, fvals);

        double maxConstraintStep = 0.0;
        for(int k=0; k<4; k++){
            maxConstraintStep = std::max(maxConstraintStep, std::abs(fvals[k]));
        }
        if(!std::isfinite(maxConstraintStep)) break;
        if(step > 0 && min4C->Status() == 0 &&
           maxConstraintStep < 1.0e-5) break;

        for(int k=0; k<4; k++){
            f4C.lambda[k] += f4C.mu * fvals[k];
        }
        f4C.mu = std::min(f4C.mu * 2.0, 1.0e8);
    }

    const double* par4C = min4C->X();

    // DIAGNOSTIC UPDATE: record the standard 4C numerical state without
    // changing its acceptance or the fitted result.
    const int status_4C = min4C->Status();
    double diagnosticConstraints4C[5];
    f4C.constraints(par4C, diagnosticConstraints4C);
    double max_constraint_4C = 0.0;
    for(int k=0; k<4; k++){
        max_constraint_4C = std::max(
            max_constraint_4C,
            std::abs(diagnosticConstraints4C[k])
        );
    }

    P4 j_4C[4];
    for(int i=0; i<4; i++) {
        double alpha = par4C[4*i+0]; double theta = clampTheta(par4C[4*i+1]);
        double phi = wrapPhi(par4C[4*i+2]); double x = par4C[4*i+3];
        double ex = std::exp(x); double m = f4C.m0[i]; if(m<0.1) m=0.1;
        j_4C[i].SetPxPyPzE(alpha*m*ex*std::sin(theta)*std::cos(phi), alpha*m*ex*std::sin(theta)*std::sin(phi), alpha*m*ex*std::cos(theta), alpha*m*std::sqrt(ex*ex+1.0));
    }
    float m1_4C = (j_4C[a] + j_4C[b]).M();
    float m2_4C = (j_4C[c] + j_4C[d]).M();

    // Report the true chi2 (without penalty/multiplier terms), i.e.
    // the actual Least-Squares estimator from Eq. (7.5) of the thesis.
    double dpar4C[16];
    for(int i=0;i<16;i++) dpar4C[i] = par4C[i] - f4C.y0[i];
    double chi2_4C = 0.0;
    for(int i=0;i<16;i++){
        double acc=0.0;
        for(int j=0;j<16;j++) acc += MyFit::Vinv[i][j]*dpar4C[j];
        chi2_4C += dpar4C[i]*acc;
    }

    // UPDATED: strict validity requires converged minimizer
    const bool standard_4C_valid =
        status_4C == 0 &&
        std::isfinite(chi2_4C) && max_constraint_4C < 1.0e-5;
    const double prob_4C = standard_4C_valid
        ? ROOT::Math::chisquared_cdf_c(chi2_4C, 4.0)
        : 0.0;

    // PULL UPDATE: standard 4C pulls are evaluated at the converged Minuit
    // solution.  The final array is replaced below only when 4C+ISR is used.
    float pull4C_standard[16];
    float pull4C_final[16];
    for(int i=0; i<16; i++){
        pull4C_standard[i] = std::numeric_limits<float>::quiet_NaN();
        pull4C_final[i] = std::numeric_limits<float>::quiet_NaN();
    }
    const int pull4C_standard_valid =
        (standard_4C_valid && computePullsStandard(f4C, par4C, pull4C_standard)) ? 1 : 0;
    for(int i=0; i<16; i++) pull4C_final[i] = pull4C_standard[i];
    int pull4C_final_valid = pull4C_standard_valid;
    float pull4C_final_yISR = std::numeric_limits<float>::quiet_NaN();

    const double MW_ISR_REF = 80.385; // [GeV], fixed W-mass reference
    const double EMAX_ISR   = isrEmax(sqrts, MW_ISR_REF);
    const double B_ISR      = isrBetaParameter(sqrts);

    double rawMissingPz = 0.0;
    for(int i=0; i<4; i++){
        rawMissingPz -= j_raw[i].Pz();
    }

    // --- 4C + ISR treatment [ISR UPDATE] ---
    // Direct application of the Sec. 7.3.2 photon model with the fifth
    // equal-mass constraint disabled. No event is filtered here.
    float m1_4C_isr = m1_4C;
    float m2_4C_isr = m2_4C;
    float chi2_4C_isr = (float)chi2_4C;
    float E_isr_4C_fitted = 0.0f;
    int isr_4C_applied = 0;

    // DIAGNOSTIC UPDATE: keep attempt, validity, minimizer, constraint and
    // fitted-y information separate. Sentinel values mean "not attempted".
    const int isr_4C_attempted =
        (standard_4C_valid && prob_4C < 0.03 &&
         EMAX_ISR > 1.0e-9) ? 1 : 0;
    int isr_4C_fit_valid = 0;
    int status_4C_isr = -999;
    double max_constraint_4C_isr = -1.0;
    double y_isr_4C_fitted = 0.0;

    if(standard_4C_valid && prob_4C < 0.03 &&
       EMAX_ISR > 1.0e-9){
        const ISRFitResult best4CISR = runBestISRFit(
            false,
            f4C.y0,
            par4C,
            f4C.m0,
            sqrts,
            EMAX_ISR,
            B_ISR,
            a, b, c, d,
            rawMissingPz,
            "4C_zero",
            "4C_pz"
        );

        // DIAGNOSTIC UPDATE: save the returned candidate even when it fails
        // the existing validity requirement. This does not select it.
        isr_4C_fit_valid = best4CISR.valid ? 1 : 0;
        status_4C_isr = best4CISR.minimizerStatus;
        max_constraint_4C_isr = best4CISR.maxConstraint;
        y_isr_4C_fitted = best4CISR.yPhoton;

        if(best4CISR.valid){
            isr_4C_applied = 1;
            m1_4C_isr = (float)best4CISR.m1;
            m2_4C_isr = (float)best4CISR.m2;
            chi2_4C_isr = (float)best4CISR.chi2;
            E_isr_4C_fitted = (float)best4CISR.EPhoton;

            // PULL UPDATE: the selected 4C+ISR solution has 17 measured
            // parameters with covariance diag(V,1).
            AugLagFunctorISR pullFit4CISR(false);
            pullFit4CISR.sqrts = sqrts;
            pullFit4CISR.Emax = EMAX_ISR;
            pullFit4CISR.bISR = B_ISR;
            pullFit4CISR.a = a; pullFit4CISR.b = b;
            pullFit4CISR.c = c; pullFit4CISR.d = d;
            for(int i=0; i<16; i++) pullFit4CISR.y0_raw[i] = f4C.y0[i];
            for(int i=0; i<4; i++) pullFit4CISR.m0[i] = f4C.m0[i];

            float pullISR[17];
            pull4C_final_valid = computePullsISR(
                pullFit4CISR, best4CISR.par, pullISR
            ) ? 1 : 0;
            for(int i=0; i<16; i++) pull4C_final[i] = pullISR[i];
            pull4C_final_yISR = pullISR[16];
        }
    }

    // DIAGNOSTIC UPDATE: recovery means that a valid ISR solution moved the
    // event above the same final 3% probability threshold.
    const int isr_4C_recovered =
        (isr_4C_applied > 0 &&
         ROOT::Math::chisquared_cdf_c(chi2_4C_isr, 4.0) > 0.03) ? 1 : 0;

    // --- 5C Fit ---
    AugLagFunctor f5C(true);
    f5C.sqrts = sqrts; f5C.a = a; f5C.b = b; f5C.c = c; f5C.d = d;
    std::copy(std::begin(f4C.y0), std::end(f4C.y0), std::begin(f5C.y0));
    std::copy(std::begin(f4C.m0), std::end(f4C.m0), std::begin(f5C.m0));

    // Same reference capture as the 4C fit; this is an implementation fix,
    // not an additional physics assumption.
    ROOT::Math::Functor functor5C([&f5C](const double *par){ return f5C(par); }, 16);
    auto min5C = std::unique_ptr<ROOT::Math::Minimizer>(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    min5C->SetFunction(functor5C);
    min5C->SetStrategy(1); min5C->SetTolerance(1e-3); min5C->SetMaxFunctionCalls(10000);

    for(int i=0; i<4; i++){
        int o = 4*i;
        min5C->SetLimitedVariable(o+0, Form("alpha%d",i), 1.0, 0.02, 0.1, 2.0);
        min5C->SetVariable(o+1, Form("theta%d",i), f5C.y0[o+1], 0.005);
        min5C->SetVariable(o+2, Form("phi%d",i), f5C.y0[o+2], 0.005);
        min5C->SetLimitedVariable(o+3, Form("x%d",i), f5C.y0[o+3], 0.01, f5C.y0[o+3]-1.5, f5C.y0[o+3]+1.5);
    }

    // UPDATED: strict convergence logic with Status()==0 requirement
    for(int step=0; step<30; step++) {
        min5C->Minimize();
        const double* parTmp = min5C->X();
        double fvals[5];
        f5C.constraints(parTmp, fvals);

        double maxConstraintStep = 0.0;
        for(int k=0; k<5; k++){
            maxConstraintStep = std::max(maxConstraintStep, std::abs(fvals[k]));
        }
        if(!std::isfinite(maxConstraintStep)) break;
        if(step > 0 && min5C->Status() == 0 &&
           maxConstraintStep < 1.0e-5) break;

        for(int k=0; k<5; k++){
            f5C.lambda[k] += f5C.mu * fvals[k];
        }
        f5C.mu = std::min(f5C.mu * 2.0, 1.0e8);
    }

    const double* par5C = min5C->X();

    // DIAGNOSTIC UPDATE: record the standard 5C numerical state without
    // changing its acceptance or the fitted result.
    const int status_5C = min5C->Status();
    double diagnosticConstraints5C[5];
    f5C.constraints(par5C, diagnosticConstraints5C);
    double max_constraint_5C = 0.0;
    for(int k=0; k<5; k++){
        max_constraint_5C = std::max(
            max_constraint_5C,
            std::abs(diagnosticConstraints5C[k])
        );
    }

    P4 j_5C[4];
    for(int i=0; i<4; i++) {
        double alpha = par5C[4*i+0]; double theta = clampTheta(par5C[4*i+1]);
        double phi = wrapPhi(par5C[4*i+2]); double x = par5C[4*i+3];
        double ex = std::exp(x); double m = f5C.m0[i]; if(m<0.1) m=0.1;
        j_5C[i].SetPxPyPzE(alpha*m*ex*std::sin(theta)*std::cos(phi), alpha*m*ex*std::sin(theta)*std::sin(phi), alpha*m*ex*std::cos(theta), alpha*m*std::sqrt(ex*ex+1.0));
    }
    float m1_5C = (j_5C[a] + j_5C[b]).M();
    float m2_5C = (j_5C[c] + j_5C[d]).M();
    float mW_5C = 0.5 * (m1_5C + m2_5C);

    double dpar5C[16];
    for(int i=0;i<16;i++) dpar5C[i] = par5C[i] - f5C.y0[i];
    double chi2_5C = 0.0;
    for(int i=0;i<16;i++){
        double acc=0.0;
        for(int j=0;j<16;j++) acc += MyFit::Vinv[i][j]*dpar5C[j];
        chi2_5C += dpar5C[i]*acc;
    }

    // UPDATED: strict validity requires converged minimizer
    const bool standard_5C_valid =
        status_5C == 0 &&
        std::isfinite(chi2_5C) &&
        max_constraint_5C < 1.0e-5;
    const double prob_5C = standard_5C_valid
        ? ROOT::Math::chisquared_cdf_c(chi2_5C, 5.0)
        : 0.0;

    // PULL UPDATE: separate standard 5C pulls from the final hybrid 5C
    // candidate.  Low-probability events use their selected 5C+ISR pulls.
    float pull5C_standard[16];
    float pull5C_final[16];
    for(int i=0; i<16; i++){
        pull5C_standard[i] = std::numeric_limits<float>::quiet_NaN();
        pull5C_final[i] = std::numeric_limits<float>::quiet_NaN();
    }
    const int pull5C_standard_valid =
        (standard_5C_valid && computePullsStandard(f5C, par5C, pull5C_standard)) ? 1 : 0;
    for(int i=0; i<16; i++) pull5C_final[i] = pull5C_standard[i];
    int pull5C_final_valid = pull5C_standard_valid;
    float pull5C_final_yISR = std::numeric_limits<float>::quiet_NaN();

    // --- 5C + ISR treatment [ISR UPDATE] ---
    // This is the reconstruction compared with/without ISR in thesis Fig. 7.9.
    // The 3% condition triggers the refit; it is not an event filter.
    float mW_5C_isr = mW_5C;
    float chi2_5C_isr = (float)chi2_5C;
    float E_isr_fitted = 0.0f;
    int isr_applied = 0;

    // DIAGNOSTIC UPDATE: explicit ISR-attempt and candidate diagnostics.
    // Existing isr_applied semantics and all final branches are unchanged.
    const int isr_5C_attempted =
        (standard_5C_valid && prob_5C < 0.03 &&
         EMAX_ISR > 1.0e-9) ? 1 : 0;
    int isr_5C_fit_valid = 0;
    int status_5C_isr = -999;
    double max_constraint_5C_isr = -1.0;
    double y_isr_5C_fitted = 0.0;

    if(standard_5C_valid && prob_5C < 0.03 &&
       EMAX_ISR > 1.0e-9){
        const ISRFitResult best5CISR = runBestISRFit(
            true,
            f5C.y0,
            par5C,
            f5C.m0,
            sqrts,
            EMAX_ISR,
            B_ISR,
            a, b, c, d,
            rawMissingPz,
            "5C_zero",
            "5C_pz"
        );

        // DIAGNOSTIC UPDATE: retain diagnostics for both valid and invalid
        // returned candidates while preserving the original selection logic.
        isr_5C_fit_valid = best5CISR.valid ? 1 : 0;
        status_5C_isr = best5CISR.minimizerStatus;
        max_constraint_5C_isr = best5CISR.maxConstraint;
        y_isr_5C_fitted = best5CISR.yPhoton;

        if(best5CISR.valid){
            isr_applied = 1;
            mW_5C_isr = 0.5f * (float)(best5CISR.m1 + best5CISR.m2);
            chi2_5C_isr = (float)best5CISR.chi2;
            E_isr_fitted = (float)best5CISR.EPhoton;

            // PULL UPDATE: compute pulls from the actually selected 5C+ISR
            // solution instead of reusing standard 5C or 4C pulls.
            AugLagFunctorISR pullFit5CISR(true);
            pullFit5CISR.sqrts = sqrts;
            pullFit5CISR.Emax = EMAX_ISR;
            pullFit5CISR.bISR = B_ISR;
            pullFit5CISR.a = a; pullFit5CISR.b = b;
            pullFit5CISR.c = c; pullFit5CISR.d = d;
            for(int i=0; i<16; i++) pullFit5CISR.y0_raw[i] = f5C.y0[i];
            for(int i=0; i<4; i++) pullFit5CISR.m0[i] = f5C.m0[i];

            float pullISR[17];
            pull5C_final_valid = computePullsISR(
                pullFit5CISR, best5CISR.par, pullISR
            ) ? 1 : 0;
            for(int i=0; i<16; i++) pull5C_final[i] = pullISR[i];
            pull5C_final_yISR = pullISR[16];
        }
    }

    // DIAGNOSTIC UPDATE: explicit final recovery flag for the 5C study.
    const int isr_5C_recovered =
        (isr_applied > 0 &&
         ROOT::Math::chisquared_cdf_c(chi2_5C_isr, 5.0) > 0.03) ? 1 : 0;

    float m_small_raw = std::min(m1_raw, m2_raw);
    float m_large_raw = std::max(m1_raw, m2_raw);

    float m_small_4C = std::min(m1_4C, m2_4C);
    float m_large_4C = std::max(m1_4C, m2_4C);

    float m_small_4C_isr = std::min(m1_4C_isr, m2_4C_isr);
    float m_large_4C_isr = std::max(m1_4C_isr, m2_4C_isr);

    // Existing output indices 0..19 are preserved for compatibility.
    // DIAGNOSTIC UPDATE: diagnostic values are appended at indices 20..35.
    return ROOT::VecOps::RVec<float>{
        m1_raw, m2_raw,                    // 0, 1
        m_small_raw, m_large_raw,          // 2, 3
        m_small_4C, m_large_4C,            // 4, 5
        (float)chi2_4C,                    // 6
        mW_5C, (float)chi2_5C,             // 7, 8
        mW_5C_isr, chi2_5C_isr,            // 9, 10
        E_isr_fitted,                      // 11
        (float)isr_applied,                // 12
        (float)EMAX_ISR,                   // 13
        (float)B_ISR,                      // 14
        m_small_4C_isr, m_large_4C_isr,    // 15, 16
        chi2_4C_isr,                       // 17
        E_isr_4C_fitted,                   // 18
        (float)isr_4C_applied,             // 19

        // DIAGNOSTIC UPDATE: standard-fit convergence diagnostics
        (float)status_4C,                  // 20
        (float)max_constraint_4C,          // 21
        (float)status_5C,                  // 22
        (float)max_constraint_5C,          // 23

        // DIAGNOSTIC UPDATE: 4C+ISR candidate diagnostics
        (float)isr_4C_attempted,           // 24
        (float)isr_4C_fit_valid,           // 25
        (float)status_4C_isr,              // 26
        (float)max_constraint_4C_isr,      // 27
        (float)y_isr_4C_fitted,            // 28
        (float)isr_4C_recovered,           // 29

        // DIAGNOSTIC UPDATE: 5C+ISR candidate diagnostics
        (float)isr_5C_attempted,           // 30
        (float)isr_5C_fit_valid,           // 31
        (float)status_5C_isr,              // 32
        (float)max_constraint_5C_isr,      // 33
        (float)y_isr_5C_fitted,            // 34
        (float)isr_5C_recovered,           // 35

        // PULL UPDATE: standard and final/hybrid 4C pulls
        pull4C_standard[0], pull4C_standard[1], pull4C_standard[2], pull4C_standard[3],
        pull4C_standard[4], pull4C_standard[5], pull4C_standard[6], pull4C_standard[7],
        pull4C_standard[8], pull4C_standard[9], pull4C_standard[10], pull4C_standard[11],
        pull4C_standard[12], pull4C_standard[13], pull4C_standard[14], pull4C_standard[15], // 36..51

        pull4C_final[0], pull4C_final[1], pull4C_final[2], pull4C_final[3],
        pull4C_final[4], pull4C_final[5], pull4C_final[6], pull4C_final[7],
        pull4C_final[8], pull4C_final[9], pull4C_final[10], pull4C_final[11],
        pull4C_final[12], pull4C_final[13], pull4C_final[14], pull4C_final[15], // 52..67

        // PULL UPDATE: standard and final/hybrid 5C pulls
        pull5C_standard[0], pull5C_standard[1], pull5C_standard[2], pull5C_standard[3],
        pull5C_standard[4], pull5C_standard[5], pull5C_standard[6], pull5C_standard[7],
        pull5C_standard[8], pull5C_standard[9], pull5C_standard[10], pull5C_standard[11],
        pull5C_standard[12], pull5C_standard[13], pull5C_standard[14], pull5C_standard[15], // 68..83

        pull5C_final[0], pull5C_final[1], pull5C_final[2], pull5C_final[3],
        pull5C_final[4], pull5C_final[5], pull5C_final[6], pull5C_final[7],
        pull5C_final[8], pull5C_final[9], pull5C_final[10], pull5C_final[11],
        pull5C_final[12], pull5C_final[13], pull5C_final[14], pull5C_final[15], // 84..99

        pull4C_final_yISR,                 // 100
        pull5C_final_yISR,                 // 101
        (float)pull4C_standard_valid,       // 102
        (float)pull4C_final_valid,          // 103
        (float)pull5C_standard_valid,       // 104
        (float)pull5C_final_valid,          // 105
        (float)standard_4C_valid,           // 106
        (float)standard_5C_valid            // 107
    };

    }

''')

# =========================================================
# 2) covariance matrix 
# =========================================================

_JSON_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    _ENERGY_CONFIG["json"],
)

def build_covariance_matrix(json_path=_JSON_PATH):
    with open(json_path) as f:
        cb = json.load(f)

    # Jet ordering: j1=index0 (highest energy) ... j4=index3 (lowest energy)
    # matching the order returned by clustering_ee_kt in FCCAnalyses.
    # Variable ordering within each jet block: alpha(0), theta(1), phi(2), x(3)
    # — same as the parametrization in Eq. (7.23) of the thesis.
    jet_vars = [
        ("alpha", 0),
        ("theta", 1),
        ("phi",   2),
        ("x",     3),
    ]

    cov = np.zeros((16, 16))
    for jet_idx in range(4):          # 0-based
        jet_label = jet_idx + 1       # JSON uses j1..j4
        for var, var_idx in jet_vars:
            key = f"delta_{var}_j{jet_label}"
            sigma = cb[key]["sigma"]
            param_idx = 4 * jet_idx + var_idx
            cov[param_idx, param_idx] = sigma ** 2

    return cov

cov = build_covariance_matrix()

# UPDATED: removed the regularization term cov += np.eye(16) * 1e-8

cov_inv = np.linalg.inv(cov)
ROOT.MyFit.setCovariance16(
    cov.flatten().tolist(),
    cov_inv.flatten().tolist(),
)

# Print diagonal sigmas for quick visual confirmation at run-time
print("=== Covariance matrix diagonal (sigma values) ===")
var_names = ["alpha", "theta", "phi", "x"]
for jet_idx in range(4):
    for var_idx, var in enumerate(var_names):
        idx = 4*jet_idx + var_idx
        print(f"  jet{jet_idx+1} {var:6s}: sigma = {np.sqrt(cov[idx,idx]):.6f}")
print("=" * 50)

# =========================================================
# 3) FCCAnalyses RDF Configuration
# =========================================================

processList = {
    _ENERGY_CONFIG["process"]: {
        "fraction": float(os.environ.get("WMASS_FRACTION", "1.0e-6")),
        "chunks": 1,
        "output": f"wmass_fit_pvalue_pulls_ISR_ecm{_ENERGY_CONFIG['tag']}"
    }
}
inputDir = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "FCCee_procDict_winter2023_IDEA.json"
outDir = "./outputs/wmass"

nCPUS = 10

doTree = True

includePaths = [
    "headers/selectQuarks.h"
]

# =========================================================

class RDFanalysis:

    def analysers(df):
        df = (
            df


            .Define("partons_all", "selectQuarks(Particle)")
            .Define("n_partons", "partons_all.size()")
            # hadronic WW truth selection
            #.Define("pass_n_partons_4", "n_partons == 4")
            .Filter("n_partons == 4", "Require only 4 partons")

            # Generator truth used later for thesis Figure 7.8.
            # With exactly four selected W-decay quarks, sqrt(s)-E_hard is
            # the energy removed from the hard WW system by ISR/beamstrahlung.
            .Define(
                "parton_e_truth",
                "FCCAnalyses::MCParticle::get_e(partons_all)"
            )
            .Define(
                "parton_pz_truth",
                "FCCAnalyses::MCParticle::get_pz(partons_all)"
            )
            .Define(
                "E_isr_true",
                f"std::max(0.0f, {float(ECM)}f - "
                "ROOT::VecOps::Sum(parton_e_truth))"
            )
            .Define(
                "pz_isr_true",
                "-ROOT::VecOps::Sum(parton_pz_truth)"
            )
            .Define("E_isr_true_collinear", "std::abs(pz_isr_true)")

            .Define("RP_px", "ReconstructedParticle::get_px(ReconstructedParticles)")
            .Define("RP_py", "ReconstructedParticle::get_py(ReconstructedParticles)")
            .Define("RP_pz", "ReconstructedParticle::get_pz(ReconstructedParticles)")
            .Define("RP_e",  "ReconstructedParticle::get_e(ReconstructedParticles)")
            .Define("RP_m",  "ReconstructedParticle::get_mass(ReconstructedParticles)")

            # --- ee_kt clustering ---
            .Define("pseudo_jets", "JetClusteringUtils::set_pseudoJets_xyzm(RP_px,RP_py,RP_pz,RP_m)")
            .Define("Jets", "JetClustering::clustering_ee_kt(2,4,0,0)(pseudo_jets)")
            .Define("jets", "JetClusteringUtils::get_pseudoJets(Jets)")

            # --- jet variables ---
            .Define("jet_px", "JetClusteringUtils::get_px(jets)")
            .Define("jet_py", "JetClusteringUtils::get_py(jets)") 
            .Define("jet_pz", "JetClusteringUtils::get_pz(jets)")
            .Define("jet_e",  "JetClusteringUtils::get_e(jets)")

            # =========================================================
            # FILTER FLAG 2: exactly 4 reconstructed jets
            # =========================================================
            .Define("pass_n_jets_4", "jet_px.size() == 4")
            .Filter("jet_px.size()==4")

         
            
            .Define("all_results", f"reconstructWW(jet_px, jet_py, jet_pz, jet_e, {float(ECM)}f)")

            .Define("mW1_raw", "all_results[0]")
            .Define("mW2_raw", "all_results[1]")
            .Define("m_small_raw", "all_results[2]")
            .Define("m_large_raw", "all_results[3]")

            # Standard 4C result
            .Define("m_small_4C", "all_results[4]")
            .Define("m_large_4C", "all_results[5]")
            .Define("chi2_4C", "all_results[6]")

            # Standard and ISR-treated 5C result
            .Define("mW_5C", "all_results[7]")
            .Define("chi2_5C", "all_results[8]")
            .Define("mW_5C_isr", "all_results[9]")
            .Define("chi2_5C_isr", "all_results[10]")
            .Define("E_isr_fitted", "all_results[11]")
            .Define("isr_applied", "all_results[12]")

            .Define("Emax_isr", "all_results[13]")
            .Define("b_isr", "all_results[14]")

            # ISR UPDATE: 4C+ISR outputs appended after legacy indices 0..14.
            .Define("m_small_4C_isr", "all_results[15]")
            .Define("m_large_4C_isr", "all_results[16]")
            .Define("chi2_4C_isr", "all_results[17]")
            .Define("E_isr_4C_fitted", "all_results[18]")
            .Define("isr_4C_applied", "all_results[19]")

            # DIAGNOSTIC UPDATE: appended outputs only; legacy indices and
            # existing reconstruction/selection branches remain unchanged.
            .Define("status_4C", "all_results[20]")
            .Define("max_constraint_4C", "all_results[21]")
            .Define("status_5C", "all_results[22]")
            .Define("max_constraint_5C", "all_results[23]")

            .Define("isr_4C_attempted", "all_results[24]")
            .Define("isr_4C_fit_valid", "all_results[25]")
            .Define("status_4C_isr", "all_results[26]")
            .Define("max_constraint_4C_isr", "all_results[27]")
            .Define("y_isr_4C_fitted", "all_results[28]")
            .Define("isr_4C_recovered", "all_results[29]")

            .Define("isr_5C_attempted", "all_results[30]")
            .Define("isr_5C_fit_valid", "all_results[31]")
            .Define("status_5C_isr", "all_results[32]")
            .Define("max_constraint_5C_isr", "all_results[33]")
            .Define("y_isr_5C_fitted", "all_results[34]")
            .Define("isr_5C_recovered", "all_results[35]")
            .Define("standard_4C_valid", "all_results[106]")
            .Define("standard_5C_valid", "all_results[107]")

            # =========================================================
            # p-value / fit goodness
            #
            # Thesis Sec. 7.3.2.2: P < 3% triggers the ISR refit.
            # No event is rejected here; final cuts are applied downstream.
            # =========================================================
            .Define(
                "prob_4C",
                "(standard_4C_valid > 0.5f) ? "
                "ROOT::Math::chisquared_cdf_c(chi2_4C, 4.0) : 0.0"
            )
            .Define(
                "prob_4C_isr",
                "(isr_4C_fit_valid > 0.5f) ? "
                "ROOT::Math::chisquared_cdf_c(chi2_4C_isr, 4.0) : 0.0"
            )
            .Define(
                "prob_final_4C",
                "(isr_4C_applied > 0.5f) ? prob_4C_isr : prob_4C"
            )
            .Define(
                "m_small_4C_final",
                "(isr_4C_applied > 0.5f) ? m_small_4C_isr : m_small_4C"
            )
            .Define(
                "m_large_4C_final",
                "(isr_4C_applied > 0.5f) ? m_large_4C_isr : m_large_4C"
            )
            .Define(
                "pass_final_4C_p03",
                "standard_4C_valid > 0.5f && prob_final_4C > 0.03"
            )

            .Define(
                "prob_5C",
                "(standard_5C_valid > 0.5f) ? "
                "ROOT::Math::chisquared_cdf_c(chi2_5C, 5.0) : 0.0"
            )
            .Define(
                "prob_5C_isr",
                "(isr_5C_fit_valid > 0.5f) ? "
                "ROOT::Math::chisquared_cdf_c(chi2_5C_isr, 5.0) : 0.0"
            )
            .Define(
                "prob_final_5C",
                "(isr_applied > 0.5f) ? prob_5C_isr : prob_5C"
            )
            .Define(
                "mW_5C_final",
                "(isr_applied > 0.5f) ? mW_5C_isr : mW_5C"
            )
            .Define(
                "pass_final_5C_p03",
                "standard_5C_valid > 0.5f && prob_final_5C > 0.03"
            )

        )

        # PULL UPDATE: the four blocks are appended after legacy/diagnostic
        # outputs 0..35. Variable order is alpha, theta, phi, x for jets 1..4.
        pull_variables = ["alpha", "theta", "phi", "x"]
        pull_blocks = [
            ("pull4C_standard", 36),
            ("pull4C_final", 52),
            ("pull5C_standard", 68),
            ("pull5C_final", 84),
        ]

        for prefix, start_index in pull_blocks:
            for jet_index in range(4):
                for variable_index, variable in enumerate(pull_variables):
                    result_index = start_index + 4 * jet_index + variable_index
                    branch = f"{prefix}_{variable}_j{jet_index + 1}"
                    df = df.Define(branch, f"all_results[{result_index}]")

        df = (
            df
            .Define("pull4C_final_yISR", "all_results[100]")
            .Define("pull5C_final_yISR", "all_results[101]")
            .Define("pull4C_standard_valid", "all_results[102]")
            .Define("pull4C_final_valid", "all_results[103]")
            .Define("pull5C_standard_valid", "all_results[104]")
            .Define("pull5C_final_valid", "all_results[105]")
        )

        return df

    def output():
        branches = [
            "mW1_raw", "mW2_raw",
            "m_small_raw", "m_large_raw",

            # Standard 4C and ISR-treated/final 4C
            "m_small_4C", "m_large_4C", "chi2_4C", "prob_4C",
            "m_small_4C_isr", "m_large_4C_isr", "chi2_4C_isr",
            "E_isr_4C_fitted", "isr_4C_applied", "prob_4C_isr",
            "m_small_4C_final", "m_large_4C_final", "prob_final_4C",

            # Standard 5C and ISR-treated/final 5C
            "mW_5C", "chi2_5C", "prob_5C",
            "mW_5C_isr", "chi2_5C_isr", "E_isr_fitted",
            "isr_applied", "prob_5C_isr",
            "mW_5C_final", "prob_final_5C",

            # ISR parametrization diagnostics
            "Emax_isr", "b_isr",

            # Generator-level radiated-energy truth for thesis Figure 7.8
            "E_isr_true", "E_isr_true_collinear", "pz_isr_true",

            # DIAGNOSTIC UPDATE: standard 4C/5C numerical diagnostics
            "standard_4C_valid", "standard_5C_valid",
            "status_4C", "max_constraint_4C",
            "status_5C", "max_constraint_5C",

            # DIAGNOSTIC UPDATE: explicit 4C+ISR attempt/result diagnostics
            "isr_4C_attempted", "isr_4C_fit_valid",
            "status_4C_isr", "max_constraint_4C_isr",
            "y_isr_4C_fitted", "isr_4C_recovered",

            # DIAGNOSTIC UPDATE: explicit 5C+ISR attempt/result diagnostics
            "isr_5C_attempted", "isr_5C_fit_valid",
            "status_5C_isr", "max_constraint_5C_isr",
            "y_isr_5C_fitted", "isr_5C_recovered",

            # Final 3% rejection flags; no .Filter() is applied in producer.
            "pass_final_4C_p03", "pass_final_5C_p03",
        ]

        pull_variables = ["alpha", "theta", "phi", "x"]
        for prefix in [
            "pull4C_standard",
            "pull4C_final",
            "pull5C_standard",
            "pull5C_final",
        ]:
            for jet_index in range(1, 5):
                for variable in pull_variables:
                    branches.append(f"{prefix}_{variable}_j{jet_index}")

        branches += [
            "pull4C_final_yISR", "pull5C_final_yISR",
            "pull4C_standard_valid", "pull4C_final_valid",
            "pull5C_standard_valid", "pull5C_final_valid",
        ]

        return branches