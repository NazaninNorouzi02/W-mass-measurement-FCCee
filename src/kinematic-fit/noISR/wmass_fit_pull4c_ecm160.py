import ROOT
import numpy as np
import json
import os

ROOT.gROOT.SetBatch(True)

# =========================================================
# 1) C++ inline code : 4C / 5C kinematic fit, ISR, and corrected 4C pulls
# =========================================================

ROOT.gInterpreter.Declare(r'''
#include <vector>
#include <cmath>
#include <algorithm>
#include <memory>
#include <iostream>
#include <limits>
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
                         const std::vector<double>& invflat) {
        for(int i=0;i<16;i++){
            for(int j=0;j<16;j++){
                V[i][j]    = vflat[i*16+j];
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

inline void buildJetsFromPars(const double *par, const double *m0, P4 j[4]) {
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
}

inline void constraints4COnly(const double *par,
                              const double *m0,
                              double sqrts,
                              double f[4]) {
    P4 j[4];
    buildJetsFromPars(par, m0, j);

    f[0] = 0.0;
    f[1] = 0.0;
    f[2] = 0.0;
    f[3] = -sqrts;

    for(int i=0; i<4; i++){
        f[0] += j[i].Px();
        f[1] += j[i].Py();
        f[2] += j[i].Pz();
        f[3] += j[i].E();
    }
}

inline bool invert4x4(const double A[4][4], double invA[4][4]) {
    double aug[4][8];

    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            aug[i][j] = A[i][j];
        }
        for(int j=0; j<4; j++){
            aug[i][j+4] = (i == j) ? 1.0 : 0.0;
        }
    }

    for(int col=0; col<4; col++){
        int pivot = col;
        double maxAbs = std::abs(aug[col][col]);

        for(int row=col+1; row<4; row++){
            double val = std::abs(aug[row][col]);
            if(val > maxAbs){
                maxAbs = val;
                pivot = row;
            }
        }

        if(maxAbs < 1.0e-18){
            return false;
        }

        if(pivot != col){
            for(int j=0; j<8; j++){
                std::swap(aug[col][j], aug[pivot][j]);
            }
        }

        double diag = aug[col][col];
        for(int j=0; j<8; j++){
            aug[col][j] /= diag;
        }

        for(int row=0; row<4; row++){
            if(row == col) continue;

            double factor = aug[row][col];
            for(int j=0; j<8; j++){
                aug[row][j] -= factor * aug[col][j];
            }
        }
    }

    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            invA[i][j] = aug[i][j+4];
        }
    }

    return true;
}

inline void jacobian4CNumeric(const double *par,
                              const double *m0,
                              double sqrts,
                              double A[4][16]) {
    double pplus[16];
    double pminus[16];
    double fplus[4];
    double fminus[4];

    for(int ip=0; ip<16; ip++){
        for(int j=0; j<16; j++){
            pplus[j] = par[j];
            pminus[j] = par[j];
        }

        double h = 1.0e-5 * std::max(1.0, std::abs(par[ip]));

        pplus[ip]  += h;
        pminus[ip] -= h;

        if(ip % 4 == 1){
            pplus[ip]  = clampTheta(pplus[ip]);
            pminus[ip] = clampTheta(pminus[ip]);
        }

        if(ip % 4 == 2){
            pplus[ip]  = wrapPhi(pplus[ip]);
            pminus[ip] = wrapPhi(pminus[ip]);
        }

        constraints4COnly(pplus,  m0, sqrts, fplus);
        constraints4COnly(pminus, m0, sqrts, fminus);

        for(int k=0; k<4; k++){
            A[k][ip] = (fplus[k] - fminus[k]) / (2.0*h);
        }
    }
}

inline void computeCorrectPulls4C(const double y0[16],
                                  const double m0[4],
                                  double sqrts,
                                  float pulls[16]) {
    for(int i=0; i<16; i++){
        pulls[i] = std::numeric_limits<float>::quiet_NaN();
    }

    double yfit[16];
    for(int i=0; i<16; i++){
        yfit[i] = y0[i];
    }

    // =========================================================
    // Iterative linearized constrained least-squares update:
    //
    // y_fit <- y_fit - V A^T (A V A^T)^(-1) f(y_fit)
    //
    // This is only for the 4C pull calculation.
    // It does not change the 5C / ISR kinematic fit.
    // =========================================================

    bool ok = true;

    for(int iter=0; iter<8; iter++){
        double f[4];
        double A[4][16];

        constraints4COnly(yfit, m0, sqrts, f);
        jacobian4CNumeric(yfit, m0, sqrts, A);

        double tmp[16][4]; // V A^T
        for(int i=0; i<16; i++){
            for(int k=0; k<4; k++){
                tmp[i][k] = 0.0;
                for(int j=0; j<16; j++){
                    tmp[i][k] += MyFit::V[i][j] * A[k][j];
                }
            }
        }

        double M[4][4]; // A V A^T
        for(int k=0; k<4; k++){
            for(int l=0; l<4; l++){
                M[k][l] = 0.0;
                for(int j=0; j<16; j++){
                    M[k][l] += A[k][j] * tmp[j][l];
                }
            }
        }

        double Minv[4][4];
        if(!invert4x4(M, Minv)){
            ok = false;
            break;
        }

        double corr[16];
        for(int i=0; i<16; i++){
            corr[i] = 0.0;
            for(int k=0; k<4; k++){
                for(int l=0; l<4; l++){
                    corr[i] += tmp[i][k] * Minv[k][l] * f[l];
                }
            }
        }

        for(int i=0; i<16; i++){
            yfit[i] -= corr[i];
        }

        double normF = std::sqrt(f[0]*f[0] + f[1]*f[1] + f[2]*f[2] + f[3]*f[3]);
        if(normF < 1.0e-7){
            break;
        }
    }

    if(!ok){
        return;
    }

    // =========================================================
    // Pull denominator:
    //
    // sigma_pull_den^2 = V_meas - V_fit
    //                  = V A^T (A V A^T)^(-1) A V
    //
    // pull_i = (y_meas_i - y_fit_i) / sqrt(sigma_pull_den_i^2)
    // =========================================================

    double A[4][16];
    jacobian4CNumeric(yfit, m0, sqrts, A);

    double tmp[16][4]; // V A^T
    for(int i=0; i<16; i++){
        for(int k=0; k<4; k++){
            tmp[i][k] = 0.0;
            for(int j=0; j<16; j++){
                tmp[i][k] += MyFit::V[i][j] * A[k][j];
            }
        }
    }

    double M[4][4]; // A V A^T
    for(int k=0; k<4; k++){
        for(int l=0; l<4; l++){
            M[k][l] = 0.0;
            for(int j=0; j<16; j++){
                M[k][l] += A[k][j] * tmp[j][l];
            }
        }
    }

    double Minv[4][4];
    if(!invert4x4(M, Minv)){
        return;
    }

    for(int i=0; i<16; i++){
        double cii = 0.0;

        for(int k=0; k<4; k++){
            for(int l=0; l<4; l++){
                cii += tmp[i][k] * Minv[k][l] * tmp[i][l];
            }
        }

        if(cii > 0.0 && std::isfinite(cii)){
            pulls[i] = (float)((y0[i] - yfit[i]) / std::sqrt(cii));
        }
    }
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

    void constraints(const double *par, double f[5]) const {
        P4 j[4];
        buildJetsFromPars(par, m0, j);

        f[0]=f[1]=f[2]=f[3]=f[4]=0.0;

        for(int i=0; i<4; i++){
            f[0] += j[i].Px();
            f[1] += j[i].Py();
            f[2] += j[i].Pz();
            f[3] += j[i].E();
        }

        f[3] -= sqrts;

        if(do5C){
            f[4] = (j[a]+j[b]).M() - (j[c]+j[d]).M();
        }
    }

    double operator()(const double *par) const {
        double chi2 = 0.0;
        double dpar[16];

        for(int i=0; i<16; i++){
            dpar[i] = par[i] - y0[i];
        }

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

// --- ISR photon: unchanged from your current code ---
inline double isrEmax(double sqrts, double mW){
    double x = 4.0 * mW * mW / (sqrts * sqrts);
    if (x >= 1.0) return 0.0;
    return 0.5 * sqrts * (1.0 - x);
}

inline double pzGamma(double y, double Emax, double b){
    double s = (y >= 0.0) ? 1.0 : -1.0;
    double val = std::erf(std::abs(y)/std::sqrt(2.0));
    return s * Emax * std::pow(val, 1.0/b);
}

struct AugLagFunctorISR {
    double y0_raw[16];
    double y0_start[17];
    double m0[4];
    float sqrts;
    double Emax;
    double bExp;
    int a, b, c, d;
    mutable double lambda[5];
    mutable double mu;

    AugLagFunctorISR() {
        for(int i=0; i<5; i++) lambda[i] = 0.0;
        mu = 10.0;
        bExp = 2.0;
    }

    void constraints(const double *par, double f[5]) const {
        P4 j[4];
        buildJetsFromPars(par, m0, j);

        double y_isr = par[16];
        double pzg = pzGamma(y_isr, Emax, bExp);
        double Eg  = std::abs(pzg);

        f[0]=f[1]=f[2]=f[3]=f[4]=0.0;

        for(int i=0; i<4; i++){
            f[0] += j[i].Px();
            f[1] += j[i].Py();
            f[2] += j[i].Pz();
            f[3] += j[i].E();
        }

        f[2] += pzg;
        f[3] += Eg;
        f[3] -= sqrts;

        f[4] = (j[a]+j[b]).M() - (j[c]+j[d]).M();
    }

    double operator()(const double *par) const {
        double chi2 = 0.0;
        double dpar[16];

        for(int i=0; i<16; i++){
            dpar[i] = par[i] - y0_raw[i];
        }

        for(int i=0; i<16; i++){
            double acc = 0.0;
            for(int j=0; j<16; j++){
                acc += MyFit::Vinv[i][j] * dpar[j];
            }
            chi2 += dpar[i] * acc;
        }

        double y_isr = par[16];
        chi2 += y_isr * y_isr;

        double f[5];
        constraints(par, f);

        double L = chi2;
        double f2 = 0.0;

        for(int k=0; k<5; k++){
            L += lambda[k] * f[k];
            f2 += f[k] * f[k];
        }

        L += 0.5 * mu * f2;
        return L;
    }
};

ROOT::VecOps::RVec<float>
reconstructWW(const ROOT::VecOps::RVec<float>& px,
              const ROOT::VecOps::RVec<float>& py,
              const ROOT::VecOps::RVec<float>& pz,
              const ROOT::VecOps::RVec<float>& E,
              float sqrts) {

    const int pairings[3][4] = { {0,1,2,3}, {0,2,1,3}, {0,3,1,2} };

    P4 j_raw[4];
    for(int i=0; i<4; i++){
        j_raw[i].SetPxPyPzE(px[i], py[i], pz[i], E[i]);
    }

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

    int a = pairings[best_ip][0];
    int b = pairings[best_ip][1];
    int c = pairings[best_ip][2];
    int d = pairings[best_ip][3];

    float m1_raw = (j_raw[a] + j_raw[b]).M();
    float m2_raw = (j_raw[c] + j_raw[d]).M();

    // =========================================================
    // 4C fit: unchanged
    // =========================================================

    AugLagFunctor f4C(false);
    f4C.sqrts = sqrts;
    f4C.a = a;
    f4C.b = b;
    f4C.c = c;
    f4C.d = d;

    for(int i=0; i<4; i++){
        double p = j_raw[i].P();
        double m = j_raw[i].M();

        double theta = (p > 0) ? std::acos(std::clamp(j_raw[i].Pz()/p, -1.0, 1.0)) : M_PI/2.0;
        double phi = std::atan2(j_raw[i].Py(), j_raw[i].Px());
        double x_init = (m > 0.1 && p > 0.1) ? std::log(std::clamp(p/m, 0.1, 100.0)) : 0.0;

        f4C.m0[i] = m;
        f4C.y0[4*i+0] = 1.0;
        f4C.y0[4*i+1] = clampTheta(theta);
        f4C.y0[4*i+2] = wrapPhi(phi);
        f4C.y0[4*i+3] = std::clamp(x_init, -4.0, 4.0);
    }

    // =========================================================
    // Corrected 4C pulls:
    // calculated separately from the old Minuit 4C fit.
    // This does not affect the kinematic outputs.
    // =========================================================

    float pull4C[16];
    computeCorrectPulls4C(f4C.y0, f4C.m0, sqrts, pull4C);

    ROOT::Math::Functor functor4C(f4C, 16);
    auto min4C = std::unique_ptr<ROOT::Math::Minimizer>(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad")
    );

    min4C->SetFunction(functor4C);
    min4C->SetStrategy(1);
    min4C->SetTolerance(1e-3);
    min4C->SetMaxFunctionCalls(5000);

    for(int i=0; i<4; i++){
        int o = 4*i;
        min4C->SetLimitedVariable(o+0, Form("alpha%d",i), 1.0, 0.02, 0.1, 2.0);
        min4C->SetVariable(o+1, Form("theta%d",i), f4C.y0[o+1], 0.005);
        min4C->SetVariable(o+2, Form("phi%d",i), f4C.y0[o+2], 0.005);
        min4C->SetLimitedVariable(o+3, Form("x%d",i), f4C.y0[o+3], 0.01, f4C.y0[o+3]-1.5, f4C.y0[o+3]+1.5);
    }

    for(int step=0; step<6; step++) {
        min4C->Minimize();
        const double* parTmp = min4C->X();

        double fvals[5];
        f4C.constraints(parTmp, fvals);

        for(int k=0; k<4; k++){
            f4C.lambda[k] += f4C.mu * fvals[k];
        }

        f4C.mu *= 1.5;
    }

    const double* par4C = min4C->X();

    P4 j_4C[4];
    buildJetsFromPars(par4C, f4C.m0, j_4C);

    float m1_4C = (j_4C[a] + j_4C[b]).M();
    float m2_4C = (j_4C[c] + j_4C[d]).M();

    double dpar4C[16];
    for(int i=0;i<16;i++){
        dpar4C[i] = par4C[i] - f4C.y0[i];
    }

    double chi2_4C = 0.0;
    for(int i=0;i<16;i++){
        double acc=0.0;
        for(int j=0;j<16;j++){
            acc += MyFit::Vinv[i][j]*dpar4C[j];
        }
        chi2_4C += dpar4C[i]*acc;
    }

    // =========================================================
    // 5C fit: unchanged
    // =========================================================

    AugLagFunctor f5C(true);
    f5C.sqrts = sqrts;
    f5C.a = a;
    f5C.b = b;
    f5C.c = c;
    f5C.d = d;

    std::copy(std::begin(f4C.y0), std::end(f4C.y0), std::begin(f5C.y0));
    std::copy(std::begin(f4C.m0), std::end(f4C.m0), std::begin(f5C.m0));

    ROOT::Math::Functor functor5C(f5C, 16);
    auto min5C = std::unique_ptr<ROOT::Math::Minimizer>(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad")
    );

    min5C->SetFunction(functor5C);
    min5C->SetStrategy(1);
    min5C->SetTolerance(1e-3);
    min5C->SetMaxFunctionCalls(5000);

    for(int i=0; i<4; i++){
        int o = 4*i;
        min5C->SetLimitedVariable(o+0, Form("alpha%d",i), 1.0, 0.02, 0.1, 2.0);
        min5C->SetVariable(o+1, Form("theta%d",i), f5C.y0[o+1], 0.005);
        min5C->SetVariable(o+2, Form("phi%d",i), f5C.y0[o+2], 0.005);
        min5C->SetLimitedVariable(o+3, Form("x%d",i), f5C.y0[o+3], 0.01, f5C.y0[o+3]-1.5, f5C.y0[o+3]+1.5);
    }

    for(int step=0; step<6; step++) {
        min5C->Minimize();
        const double* parTmp = min5C->X();

        double fvals[5];
        f5C.constraints(parTmp, fvals);

        for(int k=0; k<5; k++){
            f5C.lambda[k] += f5C.mu * fvals[k];
        }

        f5C.mu *= 1.5;
    }

    const double* par5C = min5C->X();

    P4 j_5C[4];
    buildJetsFromPars(par5C, f5C.m0, j_5C);

    float m1_5C = (j_5C[a] + j_5C[b]).M();
    float m2_5C = (j_5C[c] + j_5C[d]).M();
    float mW_5C = 0.5 * (m1_5C + m2_5C);

    double dpar5C[16];
    for(int i=0;i<16;i++){
        dpar5C[i] = par5C[i] - f5C.y0[i];
    }

    double chi2_5C = 0.0;
    for(int i=0;i<16;i++){
        double acc=0.0;
        for(int j=0;j<16;j++){
            acc += MyFit::Vinv[i][j]*dpar5C[j];
        }
        chi2_5C += dpar5C[i]*acc;
    }

    // =========================================================
    // 5C ISR treatment: unchanged
    // =========================================================

    double prob_5C = ROOT::Math::chisquared_cdf_c(chi2_5C, 5.0);

    float mW_5C_isr = mW_5C;
    float chi2_5C_isr = (float)chi2_5C;
    float E_isr_fitted = 0.0f;
    int isr_applied = 0;

    if(prob_5C < 0.03){
        isr_applied = 1;

        AugLagFunctorISR fISR;
        fISR.sqrts = sqrts;
        fISR.a = a;
        fISR.b = b;
        fISR.c = c;
        fISR.d = d;

        for(int i=0;i<16;i++){
            fISR.y0_raw[i] = f5C.y0[i];
            fISR.y0_start[i] = par5C[i];
        }

        fISR.y0_start[16] = 0.0;
        std::copy(std::begin(f5C.m0), std::end(f5C.m0), std::begin(fISR.m0));

        float mW_for_isr = (mW_5C > 40.0f && mW_5C < sqrts/2.0f) ? mW_5C : 80.385f;
        fISR.Emax = isrEmax(sqrts, mW_for_isr);

        ROOT::Math::Functor functorISR(fISR, 17);
        auto minISR = std::unique_ptr<ROOT::Math::Minimizer>(
            ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad")
        );

        minISR->SetFunction(functorISR);
        minISR->SetStrategy(1);
        minISR->SetTolerance(1e-3);
        minISR->SetMaxFunctionCalls(5000);

        for(int i=0; i<4; i++){
            int o = 4*i;
            minISR->SetLimitedVariable(o+0, Form("alpha%d",i), fISR.y0_start[o+0], 0.02, 0.1, 2.0);
            minISR->SetVariable(o+1, Form("theta%d",i), fISR.y0_start[o+1], 0.005);
            minISR->SetVariable(o+2, Form("phi%d",i), fISR.y0_start[o+2], 0.005);
            minISR->SetLimitedVariable(o+3, Form("x%d",i), fISR.y0_start[o+3], 0.01, fISR.y0_raw[o+3]-1.5, fISR.y0_raw[o+3]+1.5);
        }

        minISR->SetLimitedVariable(16, "y_isr", fISR.y0_start[16], 0.05, -5.0, 5.0);

        for(int step=0; step<6; step++) {
            minISR->Minimize();
            const double* parTmp = minISR->X();

            double fvals[5];
            fISR.constraints(parTmp, fvals);

            for(int k=0; k<5; k++){
                fISR.lambda[k] += fISR.mu * fvals[k];
            }

            fISR.mu *= 1.5;
        }

        const double* parISR = minISR->X();

        P4 j_ISR[4];
        buildJetsFromPars(parISR, fISR.m0, j_ISR);

        float m1_isr = (j_ISR[a] + j_ISR[b]).M();
        float m2_isr = (j_ISR[c] + j_ISR[d]).M();
        mW_5C_isr = 0.5f * (m1_isr + m2_isr);

        double pzg = pzGamma(parISR[16], fISR.Emax, fISR.bExp);
        E_isr_fitted = (float)std::abs(pzg);

        double dparISR[16];
        for(int i=0;i<16;i++){
            dparISR[i] = parISR[i] - fISR.y0_raw[i];
        }

        double chi2ISR = 0.0;
        for(int i=0;i<16;i++){
            double acc=0.0;
            for(int j=0;j<16;j++){
                acc += MyFit::Vinv[i][j]*dparISR[j];
            }
            chi2ISR += dparISR[i]*acc;
        }

        chi2ISR += parISR[16]*parISR[16];
        chi2_5C_isr = (float)chi2ISR;
    }

    float m_small_raw = std::min(m1_raw, m2_raw);
    float m_large_raw = std::max(m1_raw, m2_raw);
    
    float m_small_4C = std::min(m1_4C, m2_4C);
    float m_large_4C = std::max(m1_4C, m2_4C);

    ROOT::VecOps::RVec<float> out;

    out.push_back(m1_raw);              // 0
    out.push_back(m2_raw);              // 1
    out.push_back(m_small_raw);         // 2
    out.push_back(m_large_raw);         // 3
    out.push_back(m_small_4C);          // 4
    out.push_back(m_large_4C);          // 5
    out.push_back((float)chi2_4C);      // 6
    out.push_back(mW_5C);               // 7
    out.push_back((float)chi2_5C);      // 8
    out.push_back(mW_5C_isr);           // 9
    out.push_back(chi2_5C_isr);         // 10
    out.push_back(E_isr_fitted);        // 11
    out.push_back((float)isr_applied);  // 12

    for(int i=0; i<16; i++){
        out.push_back(pull4C[i]);       // 13..28
    }

    return out;
}
''')

# =========================================================
# 2) covariance matrix from all_results.json
# =========================================================

_JSON_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "all_results_ecm365.json")

def build_covariance_matrix(json_path=_JSON_PATH):
    with open(json_path) as f:
        cb = json.load(f)

    jet_vars = [
        ("alpha", 0),
        ("theta", 1),
        ("phi",   2),
        ("x",     3),
    ]

    cov = np.zeros((16, 16))

    for jet_idx in range(4):
        jet_label = jet_idx + 1

        for var, var_idx in jet_vars:
            key = f"delta_{var}_j{jet_label}"
            sigma = cb[key]["sigma"]

            param_idx = 4 * jet_idx + var_idx
            cov[param_idx, param_idx] = sigma ** 2

    return cov

cov = build_covariance_matrix()

cov += np.eye(16) * 1e-8

cov_inv = np.linalg.inv(cov)

ROOT.MyFit.setCovariance16(
    cov.flatten().tolist(),
    cov_inv.flatten().tolist()
)

print("=== Covariance matrix diagonal sigma values ===")
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
    "p8_ee_WW_ecm365": {
        "fraction": 1.0e-6,
        "chunks": 1,
        "output": "wmass_fit_pull4c_ecm365"
    }
}

inputDir = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA/"
procDict = "FCCee_procDict_winter2023_IDEA.json"
outDir = "./outputs/wmass"

ECM = 365.0

nCPUS = 10

doTree = True

includePaths = [
    "headers/selectQuarks.h"
]


MIN_JET_ENERGY_RATIO = 140.0 / 240.0
MAX_JET_PZ_RATIO      =  20.0 / 240.0

min_jet_energy = MIN_JET_ENERGY_RATIO * ECM
max_jet_pz     = MAX_JET_PZ_RATIO * ECM

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

       
            # Exactly 4 reconstructed jets
            #.Define("pass_n_jets_4", "jet_px.size() == 4")
            .Filter("jet_px.size()==4")

         
            #.Define("pass_sumJetE", f"sumJetE > {min_jet_energy}")
            #.Define("pass_sumJetPz", f"sumJetPz < {max_jet_pz}")

            # --- ISR suppression diagnostics ---
            .Define("sumJetE",  "Sum(jet_e)")
            .Define("sumJetPz", "abs(Sum(jet_pz))")
            .Filter(f"sumJetE > {min_jet_energy}")
            .Filter(f"sumJetPz < {max_jet_pz}")

            
            .Define("all_results", f"reconstructWW(jet_px, jet_py, jet_pz, jet_e, {float(ECM)}f)")

            .Define("mW1_raw", "all_results[0]")
            .Define("mW2_raw", "all_results[1]")
            .Define("m_small_raw", "all_results[2]")
            .Define("m_large_raw", "all_results[3]")
            .Define("m_small_4C", "all_results[4]")
            .Define("m_large_4C", "all_results[5]")
            .Define("chi2_4C", "all_results[6]")
            .Define("mW_5C", "all_results[7]")
            .Define("chi2_5C", "all_results[8]")
            .Define("mW_5C_isr", "all_results[9]")
            .Define("chi2_5C_isr", "all_results[10]")
            .Define("E_isr_fitted", "all_results[11]")
            .Define("isr_applied", "all_results[12]")


            # Corrected 4C pulls.
            # Order inside all_results:
            # 11 + 4*jet + variable
            # variable order: alpha, theta, phi, x
            .Define("pull4C_alpha_j1", "all_results[13]")
            .Define("pull4C_theta_j1", "all_results[14]")
            .Define("pull4C_phi_j1",   "all_results[15]")
            .Define("pull4C_x_j1",     "all_results[16]")

            .Define("pull4C_alpha_j2", "all_results[17]")
            .Define("pull4C_theta_j2", "all_results[18]")
            .Define("pull4C_phi_j2",   "all_results[19]")
            .Define("pull4C_x_j2",     "all_results[20]")

            .Define("pull4C_alpha_j3", "all_results[21]")
            .Define("pull4C_theta_j3", "all_results[22]")
            .Define("pull4C_phi_j3",   "all_results[23]")
            .Define("pull4C_x_j3",     "all_results[24]")

            .Define("pull4C_alpha_j4", "all_results[25]")
            .Define("pull4C_theta_j4", "all_results[26]")
            .Define("pull4C_phi_j4",   "all_results[27]")
            .Define("pull4C_x_j4",     "all_results[28]")

            # p-value / fit goodness
            .Define("prob_4C", "ROOT::Math::chisquared_cdf_c(chi2_4C, 4.0)")
            .Define("prob_5C", "ROOT::Math::chisquared_cdf_c(chi2_5C, 5.0)")
            .Define("prob_5C_isr", "ROOT::Math::chisquared_cdf_c(chi2_5C_isr, 5.0)")

            
            # No Event Reduction in 162.5 Gev with this filter.
            #.Define("pass_mW_5C_positive", "mW_5C > 0"
            # .Filter("mW_5C > 0")

            #.Define("pass_prob_4C_3pct", "prob_4C > 0.03f")
            #.Define("pass_prob_5C_3pct", "prob_5C > 0.03f")
            #.Define("pass_prob_5C_isr_3pct", "prob_5C_isr > 0.03f")
            #.Define("pass_prob_final_5C_3pct", "prob_final_5C > 0.03f")

            # This reproduces original ECM-dependent p-value decision.
            #.Define(
            #    "pass_selected_pvalue_cut",
            #    "prob_4C > 0.03f" if ECM < 200.0 else
            #    "((isr_applied > 0.5f && prob_5C_isr > 0.03f) || (isr_applied < 0.5f && prob_5C > 0.03f))"
            #)

            # This reproduces original ECM-dependent p-value decision.
            .Define(
                "pass_selected_pvalue_cut",
                "prob_4C > 0.03f" if ECM < 200.0 else
                "((isr_applied > 0.5f && prob_5C_isr > 0.03f) || (isr_applied < 0.5f && prob_5C > 0.03f))"
            )

            #.Filter(
            #     "prob_4C > 0.03" if ECM < 200.0 else
            #     "(isr_applied > 0.5f && prob_5C_isr > 0.03f) || (isr_applied < 0.5f && prob_5C > 0.03f)"
            # )
        )

        return df

    def output():
        return [
            "m_small_raw", "m_large_raw",
            "m_small_4C", "m_large_4C", "chi2_4C",
            "mW_5C", "chi2_5C",
            "mW_5C_isr", "chi2_5C_isr", "E_isr_fitted", "isr_applied",
            "sumJetE", "sumJetPz",
            "prob_4C", "prob_5C", "prob_5C_isr",

            "pull4C_alpha_j1", "pull4C_theta_j1", "pull4C_phi_j1", "pull4C_x_j1",
            "pull4C_alpha_j2", "pull4C_theta_j2", "pull4C_phi_j2", "pull4C_x_j2",
            "pull4C_alpha_j3", "pull4C_theta_j3", "pull4C_phi_j3", "pull4C_x_j3",
            "pull4C_alpha_j4", "pull4C_theta_j4", "pull4C_phi_j4", "pull4C_x_j4",
        ]


