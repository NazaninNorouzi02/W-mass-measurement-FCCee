#include <TBranch.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct AuditRow {
    std::string section;
    std::string item;
    std::string expression;
    std::string selection;
    Long64_t selectedEvents = 0;
    Long64_t finiteValues = 0;
    Long64_t visibleValues = 0;
    Long64_t excludedValues = 0;
    double rangeMin = std::numeric_limits<double>::quiet_NaN();
    double rangeMax = std::numeric_limits<double>::quiet_NaN();
    double meanVisible = std::numeric_limits<double>::quiet_NaN();
    double stddevVisible = std::numeric_limits<double>::quiet_NaN();
};

std::string csvQuote(const std::string& value) {
    std::string escaped;
    escaped.reserve(value.size() + 2);
    escaped.push_back('"');
    for (char character : value) {
        if (character == '"') escaped.push_back('"');
        escaped.push_back(character);
    }
    escaped.push_back('"');
    return escaped;
}

bool hasBranch(TTree& tree, const std::string& branch) {
    return tree.GetBranch(branch.c_str()) != nullptr;
}

bool requireBranches(TTree& tree, const std::vector<std::string>& branches) {
    bool complete = true;
    for (const std::string& branch : branches) {
        if (!hasBranch(tree, branch)) {
            std::cerr << "ERROR: missing branch: " << branch << '\n';
            complete = false;
        }
    }
    return complete;
}

Long64_t countEntries(TTree& tree, const std::string& selection) {
    const Long64_t count = tree.GetEntries(selection.c_str());
    if (count < 0) {
        std::cerr << "ERROR: ROOT rejected selection:\n  "
                  << selection << '\n';
    }
    return count;
}

AuditRow countRow(
    TTree& tree,
    const std::string& section,
    const std::string& item,
    const std::string& selection
) {
    AuditRow row;
    row.section = section;
    row.item = item;
    row.selection = selection;
    row.selectedEvents = countEntries(tree, selection);
    row.finiteValues = row.selectedEvents;
    row.visibleValues = row.selectedEvents;
    row.excludedValues = 0;
    return row;
}

AuditRow observableRow(
    TTree& tree,
    const std::string& section,
    const std::string& item,
    const std::string& expression,
    const std::string& selection,
    double rangeMin,
    double rangeMax
) {
    AuditRow row;
    row.section = section;
    row.item = item;
    row.expression = expression;
    row.selection = selection;
    row.rangeMin = rangeMin;
    row.rangeMax = rangeMax;
    row.selectedEvents = countEntries(tree, selection);

    const std::string finiteSelection =
        "(" + selection + ") && TMath::Finite((" + expression + "))";
    const std::string visibleSelection =
        "(" + finiteSelection + ") && ((" + expression + ") >= " +
        std::to_string(rangeMin) + ") && ((" + expression + ") < " +
        std::to_string(rangeMax) + ")";

    row.finiteValues = countEntries(tree, finiteSelection);
    const Long64_t drawn = tree.Draw(
        expression.c_str(),
        visibleSelection.c_str(),
        "goff"
    );
    if (drawn < 0) {
        std::cerr << "ERROR: ROOT rejected expression/selection:\n  "
                  << expression << "\n  " << visibleSelection << '\n';
        row.visibleValues = -1;
        row.excludedValues = -1;
        return row;
    }

    row.visibleValues = drawn;
    row.excludedValues = row.selectedEvents - row.visibleValues;
    if (drawn > 0) {
        const double* values = tree.GetV1();
        double sum = 0.0;
        double sum2 = 0.0;
        for (Long64_t index = 0; index < drawn; ++index) {
            sum += values[index];
            sum2 += values[index] * values[index];
        }
        row.meanVisible = sum / static_cast<double>(drawn);
        const double variance =
            std::max(0.0, sum2 / static_cast<double>(drawn) -
                              row.meanVisible * row.meanVisible);
        row.stddevVisible = std::sqrt(variance);
    }
    return row;
}

void printRow(const AuditRow& row) {
    std::cout << std::left << std::setw(20) << row.section
              << " | " << std::setw(43) << row.item
              << " | selected=" << std::right << std::setw(8)
              << row.selectedEvents;
    if (!row.expression.empty()) {
        std::cout << " finite=" << std::setw(8) << row.finiteValues
                  << " visible=" << std::setw(8) << row.visibleValues
                  << " excluded=" << std::setw(8) << row.excludedValues;
    }
    std::cout << '\n';
}

void writeCsv(
    const std::vector<AuditRow>& rows,
    const std::string& csvPath
) {
    std::ofstream output(csvPath.c_str());
    if (!output) {
        std::cerr << "ERROR: cannot create CSV: " << csvPath << '\n';
        return;
    }

    output
        << "section,item,expression,selection,selected_events,"
        << "finite_values,visible_values,excluded_values,range_min,"
        << "range_max,mean_visible,stddev_visible\n";
    output << std::setprecision(12);

    for (const AuditRow& row : rows) {
        output
            << csvQuote(row.section) << ','
            << csvQuote(row.item) << ','
            << csvQuote(row.expression) << ','
            << csvQuote(row.selection) << ','
            << row.selectedEvents << ','
            << row.finiteValues << ','
            << row.visibleValues << ','
            << row.excludedValues << ',';

        if (std::isfinite(row.rangeMin)) output << row.rangeMin;
        output << ',';
        if (std::isfinite(row.rangeMax)) output << row.rangeMax;
        output << ',';
        if (std::isfinite(row.meanVisible)) output << row.meanVisible;
        output << ',';
        if (std::isfinite(row.stddevVisible)) output << row.stddevVisible;
        output << '\n';
    }
}

void appendObservable(
    std::vector<AuditRow>& rows,
    TTree& tree,
    const std::string& section,
    const std::string& item,
    const std::string& expression,
    const std::string& selection,
    double rangeMin,
    double rangeMax
) {
    rows.push_back(
        observableRow(
            tree,
            section,
            item,
            expression,
            selection,
            rangeMin,
            rangeMax
        )
    );
}

}  // namespace

void audit_wmass_365(
    const char* inputPattern =
        "outputs/wmass/wmass_fit_pvalue_pulls_ISR_ecm365.root",
    const char* csvPath = "wmass_365_audit.csv"
) {
    const double probabilityCut = 0.03;
    const double massMin = 40.0;
    const double massMax = 120.0;
    const double pullMin = -5.0;
    const double pullMax = 5.0;

    TChain chain("events");
    const int filesAdded = chain.Add(inputPattern);
    const Long64_t totalEntries = chain.GetEntries();
    if (filesAdded <= 0 || totalEntries <= 0) {
        std::cerr
            << "ERROR: no non-empty 'events' trees matched:\n  "
            << inputPattern << '\n';
        return;
    }
    chain.LoadTree(0);
    chain.SetEstimate(totalEntries + 1);

    std::cout << "============================================================\n";
    std::cout << "W-mass 365 GeV audit\n";
    std::cout << "Input pattern : " << inputPattern << '\n';
    std::cout << "Files added   : " << filesAdded << '\n';
    std::cout << "Chain entries : " << totalEntries << '\n';
    std::cout << "Contributing files:\n";

    TObjArray* fileList = chain.GetListOfFiles();
    for (int index = 0; index < fileList->GetEntries(); ++index) {
        TChainElement* element =
            dynamic_cast<TChainElement*>(fileList->At(index));
        if (!element) continue;
        const char* filePath = element->GetTitle();
        TFile file(filePath, "READ");
        TTree* fileTree =
            file.IsZombie() ? nullptr : dynamic_cast<TTree*>(file.Get("events"));
        std::cout << "  [" << std::setw(2) << index + 1 << "] "
                  << filePath << " : "
                  << (fileTree ? fileTree->GetEntries() : -1)
                  << " entries\n";
    }

    std::vector<std::string> required = {
        "m_small_raw", "m_large_raw",
        "m_small_4C", "m_large_4C",
        "m_small_4C_final", "m_large_4C_final",
        "mW_5C", "mW_5C_final",
        "chi2_4C", "chi2_5C",
        "prob_4C", "prob_5C",
        "prob_4C_isr", "prob_5C_isr",
        "prob_final_4C", "prob_final_5C",
        "standard_4C_valid", "standard_5C_valid",
        "isr_4C_attempted", "isr_5C_attempted",
        "isr_4C_fit_valid", "isr_5C_fit_valid",
        "isr_4C_applied", "isr_applied",
        "isr_4C_recovered", "isr_5C_recovered",
        "E_isr_true_collinear", "E_isr_4C_fitted", "E_isr_fitted",
        "Emax_isr", "y_isr_4C_fitted", "y_isr_5C_fitted",
        "pull4C_standard_valid", "pull4C_final_valid",
        "pull5C_standard_valid", "pull5C_final_valid"
    };
    const std::vector<std::string> variables = {
        "alpha", "theta", "phi", "x"
    };
    for (const std::string& prefix : {
             "pull4C_standard", "pull4C_final",
             "pull5C_standard", "pull5C_final"
         }) {
        for (int jet = 1; jet <= 4; ++jet) {
            for (const std::string& variable : variables) {
                required.push_back(
                    prefix + "_" + variable + "_j" + std::to_string(jet)
                );
            }
        }
    }
    if (!requireBranches(chain, required)) {
        std::cerr
            << "Audit stopped. Use the ROOT file made by the current "
            << "wmass_fit_pulls_ISR_final.py producer.\n";
        return;
    }

    std::vector<AuditRow> rows;
    rows.push_back(countRow(chain, "input", "output tree entries", "1"));

    const std::string valid4 = "standard_4C_valid > 0.5";
    const std::string direct4 =
        valid4 + " && prob_4C > " + std::to_string(probabilityCut);
    const std::string triggered4 =
        valid4 + " && prob_4C < " + std::to_string(probabilityCut);
    const std::string recovered4 =
        triggered4 +
        " && isr_4C_applied > 0.5 && prob_4C_isr > " +
        std::to_string(probabilityCut);
    const std::string final4 =
        valid4 + " && prob_final_4C > " +
        std::to_string(probabilityCut);

    const std::string valid5 = "standard_5C_valid > 0.5";
    const std::string direct5 =
        valid5 + " && prob_5C > " + std::to_string(probabilityCut);
    const std::string triggered5 =
        valid5 + " && prob_5C < " + std::to_string(probabilityCut);
    const std::string recovered5 =
        triggered5 +
        " && isr_applied > 0.5 && prob_5C_isr > " +
        std::to_string(probabilityCut);
    const std::string final5 =
        valid5 + " && prob_final_5C > " +
        std::to_string(probabilityCut);

    for (const AuditRow& row : {
             countRow(chain, "cutflow 4C", "valid standard 4C", valid4),
             countRow(chain, "cutflow 4C", "direct accepted", direct4),
             countRow(chain, "cutflow 4C", "ISR-triggered standard failures", triggered4),
             countRow(chain, "cutflow 4C", "ISR attempted flag", "isr_4C_attempted > 0.5"),
             countRow(chain, "cutflow 4C", "ISR fit valid", "isr_4C_fit_valid > 0.5"),
             countRow(chain, "cutflow 4C", "ISR applied", "isr_4C_applied > 0.5"),
             countRow(chain, "cutflow 4C", "ISR recovered above 3%", recovered4),
             countRow(chain, "cutflow 4C", "ISR recovered flag", "isr_4C_recovered > 0.5"),
             countRow(chain, "cutflow 4C", "final accepted total", final4),
             countRow(chain, "cutflow 5C", "valid standard 5C", valid5),
             countRow(chain, "cutflow 5C", "direct accepted", direct5),
             countRow(chain, "cutflow 5C", "ISR-triggered standard failures", triggered5),
             countRow(chain, "cutflow 5C", "ISR attempted flag", "isr_5C_attempted > 0.5"),
             countRow(chain, "cutflow 5C", "ISR fit valid", "isr_5C_fit_valid > 0.5"),
             countRow(chain, "cutflow 5C", "ISR applied", "isr_applied > 0.5"),
             countRow(chain, "cutflow 5C", "ISR recovered above 3%", recovered5),
             countRow(chain, "cutflow 5C", "ISR recovered flag", "isr_5C_recovered > 0.5"),
             countRow(chain, "cutflow 5C", "final accepted total", final5)
         }) {
        rows.push_back(row);
    }

    rows.push_back(
        countRow(
            chain,
            "fit diagnostic",
            "both standard fits valid",
            valid4 + " && " + valid5
        )
    );
    rows.push_back(
        countRow(
            chain,
            "fit diagnostic",
            "Delta chi2 < -0.001",
            valid4 + " && " + valid5 +
            " && (chi2_5C-chi2_4C) < -0.001"
        )
    );
    rows.push_back(
        countRow(
            chain,
            "fit diagnostic",
            "Delta chi2 > 3.841",
            valid4 + " && " + valid5 +
            " && (chi2_5C-chi2_4C) > 3.841"
        )
    );
    rows.push_back(
        countRow(
            chain,
            "fit diagnostic",
            "Delta chi2 > 6.635",
            valid4 + " && " + valid5 +
            " && (chi2_5C-chi2_4C) > 6.635"
        )
    );
    rows.push_back(
        countRow(
            chain,
            "fit diagnostic",
            "probability outside [0,1]",
            "(prob_4C < 0 || prob_4C > 1 || "
            "prob_5C < 0 || prob_5C > 1 || "
            "prob_final_4C < 0 || prob_final_4C > 1 || "
            "prob_final_5C < 0 || prob_final_5C > 1)"
        )
    );

    appendObservable(
        rows, chain, "Figure 7.8", "MC effective collinear ISR",
        "E_isr_true_collinear", "isr_applied > 0.5", 0.0, 180.0
    );
    appendObservable(
        rows, chain, "Figure 7.8", "fitted 5C ISR",
        "E_isr_fitted", "isr_applied > 0.5", 0.0, 180.0
    );
    rows.push_back(
        countRow(
            chain,
            "Figure 7.8",
            "fitted within 0.5 GeV of Emax",
            "isr_applied > 0.5 && "
            "E_isr_fitted >= Emax_isr - 0.5"
        )
    );
    rows.push_back(
        countRow(
            chain,
            "Figure 7.8",
            "|fitted y_ISR| >= 4.9",
            "isr_applied > 0.5 && abs(y_isr_5C_fitted) >= 4.9"
        )
    );
    appendObservable(
        rows, chain, "Figure 7.8", "configured Emax",
        "Emax_isr", "1", 0.0, 180.0
    );

    appendObservable(
        rows, chain, "Figure 7.9", "5C mass without ISR",
        "mW_5C", valid5, massMin, massMax
    );
    appendObservable(
        rows, chain, "Figure 7.9", "5C mass with ISR treatment",
        "mW_5C_final", valid5, massMin, massMax
    );
    appendObservable(
        rows, chain, "Figure 7.9", "5C probability without ISR",
        "prob_5C", valid5, 0.0, 1.0
    );
    appendObservable(
        rows, chain, "Figure 7.9", "5C probability with ISR treatment",
        "prob_final_5C", valid5, 0.0, 1.0
    );

    appendObservable(
        rows, chain, "Figure 7.10", "raw smaller mass",
        "m_small_raw", "1", massMin, massMax
    );
    appendObservable(
        rows, chain, "Figure 7.10", "raw larger mass",
        "m_large_raw", "1", massMin, massMax
    );
    appendObservable(
        rows, chain, "Figure 7.10", "final 4C smaller mass",
        "m_small_4C_final", valid4, massMin, massMax
    );
    appendObservable(
        rows, chain, "Figure 7.10", "final 4C larger mass",
        "m_large_4C_final", valid4, massMin, massMax
    );
    appendObservable(
        rows, chain, "Figure 7.10", "final 5C mass - smaller panel",
        "mW_5C_final", valid5, massMin, massMax
    );
    appendObservable(
        rows, chain, "Figure 7.10", "final 5C mass - larger panel",
        "mW_5C_final", valid5, massMin, massMax
    );

    for (const std::string& mass : {"m_small_4C", "m_large_4C"}) {
        appendObservable(
            rows, chain, "Appendix 4C", mass + " before cut",
            mass, valid4, massMin, massMax
        );
        appendObservable(
            rows, chain, "Appendix 4C", mass + " after cut",
            mass, direct4, massMin, massMax
        );
    }
    appendObservable(
        rows, chain, "Appendix 5C", "mW_5C before cut",
        "mW_5C", valid5, massMin, massMax
    );
    appendObservable(
        rows, chain, "Appendix 5C", "mW_5C after cut",
        "mW_5C", direct5, massMin, massMax
    );

    struct PullCandidate {
        std::string fit;
        std::string candidate;
        std::string prefix;
        std::string selection;
    };
    const std::vector<PullCandidate> pullCandidates = {
        {"4C", "standard before", "pull4C_standard",
         "pull4C_standard_valid > 0.5"},
        {"4C", "standard after", "pull4C_standard",
         "pull4C_standard_valid > 0.5 && prob_4C > 0.03"},
        {"4C", "final selected", "pull4C_final",
         "pull4C_final_valid > 0.5 && prob_final_4C > 0.03"},
        {"5C", "standard before", "pull5C_standard",
         "pull5C_standard_valid > 0.5"},
        {"5C", "standard after", "pull5C_standard",
         "pull5C_standard_valid > 0.5 && prob_5C > 0.03"},
        {"5C", "final selected", "pull5C_final",
         "pull5C_final_valid > 0.5 && prob_final_5C > 0.03"}
    };
    for (const PullCandidate& candidate : pullCandidates) {
        for (int jet = 1; jet <= 4; ++jet) {
            for (const std::string& variable : variables) {
                const std::string branch =
                    candidate.prefix + "_" + variable + "_j" +
                    std::to_string(jet);
                appendObservable(
                    rows,
                    chain,
                    "Pull " + candidate.fit,
                    candidate.candidate + " jet" + std::to_string(jet) +
                        " " + variable,
                    branch,
                    candidate.selection,
                    pullMin,
                    pullMax
                );
            }
        }
    }

    std::cout << "============================================================\n";
    for (const AuditRow& row : rows) printRow(row);
    std::cout << "============================================================\n";
    writeCsv(rows, csvPath);
    std::cout << "Saved audit CSV: " << csvPath << '\n';
    std::cout << "Mass plot range: [40,120) GeV\n";
    std::cout << "Pull plot range: [-5,5)\n";
    std::cout << "Probability range: [0,1]\n";
    std::cout << "============================================================\n";
}