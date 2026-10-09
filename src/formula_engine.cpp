/**
 * HORIZON: Heavy-first Ordered Recursive Inference with Zero-loop Optimal Navigation
 * High-Performance Mass Decomposition and Isotope Modeling Engine for enviGCMS
 *
 * Algorithm Design & Key Principles:
 * 1. Heavy-element Priority Ordering: Sort search elements by descending exact mass
 *    (e.g., S > P > O > N > C > H) to minimize tree branching factor near the root.
 * 2. Mass Horizon Pruning: Precompute cumulative min/max reachable masses (horizons)
 *    to aggressively prune infeasible subtrees in O(1) at every recursion depth.
 * 3. Analytical Leaf Solving: Directly solve the lightest element (Hydrogen) via
 *    O(1) floating-point division and rounding, eliminating inner loops.
 * 4. Multi-tier Chemical Rule Filtering: Integrates Senior's valence rules (Senior 1951),
 *    Fiehn Seven Golden Rules (Kind & Fiehn 2007), DBE, and Nitrogen parity rule.
 * 5. High-Resolution Fine Isotope Modeling: Binary exponentiation convolution with
 *    instrument-aware centroid peak merging.
 * 6. High-Throughput Parallelism: Native OpenMP multithreading for batch mass sets.
 */

#include <Rcpp.h>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <cctype>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

namespace envi_iso {

// Electron rest mass in Da
const double ELECTRON_MASS = 0.000548579909;

struct Isotope {
    int offset;       // nominal mass offset from monoisotopic peak (0, 1, 2, ...)
    double mass;      // exact mass of this isotope
    double abundance; // fractional abundance (sums to 1 across isotopes of element)
};

struct ElementInfo {
    std::string symbol;
    int nominal_mass;
    double exact_mass;
    int valence;
    std::vector<Isotope> isotopes;
};

// Periodic table with IUPAC atomic masses, isotopic abundances, and valences
static const std::map<std::string, ElementInfo>& get_element_table() {
    static std::map<std::string, ElementInfo> table;
    if (table.empty()) {
        table["H"]  = {"H",  1,   1.007825032, 1, {{0, 1.007825032, 0.999885}, {1, 2.014101778, 0.000115}}};
        table["D"]  = {"D",  2,   2.014101778, 1, {{0, 2.014101778, 1.0}}};
        table["He"] = {"He", 4,   4.002603254, 0, {{0, 4.002603254, 1.0}}};
        table["Li"] = {"Li", 7,   7.016004550, 1, {{-1, 6.015122795, 0.0759}, {0, 7.016004550, 0.9241}}};
        table["Be"] = {"Be", 9,   9.012182200, 2, {{0, 9.012182200, 1.0}}};
        table["B"]  = {"B",  11, 11.009305400, 3, {{-1, 10.0129370, 0.199}, {0, 11.009305400, 0.801}}};
        table["C"]  = {"C",  12, 12.000000000, 4, {{0, 12.000000000, 0.9893}, {1, 13.003354838, 0.0107}}};
        table["N"]  = {"N",  14, 14.003074004, 3, {{0, 14.003074004, 0.99632}, {1, 15.000108898, 0.00368}}};
        table["O"]  = {"O",  16, 15.994914620, 2, {{0, 15.994914620, 0.99757}, {1, 16.999131756, 0.00038}, {2, 17.9991610, 0.00205}}};
        table["F"]  = {"F",  19, 18.998403200, 1, {{0, 18.998403200, 1.0}}};
        table["Ne"] = {"Ne", 20, 19.992440175, 0, {{0, 19.992440175, 0.9048}, {1, 20.99384668, 0.0027}, {2, 21.99138511, 0.0925}}};
        table["Na"] = {"Na", 23, 22.989769282, 1, {{0, 22.989769282, 1.0}}};
        table["Mg"] = {"Mg", 24, 23.985041700, 2, {{0, 23.985041700, 0.7899}, {1, 24.98583692, 0.1000}, {2, 25.98259293, 0.1101}}};
        table["Al"] = {"Al", 27, 26.981538630, 3, {{0, 26.981538630, 1.0}}};
        table["Si"] = {"Si", 28, 27.976926533, 4, {{0, 27.976926533, 0.92223}, {1, 28.97649470, 0.04685}, {2, 29.9737702, 0.03092}}};
        table["P"]  = {"P",  31, 30.973761998, 5, {{0, 30.973761998, 1.0}}};
        table["S"]  = {"S",  32, 31.972071000, 6, {{0, 31.972071000, 0.9499}, {1, 32.97145876, 0.0075}, {2, 33.96786690, 0.0425}, {4, 35.96708076, 0.0001}}};
        table["Cl"] = {"Cl", 35, 34.968852680, 1, {{0, 34.968852680, 0.7578}, {2, 36.96590260, 0.2422}}};
        table["Ar"] = {"Ar", 40, 39.962383122, 0, {{-4, 35.9675451, 0.003365}, {-2, 37.9627324, 0.000632}, {0, 39.962383122, 0.996003}}};
        table["K"]  = {"K",  39, 38.963706680, 1, {{0, 38.963706680, 0.932581}, {1, 39.9639984, 0.000117}, {2, 40.96182576, 0.067302}}};
        table["Ca"] = {"Ca", 40, 39.962590980, 2, {{0, 39.962590980, 0.96941}, {2, 41.9586180, 0.00647}, {3, 42.9587666, 0.00135}, {4, 43.9554818, 0.02086}}};
        table["Sc"] = {"Sc", 45, 44.955911900, 3, {{0, 44.955911900, 1.0}}};
        table["Ti"] = {"Ti", 48, 47.947946300, 4, {{-2, 45.9526316, 0.0825}, {-1, 46.9517631, 0.0744}, {0, 47.9479463, 0.7372}, {1, 48.9478700, 0.0541}, {2, 49.9447912, 0.0518}}};
        table["V"]  = {"V",  51, 50.943959500, 5, {{-1, 49.9471585, 0.0025}, {0, 50.9439595, 0.9975}}};
        table["Cr"] = {"Cr", 52, 51.940507500, 3, {{-2, 49.9460442, 0.04345}, {0, 51.9405075, 0.83789}, {1, 52.9406494, 0.09501}, {2, 53.9388804, 0.02365}}};
        table["Mn"] = {"Mn", 55, 54.938045100, 2, {{0, 54.938045100, 1.0}}};
        table["Fe"] = {"Fe", 56, 55.934937500, 3, {{-2, 53.9396105, 0.05845}, {0, 55.9349375, 0.91754}, {1, 56.9353940, 0.02119}, {2, 57.9332762, 0.00282}}};
        table["Co"] = {"Co", 59, 58.933195000, 2, {{0, 58.933195000, 1.0}}};
        table["Ni"] = {"Ni", 58, 57.935342900, 2, {{0, 57.9353429, 0.68077}, {2, 59.9307864, 0.26223}, {3, 60.9310560, 0.01140}, {4, 61.9283451, 0.03634}, {6, 63.9279660, 0.00926}}};
        table["Cu"] = {"Cu", 63, 62.929601100, 2, {{0, 62.9296011, 0.6917}, {2, 64.9277895, 0.3083}}};
        table["Zn"] = {"Zn", 64, 63.929142200, 2, {{0, 63.9291422, 0.4917}, {2, 65.9260334, 0.2773}, {3, 66.9271273, 0.0404}, {4, 67.9248442, 0.1845}, {6, 69.9253193, 0.0061}}};
        table["Ga"] = {"Ga", 69, 68.925573600, 3, {{0, 68.9255736, 0.60108}, {2, 70.9247013, 0.39892}}};
        table["Ge"] = {"Ge", 74, 73.921177800, 4, {{-4, 69.9242474, 0.2084}, {-2, 71.9220758, 0.2754}, {-1, 72.9234589, 0.0773}, {0, 73.9211778, 0.3628}, {2, 75.9214026, 0.0761}}};
        table["As"] = {"As", 75, 74.921596500, 3, {{0, 74.921596500, 1.0}}};
        table["Se"] = {"Se", 80, 79.916521300, 2, {{-6, 73.9224764, 0.0089}, {-4, 75.9192136, 0.0937}, {-3, 76.9199140, 0.0763}, {-2, 77.9173091, 0.2377}, {0, 79.9165213, 0.4961}, {2, 81.9166994, 0.0873}}};
        table["Br"] = {"Br", 79, 78.918337100, 1, {{0, 78.9183371, 0.5069}, {2, 80.9162906, 0.4931}}};
        table["Kr"] = {"Kr", 84, 83.911507000, 0, {{-6, 77.9203648, 0.0035}, {-4, 79.9163790, 0.0228}, {-2, 81.9134836, 0.1158}, {-1, 82.9141360, 0.1149}, {0, 83.911507, 0.5700}, {2, 85.9106107, 0.1730}}};
        table["Rb"] = {"Rb", 85, 84.911789738, 1, {{0, 84.911789738, 0.7217}, {2, 86.909180527, 0.2783}}};
        table["Sr"] = {"Sr", 88, 87.905612100, 2, {{-4, 83.913425, 0.0056}, {-2, 85.9092602, 0.0986}, {-1, 86.9088771, 0.0700}, {0, 87.9056121, 0.8258}}};
        table["Mo"] = {"Mo", 98, 97.905408200, 6, {{-6, 91.906811, 0.1484}, {-4, 93.9050883, 0.0925}, {-3, 94.9058421, 0.1592}, {-2, 95.9046795, 0.1668}, {-1, 96.9060216, 0.0955}, {0, 97.9054082, 0.2413}, {2, 99.907477, 0.0963}}};
        table["Ag"] = {"Ag", 107, 106.905097, 1, {{0, 106.905097, 0.51839}, {2, 108.904752, 0.48161}}};
        table["Cd"] = {"Cd", 114, 113.9033585, 2, {{-8, 105.906459, 0.0125}, {-6, 107.904184, 0.0089}, {-4, 109.903002, 0.1249}, {-3, 110.904178, 0.1280}, {-2, 111.9027578, 0.2413}, {-1, 112.9044017, 0.1222}, {0, 113.9033585, 0.2873}, {2, 115.904756, 0.0749}}};
        table["In"] = {"In", 115, 114.903878, 3, {{-2, 112.904058, 0.0429}, {0, 114.903878, 0.9571}}};
        table["Sn"] = {"Sn", 120, 119.9021947, 4, {{-8, 111.904818, 0.0097}, {-6, 113.902779, 0.0066}, {-5, 114.903342, 0.0034}, {-4, 115.901741, 0.1454}, {-3, 116.902952, 0.0768}, {-2, 117.901603, 0.2422}, {-1, 118.903308, 0.0859}, {0, 119.9021947, 0.3258}, {2, 121.904699, 0.0463}, {4, 123.9052739, 0.0579}}};
        table["Sb"] = {"Sb", 121, 120.9038157, 3, {{0, 120.9038157, 0.5721}, {2, 122.904214, 0.4279}}};
        table["Te"] = {"Te", 130, 129.9062244, 2, {{-10, 119.90402, 0.0009}, {-8, 121.9030439, 0.0255}, {-7, 122.9042700, 0.0089}, {-6, 123.9028179, 0.0474}, {-5, 124.9044307, 0.0707}, {-4, 125.9033117, 0.1884}, {-2, 127.9044631, 0.3174}, {0, 129.9062244, 0.3408}}};
        table["I"]  = {"I",  127, 126.904473, 1, {{0, 126.904473, 1.0}}};
        table["Cs"] = {"Cs", 133, 132.905451933, 1, {{0, 132.905451933, 1.0}}};
        table["Ba"] = {"Ba", 138, 137.9052472, 2, {{-8, 129.9063208, 0.00106}, {-6, 131.9050613, 0.00101}, {-4, 133.9045084, 0.02417}, {-3, 134.9056886, 0.06592}, {-2, 135.9045759, 0.07854}, {-1, 136.9058274, 0.1123}, {0, 137.9052472, 0.71698}}};
        table["Pt"] = {"Pt", 195, 194.9647911, 4, {{-5, 189.959932, 0.00012}, {-3, 191.9610380, 0.00782}, {-1, 193.9626803, 0.32967}, {0, 194.9647911, 0.33832}, {1, 195.9649514, 0.25242}, {3, 197.967893, 0.07163}}};
        table["Au"] = {"Au", 197, 196.9665687, 1, {{0, 196.9665687, 1.0}}};
        table["Hg"] = {"Hg", 202, 201.9706430, 2, {{-6, 195.965833, 0.0015}, {-4, 197.9667690, 0.0997}, {-3, 198.9682799, 0.1687}, {-2, 199.9683260, 0.2310}, {-1, 200.9703023, 0.1318}, {0, 201.9706430, 0.2986}, {2, 203.9734939, 0.0687}}};
        table["Pb"] = {"Pb", 208, 207.9766521, 2, {{-4, 203.9730436, 0.014}, {-2, 205.9744653, 0.241}, {-1, 206.9758969, 0.221}, {0, 207.9766521, 0.524}}};
        table["Bi"] = {"Bi", 209, 208.9803987, 3, {{0, 208.9803987, 1.0}}};
    }
    return table;
}

// Parse chemical formula into element counts
std::map<std::string, int> parse_formula(const std::string& formula, int& out_charge) {
    std::map<std::string, int> counts;
    out_charge = 0;
    size_t n = formula.size();
    if (n == 0) return counts;

    // Check trailing charge like "+", "-", "+1", "-1", "+2", etc.
    size_t end_idx = n;
    if (formula.back() == '+' || formula.back() == '-') {
        char sgn = formula.back();
        end_idx--;
        int mult = 1;
        out_charge = (sgn == '+') ? mult : -mult;
    } else {
        // check if trailing digit preceded by + or -
        size_t p = n - 1;
        while (p > 0 && std::isdigit(formula[p])) p--;
        if (p < n - 1 && (formula[p] == '+' || formula[p] == '-')) {
            int val = std::stoi(formula.substr(p + 1));
            out_charge = (formula[p] == '+') ? val : -val;
            end_idx = p;
        }
    }

    size_t i = 0;
    while (i < end_idx) {
        if (!std::isupper(formula[i])) {
            i++;
            continue;
        }
        std::string elem(1, formula[i]);
        i++;
        if (i < end_idx && std::islower(formula[i])) {
            elem += formula[i];
            i++;
        }
        int count = 0;
        bool has_digits = false;
        while (i < end_idx && std::isdigit(formula[i])) {
            has_digits = true;
            count = count * 10 + (formula[i] - '0');
            i++;
        }
        if (!has_digits) {
            count = 1;
        }
        counts[elem] += count;
    }
    return counts;
}

// Format map to Hill system formula
std::string format_formula(const std::map<std::string, int>& counts) {
    std::stringstream ss;
    bool has_c = (counts.count("C") && counts.at("C") > 0);

    if (has_c) {
        ss << "C";
        if (counts.at("C") > 1) ss << counts.at("C");
        if (counts.count("H") && counts.at("H") > 0) {
            ss << "H";
            if (counts.at("H") > 1) ss << counts.at("H");
        }
        for (const auto& pair : counts) {
            if (pair.first == "C" || pair.first == "H") continue;
            if (pair.second > 0) {
                ss << pair.first;
                if (pair.second > 1) ss << pair.second;
            }
        }
    } else {
        // Strict alphabetical order when no Carbon
        for (const auto& pair : counts) {
            if (pair.second > 0) {
                ss << pair.first;
                if (pair.second > 1) ss << pair.second;
            }
        }
    }
    return ss.str();
}

// Monoisotopic exact mass
double calc_exact_mass(const std::map<std::string, int>& counts, int z = 0) {
    double total = 0.0;
    const auto& table = get_element_table();
    for (const auto& pair : counts) {
        auto it = table.find(pair.first);
        if (it != table.end()) {
            total += it->second.exact_mass * pair.second;
        }
    }
    if (z != 0) {
        total = (total - z * ELECTRON_MASS) / std::abs(z);
    }
    return total;
}

// Double Bond Equivalent (Degree of Unsaturation)
double calc_dbe(const std::map<std::string, int>& counts) {
    int c = counts.count("C") ? counts.at("C") : 0;
    int si = counts.count("Si") ? counts.at("Si") : 0;
    int h = counts.count("H") ? counts.at("H") : 0;
    int f = counts.count("F") ? counts.at("F") : 0;
    int cl = counts.count("Cl") ? counts.at("Cl") : 0;
    int br = counts.count("Br") ? counts.at("Br") : 0;
    int i_cnt = counts.count("I") ? counts.at("I") : 0;
    int n = counts.count("N") ? counts.at("N") : 0;
    int p = counts.count("P") ? counts.at("P") : 0;
    return 1.0 + c + si - 0.5 * (h + f + cl + br + i_cnt) + 0.5 * (n + p);
}

// Parity & Nitrogen rule check
char calc_parity(double nominal_mass, int n_count, int charge) {
    bool masseven = (static_cast<long long>(std::round(nominal_mass)) % 2 == 0);
    bool nitrogeneven = (n_count % 2 == 0);
    bool chargeeven = (std::abs(charge) % 2 == 0);
    return (masseven ^ nitrogeneven ^ chargeeven) ? 'e' : 'o';
}

bool is_valid_nitrogen_rule(double nominal_mass, int n_count, int charge) {
    bool massodd = (static_cast<long long>(std::round(nominal_mass)) % 2 != 0);
    bool masseven = !massodd;
    bool nitrogenodd = (n_count % 2 != 0);
    bool nitrogeneven = !nitrogenodd;
    bool zodd = (std::abs(charge) % 2 != 0);
    bool zeven = !zodd;

    return (zeven && masseven && nitrogeneven) ||
           (zeven && massodd && nitrogenodd)   ||
           (zodd && masseven && nitrogenodd)   ||
           (zodd && massodd && nitrogeneven);
}

// ----------------------------------------------------------------------------
// Fiehn "Seven Golden Rules" (Kind & Fiehn, 2007) implementation
// ----------------------------------------------------------------------------

// Senior's rules check (Senior 1951)
bool check_senior_rules(const std::map<std::string, int>& counts, int charge) {
    if (counts.empty()) return false;
    int sum_valence = 0;
    int max_valence = 0;
    int total_atoms = 0;

    static const std::map<std::string, int> valences = {
        {"H", 1}, {"D", 1}, {"C", 4}, {"N", 3}, {"O", 2}, {"F", 1},
        {"P", 5}, {"S", 6}, {"Cl", 1}, {"Br", 1}, {"I", 1}, {"Si", 4},
        {"Na", 1}, {"K", 1}, {"Fe", 3}, {"B", 3}
    };

    for (const auto& pair : counts) {
        if (pair.second <= 0) continue;
        int v = 1;
        auto it = valences.find(pair.first);
        if (it != valences.end()) {
            v = it->second;
        }
        sum_valence += v * pair.second;
        if (v > max_valence) max_valence = v;
        total_atoms += pair.second;
    }

    if (total_atoms <= 1) return true;

    // Senior Rule 1: Sum of valences >= 2 * max_valence
    if (sum_valence < 2 * max_valence) return false;

    // Senior Rule 2: Parity of total valence accounting for charge
    int effective_val = sum_valence - std::abs(charge);
    if (effective_val % 2 != 0) return false;

    // Senior Rule 3: Number of atoms >= max_valence + 1
    if (total_atoms < max_valence + 1) return false;

    return true;
}

// Hydrogen to Carbon ratio check
bool check_hc_ratio(int h_count, int c_count, double min_hc = 0.1, double max_hc = 6.0) {
    if (c_count <= 0) return true;
    double hc = static_cast<double>(h_count) / c_count;
    return (hc >= min_hc && hc <= max_hc);
}

// Heteroatom to Carbon ratio checks (NOPS, halogens)
bool check_heteroatom_ratios(const std::map<std::string, int>& counts,
                             double max_nc = 4.0, double max_oc = 3.0,
                             double max_pc = 2.0, double max_sc = 3.0,
                             double max_haloc = 4.0) {
    int c = counts.count("C") ? counts.at("C") : 0;
    if (c <= 0) return true;

    int n = counts.count("N") ? counts.at("N") : 0;
    int o = counts.count("O") ? counts.at("O") : 0;
    int p = counts.count("P") ? counts.at("P") : 0;
    int s = counts.count("S") ? counts.at("S") : 0;
    int halo = (counts.count("F") ? counts.at("F") : 0) +
               (counts.count("Cl") ? counts.at("Cl") : 0) +
               (counts.count("Br") ? counts.at("Br") : 0) +
               (counts.count("I") ? counts.at("I") : 0);

    if (static_cast<double>(n) / c > max_nc) return false;
    if (static_cast<double>(o) / c > max_oc) return false;
    if (static_cast<double>(p) / c > max_pc) return false;
    if (static_cast<double>(s) / c > max_sc) return false;
    if (static_cast<double>(halo) / c > max_haloc) return false;

    return true;
}

// Complete Seven Golden Rules check
bool check_seven_golden_rules(
    const std::map<std::string, int>& counts,
    double dbe,
    int charge,
    bool senior_check = true,
    bool hc_check = true,
    bool hetero_check = true,
    double min_hc = 0.1,
    double max_hc = 6.0,
    double max_dbe = 40.0
) {
    if (dbe < -0.5 || dbe > max_dbe) return false;
    if (senior_check && !check_senior_rules(counts, charge)) return false;
    if (hc_check) {
        int h = counts.count("H") ? counts.at("H") : 0;
        int c = counts.count("C") ? counts.at("C") : 0;
        if (!check_hc_ratio(h, c, min_hc, max_hc)) return false;
    }
    if (hetero_check && !check_heteroatom_ratios(counts)) return false;
    return true;
}

// ----------------------------------------------------------------------------
// Isotopic distribution and Instrument Resolution peak merging
// ----------------------------------------------------------------------------

struct Distribution {
    std::vector<double> masses;
    std::vector<double> abundances;
};

Distribution convolve_distributions(const Distribution& d1, const Distribution& d2, size_t max_size = 25) {
    if (d1.abundances.empty()) return d2;
    if (d2.abundances.empty()) return d1;

    size_t sz1 = d1.abundances.size();
    size_t sz2 = d2.abundances.size();
    size_t out_sz = std::min(sz1 + sz2 - 1, max_size);

    Distribution out;
    out.masses.assign(out_sz, 0.0);
    out.abundances.assign(out_sz, 0.0);

    for (size_t i = 0; i < sz1; ++i) {
        if (d1.abundances[i] <= 0.0) continue;
        for (size_t j = 0; j < sz2; ++j) {
            size_t k = i + j;
            if (k >= out_sz) continue;
            double prob = d1.abundances[i] * d2.abundances[j];
            if (prob <= 0.0) continue;
            double m = d1.masses[i] + d2.masses[j];
            out.masses[k] += prob * m;
            out.abundances[k] += prob;
        }
    }

    for (size_t k = 0; k < out_sz; ++k) {
        if (out.abundances[k] > 0.0) {
            out.masses[k] /= out.abundances[k];
        }
    }
    return out;
}

Distribution apply_resolution_filter(const Distribution& raw, double resolution = 0.0, double fwhm = 0.0) {
    if (raw.masses.size() <= 1) return raw;

    double effective_fwhm = fwhm;
    if (effective_fwhm <= 0.0 && resolution > 0.0 && !raw.masses.empty()) {
        effective_fwhm = raw.masses[0] / resolution;
    }
    if (effective_fwhm <= 0.0) return raw;

    Distribution merged;
    size_t n = raw.masses.size();
    size_t i = 0;
    while (i < n) {
        double w_mass = raw.masses[i] * raw.abundances[i];
        double sum_ab = raw.abundances[i];
        size_t j = i + 1;
        while (j < n && (raw.masses[j] - raw.masses[j - 1]) <= effective_fwhm) {
            w_mass += raw.masses[j] * raw.abundances[j];
            sum_ab += raw.abundances[j];
            j++;
        }
        double centroid = (sum_ab > 0.0) ? (w_mass / sum_ab) : raw.masses[i];
        merged.masses.push_back(centroid);
        merged.abundances.push_back(sum_ab);
        i = j;
    }
    return merged;
}

Distribution calc_isotope_distribution(
    const std::map<std::string, int>& counts,
    int z = 0,
    size_t max_isotopes = 20,
    double resolution = 0.0,
    double fwhm = 0.0
) {
    const auto& table = get_element_table();
    Distribution current;
    current.masses = {0.0};
    current.abundances = {1.0};

    for (const auto& pair : counts) {
        const std::string& elem = pair.first;
        int count = pair.second;
        if (count <= 0) continue;

        auto it = table.find(elem);
        if (it == table.end()) continue;

        const auto& iso_list = it->second.isotopes;
        int min_off = 0, max_off = 0;
        for (const auto& iso : iso_list) {
            if (iso.offset < min_off) min_off = iso.offset;
            if (iso.offset > max_off) max_off = iso.offset;
        }
        int shift = -min_off;
        size_t elem_sz = max_off + shift + 1;

        Distribution elem_dist;
        elem_dist.masses.assign(elem_sz, 0.0);
        elem_dist.abundances.assign(elem_sz, 0.0);
        for (const auto& iso : iso_list) {
            size_t pos = iso.offset + shift;
            if (pos < elem_sz) {
                elem_dist.masses[pos] = iso.mass;
                elem_dist.abundances[pos] = iso.abundance;
            }
        }

        Distribution power_dist = elem_dist;
        int p = count;
        while (p > 0) {
            if (p & 1) {
                current = convolve_distributions(current, power_dist, max_isotopes);
            }
            p >>= 1;
            if (p > 0) {
                power_dist = convolve_distributions(power_dist, power_dist, max_isotopes);
            }
        }
    }

    if (z != 0) {
        for (size_t i = 0; i < current.masses.size(); ++i) {
            current.masses[i] = (current.masses[i] - z * ELECTRON_MASS) / std::abs(z);
        }
    }

    if (resolution > 0.0 || fwhm > 0.0) {
        current = apply_resolution_filter(current, resolution, fwhm);
    }

    double sum = 0.0;
    for (double a : current.abundances) sum += a;
    if (sum > 0.0) {
        for (double& a : current.abundances) a /= sum;
    }

    return current;
}

// ----------------------------------------------------------------------------
// Decomposition engine
// ----------------------------------------------------------------------------

struct DecompElement {
    std::string symbol;
    double mass;
    int min_count;
    int max_count;
};

struct DecompResult {
    std::string formula;
    double exactmass;
    double score;
    std::string valid;
    std::string parity;
    double dbe;
    double hc_ratio;
    bool senior_rules;
    std::map<std::string, int> counts;
};

void dfs_decompose(
    size_t idx,
    double rem_mass,
    const std::vector<DecompElement>& elems,
    std::vector<int>& cur_counts,
    double tol,
    std::vector<DecompResult>& results,
    const std::vector<double>& min_remain,
    const std::vector<double>& max_remain,
    int z,
    bool golden_rules = true
) {
    if (idx == elems.size() - 1) {
        const auto& el = elems[idx];
        int count = static_cast<int>(std::round(rem_mass / el.mass));
        if (count >= el.min_count && count <= el.max_count) {
            double final_rem = rem_mass - count * el.mass;
            if (std::abs(final_rem) <= tol) {
                cur_counts[idx] = count;
                std::map<std::string, int> f_map;
                double nominal_mass = 0.0;
                int n_count = 0;
                const auto& table = get_element_table();
                for (size_t i = 0; i < elems.size(); ++i) {
                    if (cur_counts[i] > 0) {
                        f_map[elems[i].symbol] = cur_counts[i];
                    }
                    auto it = table.find(elems[i].symbol);
                    if (it != table.end()) {
                        nominal_mass += it->second.nominal_mass * cur_counts[i];
                    }
                    if (elems[i].symbol == "N") n_count = cur_counts[i];
                }

                double exact = calc_exact_mass(f_map, z);
                double dbe = calc_dbe(f_map);
                bool valid_n = is_valid_nitrogen_rule(nominal_mass, n_count, z);
                char par = calc_parity(nominal_mass, n_count, z);

                int h_cnt = f_map.count("H") ? f_map["H"] : 0;
                int c_cnt = f_map.count("C") ? f_map["C"] : 0;
                double hc = (c_cnt > 0) ? (static_cast<double>(h_cnt) / c_cnt) : 0.0;
                bool senior_pass = check_senior_rules(f_map, z);

                bool is_valid = false;
                if (golden_rules) {
                    is_valid = valid_n && check_seven_golden_rules(f_map, dbe, z);
                } else {
                    is_valid = valid_n && (dbe >= -0.5);
                }

                DecompResult res;
                res.formula = format_formula(f_map);
                res.exactmass = exact;
                res.score = 1.0;
                res.dbe = dbe;
                res.hc_ratio = hc;
                res.senior_rules = senior_pass;
                res.parity = std::string(1, par);
                res.valid = is_valid ? "Valid" : "Invalid";
                res.counts = f_map;
                results.push_back(res);
            }
        }
        return;
    }

    const auto& el = elems[idx];
    int start = el.max_count;
    int end = el.min_count;

    for (int c = start; c >= end; --c) {
        double new_rem = rem_mass - c * el.mass;
        if (new_rem < min_remain[idx + 1] - tol) continue;
        if (new_rem > max_remain[idx + 1] + tol) break;

        cur_counts[idx] = c;
        dfs_decompose(idx + 1, new_rem, elems, cur_counts, tol, results, min_remain, max_remain, z, golden_rules);
    }
}

// Convert native C++ decomposition results to Rcpp::List
List format_decomp_result_list(
    const std::vector<DecompResult>& results,
    int z,
    int maxisotopes,
    double resolution = 0.0,
    double fwhm = 0.0
) {
    size_t nr = results.size();
    CharacterVector formulas(nr);
    NumericVector scores(nr);
    NumericVector exactmasses(nr);
    IntegerVector charges(nr, z);
    CharacterVector parities(nr);
    CharacterVector valids(nr);
    NumericVector dbes(nr);
    NumericVector hc_ratios(nr);
    LogicalVector senior_flags(nr);
    List isotopes(nr);

    for (size_t i = 0; i < nr; ++i) {
        formulas[i] = results[i].formula;
        scores[i] = results[i].score;
        exactmasses[i] = results[i].exactmass;
        parities[i] = results[i].parity;
        valids[i] = results[i].valid;
        dbes[i] = results[i].dbe;
        hc_ratios[i] = results[i].hc_ratio;
        senior_flags[i] = results[i].senior_rules;

        auto dist = calc_isotope_distribution(results[i].counts, z, maxisotopes, resolution, fwhm);
        size_t k = std::min((size_t)maxisotopes, dist.masses.size());
        NumericMatrix iso_mat(2, k);
        for (size_t j = 0; j < k; ++j) {
            iso_mat(0, j) = dist.masses[j];
            iso_mat(1, j) = dist.abundances[j];
        }
        isotopes[i] = iso_mat;
    }

    return List::create(
        Named("formula") = formulas,
        Named("score") = scores,
        Named("exactmass") = exactmasses,
        Named("charge") = charges,
        Named("parity") = parities,
        Named("valid") = valids,
        Named("DBE") = dbes,
        Named("isotopes") = isotopes,
        Named("hc_ratio") = hc_ratios,
        Named("senior_rules") = senior_flags
    );
}

} // namespace envi_iso

// ----------------------------------------------------------------------------
// Rcpp Exported Functions
// ----------------------------------------------------------------------------

//' @title Calculate Exact Mass and Isotope Distribution using the HORIZON Engine
//' @description Computes exact monoisotopic mass, DBE, and theoretical isotopic pattern using the HORIZON engine.
//' @param formula Character string of the molecular formula (e.g. "C6H12O6", "H2O").
//' @param z Integer charge state (default: 0).
//' @param maxisotopes Maximum number of isotopic peaks to return (default: 10).
//' @param resolution Mass resolution (m / FWHM) for merging fine isotope peaks (default: 0).
//' @param fwhm Full width at half maximum in Da for isotope peak merging (default: 0).
//' @return A list compatible with Rdisop::getMolecule.
//' @keywords internal
// [[Rcpp::export]]
List rcpp_get_molecule(
    std::string formula,
    int z = 0,
    int maxisotopes = 10,
    double resolution = 0.0,
    double fwhm = 0.0
) {
    int parsed_z = 0;
    auto counts = envi_iso::parse_formula(formula, parsed_z);
    int final_z = (z != 0) ? z : parsed_z;

    double exact = envi_iso::calc_exact_mass(counts, final_z);
    double dbe = envi_iso::calc_dbe(counts);

    double nominal_mass = 0.0;
    int n_count = counts.count("N") ? counts["N"] : 0;
    const auto& table = envi_iso::get_element_table();
    for (const auto& pair : counts) {
        auto it = table.find(pair.first);
        if (it != table.end()) nominal_mass += it->second.nominal_mass * pair.second;
    }

    char par = envi_iso::calc_parity(nominal_mass, n_count, final_z);
    bool valid_n = envi_iso::is_valid_nitrogen_rule(nominal_mass, n_count, final_z);
    bool golden_valid = valid_n && envi_iso::check_seven_golden_rules(counts, dbe, final_z);

    auto dist = envi_iso::calc_isotope_distribution(counts, final_z, maxisotopes, resolution, fwhm);
    size_t k = std::min((size_t)maxisotopes, dist.masses.size());
    NumericMatrix iso_mat(2, k);
    for (size_t i = 0; i < k; ++i) {
        iso_mat(0, i) = dist.masses[i];
        iso_mat(1, i) = dist.abundances[i];
    }

    return List::create(
        Named("formula") = envi_iso::format_formula(counts),
        Named("score") = 1.0,
        Named("exactmass") = exact,
        Named("charge") = final_z,
        Named("parity") = std::string(1, par),
        Named("valid") = golden_valid ? "Valid" : "Invalid",
        Named("DBE") = dbe,
        Named("isotopes") = List::create(iso_mat)
    );
}

//' @title Check Formula against Fiehn Seven Golden Rules
//' @description Tests molecular formula against Senior's valences, H/C ratio, heteroatom ratios, and DBE.
//' @param formula Character molecular formula string.
//' @param z Charge state (default: 0).
//' @param min_hc Minimum H/C ratio (default: 0.1).
//' @param max_hc Maximum H/C ratio (default: 6.0).
//' @param max_dbe Maximum DBE (default: 40.0).
//' @return A list with rule pass/fail booleans and summary details.
//' @keywords internal
// [[Rcpp::export]]
List rcpp_check_golden_rules(
    std::string formula,
    int z = 0,
    double min_hc = 0.1,
    double max_hc = 6.0,
    double max_dbe = 40.0
) {
    int parsed_z = 0;
    auto counts = envi_iso::parse_formula(formula, parsed_z);
    int final_z = (z != 0) ? z : parsed_z;

    double dbe = envi_iso::calc_dbe(counts);
    int h = counts.count("H") ? counts.at("H") : 0;
    int c = counts.count("C") ? counts.at("C") : 0;
    double hc = (c > 0) ? (static_cast<double>(h) / c) : 0.0;

    bool senior_pass = envi_iso::check_senior_rules(counts, final_z);
    bool hc_pass = (c <= 0) || (hc >= min_hc && hc <= max_hc);
    bool hetero_pass = envi_iso::check_heteroatom_ratios(counts);
    bool dbe_pass = (dbe >= -0.5 && dbe <= max_dbe);

    double nominal_mass = 0.0;
    int n_count = counts.count("N") ? counts.at("N") : 0;
    const auto& table = envi_iso::get_element_table();
    for (const auto& pair : counts) {
        auto it = table.find(pair.first);
        if (it != table.end()) nominal_mass += it->second.nominal_mass * pair.second;
    }
    bool nitrogen_pass = envi_iso::is_valid_nitrogen_rule(nominal_mass, n_count, final_z);

    bool overall = senior_pass && hc_pass && hetero_pass && dbe_pass && nitrogen_pass;

    return List::create(
        Named("valid") = overall,
        Named("senior_rules") = senior_pass,
        Named("hc_ratio") = hc,
        Named("hc_pass") = hc_pass,
        Named("heteroatom_pass") = hetero_pass,
        Named("dbe") = dbe,
        Named("dbe_pass") = dbe_pass,
        Named("nitrogen_rule") = nitrogen_pass
    );
}

//' @title Calculate Isotopic Pattern Similarity Scores
//' @description Computes cosine similarity, weighted dot product, and log-likelihood
//'   between experimental peaks and theoretical isotopic distribution.
//' @param obs_mz Numeric vector of observed m/z values.
//' @param obs_int Numeric vector of observed peak intensities.
//' @param theo_mz Numeric vector of theoretical isotopic m/z values.
//' @param theo_int Numeric vector of theoretical isotopic abundances.
//' @param tolerance Mass matching tolerance (default: 0.005 Da or 10 ppm).
//' @param ppm Logical indicating if tolerance is in ppm (default: FALSE).
//' @return A list containing cosine score, weighted cosine score, log-likelihood,
//'   and matching details.
//' @keywords internal
// [[Rcpp::export]]
List rcpp_score_isotopes(
    NumericVector obs_mz,
    NumericVector obs_int,
    NumericVector theo_mz,
    NumericVector theo_int,
    double tolerance = 0.005,
    bool ppm = false
) {
    size_t n_obs = obs_mz.size();
    size_t n_theo = theo_mz.size();
    if (n_obs == 0 || n_theo == 0) {
        return List::create(
            Named("cosine") = 0.0,
            Named("weighted_cosine") = 0.0,
            Named("log_likelihood") = -1e6,
            Named("matched_peaks") = 0
        );
    }

    std::vector<double> norm_theo(n_theo, 0.0);
    double sum_t = 0.0;
    for (size_t i = 0; i < n_theo; ++i) sum_t += theo_int[i];
    if (sum_t > 0) {
        for (size_t i = 0; i < n_theo; ++i) norm_theo[i] = theo_int[i] / sum_t;
    }

    std::vector<double> norm_obs(n_obs, 0.0);
    double sum_o = 0.0;
    for (size_t i = 0; i < n_obs; ++i) sum_o += obs_int[i];
    if (sum_o > 0) {
        for (size_t i = 0; i < n_obs; ++i) norm_obs[i] = obs_int[i] / sum_o;
    }

    std::vector<double> aligned_obs(n_theo, 0.0);
    int matched_count = 0;

    for (size_t j = 0; j < n_theo; ++j) {
        double tm = theo_mz[j];
        double tol = ppm ? (tm * tolerance * 1e-6) : tolerance;
        double best_diff = tol + 1.0;
        int best_idx = -1;

        for (size_t i = 0; i < n_obs; ++i) {
            double diff = std::abs(obs_mz[i] - tm);
            if (diff <= tol && diff < best_diff) {
                best_diff = diff;
                best_idx = static_cast<int>(i);
            }
        }

        if (best_idx >= 0) {
            aligned_obs[j] = norm_obs[best_idx];
            matched_count++;
        }
    }

    double dot_prod = 0.0, norm_o_sq = 0.0, norm_t_sq = 0.0;
    double w_dot_prod = 0.0, w_norm_o_sq = 0.0, w_norm_t_sq = 0.0;
    double log_lik = 0.0;

    for (size_t j = 0; j < n_theo; ++j) {
        double o = aligned_obs[j];
        double t = norm_theo[j];
        double m = theo_mz[j];

        dot_prod += o * t;
        norm_o_sq += o * o;
        norm_t_sq += t * t;

        double wo = std::sqrt(std::max(0.0, o)) * std::sqrt(m);
        double wt = std::sqrt(std::max(0.0, t)) * std::sqrt(m);
        w_dot_prod += wo * wt;
        w_norm_o_sq += wo * wo;
        w_norm_t_sq += wt * wt;

        if (o > 0) {
            log_lik += o * std::log(std::max(1e-10, t));
        }
    }

    double cosine = (norm_o_sq > 0 && norm_t_sq > 0) ? (dot_prod / (std::sqrt(norm_o_sq) * std::sqrt(norm_t_sq))) : 0.0;
    double w_cosine = (w_norm_o_sq > 0 && w_norm_t_sq > 0) ? (w_dot_prod / (std::sqrt(w_norm_o_sq) * std::sqrt(w_norm_t_sq))) : 0.0;

    return List::create(
        Named("cosine") = cosine,
        Named("weighted_cosine") = w_cosine,
        Named("log_likelihood") = log_lik,
        Named("matched_peaks") = matched_count
    );
}

// Internal helper to set up elements and search boundaries
static std::vector<envi_iso::DecompElement> setup_decomp_elements(
    SEXP elements,
    const std::string& minElements,
    const std::string& maxElements
) {
    std::map<std::string, int> min_bounds = {{"C", 0}, {"H", 0}, {"N", 0}, {"O", 0}, {"P", 0}, {"S", 0}};
    std::map<std::string, int> max_bounds = {{"C", 100}, {"H", 200}, {"N", 50}, {"O", 50}, {"P", 10}, {"S", 10}};

    if (!minElements.empty()) {
        int dummy_z = 0;
        auto parsed_min = envi_iso::parse_formula(minElements, dummy_z);
        for (const auto& p : parsed_min) min_bounds[p.first] = p.second;
    }
    if (!maxElements.empty()) {
        int dummy_z = 0;
        auto parsed_max = envi_iso::parse_formula(maxElements, dummy_z);
        for (const auto& p : parsed_max) max_bounds[p.first] = p.second;
    }

    std::vector<std::string> elem_symbols;
    if (Rf_isString(elements) && Rf_length(elements) > 0) {
        CharacterVector cv(elements);
        for (int i = 0; i < cv.size(); ++i) {
            std::string s = as<std::string>(cv[i]);
            int dummy_z = 0;
            auto parsed = envi_iso::parse_formula(s, dummy_z);
            for (const auto& p : parsed) {
                if (std::find(elem_symbols.begin(), elem_symbols.end(), p.first) == elem_symbols.end()) {
                    elem_symbols.push_back(p.first);
                }
            }
        }
    } else if (Rf_isNewList(elements) && Rf_length(elements) > 0) {
        List lst(elements);
        CharacterVector nms = lst.names();
        if (nms.size() > 0) {
            for (int i = 0; i < nms.size(); ++i) {
                std::string s = as<std::string>(nms[i]);
                if (!s.empty()) elem_symbols.push_back(s);
            }
        } else {
            for (int i = 0; i < lst.size(); ++i) {
                SEXP item = lst[i];
                if (Rf_isNewList(item)) {
                    List sub(item);
                    if (sub.containsElementNamed("name")) {
                        elem_symbols.push_back(as<std::string>(sub["name"]));
                    }
                }
            }
        }
    }

    if (elem_symbols.empty()) {
        elem_symbols = {"C", "H", "N", "O", "P", "S"};
    }

    const auto& table = envi_iso::get_element_table();
    std::vector<envi_iso::DecompElement> decomp_elems;
    for (const auto& sym : elem_symbols) {
        auto it = table.find(sym);
        if (it != table.end()) {
            int mn = min_bounds.count(sym) ? min_bounds[sym] : 0;
            int mx = max_bounds.count(sym) ? max_bounds[sym] : 50;
            decomp_elems.push_back({sym, it->second.exact_mass, mn, mx});
        }
    }

    std::sort(decomp_elems.begin(), decomp_elems.end(), [](const envi_iso::DecompElement& a, const envi_iso::DecompElement& b) {
        if (a.symbol == "H") return false;
        if (b.symbol == "H") return true;
        return a.mass > b.mass;
    });

    return decomp_elems;
}

//' @title Decompose Mass into Candidate Chemical Formulas using the HORIZON Algorithm
//' @description Finds all combinations of elements matching a mass within given error bounds using the HORIZON (Heavy-first Ordered Recursive Inference with Zero-loop Optimal Navigation) algorithm.
//' @param mass Target mass or m/z value.
//' @param ppm Tolerance in ppm (default: 2.0).
//' @param mzabs Absolute mass deviation in Da (default: 0.0001).
//' @param elements Allowed chemical elements (string or vector).
//' @param minElements Lower bounds on elements (e.g. "C0H0", "C1").
//' @param maxElements Upper bounds on elements (e.g. "C50H50").
//' @param z Charge state (default: 0).
//' @param maxisotopes Maximum number of isotopes to compute (default: 10).
//' @param golden_rules Whether to apply Fiehn Seven Golden Rules (default: TRUE).
//' @param resolution Mass resolution for isotope merging (default: 0).
//' @param fwhm Full width at half maximum in Da for isotope peak merging (default: 0).
//' @return A list compatible with Rdisop::decomposeMass.
//' @keywords internal
// [[Rcpp::export]]
List rcpp_decompose_mass(
    double mass,
    double ppm = 2.0,
    double mzabs = 0.0001,
    SEXP elements = R_NilValue,
    std::string minElements = "",
    std::string maxElements = "",
    int z = 0,
    int maxisotopes = 10,
    bool golden_rules = true,
    double resolution = 0.0,
    double fwhm = 0.0
) {
    double tol = mzabs + (ppm * 1e-6 * mass);
    double target_neutral_mass = mass;
    if (z != 0) {
        target_neutral_mass = mass * std::abs(z) + z * envi_iso::ELECTRON_MASS;
    }

    auto decomp_elems = setup_decomp_elements(elements, minElements, maxElements);
    size_t num_elems = decomp_elems.size();
    if (num_elems == 0) {
        return List::create(
            Named("formula") = CharacterVector(),
            Named("score") = NumericVector(),
            Named("exactmass") = NumericVector(),
            Named("charge") = IntegerVector(),
            Named("parity") = CharacterVector(),
            Named("valid") = CharacterVector(),
            Named("DBE") = NumericVector(),
            Named("isotopes") = List(),
            Named("hc_ratio") = NumericVector(),
            Named("senior_rules") = LogicalVector()
        );
    }

    std::vector<double> min_remain(num_elems, 0.0);
    std::vector<double> max_remain(num_elems, 0.0);
    for (int i = (int)num_elems - 1; i >= 0; --i) {
        double mn = decomp_elems[i].min_count * decomp_elems[i].mass;
        double mx = decomp_elems[i].max_count * decomp_elems[i].mass;
        if (i == (int)num_elems - 1) {
            min_remain[i] = mn;
            max_remain[i] = mx;
        } else {
            min_remain[i] = min_remain[i + 1] + mn;
            max_remain[i] = max_remain[i + 1] + mx;
        }
    }

    std::vector<int> cur_counts(num_elems, 0);
    std::vector<envi_iso::DecompResult> results;
    envi_iso::dfs_decompose(
        0, target_neutral_mass, decomp_elems, cur_counts, tol, results,
        min_remain, max_remain, z, golden_rules
    );

    std::sort(results.begin(), results.end(), [mass](const envi_iso::DecompResult& a, const envi_iso::DecompResult& b) {
        return std::abs(a.exactmass - mass) < std::abs(b.exactmass - mass);
    });

    return envi_iso::format_decomp_result_list(results, z, maxisotopes, resolution, fwhm);
}

//' @title Decompose a Vector of Masses in Parallel using the HORIZON Algorithm
//' @description Vectorized, multi-threaded chemical formula decomposition across thousands of masses using the HORIZON (Heavy-first Ordered Recursive Inference with Zero-loop Optimal Navigation) algorithm.
//' @param masses Numeric vector of target masses.
//' @param ppm Tolerance in ppm (default: 2.0).
//' @param mzabs Absolute mass deviation in Da (default: 0.0001).
//' @param elements Allowed chemical elements.
//' @param minElements Lower bounds on elements.
//' @param maxElements Upper bounds on elements.
//' @param z Charge state (default: 0).
//' @param maxisotopes Maximum number of isotopes to compute (default: 10).
//' @param golden_rules Whether to apply Fiehn Seven Golden Rules (default: TRUE).
//' @param resolution Mass resolution for isotope merging (default: 0).
//' @param fwhm Full width at half maximum in Da for isotope peak merging (default: 0).
//' @param nthreads Number of OpenMP threads to use (default: 1).
//' @return A list of decomposition results corresponding to each input mass.
//' @keywords internal
// [[Rcpp::export]]
List rcpp_decompose_masses(
    NumericVector masses,
    double ppm = 2.0,
    double mzabs = 0.0001,
    SEXP elements = R_NilValue,
    std::string minElements = "",
    std::string maxElements = "",
    int z = 0,
    int maxisotopes = 10,
    bool golden_rules = true,
    double resolution = 0.0,
    double fwhm = 0.0,
    int nthreads = 1
) {
    size_t n_masses = masses.size();
    if (n_masses == 0) return List();

    auto decomp_elems = setup_decomp_elements(elements, minElements, maxElements);
    size_t num_elems = decomp_elems.size();
    if (num_elems == 0) {
        List empty_out(n_masses);
        for (size_t i = 0; i < n_masses; ++i) {
            empty_out[i] = List::create(
                Named("formula") = CharacterVector(),
                Named("score") = NumericVector(),
                Named("exactmass") = NumericVector(),
                Named("charge") = IntegerVector(),
                Named("parity") = CharacterVector(),
                Named("valid") = CharacterVector(),
                Named("DBE") = NumericVector(),
                Named("isotopes") = List(),
                Named("hc_ratio") = NumericVector(),
                Named("senior_rules") = LogicalVector()
            );
        }
        return empty_out;
    }

    std::vector<double> min_remain(num_elems, 0.0);
    std::vector<double> max_remain(num_elems, 0.0);
    for (int i = (int)num_elems - 1; i >= 0; --i) {
        double mn = decomp_elems[i].min_count * decomp_elems[i].mass;
        double mx = decomp_elems[i].max_count * decomp_elems[i].mass;
        if (i == (int)num_elems - 1) {
            min_remain[i] = mn;
            max_remain[i] = mx;
        } else {
            min_remain[i] = min_remain[i + 1] + mn;
            max_remain[i] = max_remain[i + 1] + mx;
        }
    }

    std::vector<std::vector<envi_iso::DecompResult>> all_results(n_masses);

#ifdef _OPENMP
    int max_threads = omp_get_num_procs();
    int actual_threads = (nthreads > 0) ? std::min(nthreads, max_threads) : 1;
#else
    int actual_threads = 1;
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(actual_threads) if(actual_threads > 1)
#endif
    for (size_t i = 0; i < n_masses; ++i) {
        double m = masses[i];
        double tol = mzabs + (ppm * 1e-6 * m);
        double target_neutral_mass = m;
        if (z != 0) {
            target_neutral_mass = m * std::abs(z) + z * envi_iso::ELECTRON_MASS;
        }

        std::vector<int> cur_counts(num_elems, 0);
        std::vector<envi_iso::DecompResult> local_res;
        envi_iso::dfs_decompose(
            0, target_neutral_mass, decomp_elems, cur_counts, tol, local_res,
            min_remain, max_remain, z, golden_rules
        );

        std::sort(local_res.begin(), local_res.end(), [m](const envi_iso::DecompResult& a, const envi_iso::DecompResult& b) {
            return std::abs(a.exactmass - m) < std::abs(b.exactmass - m);
        });

        all_results[i] = std::move(local_res);
    }

    List final_out(n_masses);
    for (size_t i = 0; i < n_masses; ++i) {
        final_out[i] = envi_iso::format_decomp_result_list(all_results[i], z, maxisotopes, resolution, fwhm);
    }

    return final_out;
}
