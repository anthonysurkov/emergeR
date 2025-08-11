// NOTE - This version is nonfunctional on this machine. It is here for consistency's sake.
//        Originally built on Anthony's local machine.

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <unordered_map>
#include <tuple>

extern "C" {
    #include <ViennaRNA/fold.h>
    #include <ViennaRNA/utils.h>
}

// CONSTS
std::string TARGET_FIVE =
    "GGGAAUGGCGCUUUGG"
    "AGAACAACUCUAGGCA"
    "GAGGUCUCAAAUUUCU"
    "UUCUUUCUUUGAGACCU";
std::string TARGET_THREE =
    "UUGUUCUCCAUUGUACCCUC";

const char* HOME = std::getenv("HOME");
const std::string INFILE = std::string(HOME) +
    "/Documents/ucd/beal/R255XE488Q/data/top_editors_by_cluster.csv";
const std::string OUTFILE = std::string(HOME) +
    "/Documents/ucd/beal/R255XE488Q/data/mfe.csv";

const int NUM_GRPS = 7;

// STRUCTS
struct Result {
    std::string n10;
    std::string edit;
    float       mfe;
    std::vector<std::string> partner;
};


/* HELPER FUNCTIONS * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Print sequence, structure prediction, and free energy prediction.
void print_mfe(
    const std::string sequence,
    const char *structure_data,
    float mfe
) {
    std::cout << "Seq:  " << sequence << std::endl;
    std::cout << "Fold: " << structure_data << std::endl;
    std::cout << "MFE:  " << mfe << " kcal/mol" << std::endl;

}

// Compute MFEs for a guide.
Result do_this_guide(
    const std::string &n10,
    const std::string edit
) {
    std::string sequence = TARGET_FIVE + n10 + TARGET_THREE;
    int L = (int)sequence.length();

    // Set up containers
    std::vector<char> buf(L + 1);

    // Compute mfe
    float mfe = fold(sequence.c_str(), buf.data());
    std::string structure(buf.data());

    // Find base pairing partners
    short *ptable = vrna_ptable(structure.c_str());

    std::vector<std::string> partner(L);
    for (int i = 1; i <= L; ++i) {
        int j = ptable[i];
        if (j > 0) {
            partner[i-1] = std::string(1, sequence[j-1])
                         + std::to_string(j);
        }
    }

    // Clean up
    free_arrays();
    free(ptable);

    return { n10, edit, mfe, std::move(partner) };
}

// Insert commas into Vienna's structure output for easier processing later.
std::string commaify(const std::string &s)
{
    std::string out;
    out.reserve(s.size() * 2);
    for (size_t i = 0; i < s.size(); ++i) {
        out.push_back(s[i]);
        if (i + 1 < s.size()) {
            out.push_back(',');
        }
    }

    return out;
}


/* MAIN * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

int main() {
    // Configure io handlers
    std::ifstream in(INFILE);
    std::ofstream out(OUTFILE);
    if (!in) {
        std::cerr << "Failed to open " << INFILE << '\n';
        return 1;
    }
    if (!out) {
        std::cerr << "Failed to open " << OUTFILE << '\n';
    }

    // Discard infile header, write outfile header
    std::string header =
        "group,n10,edit,mfe,"
        "G1,G2,G3,A4,A5,U6,G7,G8,C9,G10,C11,U12,"
        "U13,U14,G15,G16,A17,G18,A19,A20,C21,A22,A23,C24,"
        "U25,C26,U27,A28,G29,G30,C31,A32,G33,A34,G35,G36,"
        "U37,C38,U39,C40,A41,A42,A43,U44,U45,U46,C47,U48,"
        "U49,U50,C51,U52,U53,U54,C55,U56,U57,U58,G59,A60,"
        "G61,A62,C63,C64,U65,N1,N2,N3,N4,N5,N6,N7,"
        "N8,N9,N10,U76,U77,G78,U79,U80,C81,U82,C83,C84,"
        "A85,U86,U87,G88,U89,A90,C91,C92,C93,U94,C95";

    std::string line;
    std::getline(in, line);
    out << header << '\n';

    // map: groupId -> list of n10s
    std::unordered_map<int, std::vector<Result>> groups;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        size_t p1 = line.find(',');
        if (p1 == std::string::npos) continue;
        size_t p2 = line.find(',', p1+1);
        if (p2 == std::string::npos) continue;

        std::string n10 = line.substr(0, p1);
        std::string edit = line.substr(p1+1, p2-p1-1);
        int         grp = std::stoi(line.substr(p2+1));

        Result res = do_this_guide(n10, edit);
        groups[grp].push_back(std::move(res));
    }

    for (auto const& [grpID, vec] : groups) {
        for (auto const& r : vec) {
            out << grpID << ','
                << r.n10 << ','
                << r.edit << ','
                << r.mfe;
            for (auto const& p : r.partner) {
                out << ',' << (p.empty() ? "" : p);
            }
            out << '\n';
        }
    }

    return 0;
}