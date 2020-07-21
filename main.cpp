#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <string>
#include <fstream>
#include <regex>
#include <cctype>
#include <gmp.h>

using namespace std;
// data unit
struct node {
    string id;
    double val;
};
// gmp mpq data unit
struct mpq_node {
    string id;
    mpq_t val;
};

// reference data (const)
unordered_set<string> ref_set;
// mapping: string&int (clear each time)
unordered_map<string, int> trans;
unordered_map<int, string> reverse_trans;
// Adjacency Matrix (clear each time)
vector<vector<int> > matrix, path;
// Adjacency List (const)
unordered_map<string, unordered_map<string, bool> > m;
// Answer sheet (clear each time)
vector<node> ans;

vector<int> path_count;
set<string> methods;

string ref_data, input_data;
// step of comparison result
int step = 100;

// check and input data
void init();

// trans list to matrix
void trans_matrix(int initial_val);

// comparator
bool cmp(const node &a, const node &b);

bool mpq_unit_cmp(const mpq_node &a, const mpq_node &b);

// output results
// vector<node> &cur: sorted data
// vector<double> rate: comparison results
// const string &method: algorithm name
void output(vector<node> &cur, vector<int> count, const string &method);

// compare with reference data
vector<int> compare_ref();

// BC: count pass through between "start" and "end"
void get_path_count(const int &start, const int &end);

// get mpz_t nu&de of a Combination Number
// nCr(b,a)
void mpz_comb_nu_de_ui(mpz_ptr nu, mpz_ptr de, unsigned int a, unsigned int b);

// get result of a Combination Number
// nCr(b,a)
void mpz_comb_ui(mpz_ptr cur, unsigned int a, unsigned int b);

// Algorithm: Degree Centrality with p-value
void DC_P();

// Algorithm: Eigenvector Centrality
void EC();

// Algorithm: Degree Centrality
void DC();

// Algorithm: Closeness Centrality
void CC();

// Algorithm: Betweenness Centrality
void BC();

// method switcher
void check(const string &method);

// output help information
void help();

int main(int argc, char **argv) {
    if (argc == 1) {
        help();
        system("pause");
        return 0;
    }
    for (int i = 1; i < argc; i++) {
        string cur = argv[i];
        if (cur == "-h") {
            help();
            return 0;
        } else if (cur == "-p") {
            methods.insert("DC_P");
        } else if (cur == "-e") {
            methods.insert("EC");
        } else if (cur == "-d") {
            methods.insert("DC");
        } else if (cur == "-c") {
            methods.insert("CC");
        } else if (cur == "-b") {
            methods.insert("BC");
        } else if (cur == "-a") {
            methods.insert("ALL");
        } else if (cur == "-r") {
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                ref_data = argv[i + 1];
                i++;
            }
        } else if (cur == "-i") {
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                input_data = argv[i + 1];
                i++;
            }
        } else if (cur == "-s") {
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                bool flag = true;
                for (auto &it : string(argv[i + 1])) {
                    if (!isdigit(it)) {
                        flag = false;
                    }
                }
                if (flag) {
                    step = stoi(argv[i + 1]);
                    i++;
                }
            }
        } else {
            cout << "Unknown arg: " + cur << endl;
            return 1;
        }
    }
    if (input_data.empty() || ref_data.empty() || methods.empty()) {
        printf("Invalid args!\nUse -h for help.");
        return 1;
    }
    init();
    if (*methods.begin() == "ALL") {
        methods = {"BC", "CC", "DC", "EC", "DC_P"};
    }
    for (auto &it : methods) {
        check(it);
    }
    return 0;
}

void help() {
    printf("===== Center Protein Calc =====\n\n");
    printf("Usage: [option] [value]\n");
    printf("-h See this.\n");
    printf("-i Specific input data path.\n");
    printf("-r Specific reference data path.\n");
    printf("-s (optional) Specific step (default 100).\n");
    printf("-a Use 5 centrality algorithms together (BC, CC, DC, EC, DC_P).\n");
    printf("-b Use algorithm Betweenness Centrality (BC).\n");
    printf("-c Use algorithm Closeness Centrality (CC).\n");
    printf("-d Use algorithm Degree Centrality (DC).\n");
    printf("-e Use algorithm Eigenvector Centrality (EC).\n");
    printf("-p Use algorithm Degree Centrality with p-value (DC_P).\n\n");
    printf("=========== Caution ===========\n\n");
    printf("Must have at least one of ['-a', '-b', '-c', '-d', '-e', '-p'].\n");
    printf("Must have '-r' and 'your refer data's path'.\n");
    printf("Must have '-i' and  'your input data's path'.\n\n");
    printf("============ Tips ============\n\n");
    printf("Algorithm BC or CC will take a long trip to run (like O(N^3)), but it works!\n");
    printf("Use '-b -c' together (save you 50%% time)\n\n");
    printf("============ About ============\n\n");
    printf("Author: bipy@GitHub\n");
    printf("Version: 20200721.1\n\n");
}

bool cmp(const node &a, const node &b) {
    return a.val > b.val;
}

bool mpq_unit_cmp(const mpq_node &a, const mpq_node &b) {
    return mpq_cmp(a.val, b.val) < 0;
}

void init() {
    // init input
    ifstream data_in(input_data, ios::in);
    if (!data_in) {
        cout << input_data + " NOT FOUND!" << endl;
        return;
    }
    string node_a, node_b, s;
    getline(data_in, node_a);
    if (!regex_match(node_a, regex("^node.*$"))) {
        data_in.seekg(ios::beg);
    }
    while (getline(data_in, node_a, '\t') && getline(data_in, node_b)) {
        m[node_a][node_b] = m[node_b][node_a] = true;
    }
    data_in.close();
    // init reference
    ifstream ref_in(ref_data, ios::in);
    if (!ref_in) {
        cout << ref_data + " NOT FOUND!" << endl;
        return;
    }
    while (getline(ref_in, s)) {
        ref_set.insert(s);
    }
    ref_in.close();
}

void trans_matrix(int initial_val) {
    matrix.clear();
    trans.clear();
    reverse_trans.clear();
    matrix.resize(m.size(), vector<int>(m.size(), initial_val));
    int index = 0;
    for (auto &it : m) {
        trans[it.first] = index;
        reverse_trans[index] = it.first;
        index++;
    }
    for (auto &i : m) {
        for (auto &j : i.second) {
            matrix[trans[i.first]][trans[j.first]] = 1;
        }
    }
}

void output(vector<node> &cur, vector<int> count, const string &method) {
    // output filename
    string dest = input_data.substr(0, input_data.size() - 4) + " output_" + method + ".txt";
    ofstream ans_out(dest, ios::out);
    // count | correct
    for (int i = 0; i < count.size(); i++) {
        ans_out << (i + 1) * step << " " << count[i] << endl;
    }
    ans_out << endl;
    for (auto &it : cur) {
        ans_out << it.id << endl;
    }
    ans_out.close();
}

vector<int> compare_ref() {
    vector<int> rt;
    int count = 0;
    for (int i = 0; i < ans.size(); i++) {
        if (ref_set.find(ans[i].id) != ref_set.end()) {
            count++;
        }
        if (i % step == 0) {
            rt.emplace_back(count);
        }
    }
    return rt;
}

void get_path_count(const int &start, const int &end) {
    if (path[start][end] < 0) {
        return;
    }
    path_count[path[start][end]]++;
    get_path_count(start, path[start][end]);
    get_path_count(path[start][end], end);
}

void check(const string &method) {
    printf("%s start!\n", method.c_str());
    if (method == "DC") {
        DC();
    } else if (method == "CC") {
        CC();
    } else if (method == "BC") {
        BC();
    } else if (method == "EC") {
        EC();
    } else if (method == "DC_P") {
        DC_P();
    }
    sort(ans.begin(), ans.end(), cmp);
    output(ans, compare_ref(), method);
    ans.clear();
    printf("%s complete!\n\n", method.c_str());
}

void floyd(bool count_path) {
    trans_matrix(INT16_MAX);
    int size = matrix.size();
    if (count_path) {
        path.resize(size, vector<int>(size, -1));
    }
    for (int k = 0; k < size; k++) {
        // progress interface
        if (k % 50 == 0) {
            printf("Floyd: %.1f %%\n", 100 * static_cast<double>(k) / size);
        }
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (matrix[i][j] > matrix[i][k] + matrix[k][j]) {
                    matrix[i][j] = matrix[i][k] + matrix[k][j];
                    if (count_path) {
                        path[i][j] = k;
                    }
                }
            }
        }
    }
}

void mpz_comb_nu_de_ui(mpz_ptr nu, mpz_ptr de, unsigned int a, unsigned int b) {
    mpz_init_set_ui(nu, b);
    mpz_init_set_ui(de, a);
    for (unsigned int i = 1; i < a; i++) {
        mpz_mul_ui(nu, nu, i + b - a);
        mpz_mul_ui(de, de, i);
    }
}

void mpz_comb_ui(mpz_ptr cur, unsigned int a, unsigned int b) {
    if (a > b) {
        mpz_comb_ui(cur, b, a);
        return;
    }
    mpz_t de, nu;
    mpz_comb_nu_de_ui(nu, de, a, b);
    mpz_init(cur);
    mpz_div(cur, nu, de);
}

// unused method
//void mpz_comb(mpz_t &ptr, mpz_t &a, mpz_t &b) {
//    if (mpz_cmp(a, b) > 0) {
//        mpz_comb(ptr, b, a);
//        return;
//    }
//    mpz_t de, nu, i, b_a;
//    mpz_init_set(nu, b);
//    mpz_init_set(de, a);
//    mpz_sub(b_a, b, a);
//    for (mpz_init_set_ui(i, 2); mpz_cmp(i, a) < 0; mpz_add_ui(i, i, 1)) {
//        mpz_mul(de, de, i);
//    }
//    for (mpz_add_ui(i, b_a, 1); mpz_cmp(i, b) < 0; mpz_add_ui(i, i, 1)) {
//        mpz_mul(nu, nu, i);
//    }
//    mpz_div(ptr, nu, de);
//}

void EC() {
    trans_matrix(-1);
    int size = matrix.size();
    for (int i = 0; i < size; i++) {
        double l = 0.0;
        for (auto &j:matrix[i]) {
            if (j != -1) {
                l += static_cast<double>(m[reverse_trans[j]].size()) * j;
            }
        }
        // push into ans
        ans.emplace_back(node{reverse_trans[i], l});
    }
}

void DC() {
    // push indegree into ans
    for (auto &it : m) {
        ans.emplace_back(node{it.first, static_cast<double>(it.second.size())});
    }
}

void CC() {
    // set ensures Algo:BC is ahead of Algo:CC
    // if floyd has processed before, skip
    if (path.empty()) {
        floyd(false);
    }
    int size = matrix.size();
    for (int i = 0; i < size; i++) {
        int sum = 0;
        for (int j = 0; j < size; j++) {
            sum += matrix[i][j];
        }
        // use multiplicative inverse * 10,000 as value push into ans
        ans.emplace_back(node{reverse_trans[i], 10000.0 / sum});
    }
}

void BC() {
    floyd(true);
    int size = path.size();
    path_count.resize(size);
    for (int i = 0; i < size; i++) {
        // progress interface
        if (i % 200 == 0) {
            printf("Counting: %.1f %%\n", 100 * static_cast<double>(i) / size);
        }
        // count passed through vertexes
        for (int j = i + 1; j < size; j++) {
            get_path_count(i, j);
        }
    }
    // push into ans
    for (int i = 0; i < size; i++) {
        ans.emplace_back(node{reverse_trans[i], static_cast<double>(path_count[i])});
    }
}

void DC_P() {
    /**
     *  To ensure C_N_2 is less than 0xffffffffU
     *  N must be less than about 92,000
     *  requirement: libgmp 6.20
     */

    // N: the number of vertexes
    // M: the number of edges
    unsigned int N = m.size(), M = 0;
    // calc M
    for (auto &it : m) {
        M += it.second.size();
    }
    M /= 2;
    // limit: min(N - 1, M) (> 0)
    unsigned int limit = min(N - 1, M);
    // C_N_2: nCr(N,2)
    // de: denominator (pi)
    mpz_t C_N_2, de;
    mpz_comb_ui(C_N_2, 2, N);
    mpz_comb_ui(de, M, mpz_get_ui(C_N_2));
    // dp array
    vector<mpz_t> mpz_nu(limit + 1);
    // sorted list
    vector<mpq_node> mpq_v(N);
    // All about pi's numerator
    mpz_t de_0, nu_0, temp_0, de_1, nu_1, temp_1;
    mpz_init_set_ui(nu_0, 1);
    mpz_init_set_ui(de_0, 1);

    // foreach d: Calc pd O(N)
    for (unsigned int d = 0; d <= limit; d++) {
        // progress interface
        if (d % 100 == 0) {
            printf("p-value: %.1f %%\n", 100 * static_cast<double>(d) / limit);
        }
        // Let nCr(N,0) and nCr(N,N) be 1
        if (d != 0 && d != limit) {
            // trans to next form
            mpz_mul_ui(nu_0, nu_0, (N - 1) - d + 1);
            mpz_mul_ui(de_0, de_0, d);
            // (pi's numerator)'s first part
            mpz_divexact(temp_0, nu_0, de_0);
            // init (pi's numerator)'s second part
            if (d == 1) {
                mpz_comb_nu_de_ui(nu_1, de_1, M - 1, mpz_get_ui(C_N_2) - (N - 1));
            } else {
                // trans to next form
                // use multiplicative inverse instead of divsion
                mpz_mul_ui(nu_1, nu_1, M - d + 1);
                mpz_mul_ui(de_1, de_1, mpz_get_ui(C_N_2) - (N - 1) - (M - d));
            }
            mpz_divexact(temp_1, nu_1, de_1);
        } else {
            mpz_init_set_ui(temp_0, 1);
            mpz_init_set_ui(temp_1, 1);
        }
        // put first part cross second part into array
        mpz_init(mpz_nu[d]);
        mpz_mul(mpz_nu[d], temp_0, temp_1);
        // recursion
        if (d > 0) {
            mpz_add(mpz_nu[d], mpz_nu[d], mpz_nu[d - 1]);
        }
    }

    // foreach vertex i, calculate pi
    // init
    mpz_t pi_nu;
    mpz_init(pi_nu);
    mpq_t pi, pi_f_nu, pi_f_de;
    mpq_inits(pi, pi_f_nu, pi_f_de, NULL);
    mpq_set_z(pi_f_de, de);

    // calculate pi
    unsigned int i = 0;
    for (auto &it : m) {
        mpq_v[i].id = it.first;
        mpq_init(mpq_v[i].val);
        int d = it.second.size();
        if (d == 0) {
            mpz_set(pi_nu, mpz_nu.back());
        } else {
            mpz_sub(pi_nu, mpz_nu.back(), mpz_nu[d - 1]);
        }
        // trans to mpq_t
        mpq_set_z(pi_f_nu, pi_nu);
        mpq_div(mpq_v[i].val, pi_f_nu, pi_f_de);
        i++;
    }
    // ascending sort by pi
    sort(mpq_v.begin(), mpq_v.end(), mpq_unit_cmp);
    // push into ans according to the previous sort result
    double index = N + 1.0;
    for (auto &it : mpq_v) {
        ans.emplace_back(node{it.id, index--});
    }
}