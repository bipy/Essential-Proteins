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
#include <cmath>

using namespace std;
// data unit
struct node {
    string id;
    double val;
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
// Log(factorial) record
vector<double> log_fac;

vector<int> path_count;
set<string> methods;

string ref_data, input_data;
// step of comparison result
int step = 100;

// check and input data
bool init();

// trans list to matrix
void trans_matrix(int initial_val);

// comparator
bool cmp_less(const node &a, const node &b);

bool cmp_greater(const node &a, const node &b);

// output results
// vector<node> &cur: sorted data
// vector<double> rate: comparison results
// const string &method: algorithm name
void output(vector<node> &cur, vector<int> count, const string &method);

// compare with reference data
vector<int> compare_ref();

// BC: count pass through between "start" and "end"
void get_path_count(const int &start, const int &end);

// return log(x!)
double get_log_fac(unsigned int x);

// return log(nCr(N,r))
double get_log_comb(unsigned int N, unsigned int r);

// return nCr(N,r)
unsigned int get_exact_comb(unsigned int N, unsigned int r);

// vector normalization
double normalize(vector<double> &v);

// Algorithm: Significance-Based Essential Protein Discovery
void SigEP();

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
            methods.insert("SigEP");
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
        return 2;
    }
    if (!init()) {
        return 3;
    }
    if (*methods.begin() == "ALL") {
        methods = {"BC", "CC", "DC", "EC", "SigEP"};
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
    printf("-a Use 5 centrality algorithms together (BC, CC, DC, EC, SigEP).\n");
    printf("-b Use algorithm Betweenness Centrality (BC).\n");
    printf("-c Use algorithm Closeness Centrality (CC).\n");
    printf("-d Use algorithm Degree Centrality (DC).\n");
    printf("-e Use algorithm Eigenvector Centrality (EC).\n");
    printf("-p Use algorithm Significance-Based Essential Protein Discovery (SigEP).\n\n");
    printf("=========== Caution ===========\n\n");
    printf("Must have at least one of ['-a', '-b', '-c', '-d', '-e', '-p'].\n");
    printf("Must have '-r' and 'your refer data's path'.\n");
    printf("Must have '-i' and  'your input data's path'.\n\n");
    printf("============ Tips ============\n\n");
    printf("Algorithm BC or CC will take a long trip to run (like O(N^3)), but it works!\n");
    printf("Use '-b -c' together (save you 50%% time)\n\n");
    printf("============ About ============\n\n");
    printf("Author: bipy@GitHub\n");
    printf("Version: 20200723.1\n\n");
}

bool cmp_less(const node &a, const node &b) {
    return a.val < b.val;
}

bool cmp_greater(const node &a, const node &b) {
    return a.val > b.val;
}


bool init() {
    // init input
    ifstream data_in(input_data, ios::in);
    if (!data_in) {
        cout << input_data + " NOT FOUND!" << endl;
        return false;
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
        return false;
    }
    while (getline(ref_in, s)) {
        ref_set.insert(s);
    }
    ref_in.close();
    return true;
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
    } else if (method == "SigEP") {
        SigEP();
    }
    output(ans, compare_ref(), method);
    ans.clear();
    printf("%s complete!\n\n", method.c_str());
}

void floyd(bool count_path) {
    // trans adjacent list to matrix
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

double get_log_fac(unsigned int x) {
    for (unsigned int i = log_fac.size(); i <= x; i++) {
        log_fac.push_back(log_fac.back());
        log_fac[i] += log2(i);
    }
    return log_fac[x];
}

double get_log_comb(unsigned int N, unsigned int r) {
    if (N < r) {
        return get_log_comb(r, N);
    }
    // nCr(N,r) = N! / r! / (N - r)!
    return get_log_fac(N) - get_log_fac(r) - get_log_fac(N - r);
}

unsigned int get_exact_comb(unsigned int N, unsigned int r) {
    return lround(exp2(get_log_comb(N, r)));
}

double normalize(vector<double> &v) {
    double sum = 0.0;
    for (auto &it:v) {
        sum += it;
    }
    sum = sqrt(sum);
    for (auto &it:v) {
        it /= sum;
    }
    return sum;
}

void EC() {
    // trans adjacent list to matrix
    trans_matrix(-1);
    int size = matrix.size();
    // Eigenvector
    vector<double> cur(size, -1), last(size, -1);
    for (int i = 0; i < size; i++) {
        cur[i] = double(m[reverse_trans[i]].size());
    }
    double V_cur = 0.1, V_last = 0.0;
    // travel until V_cur == V_last
    while (fabs(V_cur - V_last) > 0.000001) {
        last = cur;
        V_last = V_cur;
        // foreach row
        for (int i = 0; i < size; i++) {
            cur[i] = 0.0;
            // foreach col of row i
            for (int j = 0; j < size; j++) {
                if (matrix[i][j] != -1) {
                    cur[i] += last[j] * matrix[i][j];
                }
            }
        }
        // normalize cur
        V_cur = normalize(cur);
    }
    // push into ans
    for (int i = 0; i < size; i++) {
        ans.emplace_back(node{reverse_trans[i], cur[i]});
        //printf("%s %f\n", reverse_trans[i].c_str(), cur[i]);
    }
    // descending sort
    sort(ans.begin(), ans.end(), cmp_greater);
}

void DC() {
    // push indegree into ans
    for (auto &it : m) {
        ans.emplace_back(node{it.first, static_cast<double>(it.second.size())});
    }
    // descending sort
    sort(ans.begin(), ans.end(), cmp_greater);
}

void CC() {
    // STL "set" ensures Algo:BC is ahead of Algo:CC
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
        // use distance sum as value push into ans
        ans.emplace_back(node{reverse_trans[i], static_cast<double>(sum)});
    }
    // ascending sort
    sort(ans.begin(), ans.end(), cmp_less);
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
    // descending sort
    sort(ans.begin(), ans.end(), cmp_greater);
}

void SigEP() {
    /**
     *  To ensure C_N_2 is less than 0xffffffffU
     *  N must be less than about 92,000
     */
    // N: the number of vertexes
    // M: the number of edges
    unsigned int N = m.size(), M = 0;
    // calc M
    for (auto &it : m) {
        M += it.second.size();
    }
    M /= 2;
    // init log_fac: log(0!) = 0.0
    log_fac.push_back(0.0);
    // C_N_2: nCr(N,2)
    unsigned int C_N_2 = get_exact_comb(N, 2);
    // de: log(denominator (pi))
    double de = get_log_comb(C_N_2, M);
    // foreach vertex
    for (auto &it : m) {
        // p: log(pi)
        double p = 0.0;
        // e: neighbour's edges
        // cur_d: current vertex's degree
        int e = 0, cur_d = it.second.size();
        // count e
        for (auto i = it.second.begin(); i != it.second.end(); i++) {
            for (auto j = i; j != it.second.end(); j++) {
                if (i == j) continue;
                if (m[i->first].find(j->first) != m[i->first].end()) {
                    e++;
                }
            }
        }
        // when d == 1, specific c = 0, then ignore est
        // c: local clustering coefficient
        double c = cur_d == 1 ? 0.0 : e / static_cast<double>(cur_d * (cur_d - 1));
        // est: the third arg of min()
        unsigned int est = c == 0.0 ? M : static_cast<unsigned int>(
                floor((c - 2 + sqrt(pow(2 - c, 2) + 8 * c * M)) / (2 * c))
        );
        // beta: the upper limit of degree
        unsigned int beta = min(min(est, M), N - 1);

        // foreach degree: from current vertex's degree to beta's
        for (unsigned int d = cur_d; d < beta; d++) {
            // C_d_2: nCr(d,2)
            unsigned int C_d_2 = get_exact_comb(d, 2);
            // cc: c cross C_d_2
            auto cc = static_cast<unsigned int>(ceil(c * C_d_2));
            // nu_0: the first term of numerator
            double nu_0 = get_log_comb(N - 1, d);
            // nu_1: the second term of numerator
            double nu_1 = get_log_comb(C_d_2, cc);
            // nu_2: the third term of numerator
            double nu_2 = get_log_comb(C_N_2 - (N - 1) - cc, M - d - cc);

            p += nu_0 + nu_1 + nu_2;
        }
        p += log(N) - de;
        ans.emplace_back(node{it.first, p});
    }
    // ascending sort
    sort(ans.begin(), ans.end(), cmp_less);
}