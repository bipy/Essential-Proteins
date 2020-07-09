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

vector<int> path_count;
set<string> methods;

string ref_data, input_data;
int step = 100;

// check and input data
void init();

// trans list to matrix
void trans_matrix(int initial_val);

// comparator
bool cmp(const node &a, const node &b);

// output results
// vector<node> &cur: sorted data
// vector<double> rate: comparison results
// const string &method: algorithm name
void output(vector<node> &cur, vector<double> rate, const string &method);

// compare with reference data
// int step: step of comparison result
vector<double> compare_ref();

// BC: count pass through between "start" and "end"
void get_path_count(const int &start, const int &end);

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
                    step = stoi(argv[i+1]);
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
        methods = {"BC", "CC", "DC", "EC"};
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
    printf("-a Use 4 centrality algorithms together (BC, CC, DC, EC).\n");
    printf("-b Use algorithm Betweenness Centrality (BC).\n");
    printf("-c Use algorithm Closeness Centrality (CC).\n");
    printf("-d Use algorithm Degree Centrality (DC).\n");
    printf("-e Use algorithm Eigenvector Centrality (EC).\n\n");
    printf("=========== Caution ===========\n\n");
    printf("Must have at least one of ['-a', '-b', '-c', '-d', '-e'].\n");
    printf("Must have '-r' and 'your refer data's path'.\n");
    printf("Must have '-i' and  'your input data's path'.\n\n");
    printf("============ Tips ============\n\n");
    printf("Algorithm BC or CC will take a long trip to run (like O(N^3)), but it works!\n");
    printf("Use '-b -c' together (save you 50%% time)\n\n");
    printf("============ About ============\n\n");
    printf("Author: bipy@GitHub\n");
    printf("Version: 20200709.1\n\n");
}

bool cmp(const node &a, const node &b) {
    return a.val > b.val;
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

void output(vector<node> &cur, vector<double> rate, const string &method) {
    string dest = input_data.substr(0, input_data.size() - 4) + " output_" + method + ".txt";
    ofstream ans_out(dest, ios::out);
    for (int i = 0; i < rate.size(); i++) {
        ans_out << (i + 1) * step << " " << rate[i] << endl;
    }
    ans_out << endl;
    for (auto &it : cur) {
        ans_out << it.id << endl;
    }
    ans_out.close();
}

vector<double> compare_ref() {
    vector<double> rt;
    int count = 0;
    for (int i = 0; i < ans.size(); i++) {
        if (ref_set.find(ans[i].id) != ref_set.end()) {
            count++;
        }
        if (i != 0 && i % step == 0) {
            rt.emplace_back(count / (double) i);
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
        ans.emplace_back(node{reverse_trans[i], l});
    }
}

void DC() {
    for (auto &it : m) {
        ans.emplace_back(node{it.first, static_cast<double>(it.second.size())});
    }
}

void CC() {
    if (path.empty()) {
        floyd(false);
    }
    int size = matrix.size();
    for (int i = 0; i < size; i++) {
        int sum = 0;
        for (int j = 0; j < size; j++) {
            sum += matrix[i][j];
        }
        ans.emplace_back(node{reverse_trans[i], 10000.0 / sum});
    }
}

void BC() {
    floyd(true);
    int size = path.size();
    path_count.resize(size);
    for (int i = 0; i < size; i++) {
        if (i % 200 == 0) {
            printf("Counting: %.1f %%\n", 100 * static_cast<double>(i) / size);
        }
        for (int j = i + 1; j < size; j++) {
            get_path_count(i, j);
        }
    }
    for (int i = 0; i < size; i++) {
        ans.emplace_back(node{reverse_trans[i], static_cast<double>(path_count[i])});
    }
}