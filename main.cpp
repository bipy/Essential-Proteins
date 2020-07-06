#include <iostream>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <fstream>
#include <regex>

using namespace std;
struct node {
    string id;
    double val;
};

unordered_set<string> ref_set;
unordered_map<string, int> trans;
unordered_map<int, string> reverse_trans;
vector<vector<int> > matrix, path;
unordered_map<string, unordered_map<string, bool> > m;
vector<node> ans;

string ref_data, input_data;

void init();

void trans_matrix();

bool cmp(const node &a, const node &b);

void output(vector<node> &cur, vector<double> rate, const string &method);

vector<double> compare_ref(const int step);

void DC();

void CC();

void BC();

void check(const string &method);

void help();

int main(int argc, char **argv) {
    string method;
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
        } else if (cur == "-d") {
            method = "DC";
        } else if (cur == "-c") {
            method = "CC";
        } else if (cur == "-b") {
            method = "BC";
        } else if (cur == "-r") {
            if (i + 1 < argc) {
                ref_data = argv[i + 1];
                i++;
            }
        } else if (cur == "-i") {
            if (i + 1 < argc) {
                input_data = argv[i + 1];
                i++;
            }
        }
    }
    if (input_data.empty() || ref_data.empty() || method.empty()) {
        printf("Invalid args!\nUse -h for help.");
        return 1;
    }
    init();
    check(method);
    return 0;
}

bool cmp(const node &a, const node &b) {
    return a.val > b.val;
}

void init() {
    ifstream fin(input_data, ios::in);
    if (!fin) {
        cout << input_data + " NOT FOUND!" << endl;
        return;
    }
    string node_a, node_b;
    getline(fin, node_a);
    if (!regex_match(node_a, regex("^node.*$"))) {
        fin.seekg(ios::beg);
    }
    while (getline(fin, node_a, '\t') && getline(fin, node_b)) {
        m[node_a][node_b] = m[node_b][node_a] = true;
    }
    fin.close();
}

void trans_matrix() {
    matrix.resize(m.size(), vector<int>(m.size(), INT16_MAX));
    int index = 0;
    for (auto it = m.begin(); it != m.end(); it++, index++) {
        trans[it->first] = index;
        reverse_trans[index] = it->first;
    }
    for (auto it = m.begin(); it != m.end(); it++) {
        for (auto i = it->second.begin(); i != it->second.end(); i++) {
            int a = trans[it->first], b = trans[i->first];
            matrix[a][b] = 1;
        }
    }
}

void output(vector<node> &cur, vector<double> rate, const string &method) {
    string dest = input_data.substr(0, input_data.size() - 4) + " output_" + method + ".txt";
    ofstream fout(dest, ios::out);
    for (int i = 0; i < rate.size(); i++) {
        fout << (i + 1) * 100 << " " << rate[i] << endl;
    }
    for (auto it = cur.begin(); it != cur.end(); it++) {
        fout << it->id << endl;
    }
    fout.close();
}

vector<double> compare_ref(const int step) {
    vector<double> rt;
    ifstream fin(ref_data, ios::in);
    if (!fin) {
        cout << ref_data + " NOT FOUND!" << endl;
        return rt;
    }
    string s;
    while (getline(fin, s)) {
        ref_set.insert(s);
    }
    int count = 0;
    for (int i = 0; i < ans.size(); i++) {
        if (ref_set.find(ans[i].id) != ref_set.end()) {
            count++;
        }
        if (i != 0 && i % step == 0) {
            rt.push_back(count / (double) i);
        }
    }
    return rt;
}

void check(const string &method) {
    if (method == "DC") {
        DC();
    } else if (method == "CC") {
        CC();
    } else if (method == "BC") {
        BC();
    }
    sort(ans.begin(), ans.end(), cmp);
    output(ans, compare_ref(100), method);
}

void help() {
    printf("==== Center Protein Calc ====\n");
    printf("Usage: [option] [value]\n");
    printf("-h See this.\n");
    printf("-i Specific input data path.\n");
    printf("-r Specific reference data path.\n");
    printf("-d Use algorithm DC.\n");
    printf("-c Use algorithm CC.\n");
    printf("-b Use algorithm BC.\n");
    printf("=========================\n");
    printf("Author: bipy@GitHub\n");
    printf("Version: 20200707.1\n");
}

void floyd(bool count_path) {
    trans_matrix();
    int size = matrix.size();
    if (count_path) {
        path.resize(size, vector<int>(size, 0));
    }
    for (int k = 0; k < size; k++) {
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

void DC() {
    for (auto it = m.begin(); it != m.end(); it++) {
        ans.push_back(node{it->first, static_cast<double>(it->second.size())});
    }
}

void CC() {
    floyd(false);
    int size = matrix.size();
    for (int i = 0; i < size; i++) {
        int sum = 0;
        for (int j = 0; j < size; j++) {
            sum += matrix[i][j];
        }
        ans.push_back(node{reverse_trans[i], 10000.0 / sum});
    }
}

void BC() {
    floyd(true);
    int size = matrix.size();
    vector<int> pass_count(size, 0);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            pass_count[path[i][j]]++;
        }
    }
    for (int i = 0; i < size; i++) {
        ans.push_back(node{reverse_trans[i], static_cast<double>(pass_count[i])});
    }
}