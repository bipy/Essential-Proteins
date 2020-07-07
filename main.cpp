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
vector<vector<int> > matrix;
unordered_map<string, unordered_map<string, bool> > m;
vector<node> ans;
vector<int> path_count;
vector<string> methods;

string ref_data, input_data;

void init();

void trans_matrix();

bool cmp(const node &a, const node &b);

void output(vector<node> &cur, vector<double> rate, const string &method);

vector<double> compare_ref(int step);

void DC();

void CC();

void BC();

void check(const string &method);

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
        } else if (cur == "-d") {
            methods.emplace_back("DC");
        } else if (cur == "-c") {
            methods.emplace_back("CC");
        } else if (cur == "-b") {
            methods.emplace_back("BC");
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
    if (input_data.empty() || ref_data.empty() || methods.empty()) {
        printf("Invalid args!\nUse -h for help.");
        return 1;
    }
    init();
    sort(methods.begin(), methods.end());
    for (auto &it:methods) {
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
    printf("-d Use algorithm DC.\n");
    printf("-c Use algorithm CC.\n");
    printf("-b Use algorithm BC.\n\n");
    printf("=========== Caution ===========\n\n");
    printf("Must have at least one of '-b', '-c' or '-d'.\n");
    printf("Must have '-r' and 'your refer data's path'.\n");
    printf("Must have '-i' and  'your input data's path'.\n");
    printf("Algorithm BC or CC will take a long trip to run (like O(N^3)), but it works!\n\n");
    printf("============ Tips ============\n\n");
    printf("Use '-b -c' together (save 50%% time)\n\n");
    printf("============ About ============\n\n");
    printf("Author: bipy@GitHub\n");
    printf("Version: 20200707.2\n\n");
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

void trans_matrix() {
    matrix.resize(m.size(), vector<int>(m.size(), INT16_MAX));
    int index = 0;
    for (auto &it : m) {
        trans[it.first] = index;
        reverse_trans[index] = it.first;
        index++;
    }
    for (auto &it : m) {
        for (auto &i:it.second) {
            int a = trans[it.first], b = trans[i.first];
            matrix[a][b] = 1;
        }
    }
}

void output(vector<node> &cur, vector<double> rate, const string &method) {
    string dest = input_data.substr(0, input_data.size() - 4) + " output_" + method + ".txt";
    ofstream ans_out(dest, ios::out);
    for (int i = 0; i < rate.size(); i++) {
        ans_out << (i + 1) * 100 << " " << rate[i] << endl;
    }
    ans_out << endl;
    for (auto &it : cur) {
        ans_out << it.id << endl;
    }
    ans_out.close();
}

vector<double> compare_ref(const int step) {
    vector<double> rt;
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
    cout << method + " start!" << endl;
    if (method == "DC") {
        DC();
    } else if (method == "CC") {
        CC();
    } else if (method == "BC") {
        BC();
    }
    sort(ans.begin(), ans.end(), cmp);
    output(ans, compare_ref(100), method);
    cout << method + " complete!" << endl;
}

void floyd(bool count_path) {
    trans_matrix();
    int size = matrix.size();
    if (count_path) {
        path_count.resize(size, 0);
    }
    for (int k = 0; k < size; k++) {
        if (k % 100 == 0) {
            printf("Progress: %.1f %%\n", 100 * static_cast<double>(k) / size);
        }
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (matrix[i][j] > matrix[i][k] + matrix[k][j]) {
                    matrix[i][j] = matrix[i][k] + matrix[k][j];
                    if (count_path) {
                        path_count[k]++;
                    }
                }
            }
        }
    }
}

void DC() {
    for (auto &it : m) {
        ans.push_back(node{it.first, static_cast<double>(it.second.size())});
    }
}

void CC() {
    if (path_count.empty()) {
        floyd(false);
    }
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
    int size = path_count.size();
    for (int i = 0; i < size; i++) {
        ans.push_back(node{reverse_trans[i], static_cast<double>(path_count[i])});
    }
}