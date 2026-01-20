#include <algorithm>
#include <bits/stdc++.h>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdint>
#include <chrono>
#include <stdexcept>

using namespace std;

enum class Initial_Place {base, fanout};

struct Instance{
    int n = 0;
    vector<int> block_num;
    unordered_map<int, int> num_to_index;
    vector<vector<int>>nets;
    vector<vector<int>> block_nets;
    vector<pair<int, int>> communities;
    vector<vector<int>>partners;
};

struct SearchState{
    vector<int> side; //-1 is assigned, 0 is left, 1 is right
    int left_count = 0, right_count = 0, target = 0;
    vector<int> net_left, net_right;
    vector<char> net_crossed;
    int crossing_count = 0;
    int comm_mismatches = 0;
};

void usage(const char* prog){
    cerr << "Usage: " << prog << " -file_name <path> -init_place <base|fanout>\n";
    exit(1);
}

struct VisTree {
  struct Node {
    int parent;  // -1 for root
    int block;
    int side;    // 0 = left and 1 = right
    int depth;
  };
  std::vector<Node> nodes;

  int add(int parent, int block, int side, int depth) {
    nodes.push_back({parent, block, side, depth});
    return (int)nodes.size() - 1;
  }
};

Instance read_file(const string &path){
    ifstream fin(path);
    if(!fin){
        cerr << "ERROR: cannot open file: " << path << "\n";
        exit(1);
    }
    vector<pair<int, vector<int>>> tmp;
    while(true){
        int b;
        if(!(fin >> b)) break;
        if(b == -1) break;

        vector<int> nets_of_b;
        while(true){
            int net;
            if(!(fin >> net)) throw runtime_error("Invalid net list in blocks section");
            if(net == -1) break;
            nets_of_b.push_back(net);
        }
        tmp.push_back({b, move(nets_of_b)});
    }
    vector<pair<int, int>> tmp_comm;
    while(true){
        int u;
        if(!(fin >> u)) break;
        if(u == -1) break;
        int v;
        if (!(fin >> v)) throw runtime_error("Incorrect community format (expected a pair).");
        tmp_comm.push_back({u, v});
    }
    //Map block numbers to indices
    Instance inst;
    {
        vector<int> ids;
        ids.reserve(tmp.size());
        for(auto &p: tmp) ids.push_back(p.first);
        sort(ids.begin(), ids.end());
        ids.erase(unique(ids.begin(), ids.end()), ids.end());

        inst.n = (int)ids.size();
        if (inst.n == 0) throw runtime_error("No blocks found.");
        if (inst.n % 2 != 0) throw runtime_error("Number of blocks must be even for equal-sized bi-partition.");
        inst.block_num = ids;
        for (int i = 0; i < inst.n; ++i) inst.num_to_index[ids[i]] = i;
    }
    unordered_map<int, int> net_to_index;
    auto get_net_index = [&](int net_id)->int{
        auto it = net_to_index.find(net_id);
        if(it != net_to_index.end()) return it->second;
        int k = (int)inst.nets.size();
        net_to_index[net_id] = k;
        inst.nets.push_back({});
        return k;
    };
    inst.block_nets.assign(inst.n, {});
    for(auto &p : tmp){
        int bi = inst.num_to_index[p.first];
        for(int net_id: p.second){
            int k = get_net_index(net_id);
            inst.nets[k].push_back(bi);
            inst.block_nets[bi].push_back(k);
        }
    }
    for(auto &pr: tmp_comm){
        auto it1 = inst.num_to_index.find(pr.first);
        auto it2 = inst.num_to_index.find(pr.second);
        if(it1 != inst.num_to_index.end() && it2 != inst.num_to_index.end()){
            int u = it1->second, v = it2->second;
            if(u != v) inst.communities.push_back({u, v});
        }
    }
    inst.partners.assign(inst.n, {});
    for (auto &pr : inst.communities) {
      inst.partners[pr.first].push_back(pr.second);
      inst.partners[pr.second].push_back(pr.first);
    }
    return inst;
}

int full_cost(const Instance &I, const vector<int> &side){ //cpmplete it
    int cross = 0;
    for(const auto &net: I.nets){
        bool has_L = false, has_R = false;
        for(int b: net){
            if(side[b] == 0) has_L = true;
            else if(side[b] == 1) has_R = true;
            if(has_L && has_R){
                cross++;
                break;
            }
        }
    }
    int comm = 0;
    for(const auto &pr: I.communities){
        int u = pr.first;
        int v = pr.second;
        if (side[u] != -1 && side[v] != -1 && side[u] != side[v]) comm++; //count communities
    }
    return cross + comm;
}

vector<int> degree_fanout(const Instance &I){
    vector<int> deg(I.n, 0);
    for(int b = 0; b < I.n; ++b) deg[b] = (int)I.block_nets[b].size();
    return deg;
}

vector<int> greedy_initial_assignment(const Instance &I, const vector<int> &order){
    vector<int> side(I.n, -1);
    int target = I.n/2;
    int left_count = 0, right_count = 0;

    vector<int> net_L(I.nets.size(),0), net_R(I.nets.size(),0);
    vector<char>crossed(I.nets.size(),0);
    int cross = 0;

    auto delta_cross = [&](int b, int s)->int{
        int add = 0;
        for(int k : I.block_nets[b]){
            int l = net_L[k], r = net_R[k];
            bool already = crossed[k];
            if(s == 0){
                if(!already && r > 0) add++;
            }else{
                if(!already && l > 0) add++;
            }
        }
        return add;
    };

    auto delta_comm = [&](int b, int s)->int {
      int add = 0;
      for (int v : I.partners[b]) {
        if (side[v] != -1 && side[v] != s) add++;
      }
      return add;
    };

    for(int b: order){
        int bestS = -1, bestDelta = INT_MAX;
        for(int s = 0; s <= 1; ++s){
            if(s == 0 && left_count >= target) continue;
            if(s == 1 && right_count >= target) continue;
            int dc = delta_cross(b, s) + delta_comm(b, s);
            if (dc < bestDelta) {
                bestDelta = dc;
                bestS = s;
            }
        }
        if(bestS == -1){
            bestS = (left_count < target) ? 0 : 1;
        }
        side[b] = bestS;
        if (bestS == 0) left_count++; else right_count++;
        for (int k : I.block_nets[b]) {

            if (bestS == 0) net_L[k]++; else net_R[k]++;
            if (!crossed[k] && net_L[k] > 0 && net_R[k] > 0) {
                crossed[k] = 1; cross++;
            }
        }
    }
    return side;
}

void build_order(const Instance &I, Initial_Place mode, vector<int> &order, vector<int> &seed_assign){
    vector<int> index(I.n);
    iota(index.begin(), index.end(), 0);

    if(mode == Initial_Place::base){
        sort(index.begin(), index.end(), [&](int a, int b){
            return I.block_num[a] < I.block_num[b];
        });
        order = index;
    }else{ //fanout in descending degree
        auto deg = degree_fanout(I);
        sort(index.begin(), index.end(), [&](int a, int b){
            if(deg[a] != deg[b]) return deg[a] > deg[b];
            return I.block_num[a] < I.block_num[b];
        });
        order = index;
    }
    seed_assign = greedy_initial_assignment(I, order);
}

VisTree G_VIS_TREE;

struct Branch_n_Bound{
    int nodes_visited = 0;
    const Instance &I;
    vector<int> order;
    vector<int> pref_side;
    int target;
    int cur_parent = -1;

    int best_cost = INT_MAX;
    vector<int> best_side;

    SearchState S;

    Branch_n_Bound(const Instance &inst, const vector<int> &ord, const vector<int> &seed)
        : I(inst), order(ord), pref_side(seed) {
        target = I.n/2;
        S.side.assign(I.n, -1);
        S.left_count = S.right_count = 0;
        S.target = target;
        S.net_left.assign(I.nets.size(), 0);
        S.net_right.assign(I.nets.size(), 0);
        S.net_crossed.assign(I.nets.size(), 0);
        S.crossing_count = 0;
        S.comm_mismatches = 0;
    }

    bool balance(int depth) const{
        int rem = I.n - depth;
        if (S.left_count > target || S.right_count > target) return false;
        if (S.left_count + rem < target) return false;
        if (S.right_count + rem < target) return false;
        return true;
    }

    int LB_forced() const {
        int lb = S.crossing_count + S.comm_mismatches;

        int Lrem = target - S.left_count;
        int Rrem = target - S.right_count;

        for (int k = 0; k < (int)I.nets.size(); ++k) {
            if (S.net_crossed[k]) continue;

            int aL = S.net_left[k];
            int aR = S.net_right[k];
            int sz = (int)I.nets[k].size();
            int r  = sz - aL - aR;
            if (r <= 0) continue;

            if (aL > 0 && aR == 0) {
                if (Lrem < r) lb += 1;
            }
            else if (aR > 0 && aL == 0) {
                if (Rrem < r) lb += 1;
            }
            else if (aL == 0 && aR == 0) {
                if (Lrem < sz && Rrem < sz) lb += 1;
            }
        }

        for (auto &pr : I.communities) {
            int u = pr.first, v = pr.second;
            int su = S.side[u];
            int sv = S.side[v];

            if (su != -1 && sv != -1) continue;

            if (su != -1 && sv == -1) {
                if ((su == 0 && Lrem == 0) || (su == 1 && Rrem == 0)) lb += 1;
            } else if (su == -1 && sv != -1) {
                if ((sv == 0 && Lrem == 0) || (sv == 1 && Rrem == 0)) lb += 1;
            } else {
                if (Lrem < 2 && Rrem < 2) lb += 1;
            }
        }

        return lb;
    }

    int LB() const{
        return S.crossing_count + S.comm_mismatches;
    }

    struct Delta{
        int b, s;
        vector<int> nets_newly_crossed;
        int comm_added = 0;
    };

    Delta apply(int b, int s){
        Delta d{b, s, {}, 0};
        for(int v: I.partners[b]){
            if(S.side[v] != -1 && S.side[v] != s){
                S.comm_mismatches++;
                d.comm_added++;
            }
        }
        for(int k: I.block_nets[b]){
            if(s == 0) S.net_left[k]++;
            else S.net_right[k]++;
            if(!S.net_crossed[k] && S.net_left[k] > 0 && S.net_right[k] > 0){
                S.net_crossed[k] = 1;
                S.crossing_count++;
                d.nets_newly_crossed.push_back(k);
            }
        }
        S.side[b] = s;
        if(s == 0) S.left_count++;
        else S.right_count++;
        return d;
    }

    void undo(const Delta &d){
        if (d.s == 0) S.left_count--; else S.right_count--;
        S.side[d.b] = -1;
        for (int k : I.block_nets[d.b]) {
          if (d.s == 0) S.net_left[k]--; else S.net_right[k]--;
        }
        for (int k : d.nets_newly_crossed) {
          S.net_crossed[k] = 0;
          S.crossing_count--;
        }
        S.comm_mismatches -= d.comm_added;
    }

    //dfs
    void dfs(int depth){
        nodes_visited++;
        if(depth == I.n){
            int cost = S.crossing_count + S.comm_mismatches;
            if(cost < best_cost){
                best_cost = cost;
                best_side = S.side;
            }
            return;
        }
        if(!balance(depth))return;

        int lb = LB_forced();
        if (lb >= best_cost) return;
        int b = order[depth];
        int first  = pref_side[b];
        int second = 1 - first;
        auto try_place = [&](int s) {
            if (s == 0 && S.left_count  >= target) return;
            if (s == 1 && S.right_count >= target) return;
            Delta d = apply(b, s);
            int my_index = G_VIS_TREE.add(cur_parent, I.block_num[b], s, depth);
            int saved_parent = cur_parent;
            cur_parent = my_index;
            if (balance(depth + 1)) {
                int lb2 = LB_forced();
                if (lb2 < best_cost) {
                    dfs(depth+1);
                }
            }
            cur_parent = saved_parent;
            undo(d);
        };
        try_place(first);
        if (depth > 0) {
            try_place(second);
        }
    }
};

int main(int argc, char **argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string file_path;
    Initial_Place init_place = Initial_Place::base;

    for(int i = 1; i < argc; i++){
        string a = argv[i];
        if(a == "-file_name"){
            if(i+1 >= argc) usage(argv[0]);
            file_path = argv[++i];
        }else if(a == "-init_place"){
            if(i+1 >= argc) usage(argv[0]);
            string v = argv[++i];
            if(v =="base") init_place = Initial_Place::base;
            else if(v =="fanout") init_place = Initial_Place::fanout;
            else usage(argv[0]);
        }else{
            usage(argv[0]);
        }
    }
    if(file_path.empty()) usage(argv[0]);

    Instance I = read_file(file_path);

    vector<int>order, seed;
    build_order(I, init_place, order, seed);
    int seed_cost = full_cost(I, seed);

    Branch_n_Bound solver(I, order, seed);
    solver.best_cost = seed_cost;
    solver.best_side = seed;

    G_VIS_TREE.nodes.clear();
    solver.cur_parent = -1;
    auto start = std::chrono::high_resolution_clock::now();
    solver.dfs(0);

    auto end = std::chrono::high_resolution_clock::now();
    double runtime_ms = std::chrono::duration<double, std::milli>(end - start).count();

    vector<int> left_ids, right_ids;
    for (int i = 0; i < I.n; ++i){
        if (solver.best_side[i] == 0) left_ids.push_back(I.block_num[i]);
        else right_ids.push_back(I.block_num[i]);
    }
    sort(left_ids.begin(), left_ids.end());
    sort(right_ids.begin(), right_ids.end());
    cout << "Optimal cost = " << solver.best_cost << "\n";
    cout << "Left  (" << left_ids.size()  << "): ";
    for (size_t i = 0; i < left_ids.size(); ++i) {
        if (i) cout << " ";
        cout << left_ids[i];
    }
    cout << "\nRight (" << right_ids.size() << "): ";
    for (size_t i = 0; i < right_ids.size(); ++i) {
        if (i) cout << " ";
        cout << right_ids[i];
    }
    cout << "\n";
    cout << "Nodes traversed = " << solver.nodes_visited << "\n";
    cout << "Runtime (ms) = " << runtime_ms << "\n";
    return 0;
}
