#include <bits/stdc++.h>
using namespace std;
using ll = long long;
#define all(a) (a).begin(), (a).end()
#define rall(a) (a).rbegin(), (a).rend()
#define gen(a, b) uniform_int_distribution<int>(a, b)(rnd)
mt19937 rnd(chrono::steady_clock::now().time_since_epoch().count());
#ifndef LOCAL
#define endl '\n'
#endif

signed main() {
    ios::sync_with_stdio(0);
    cin.tie(0);
#ifdef LOCAL
    freopen("test.txt", "r", stdin);
#endif
    $END$
}


#include <bits/stdc++.h>

using namespace std;

template <typename A, typename B>
string to_string(pair<A, B> p);

template <typename A, typename B, typename C>
string to_string(tuple<A, B, C> p);

template <typename A, typename B, typename C, typename D>
string to_string(tuple<A, B, C, D> p);

string to_string(const string& s) {
    return '"' + s + '"';
}

string to_string(const char* s) {
    return to_string((string) s);
}

string to_string(bool b) {
    return (b ? "true" : "false");
}

string to_string(vector<bool> v) {
    bool first = true;
    string res = "{";
    for (int i = 0; i < static_cast<int>(v.size()); i++) {
        if (!first) {
            res += ", ";
        }
        first = false;
        res += to_string(v[i]);
    }
    res += "}";
    return res;
}

template <size_t N>
string to_string(bitset<N> v) {
    string res = "";
    for (size_t i = 0; i < N; i++) {
        res += static_cast<char>('0' + v[i]);
    }
    return res;
}

template <typename A>
string to_string(A v) {
    bool first = true;
    string res = "{";
    for (const auto &x : v) {
        if (!first) {
            res += ", ";
        }
        first = false;
        res += to_string(x);
    }
    res += "}";
    return res;
}

template <typename A, typename B>
string to_string(pair<A, B> p) {
    return "(" + to_string(p.first) + ", " + to_string(p.second) + ")";
}

template <typename A, typename B, typename C>
string to_string(tuple<A, B, C> p) {
    return "(" + to_string(get<0>(p)) + ", " + to_string(get<1>(p)) + ", " + to_string(get<2>(p)) + ")";
}

template <typename A, typename B, typename C, typename D>
string to_string(tuple<A, B, C, D> p) {
    return "(" + to_string(get<0>(p)) + ", " + to_string(get<1>(p)) + ", " + to_string(get<2>(p)) + ", " + to_string(get<3>(p)) + ")";
}

void debug_out() { cerr << endl; }

template <typename Head, typename... Tail>
void debug_out(Head H, Tail... T) {
    cerr << " " << to_string(H);
    debug_out(T...);
}

#define debug(...) cerr << "[" << #__VA_ARGS__ << "]:", debug_out(__VA_ARGS__)


struct SegmentTreeMass
{
    struct Node
    {
        int sum = 0;
        int add = 0;

        Node(int Sum = 0, int Add = 0)
        {
            sum = Sum;
            add = Add;
        }
    };

    Node neutral = Node();
    int n;
    vector<Node> tree;

    SegmentTreeMass(vector<int>& start)
    {
        n = start.size();
        tree.resize(4 * n + 228);
        init(start);
    }

    SegmentTreeMass(int N)
    {
        n = N;
        tree.resize(4 * n + 228);
        vector<int> start(n, 0);
        init(start);
    }

    Node merge(Node n1, Node n2)
    {
        return Node(n1.sum + n2.sum);
    }

    void fix(int v, int l, int r)
    {
        tree[v] = merge(tree[v * 2 + 1], tree[v * 2 + 2]);
    }

    void apply(int v, int l, int r, int val)
    {
        tree[v].add += val;
        tree[v].sum += val * (r - l);
    }

    void init(int v, int l, int r, vector<int>& start)
    {
        if (l + 1 == r)
        {
            tree[v] = Node(start[l], 0);
            return;
        }

        int m = (r + l) / 2;
        init(v * 2 + 1, l, m, start);
        init(v * 2 + 2, m, r, start);
        fix(v, l, r);
    }

    void init(vector<int>& start)
    {
        init(0, 0, n, start);
    }

    void push(int v, int l, int r)
    {
        int m = (r + l) / 2;

        apply(v * 2 + 1, l, m, tree[v].add);
        apply(v * 2 + 2, m, r, tree[v].add);
        tree[v].add = 0;
    }

    // [l: r)
    Node query(int v, int l, int r, int ql, int qr)
    {
        if (ql <= l && r <= qr)
        {
            return tree[v];
        }
        if (r <= ql || qr <= l)
        {
            return neutral;
        }

        int m = (r + l) / 2;
        push(v, l, r);
        return merge(
                query(v * 2 + 1, l, m, ql, qr),
                query(v * 2 + 2, m, r, ql, qr));

    }

    Node query(int ql, int qr)
    {
        return query(0, 0, n, ql, qr);
    }

    void add(int v, int l, int r, int ql, int qr, int val)
    {
        if (ql <= l && r <= qr)
        {
            apply(v, l, r, val);
            return;
        }
        if (r <= ql || qr <= l)
        {
            return;
        }

        int m = (r + l) / 2;
        push(v, l, r);
        add(v * 2 + 1, l, m, ql, qr, val);
        add(v * 2 + 2, m, r, ql, qr, val);
        fix(v, l, r);
    }

    void add(int ql, int qr, int val)
    {
        add(0, 0, n, ql, qr, val);
    }

    void add(int idx, int val)
    {
        add(0, 0, n, idx, idx + 1, val);
    }
};

struct sufmass {
    string s;
    vector<int> P;
    vector<int> lcp;
    void SORT(vector<int>& p, vector<int>& c) {
        int n = c.size();
        vector<int> cnt(n), pos(n);
        for (auto x : c)
            cnt[x]++;
        for (int i = 1; i < n; i++)
            pos[i] = pos[i - 1] + cnt[i - 1];
        vector<int> np(n);
        for (auto i : p) {
            auto& j = pos[c[i]];
            np[j++] = i;
        }
        p = np;
    }
    void get_lcp() {
        int n = s.size();
        lcp.assign(n, {});
        vector<int> pos(n);
        for (int i = 0; i < n; ++i) {
            pos[P[i]] = i;
        }
        int l = 0;
        for (int ii = 0; ii < n - 1; ++ii) {
            if (pos[ii] == 0) {
                l--;
                l = max<int>(0, l);
                continue;
            }
            int I = ii;
            int J = P[pos[ii] - 1];
            while (s[(I + l) % n] == s[(J + l) % n]) {
                l++;
            }
            lcp[pos[ii]] = l;
            l--;
            l = max<int>(0, l);
        }
    }
    void build(string S) {
        s = S;
        s += "$";
        int n = s.size();
        vector<pair<int, int>> sp;
        for (int i = 0; i < n; i++)
            sp.push_back({s[i] - 'a', i});
        sort(all(sp));
        vector<int> p(n), c(n);
        for (int i = 0; i < n; i++)
            p[i] = sp[i].second;
        for (int i = 1; i < n; i++)
            c[p[i]] = c[p[i - 1]] + (sp[i].first != sp[i - 1].first);
        for (int k = 1; k < n; k <<= 1) {
            for (int i = 0; i < n; i++)
                p[i] = (p[i] - k + n) % n;
            SORT(p, c);
            vector<int> C(n);
            pair<int, int> prev = {c[p[0]], c[(p[0] + k) % n]};
            for (int i = 1; i < n; i++) {
                pair<int, int> cur = {c[p[i]], c[(p[i] + k) % n]};
                C[p[i]] = C[p[i - 1]] + (cur != prev);
                prev = cur;
            }
            c = C;
        }
        P = p;
        get_lcp();
    }
};
from os import system
system("g++ stupid.cpp -o stupid") 
system("g++ sol.cpp -o sol") 
system("g++ gen.cpp -o gen") 
while True:
    system("gen.exe >test.txt")
    system("stupid.exe < test.txt > st.out")
    system("sol.exe < test.txt > sol.out")
    if open("st.out").read() == open("sol.out").read():
        print("OK")
    else:
        print("WA")
        print(open("test.txt").read())
        print("sol: ")
        print(open("sol.out").read(), end="")
        print("stupid: ")
        print(open("st.out").read(), end="")
        break

int gcd(int a, int b, int &x, int &y) {
    if (a == 0) {
        x = 0;
        y = 1;
        return b;
    }
    int x1, y1;
    int d = gcd(b % a, a, x1, y1);
    x = y1 - (b / a) * x1;
    y = x1;
    return d;
}
