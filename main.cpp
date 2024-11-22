//#define _GLIBCXX_DEBUG
#include "bits/stdc++.h"
#include "optimization.h" // http://acm.math.spbu.ru/~sk1/algo/lib

using namespace std;
#define sz(a)(int)a.size()
#define all(a) a.begin(),a.end()
#define err(...) fprintf(stderr, "%.2f : ", 1. * clock() / CLOCKS_PER_SEC), fprintf(stderr, __VA_ARGS__), fflush(stderr)
namespace BwT {
    namespace Encode {
        vector<int> ans;

        vector<int> sufmas(const vector<int> &text) {
            int n = sz(text);
            vector<int> cnt(max(400, n)), pn(n), cn(n);
            vector<int> p(n), c(n);
            vector<pair<char, int>> cr;
            for (int i = 0; i < n; ++i)
                ++cnt[text[i]];
            for (int i = 1; i < 400; ++i)
                cnt[i] += cnt[i - 1];
            for (int i = 0; i < n; ++i)
                p[--cnt[text[i]]] = i;
            c[p[0]] = 0;
            int cl = 1;
            for (int i = 1; i < n; ++i) {
                if (text[p[i]] != text[p[i - 1]]) ++cl;
                c[p[i]] = cl - 1;
            }
            ++cl;

            for (int step = 0; (1 << step) < n; ++step) {
                cnt.assign(cl, 0);
                for (int i = 0; i < n; ++i) {
                    pn[i] = p[i] - (1 << step);
                    if (pn[i] < 0)
                        pn[i] += n;
                }
                for (int i = 0; i < n; ++i)
                    ++cnt[c[i]];
                for (int i = 1; i < cl; ++i)
                    cnt[i] += cnt[i - 1];
                for (int i = n - 1; i >= 0; --i)
                    p[--cnt[c[pn[i]]]] = pn[i];
                cn[p[0]] = 0;
                cl = 1;
                for (int i = 1; i < n; ++i) {
                    int m1 = p[i - 1] + (1 << step);
                    if (m1 >= n)
                        m1 -= n;
                    int m2 = p[i] + (1 << step);
                    if (m2 >= n)
                        m2 -= n;
                    if (c[p[i]] != c[p[i - 1]] || c[m1] != c[m2])
                        ++cl;
                    cn[p[i]] = cl - 1;
                }
                ++cl;
                c.swap(cn);
            }
            return p;
        }

        vector<int> encode(const vector<int> &text) {
            vector<int> suf = sufmas(text);
            int x = suf[0];
            while (x > 0) {
                ans.push_back(x % 10);
                x /= 10;
            }
            while (sz(ans) < 8)
                ans.push_back(0);
            reverse(all(ans));
            for (auto &y:suf)
                ans.push_back(text[(y + sz(text) - 1) % sz(text)]);
            return ans;
        }
    }
    namespace Decode {
        vector<int> ans;

        vector<int> decode(const vector<int> &text) {
            int p0 = 0;
            for (int i = 0; i < 8; ++i) {
                p0 *= 10;
                p0 += text[i];
            }
            vector<int> r;
            vector<pair<int, int>> l;
            for (int i = 8; i < sz(text); ++i) {
                r.push_back(text[i]);
//                cout << text[i] << " ";
                l.emplace_back(text[i], i - 8);
            }
            sort(all(l));
            int cur = 0;
            vector<int> t;
            while (sz(t) < sz(r)) {
                t.push_back(l[cur].first);
                cur = l[cur].second;
            }
            cur = (sz(t) - p0) % sz(t);
            while (sz(ans) < sz(t)) {
                ans.push_back(t[cur]);
                (++cur) %= sz(t);
            }
            return ans;
        }
    }
}
namespace MTF {
    namespace Encode {
        vector<int> ans;

        vector<int> encode(const vector<int> &text) {

            return ans;
        }
    }
    namespace Decode {
        vector<int> ans;

        vector<int> decode(const vector<int> &text) {

            return ans;
        }
    }
}
namespace RLE {
    namespace Encode {
        vector<int> ans;

        vector<int> encode(const vector<int> &text) {
            for (int i = 0; i < sz(text);) {
                ans.push_back((text[i] << 1) + isdigit(text[i]));
                ++i;
                int cnt = 0;
                while (i < sz(text) && text[i] == ans.back() >> 1)
                    ++i, ++cnt;
                if (cnt > 0) {
                    string t = to_string(cnt);
                    for (auto &x:t)
                        ans.push_back(x << 1);
                }
            }

            return ans;
        }
    }
    namespace Decode {
        vector<int> ans;

        vector<int> decode(const vector<int> &text) {
            for (int i = 0; i < sz(text);) {
                ans.push_back(text[i] >> 1);
                ++i;
                string cnt;
                while (i < sz(text) && isdigit(text[i] >> 1) && !(text[i] & 1)) {
                    cnt.push_back(text[i] >> 1);
                    ++i;
                }
                if (cnt.empty())
                    continue;
                int e = stoi(cnt);
                for (int j = 0; j < e; ++j)
                    ans.push_back(ans.back());
            }
            return ans;
        }
    }
}

namespace Huffman {
    namespace Encode {
        vector<int> a;
        vector<int> ans;

        char digit16(int x) {
            return x < 10 ? '0' + x : 'A' + x - 10;
        }

        void addBits(int x, int k) {
            assert(x < (1 << k));
            for (int i = 0; i < k; ++i)
                a.push_back((x >> i) & 1);
        }

        void flush() {
            for (int i = 0; i < sz(a); i += 4) {
                int k = 0;
                for (int j = i; j < i + 4; ++j)
                    k += a[j] * (1 << (j - i));
                ans.push_back(digit16(k));
            }
        }


        struct Trie {
            struct item {
                int i, s, real, l, r;

                item(int i = -1, int s = 0, int real = -1, int l = -1, int r = -1) : i(i), s(s),
                                                                                     real(real), l(l), r(r) {};
            };

            struct cmp {
                bool operator()(item u, item v) const {
                    return make_pair(u.s, u.i) < make_pair(v.s, v.i);
                }
            };

            map<int, vector<int>> mp;
            vector<item> t;
            set<item, cmp> st;
            map<int, int> pos;
            int cur = 0;

            Trie() {}


            void add(int x, int s) {
                pos[x] = cur;
                t.emplace_back(pos[x], s, x);
                st.insert(t.back());
                ++cur;
            }

            void algo() {
                while (sz(st) > 1) {
                    item l = *st.begin();
                    st.erase(st.begin());
                    item r = *st.begin();
                    st.erase(st.begin());
                    t.emplace_back(item(cur, l.s + r.s, -1, l.i, r.i));
                    st.insert(t.back());
                    ++cur;
                }
            }

            void dfs(int v, vector<int> &vec) {
                if (t[v].real != -1) {
                    a.push_back(1);
                    addBits(t[v].real, 8);
                    mp[t[v].real] = vec;
                } else {
                    a.push_back(0);
                    vec.push_back(0);
                    dfs(t[v].l, vec);
                    vec.back() = 1;
                    dfs(t[v].r, vec);
                    vec.pop_back();
                }
            }

            void encode(const vector<int> &text) {
                for (auto x:text)
                    for (auto y:mp[x])
                        a.push_back(y);
            }
        };

        void build(const vector<int> &text) {
            Trie trie;
            map<int, int> cnt;
            for (auto x:text)
                cnt[x]++;
            for (auto[x, c]:cnt)
                trie.add(x, c);
            trie.algo();
            vector<int> empty;
            trie.dfs(trie.cur - 1, empty);
            trie.encode(text);
        }

        vector<int> encode(const vector<int> &text) {
            addBits(sz(text), 28);
            build(text);
            flush();
            return ans;
        }
    }

    namespace Decode {
        vector<int> a;
        vector<int> ans;
        int cur = 0;

        int digit16(int x) {
            return isdigit(x) ? x - '0' : x - 'A' + 10;
        }

        void prepare(const vector<int> &text) {
            for (auto &c : text) {
                int u = digit16(c);
                for (int i = 0; i < 4; ++i)
                    a.push_back((u >> i) & 1);
            }
        }

        int read(int k) {
            int res = 0;
            for (int i = 0; i < k; ++i, ++cur)
                res += (1 << i) * a[cur];
            return res;
        }

        struct Trie {
            int last = -1;

            struct item {
                int i, c, l, r;

                item(int i = -1, int c = -1, int l = -1, int r = -1) : i(i), c(c), l(l), r(r) {};
            };

            vector<item> t;

            Trie() {};

            int add(int c = -1, int l = -1, int r = -1) {
                t.emplace_back(++last, c, l, r);
                return last;
            }

            void dfs(int v) {
                if (a[cur] == 1) {
                    ++cur;
                    t[v].c = read(8);
                } else {
                    ++cur;
                    t[v].l = add();
                    dfs(t[v].l);
                    t[v].r = add();
                    dfs(t[v].r);
                }
            }

            int root;

            void algo() {
                root = add();
                dfs(root);
            }

            void recover() {
                int v = root;
                while (cur < sz(a) && t[v].c == -1) {
                    if (a[cur] == 0)
                        v = t[v].l;
                    else
                        v = t[v].r;
                    ++cur;
                }
                if (t[v].c != -1)
                    ans.push_back(t[v].c);
            }
        };


        vector<int> decode(const vector<int> &text) {
            err("decode len=%ld\n", text.size());
            prepare(text);
            int n = read(28);
            Trie trie;
            trie.algo();
            while (cur < sz(a) && sz(ans) < n)
                trie.recover();
            return ans;
        }
    }
}

namespace Coder {
    void flush(const vector<int> &res) {
        for (auto &x:res)
            writeChar(x);
    }

    void code(vector<int> &text) {
        flush(Huffman::Encode::encode(RLE::Encode::encode(BwT::Encode::encode(text))));
//        flush(Huffman::Encode::encode(text));
    }
};
namespace Decoder {
    void flush(const vector<int> &res) {
        for (auto &x:res)
            writeChar(x);
    }

    void decode(vector<int> &text) {
        err("decode len=%ld\n", text.size());
        flush(BwT::Decode::decode(RLE::Decode::decode(Huffman::Decode::decode(text))));
//        flush(BwT::Decode::decode(RLE::Decode::decode(Huffman::Decode::decode(text))));
    }
};

int main() {
    freopen("output.txt", "w", stdout);
    vector<int> text;
    while (!isEof())
        text.push_back(getChar());

    auto isEncoded = [&]() {
        for (auto c : text)
            if (!('0' <= c && c <= '9') && !('A' <= c && c <= 'F'))
                return 0;
        return 1;
    };
    (!isEncoded() ? Coder::code : Decoder::decode)(text);
}