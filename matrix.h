#ifndef MATRIX_H
#define MATRIX_H
#pragma GCC optimize(2)

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

constexpr char error[] = "計算被中止，係由除數非預期為零。";
constexpr char no_solution[] = "計算被中止，係由該聯立一次方程式無解。";

template <typename type> type get_abs(type x) {return x > type(0) ? x : -x;}
template <typename type> inline void Swap(type &a, type &b) {type c = a; a = b, b = c;}
template <typename type> type find_min(type x, type y) {return x < y ? move(x) : move(y);}
template <typename type> type find_max(type x, type y) {return x > y ? move(x) : move(y);}
template <typename type> type find_gcd(type x, type y) {return !y ? move(x) : move(find_gcd(y, x % y));}
template <typename type> type find_lcm(type x, type y) {return x / move(find_gcd(x, get_abs(y))) * move(y);}

class Super;  // 任意位數整數之類別
class frac;   // 任意精度有理數之類別
template <typename type> class VectorX;         // 任意維向量空間之模板類別
template <typename type> class Matrix;          // 任意維矩陣空間之模板類別
template <typename type> class Square_Matrix;  // 任意維正方矩陣空間之模板類別

class Super {
private:
    bool is_nega;
    vector<int> dgts;  // 儲存每一位數字，最低位數值儲存於dgts[0]，類推之。
    void remove0s(void) {
        while (dgts.size() > 1 and !dgts.back()) dgts.pop_back();
    }
    void half(void) {
        for (int i = dgts.size() - 1, remain = 0; i >= 0; --i, remain *= 10)
            dgts[i] += remain, remain = dgts[i] % 2, dgts[i] >>= 1;
        remove0s();
    }

public:
    Super(void) : is_nega(false) {dgts.push_back(0);}
    Super(int num) {*this = to_string(num), is_nega = num < 0; remove0s();}
    Super(const string &num) {
        is_nega = (num[0] == '-');
        dgts.resize(num.size() - is_nega);
        for (int i = num.size() - 1; i >= is_nega; --i)
            dgts[num.size() - 1 - i] = num[i] - '0';
    }
    Super(const Super &num) : is_nega(num.is_nega), dgts(num.dgts) {}  
    Super(Super &&num) noexcept : is_nega(num.is_nega), dgts(move(num.dgts)) {}  

    bool operator ! (void) const {return dgts.size() == 1 and !dgts[0];}
    Super operator - (void) const {
        if (!*this) return *this;
        Super x = *this; x.is_nega = !x.is_nega; return move(x);
    }
    Super &operator = (const Super &num) {dgts = num.dgts, is_nega = num.is_nega; return *this;}
    Super &operator = (Super &&num) {
        is_nega = num.is_nega, dgts = move(num.dgts);
        return *this;
    }

    Super operator + (const Super &other) const {
        if (is_nega and !other.is_nega) return other - (-*this);
        if (!is_nega and other.is_nega) return *this - (-other);
        if (is_nega and other.is_nega) return -((-*this) + (-other));
        Super result;
        int maxLength = find_max(dgts.size(), other.dgts.size());
        result.dgts.resize(maxLength);
        int carry = 0;
        for (int i = 0; i < maxLength; ++i) {
            int sum = carry;
            if (i < dgts.size()) sum += dgts[i];
            if (i < other.dgts.size()) sum += other.dgts[i];
            result.dgts[i] = sum % 10;
            carry = sum / 10;
        }
        if (carry > 0) result.dgts.push_back(carry);
        return move(result);
    }
    Super operator - (const Super &other) const {
        if (is_nega and !other.is_nega) return -(-*this + other);
        if (!is_nega and other.is_nega) return *this + -other;
        if (is_nega and other.is_nega) return -other - (-*this);
        if (*this < other) return -(other - *this);
        Super result;
        result.dgts.resize(dgts.size());
        int borrow = 0;
        for (int i = 0; i < dgts.size(); ++i) {
            int diff = dgts[i] - borrow;
            if (i < other.dgts.size()) diff -= other.dgts[i];
            if (diff < 0)  diff += 10, borrow = 1;
            else borrow = 0;
            result.dgts[i] = diff;
        }
        result.remove0s();
        return move(result);
    }
    Super operator * (const Super &other) const {
        Super result;
        result.dgts.resize(dgts.size() + other.dgts.size(), 0);
        for (int i = 0; i < dgts.size(); ++i) {
            int carry = 0;
            for (int j = 0; j < other.dgts.size(); ++j) {
                int product = dgts[i] * other.dgts[j] + carry + result.dgts[i + j];
                result.dgts[i + j] = product % 10;
                carry = product / 10;
            }
            if (carry > 0) result.dgts[i + other.dgts.size()] += carry;
        }
        result.remove0s();
        result.is_nega = is_nega ^ other.is_nega;
        if (!result) result.is_nega = false;
        return move(result);
    }
    Super operator / (const Super &other) const {
        if (!other) {cerr << error << endl; return *this;}
        if (is_nega and !other.is_nega) return -(-*this / other);
        if (!is_nega and other.is_nega) return -(*this / -other);
        if (is_nega and other.is_nega) return -*this / -other;
        Super l = 0, r = *this, mid, remain;
        while (l <= r) {
            mid = l + r; mid.half();
            remain = *this - mid * other;
            if (remain >= 0 and remain < other) break;
            if (remain < 0) r = mid - 1;
            else l = mid + 1;
        }
        return move(mid);
    }
    Super operator % (const Super &other) const {
        if (!other) {cerr << error << endl; return -1;}
        return move(*this - *this / other * other);
    }

    Super &operator += (const Super &b) {*this = *this + b; return *this;}
    Super &operator *= (const Super &b) {*this = *this * b; return *this;}
    Super &operator -= (const Super &b) {*this = *this - b; return *this;}
    Super &operator /= (const Super &b) {*this = *this / b; return *this;}
    Super &operator %= (const Super &b) {*this = *this % b; return *this;}

    bool operator < (const Super &other) const {
        if (is_nega and !other.is_nega) return true;
        if (!is_nega and other.is_nega) return false;
        if (is_nega and other.is_nega) return -*this > -other;
        if (dgts.size() != other.dgts.size())
            return dgts.size() < other.dgts.size();
        for (int i = dgts.size() - 1; i >= 0; --i)
            if (dgts[i] != other.dgts[i])
                return dgts[i] < other.dgts[i];
        return false;
    }
    bool operator > (const Super &other) const {return other < *this;}

    bool operator <= (const Super &other) const {return !(other < *this);}
    bool operator >= (const Super &other) const {return !(*this < other);}
    bool operator != (const Super &other) const {return other < *this || *this < other;}
    bool operator == (const Super &other) const {return !(other < *this) && !(other > *this);}

    friend istream &operator >> (istream &in, Super &num) {
        string s; in >> s;
        num = s; return in;
    }
    friend ostream &operator << (ostream &os, const Super &num) {
        if (num.is_nega) os << '-';
        for (int i = num.dgts.size() - 1; i >= 0; --i)
            os << num.dgts[i];
        return os;
    }
};

#define type Super
class frac {
private:
    type nm, dm;
    inline void self_simplify(void) {
        type gcd = find_gcd(nm, get_abs(dm));  // 求最大公約數時第二引數不得為負，原因未明。
        nm /= gcd, dm /= gcd;
        if (dm < 0) dm = -dm, nm = -nm;
    }

public:
    frac(void) : nm(0), dm(1) {}
    frac(int num) : nm(num), dm(1) {}
    frac(type num) : nm(num), dm(1) {}
    frac(type nm, type dm) : nm(nm), dm(dm) {self_simplify();}
    operator type() {return nm / dm;}
    operator bool() const {return nm != 0;}
    frac(const frac &num) : nm(num.nm), dm(num.dm) {}
    frac(frac &&num) noexcept : nm(move(num.nm)), dm(move(num.dm)) {}

    frac operator + (const frac &other) {
        type lcm = find_lcm(this->dm, other.dm);
        return move(frac(lcm / this->dm * this->nm + lcm / other.dm * other.nm, lcm));
    }
    frac operator - (const frac &other) {
        type lcm = find_lcm(this->dm, other.dm);
        return move(frac(lcm / this->dm * this->nm - lcm / other.dm * other.nm, lcm));
    }
    frac operator * (const frac &other) {return move(frac(this->nm * other.nm, this->dm * other.dm));}
    frac operator / (const frac &other) {
        if (!other) {cerr << error << endl; return -1;}
        return move(frac(this->nm * other.dm, this->dm * other.nm));
    }

    bool operator ! (void) const {return !nm;}
    frac operator - (void) const {return move(frac(-nm, dm));}
    frac &operator = (type num) {this->nm = num, this->dm = 1; return *this;}
    frac &operator = (const frac &num) {this->nm = num.nm, this->dm = num.dm; return *this;}
    frac &operator = (frac &&num) noexcept {this->nm = move(num.nm), this->dm = move(num.dm); return *this;}

    frac &operator += (const frac &other) {
        type lcm = find_lcm(this->dm, other.dm);
        this->nm = lcm / this->dm * this->nm + lcm / other.dm * other.nm;
        this->dm = lcm; this->self_simplify();
        return *this;
    }
    frac &operator -= (const frac &other) {
        type lcm = find_lcm(this->dm, other.dm);
        this->nm = lcm / this->dm * this->nm - lcm / other.dm * other.nm;
        this->dm = lcm; this->self_simplify();
        return *this;
    }
    frac &operator *= (const frac &other) {
        this->nm *= other.nm, this->dm *= other.dm;
        this->self_simplify(); return *this;
    }
    frac &operator /= (const frac &other) {
        if (!other) {cerr << error << endl; return *this;}
        this->nm *= other.dm, this->dm *= other.nm;
        this->self_simplify(); return *this;
    }

    bool operator < (const frac &other) const {
        type lcm = find_lcm(this->dm, other.dm);
        return lcm / this->dm * this->nm < lcm / other.dm * other.nm;
    }
    bool operator > (const frac &other) const {return other < *this;}
    bool operator <= (const frac &other) const {return !(other < *this);}
    bool operator >= (const frac &other) const {return !(*this < other);}
    bool operator != (const frac &other) const {return other < *this || *this < other;}
    bool operator == (const frac &other) const {return !(other < *this) && !(other > *this);}

    friend istream &operator >> (istream &in, frac &num) {
        in >> num.nm; num.dm = 1;  // 目前僅支援以整數型態輸入數值
        return in;
    }
    friend ostream &operator << (ostream &os, const frac &num) {
        os << num.nm;
        if (num.dm != 1 and num) os << '/' << num.dm;
        return os;
    }
};
#undef type

template <typename type>
class VectorX {  // 由本類別所宣告變數參與之運算，若有維度不一之情事，槪以維度較少者為準，並捨棄維度較多者之多餘維度。
protected:
    size_t size;
    vector<type> vec;

public:
    VectorX(void) : size(0) {}
    VectorX(size_t size) : size(size) {vec.resize(size);}
    VectorX(vector<int> &num) : size(num.size()) {
        vec.resize(size);
        for (int i = 0; i < size; ++i) vec[i] = move(num[i]);
    }
    VectorX(vector<type> &num) : size(num.size()) {this->vec = num;}
    VectorX(VectorX<type> &num) : size(num.size) {this->vec = num.vec;}
    VectorX(VectorX<type> &&num) noexcept : size(num.size) {this->vec = move(num.vec);}

    size_t dimension(void) {return size;}
    type &operator [] (int index) {return vec[index - 1];}
    inline void set_size(size_t size) {this->size = size; vec.resize(size);}
    VectorX &operator = (const VectorX &other) {
        this->size = other.size;
        vec = other.vec; return *this;
    }
    VectorX &operator = (VectorX &&other) noexcept {
        this->size = other.size;
        vec = move(other.vec); return *this;
    }
    inline void exchange(int a, int b) {type c = (*this)[a]; (*this)[a] = (*this)[b]; (*this)[b] = c;}

    VectorX operator + (const VectorX &other) {
        VectorX<type> d(find_min(this->size, other.size));
        for (int i = 0; i < d.size; ++i)
            d.vec[i] = this->vec[i] + other.vec[i];
        return move(d);
    }
    VectorX operator - (const VectorX &other) {
        VectorX<type> d(find_min(this->size, other.size));
        for (int i = 0; i < d.size; ++i)
            d.vec[i] = this->vec[i] - other.vec[i];
        return move(d);
    }
    type operator * (const VectorX &other) {
        type result = 0;
        int limit = find_min(this->size, other.size);
        for (int i = 0; i < limit; ++i) result += vec[i] * other.vec[i];
        return move(result);
    }
    Matrix<type> operator * (const Matrix<type> &other) {
        return move(Matrix<type>(*this).transposed_matrix() * other);
    }
    VectorX operator ^ (const VectorX &other) {  // 過載為三維向量之外積
        VectorX<type> d(3);
        d[0] = vec[1] * other.vec[2] - vec[2] * other.vec[1];
        d[1] = vec[2] * other.vec[0] - vec[0] * other.vec[2];
        d[2] = vec[0] * other.vec[1] - vec[1] * other.vec[0];
        return move(d);
    }

    VectorX &operator += (const VectorX &other) {
        int limit = find_min(this->size, other.size);
        for (int i = 0; i < limit; ++i)
            vec[i] = this->vec[i] + other.vec[i];
        while (size > limit) vec.pop_back(), --size;
        return *this;
    }
    VectorX &operator -= (const VectorX &other) {
        int limit = find_min(this->size, other.size);
        for (int i = 0; i < limit; ++i)
            vec[i] = this->vec[i] - other.vec[i];
        while (size > limit) vec.pop_back(), --size;
        return *this;
    }
    VectorX &operator ^= (const VectorX &other) {*this = *this ^ other; return *this;}

    friend istream &operator >> (istream &in, VectorX &num) {
        for (int i = 0; i < num.size; ++i) in >> num.vec[i];
        return in;
    }
    friend ostream &operator << (ostream &os, const VectorX &num) {
        for (auto i : num.vec) os << i << " ";
        return os;
    }
};

#define SIZE 23  // 短絀時，視實際之需要調整陣列之規模
template <typename type>  // 宣告此類別時，原則上將模版型別設定為frac，以最大限度保留計算精度
class Matrix {
private:
    inline void line_copy(int org, int tar) {
        for (int i = 1; i <= m; ++i) mt[tar][i] = mt[org][i];
    }
    inline void row_coy(int org, int tar) {
        for (int i = 1; i <= n; ++i) mt[i][tar] = mt[i][org];
    }
    inline void to_exit(void) {cout << no_solution << endl; return Matrix(0);}

protected:
    size_t n, m;
    type mt[SIZE][SIZE];

    inline void display(void) {  // 此函式僅於調試程式時用之
        cerr << "n = " << n << ", m = " << m << endl;
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= m; ++j)
                cerr << mt[i][j] << (j < m ? " " : "");
            cerr << "\n";
        }
        cerr << "\n";
    }
    inline void combine(const Matrix &other) {
        for (int i = 1; i <= this->n; ++i)
            for (int j = 1; j <= other.m; ++j)
                this->mt[i][this->m + j] = other.mt[i][j];
        this->m += other.m;
    }
    inline void separate(size_t k) {
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= n; ++j)
                mt[i][j] = mt[i][j + n];
        this->m -= k;
    }
    inline bool solvable(void) {
        for (int i = 1; i <= n; ++i)
            for (int j = 1; j <= n; ++j)
                if (mt[i][j] != (type)(i == j)) return false;
        return true;
    }

    inline void line_multiply(int k, type p) {
        for (int i = 1; i <= m; ++i) mt[k][i] *= p;
    }
    inline void line_exchange(int a, int b) {
        line_copy(a, 0), line_copy(b, a), line_copy(0, b);
    }
    inline void multi_sum(int a, int b, type p) {
        for (int i = 1; i <= m; ++i) mt[a][i] += mt[b][i] * p;
    }

public:
    Matrix(void) : n(0), m(0) {}
    Matrix(size_t n, size_t m) : n(n), m(m) {
        if (typeid(type) != typeid(frac)) memset(mt, 0, sizeof(mt));
    }
    Matrix(const Matrix &num) : n(num.n), m(0) {combine(num);}
    static Matrix<type> identity(size_t n) {  // 使用此建構子於生成n階單位矩陣
        Matrix<type> p(n, n);
        if (typeid(type) != typeid(frac)) memset(p.mt, 0, sizeof(p.mt));
        for (int i = 1; i <= n; ++i) p.mt[i][i] = type(1);
        return p;
    }
    Matrix(const Matrix &a, const Matrix &b) {
        for (int i = 1; i <= a.n; ++i) {
            for (int j = 1; j <= a.m; ++j) this->mt[i][j] = a.mt[i][j];
            for (int j = 1; j <= b.m; ++j) this->mt[i][a.m + j] = b.mt[i][j];
        }
        this->n = a.n, this->m = a.m + b.m;
    }
    Matrix(VectorX<type> &other) : n(other.dimension()), m(1) {  // 從n維向量轉換為n維矩陣時，默認為豎置方向
        for (int i = 1; i <= n; ++i) this->mt[i][1] = other[i - 1];
    }
    inline void set_nm(size_t n, size_t m) {this->n = n, this->m = m;}

    Matrix transposed_matrix(void) {
        Matrix<type> d(m, n);
        for (int i = 1; i <= d.n; ++i)
            for (int j = 1; j <= d.m; ++j)
                d.mt[i][j] = this->mt[j][i];
        return move(d);
    }
    inline void matrix_simplify(void) {
        pair<int, int> ptr = make_pair(1, 1);
        do {
            int ori_pos = ptr.first;
            while (!mt[ptr.first][ptr.second] and ptr.first <= n) ++ptr.first;
            if (ptr.first > n) {
                ptr.second++, ptr.first = ori_pos;
                continue;
            }
            if (ori_pos != ptr.first)
                line_exchange(ori_pos, ptr.first), ptr.first = ori_pos;
            for (int i = 1; i <= n; ++i) {
                if (i == ptr.first) continue;
                multi_sum(i, ptr.first, -mt[i][ptr.second] / mt[ptr.first][ptr.second]);
            }
            line_multiply(ptr.first, type(1) / mt[ptr.first][ptr.second]);
            ptr.first++, ptr.second++;
        } while (ptr.first <= n && ptr.second <= m);
    }

    Matrix solve_liner_equations(void) {
        if (n + 1 != m) to_exit();
        Matrix<type> d = *this;
        d.matrix_simplify();
        if (!d.solvable()) to_exit();
        d.separate(n); return move(d);
    }
    Matrix solve_liner_equations(const Matrix &other) {
        if (other.m != 1) to_exit();
        Matrix<type> d = *this;
        d.combine(other), d.matrix_simplify();
        if (!d.solvable()) to_exit();
        d.separate(n); return move(d);
    }

    Matrix &operator = (const Matrix &other) {
        for (int i = 1; i <= other.n; i++)
            for (int j = 1; j <= other.m; j++)
                this->mt[i][j] = other.mt[i][j];
        this->n = other.n, this->m = other.m;
        return *this;
    }
    Matrix operator + (const Matrix &other) const {
        Matrix<type> d(n, m);
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= m; j++)
                d.mt[i][j] += other.mt[i][j];
        return move(d);
    }
    Matrix operator - (const Matrix &other) const {
        Matrix<type> d(n, m);
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= m; j++)
                d.mt[i][j] -= other.mt[i][j];
        return move(d);
    }
    Matrix operator * (const Matrix &other) {
        Matrix<type> d(this->n, other.m);
        for (int i = 1; i <= this->n; i++)
            for (int j = 1; j <= other.m; j++)
                for (int k = 1; k <= this->m; k++)
                    d.mt[i][j] += mt[i][k] * other.mt[k][j];
        return move(d);
    }
    Matrix operator / (const Square_Matrix<type> &other) {return *this * other.inversed_matrix();}

    Matrix &operator += (const Matrix &other) {
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= m; j++)
                mt[i][j] += other.mt[i][j];
        return *this;
    }
    Matrix &operator -= (const Matrix &other) {
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= m; j++)
                mt[i][j] -= other.mt[i][j];
        return *this;
    }
    Matrix &operator *= (const Matrix &other) {*this = *this * other; return *this;}
    Matrix &operator /= (const Matrix &other) {*this = *this / other; return *this;}

    friend istream &operator >> (istream &in, Matrix &num) {
        for (int i = 1; i <= num.n; ++i)
            for (int j = 1; j <= num.m; ++j)
                in >> num.mt[i][j];
        return in;
    }
    friend ostream &operator << (ostream &os, const Matrix &num) {
        for (int i = 1; i <= num.n; ++i) {
            for (int j = 1; j <= num.m; ++j) {
                os << num.mt[i][j];
                if (j < num.m) os << " ";
            }
            os << "\n";
        }
        return os;
    }
};

template <typename type>
class Square_Matrix : public Matrix<type> {
private:
    inline void inverse(void) {
        this->combine(move(Matrix<type>(this->n)));
        this->matrix_simplify(), this->separate(this->n);
    }
    type cofactor(int pos_i, int pos_j) {
        Square_Matrix<type> p(this->n - 1);
        for (int i = 1; i < this->n; ++i)
            for (int j = 1; j < this->m; ++j)
                p.mt[i][j] = this->mt[i < pos_i ? i : i + 1][j < pos_j ? j : j + 1];
        return p.det_calc(p.n) * (type)((pos_i + pos_j) % 2 ? -1 : 1);
    }
    type det_calc(size_t k) {
        if (k == 1) return this->mt[this->n][this->n];
        int minus = 1, pos = this->n - k + 1, find = pos;
        while (!this->mt[find][pos] and find <= this->n) ++find, minus = -1;
        if (find > this->n) return (type)0;
        if (pos != find) this->line_exchange(find, pos);
        for (int i = pos + 1; i <= this->n; ++i) 
            this->multi_sum(i, pos, -this->mt[i][pos] / this->mt[pos][pos]);
        return type(minus) * this->mt[pos][pos] * det_calc(k - 1);
    }
    inline void calc(type k) {
        for (int i = 1; i <= this->n; ++i)
            for (int j = 1; j <= this->m; ++j)
                this->mt[i][j] /= k;
        return;
    }

public:
    Square_Matrix(size_t n) {
        this->n = this->m = n;
        if (typeid(type) != typeid(frac)) memset(this->mt, 0, sizeof(this->mt));
    }
    static Square_Matrix<type> identity(size_t n) {  // 使用此建構子於生成n階單位矩陣
        Square_Matrix<type> p(n);
        if (typeid(type) != typeid(frac)) memset(p.mt, 0, sizeof(p.mt));
        for (int i = 1; i <= n; ++i) p.mt[i][i] = type(1);
        return p;
    }
    Square_Matrix operator * (const Square_Matrix &other) {
        Square_Matrix<type> d(this->n);
        for (int i = 1; i <= this->n; i++)
            for (int j = 1; j <= other.m; j++)
                for (int k = 1; k <= this->m; k++)
                    d.mt[i][j] += this->mt[i][k] * other.mt[k][j];
        return d;
    }

    type determinant(void) {
        Square_Matrix<type> p = *this;
        return move(p.det_calc(this->n));
    }

    Square_Matrix inversed_matrix(void) const {
        Square_Matrix<type> p = *this;
        p.inverse(); return move(p);
    }
    Square_Matrix cofator_matrix(void) {
        Square_Matrix<type> p(this->n);
        for (int i = 1; i <= this->n; ++i)
            for (int j = 1; j <= this->m; ++j)
                p.mt[i][j] = move(cofactor(j, i));
        return move(p);
    }
    Square_Matrix inv_mtrx_via_cfctr(void) {
        Square_Matrix<type> p(cofator_matrix());
        type det = determinant();
        if (!det) {cerr << error << endl; return 0;}
        p.calc(det); return move(p);
    }

    #define n this->n
    inline void LU_decomposition(VectorX<type> &B) {
        for (int k = 1; k <= n; ++k) {
            int tmp = k;
            for (int i = k + 1; i <= n; ++i)
                if (get_abs(this->mt[i][k]) > get_abs(this->mt[tmp][k])) tmp = i;
            if (!this->mt[tmp][k]) {cout << no_solution << endl; return;}
            this->line_exchange(k, tmp);
            Swap(B[k], B[tmp]);

            for (int i = k + 1; i <= n; ++i) {
                this->mt[i][k] /= this->mt[k][k];
                for (int j = k + 1; j <= n; ++j)
                    this->mt[i][j] -= this->mt[i][k] * this->mt[k][j];
            }
        }
    }
    VectorX<type> LU_solver(VectorX<type> B) {
        Square_Matrix<type> o = *this;
        o.LU_decomposition(B);
        cout << o << endl;

        for (int i = 1; i <= n; ++i)
            for (int j = 1; j < i; ++j)
                B[i] -= o.mt[i][j] * B[j];
        for (int i = n; i; --i) {
            for (int j = i + 1; j <= n; ++j)
                B[i] -= o.mt[i][j] * B[j];
            B[i] /= o.mt[i][i];
        }

        return move(B);
    }

    #define vctr for (int i = 1; i <= n; ++i)
    pair <Square_Matrix, Square_Matrix> GramQRdecomposition(void) {
        Square_Matrix<type> Q(n), R(n);
        for (int k = 1; k <= n; ++k) {
            vctr R.mt[k][k] += this->mt[i][k] * this->mt[i][k];
            R.mt[k][k] = sqrt(R.mt[k][k]);
            vctr Q.mt[i][k] = this->mt[i][k] / R.mt[k][k];

            for (int j = k + 1; j <= n; ++j) {
                vctr R.mt[k][j] += Q.mt[i][k] * this->mt[i][j];
                vctr this->mt[i][j] -= Q.mt[i][k] * R.mt[k][j];
            }
        }
        return make_pair(Q, R);
    }
    #undef vctr
    pair<Square_Matrix, Square_Matrix> HouseQRdecomposition() {
        Square_Matrix<type> Q = Square_Matrix<type>::identity(n);
        Square_Matrix<type> R = *this;

        for (int k = 1; k <= n; ++k) {
            vector<type> v(n + 1, 0);
            type norm_x = 0;
            for (int i = k; i <= n; ++i) norm_x += R.mt[i][k] * R.mt[i][k];
            norm_x = (R.mt[k][k] >= 0) ? -sqrt(norm_x) : sqrt(norm_x);
            v[k] = R.mt[k][k] - norm_x;
            for (int i = k + 1; i <= n; ++i) v[i] = R.mt[i][k];

            type v_norm = 0;
            for (int i = k; i <= n; ++i) v_norm += v[i] * v[i];
            v_norm = sqrt(v_norm);

            if (v_norm < 1e-10) continue;  // 無需向上翻轉的情形
            for (int i = k; i <= n; ++i) v[i] /= v_norm;

            Square_Matrix<type> H = Square_Matrix<type>::identity(n);
            for (int i = k; i <= n; ++i)
                for (int j = k; j <= n; ++j)
                    H.mt[i][j] -= 2 * v[i] * v[j];
 
            Q = Q * H;  // H矩陣本身軸對稱，無需再倒置
            R = H * R;
        }

        return make_pair(Q, R);
    }

    inline void LU_print(void) {
        cout << "L:" << endl;
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j < i; ++j) cout << this->mt[i][j] << " ";
            cout << "1 ";
            for (int j = i + 1; j <= n; ++j) cout << "0 ";
            cout << endl;
        }
        cout << "R:" << endl;
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j < i; ++j) cout << "0 ";
            for (int j = i; j <= n; ++j) cout << this->mt[i][j] << " ";
            cout << endl;
        }
        cout << endl;
    }
    inline friend void RU_print(pair<Square_Matrix<type>, Square_Matrix<type>> &p) {
        cout << "Q:" << endl << p.first;
        cout << "R:" << endl << p.second << endl;
    }
    #undef n
};
#undef SIZE

#endif MATRIX_H