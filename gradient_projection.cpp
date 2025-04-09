#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

vector<double> matvec(const vector<vector<double>>& A, const vector<double>& x) {
    vector<double> res(A.size(), 0.0);
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < x.size(); ++j) {
            res[i] += A[i][j] * x[j];
        }
    }
    return res;
}

vector<vector<double>> matmul(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    vector<vector<double>> res(A.size(), vector<double>(B[0].size(), 0.0));
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t k = 0; k < B.size(); ++k) {
            for (size_t j = 0; j < B[0].size(); ++j) {
                res[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return res;
}

vector<vector<double>> transpose(const vector<vector<double>>& A) {
    vector<vector<double>> res(A[0].size(), vector<double>(A.size()));
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[0].size(); ++j) {
            res[j][i] = A[i][j];
        }
    }
    return res;
}

vector<vector<double>> inverse(vector<vector<double>> A) {
    int n = A.size();
    vector<vector<double>> res(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) res[i][i] = 1.0;

    for (int col = 0; col < n; ++col) {
        int pivot = col;
        for (int row = col + 1; row < n; ++row) {
            if (abs(A[row][col]) > abs(A[pivot][col])) pivot = row;
        }
        swap(A[col], A[pivot]);
        swap(res[col], res[pivot]);

        double div = A[col][col];
        for (int j = 0; j < n; ++j) {
            A[col][j] /= div;
            res[col][j] /= div;
        }

        for (int row = 0; row < n; ++row) {
            if (row != col && A[row][col] != 0.0) {
                double factor = A[row][col];
                for (int j = 0; j < n; ++j) {
                    A[row][j] -= A[col][j] * factor;
                    res[row][j] -= res[col][j] * factor;
                }
            }
        }
    }
    return res;
}

vector<double> compute_lagrange_multipliers(
    const vector<vector<double>>& A,
    const vector<double>& grad
) {
    auto AT = transpose(A);
    auto AAT = matmul(A, AT);
    auto AAT_inv = inverse(AAT);
    auto neggrad = grad;
    for (auto& si : neggrad) si = -si;
    return matvec(AAT_inv, matvec(A, neggrad));
}

vector<double> rosen_modified(
    const vector<double>& x0,
    vector<vector<double>>& A,
    vector<double>& b,
    const vector<vector<double>>& F2,
    const vector<double>& g2,
    double tol = 1e-6,
    int max_iter = 1000
) {
    vector<double> x = x0;

    for (int iter = 0; iter < max_iter; ++iter) {
        vector<double> grad(x.size());
        grad[0] = 2 * x[0];// + x[1] + x[3];
        grad[1] = 2 * x[1];// + x[0] + x[2];
        grad[2] = 2 * x[2];// + x[1] + x[3];
        grad[3] = 2 * x[3];// + x[2] + x[0];

        auto AT = transpose(A);
        auto AAT_inv = inverse(matmul(A, AT));
        auto P_term = matmul(AT, matmul(AAT_inv, A));
        vector<vector<double>> P(A[0].size(), vector<double>(A[0].size(), 0.0));
        for (size_t i = 0; i < P.size(); ++i) {
            P[i][i] = 1.0;
            for (size_t j = 0; j < P[0].size(); ++j) {
                P[i][j] -= P_term[i][j];
            }
        }

        vector<double> S = matvec(P, grad);
        for (auto& si : S) si = -si;

        bool all_zero = true;
        for (double si : S) {
            if (abs(si) > tol) {
                all_zero = false;
                break;
            }
        }

        if (!all_zero) {
            double alpha = 0.01;
            for (size_t i = 0; i < F2.size(); ++i) {
                double FiS = 0.0;
                for (size_t j = 0; j < S.size(); ++j) {
                    FiS += F2[i][j] * S[j];
                }
                if (FiS > 0) {
                    double Fix = 0.0;
                    for (size_t j = 0; j < x.size(); ++j) {
                        Fix += F2[i][j] * x[j];
                    }
                    alpha = min(alpha, (g2[i] - Fix) / FiS);
                }
            }

            for (size_t i = 0; i < x.size(); ++i) {
                x[i] += alpha * S[i];
            }
        } else {
            bool has_negative = false;
            vector<double> W = compute_lagrange_multipliers(A, grad);
            for (auto& wi : W) wi = -wi;
            for (size_t i = 0; i < W.size(); ++i) {
                if (W[i] < -tol) {
                    A.erase(A.begin() + i);
                    b.erase(b.begin() + i);
                    has_negative = true;
                    break;
                }
            }

            if (!has_negative) {
                cout << "Solution on iter: " << iter << endl;
                return x;
            }
        }
    }
    return x;
}

int main() {
    // x1^2 + x2^2 + x3^2 + x4^2
    vector<vector<double>> A = {
        {2, 3, 1, 4},  // 2*x1 + 3*x2 + 1*x3 + 4*x4 = 10
        {1, -2, 1, -1},   // x1 - 2 * x2 + x3 - x4 = 1
        {1, 1, 2, 1},   // x1 + x2 + 2*x3 + x4 <= 5
        {2, 3, 1, 2},   // 2*x1 + 3*x2 + 1*x3 + 2*x4 <= 8
    };
    vector<double> b = {10, 1, 5, 8};

    vector<vector<double>> F2 = { 
        {1, 2, 3, 1},   // x1 + 2*x2 + 3*x3 + x4 <= 7
        {3, 1, 2, 2}    // 3*x1 + 1*x2 + 2*x3 + 2*x4 <= 9
    };
    vector<double> g2 = { 7, 9};

    vector<double> x0 = {1, 1, 1, 1};

    auto x_opt = rosen_modified(x0, A, b, F2, g2);

    cout << "Opt x: ";
    for (double xi : x_opt) cout << xi << " ";
    cout << endl;

    return 0;
}