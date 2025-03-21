#include <iostream>
#include <cmath>
using namespace std;


double f(double x) {
    return x * x - 2 * x - 2 * cos(x);
}


int f_calls = 0;

double uniform_search(double a, double b, int N, int& counter) {
    double step = (b - a) / N;
    double min_x = a;
    double min_val = f(a);
    counter = 1; 

    for (int i = 1; i <= N; ++i) {
        double x = a + i * step;
        double current_val = f(x);
        counter++;
        if (current_val < min_val) {
            min_val = current_val;
            min_x = x;
        }
    }
    return min_x;
}


double dichotomy(double a, double b, double epsilon, double tol, int& counter) {
    counter = 0;
    while (b - a > tol) {
        double mid = (a + b) / 2;
        double x1 = mid - epsilon / 2;
        double x2 = mid + epsilon / 2;

        x1 = max(a, x1);
        x2 = min(b, x2);

        double f1 = f(x1);
        double f2 = f(x2);
        counter += 2;

        if (f1 < f2) {
            b = mid;
        } else {
            a = mid;
        }
    }
    return (a + b) / 2;
}


void print_theoretical_estimates(double epsilon) {
    double uniform_est = 1.0 / epsilon;
    double dichotomy_est = ceil(2 * log2(1.0 / epsilon) / log2(2));
    cout << "Теоретические оценки:\n";
    cout << "  Равномерный поиск: ~" << uniform_est << " вызовов\n";
    cout << "  Дихотомия: ~" << dichotomy_est << " вызовов\n\n";
}

int main() {
    double a = 0.0;
    double b = 1.0;
    double epsilons[] = {1e-1, 1e-2, 1e-3};

    for (double tol : epsilons) {
        cout << "Точность: " << tol << "\n";
        print_theoretical_estimates(tol);

 
        int uniform_calls;
        int N = (b - a) / tol;
        double uni_min = uniform_search(a, b, N, uniform_calls);
        cout << "Равномерный поиск\n";
        cout << "  Вызовов функции: " << uniform_calls << "\n";
        cout << "  x_min = " << uni_min << "\n  f(x_min) = " << f(uni_min) << "\n\n";

        int dichotomy_calls;
        double dich_min = dichotomy(a, b, tol/2, tol, dichotomy_calls);
        cout << "Дихотомия\n";
        cout << "  Вызовов функции: " << dichotomy_calls << "\n";
        cout << "  x_min = " << dich_min << "\n  f(x_min) = " << f(dich_min) << "\n\n";
    }

    return 0;
}