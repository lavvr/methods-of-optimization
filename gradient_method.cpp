#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <functional>

double f( double X, double Y )
{
  return X * X + Y * Y - log(X + Y);
}

std::vector<double> grad( double X, double Y )
{
  return {2 * X - 1.0 / (X + Y), 2 * Y - 1.0 / (X + Y)};
}

std::vector<std::vector<double>> hessian( double X, double Y )
{
  std::vector<double> r1(2), r2(2);

  r1[0] = 2 + 1.0 / ((X + Y) * (X + Y));
  r1[1] = 1.0 / ((X + Y) * (X + Y));

  r2[0] = r1[1];
  r2[1] = r1[0];
  return {r1, r2};
}

double norm( const std::vector<double> &V )
{
  double n = 0;

  for (double v : V)
    n += v * v;
  return sqrt(n);
}

std::vector<double> gradStep( double X0, double Y0, double a0, double eps, double gamma, double delta )
{
  double X = X0, Y = Y0;
  double a = a0;
  int k = 0;
  std::ofstream log_file("log_convergence.txt");

  while (1)
  {
    std::vector<double> gr = grad(X, Y);
    double gradNorm = norm(gr);
    if (gradNorm <= delta)
    {
      printf("%d\n", k);
      return {X, Y};
    }

    log_file << k++ << " " << log(gradNorm) << '\n';

    while (f(X - a * gr[0], Y - a * gr[1]) - f(X, Y) > -eps * a * gradNorm)
      a *= gamma;

    X = X + a * (-gr[0]);
    Y = Y + a * (-gr[1]);
  }
}

std::vector<double> solveSystem( std::vector<std::vector<double>> A, const std::vector<double> &b )
{
  int n = A.size();
  for (int i = 0; i < n; i++)
    A[i].push_back(b[i]);

  for (int i = 0; i < n; i++)
  {
    int pivot = i;
    for (int j = i + 1; j < n; j++)
      if (std::fabs(A[j][i]) > std::fabs(A[pivot][i]))
        pivot = j;
    /// !!
    if (fabs(A[pivot][i]) < 0.0001)
      return {};

    std::swap(A[i], A[pivot]);

    for (int j = i + 1; j <= n; j++)
      A[i][j] /= A[i][i];

    for (int j = 0; j < n; j++)
      if (j != i)
      {
        double factor = A[j][i];

        for (int k = i; k <= n; k++)
          A[j][k] -= factor * A[i][k];
      }
  }

  std::vector<double> res(n);

  for (int i = 0; i < n; i++)
    res[i] = A[i][n];

  return res;
}

double dichotomy( std::function<double(double)> f, double a, double b, double precision)
{
  double eps = 0.001 * (b - a);
  double mid = (a + b) / 2.0;
  double f_left, f_right;

  if (b - a < precision)
    return (a + b) / 2;

  f_left = f(mid - eps);
  f_right = f(mid + eps);

  if (f_left <= f_right)
    return dichotomy(f, a, mid + eps, precision);
  return dichotomy(f, mid - eps, b, precision);
}

std::vector<double> gradSecondOrderDescent( double X0, double Y0, double a0, double delta )
{
  double X = X0, Y = Y0;
  double a = a0;
  int k = 0;

  while (1)
  {
    std::vector<double> gr = grad(X, Y);
    double gradNorm = norm(gr);
    if (gradNorm <= delta)
    {
      printf("%d\n", k);
      return {X, Y};
    }
    k++;

    std::vector<std::vector<double>> H = hessian(X, Y);
    std::vector<double> p = solveSystem(H, gr);

    auto phi = [&](double alpha) {
      std::vector<double> x_new(2);
      x_new[0] = X, x_new[1] = Y;

      for (size_t i = 0; i < x_new.size(); ++i)
        x_new[i] -= alpha * p[i];
      return f(x_new[0], x_new[1]);
    };

    double alpha = dichotomy(phi, 0.0, 1.0, delta);

    X = X - a * (p[0]);
    Y = Y - a * (p[1]);
  }
}

std::vector<double> gradCoordDescent( double X0, double Y0, double delta )
{
  double X = X0, Y = Y0;
  double change = 1e100;
  int k = 0;

  while (1)
  {
    if (change <= delta)
    {
      printf("%d\n", k);
      return {X, Y};
    }
    k++;

    change = 0;
    for (int i = 0; i < 2; i++)
    {
      auto phi = [&](double alpha)
      {
        std::vector<double> x_new(2);
        x_new[0] = X, x_new[1] = Y;

        x_new[i] += alpha;
        return f(x_new[0], x_new[1]);
      };

      double alpha_opt = dichotomy(phi, -1.0, 1.0, delta);

      if (i == 0)
        X += alpha_opt;
      else
        Y += alpha_opt;
      change += abs(alpha_opt);
    }
  }
}

int main( void )
{
  std::vector<double> mi = gradStep(1, 1, 0.3, 1e-7, 0.5, 0.00001);

  printf("%.5f %.5f %.5f\n", mi[0], mi[1], f(mi[0], mi[1]));

  mi = gradSecondOrderDescent(1, 1, 0.5, 0.00001);
  printf("%.5f %.5f %.5f\n", mi[0], mi[1], f(mi[0], mi[1]));

  mi = gradCoordDescent(1, 1, 0.00001);
  printf("%.5f %.5f %.5f\n", mi[0], mi[1], f(mi[0], mi[1]));
}