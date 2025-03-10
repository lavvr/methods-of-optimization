#include <iostream>
#include <algorithm>
#include <climits>
#include <vector>
#include <utility>
#include <locale.h>

using namespace std;

const int M = 4;
const int N = 5;

struct cell {
    int val;   
    int x;       
    bool basis; 
};

void printPlan(cell cost[M][N]) {
    cout << "Матрица распределения:" << endl;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            cout << cost[i][j].x << "\t";
        }
        cout << endl;
    }
}

void computePotentials(cell cost[M][N], int u[M], int v[N]) {
    bool uKnown[M] = { false };
    bool vKnown[N] = { false };

    u[0] = 0;
    uKnown[0] = true;

    bool updated = true;
    while (updated) {
        updated = false;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if (cost[i][j].basis) {
                    if (uKnown[i] && !vKnown[j]) {
                        v[j] = cost[i][j].val - u[i];
                        vKnown[j] = true;
                        updated = true;
                    }
                    else if (vKnown[j] && !uKnown[i]) {
                        u[i] = cost[i][j].val - v[j];
                        uKnown[i] = true;
                        updated = true;
                    }
                }
            }
        }
    }
}

bool findCycle(int si, int sj, int ci, int cj, int lastDir, vector<pair<int, int>>& path, cell cost[M][N]) {
    if (path.size() >= 4 && ci == si && cj == sj) {
        return true;
    }
    if (lastDir == 0 || lastDir == 2) {
        for (int j = 0; j < N; j++) {
            if (j == cj) continue;
            if ((ci == si && j == sj) || cost[ci][j].basis) {
                bool inPath = false;
                for (auto& p : path) {
                    if (p.first == ci && p.second == j) {
                        inPath = true;
                        break;
                    }
                }
                if (inPath && !(ci == si && j == sj))
                    continue;
                path.push_back({ ci, j });
                if (findCycle(si, sj, ci, j, 1, path, cost))
                    return true;
                path.pop_back();
            }
        }
    }
    if (lastDir == 0 || lastDir == 1) {
        for (int i = 0; i < M; i++) {
            if (i == ci) continue;
            if ((i == si && cj == sj) || cost[i][cj].basis) {
                bool inPath = false;
                for (auto& p : path) {
                    if (p.first == i && p.second == cj) {
                        inPath = true;
                        break;
                    }
                }
                if (inPath && !(i == si && cj == sj))
                    continue;
                path.push_back({ i, cj });
                if (findCycle(si, sj, i, cj, 2, path, cost))
                    return true;
                path.pop_back();
            }
        }
    }
    return false;
}

void initialPlanDoublePreference(cell cost[M][N], int supply[], int demand[]) {
    bool rowMark[M][N]; 
    bool colMark[M][N]; 
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            rowMark[i][j] = false;
            colMark[i][j] = false;
        }
    }

    for (int i = 0; i < M; i++) {
        int minVal = INT_MAX;
        for (int j = 0; j < N; j++) {
            if (cost[i][j].val < minVal) {
                minVal = cost[i][j].val;
            }
        }
        for (int j = 0; j < N; j++) {
            if (cost[i][j].val == minVal) {
                rowMark[i][j] = true;
            }
        }
    }

    for (int j = 0; j < N; j++) {
        int minVal = INT_MAX;
        for (int i = 0; i < M; i++) {
            if (cost[i][j].val < minVal) {
                minVal = cost[i][j].val;
            }
        }
        for (int i = 0; i < M; i++) {
            if (cost[i][j].val == minVal) {
                colMark[i][j] = true;
            }
        }
    }

    struct PrefCell {
        int i, j, priority;
    };
    vector<PrefCell> cells;
    cells.reserve(M * N);

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            int p = 0;
            if (rowMark[i][j]) p++;
            if (colMark[i][j]) p++;
            cells.push_back({ i, j, p });
        }
    }

    sort(cells.begin(), cells.end(), [](const PrefCell& a, const PrefCell& b) {
        return a.priority > b.priority;
        });
    for (auto& c : cells) {
        int i = c.i;
        int j = c.j;
        if (supply[i] > 0 && demand[j] > 0) {
            int allocation = min(supply[i], demand[j]);
            cost[i][j].x = allocation;
            cost[i][j].basis = true;
            supply[i] -= allocation;
            demand[j] -= allocation;
        }
    }

    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            if (supply[i] > 0 && demand[j] > 0) {
                int allocation = min(supply[i], demand[j]);
                cost[i][j].x += allocation; 
                cost[i][j].basis = true;
                supply[i] -= allocation;
                demand[j] -= allocation;
            }
        }
    }
}


int main() {
    setlocale(LC_CTYPE, "Russian");


    int a[M] = { 14, 16, 11, 16 }; 
    int b[N] = { 3, 15, 18, 2, 19 }; 

    cell cost[M][N] = {
        { {1, 0, false}, {4, 0, false}, {10, 0, false}, {7, 0, false}, {9, 0, false} },
        { {9, 0, false}, {5, 0, false}, {17, 0, false}, {16, 0, false}, {12, 0, false} },
        { {8, 0, false}, {6, 0, false}, {14, 0, false}, {10, 0, false}, {17, 0, false} },
        { {11, 0, false}, {8, 0, false}, {20, 0, false}, {10, 0, false}, {15, 0, false} }
    };

    int supplyDP[M], demandDP[N];
    for (int i = 0; i < M; i++) {
        supplyDP[i] = a[i];
    }
    for (int j = 0; j < N; j++) {
        demandDP[j] = b[j];
    }

    initialPlanDoublePreference(cost, supplyDP, demandDP);

    cout << "метод двойного предпочтения:" << endl;
    printPlan(cost);

    while (true) {
        int u[M], v[N];
        computePotentials(cost, u, v);

        int maxDelta = 0;
        int enterI = -1, enterJ = -1;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if (!cost[i][j].basis) {
                    int delta = u[i] + v[j] - cost[i][j].val;
                    if (delta > maxDelta) {
                        maxDelta = delta;
                        enterI = i;
                        enterJ = j;
                    }
                }
            }
        }

        if (maxDelta <= 0)
            break;

        vector<pair<int, int>> cycle;
        cycle.push_back({ enterI, enterJ });
        bool cycleFound = findCycle(enterI, enterJ, enterI, enterJ, 0, cycle, cost);
        if (!cycleFound) {
            cout << "Цикл не найден." << endl;
            break;
        }

        int theta = INT_MAX;
        for (size_t k = 1; k < cycle.size() - 1; k += 2) {
            int i = cycle[k].first;
            int j = cycle[k].second;
            theta = min(theta, cost[i][j].x);
        }
        for (size_t k = 0; k < cycle.size() - 1; k++) {
            int i = cycle[k].first;
            int j = cycle[k].second;
            if (k % 2 == 0)
                cost[i][j].x += theta;
            else
                cost[i][j].x -= theta;
        }
        cost[enterI][enterJ].basis = true;
        for (size_t k = 1; k < cycle.size() - 1; k += 2) {
            int i = cycle[k].first;
            int j = cycle[k].second;
            if (cost[i][j].x == 0)
                cost[i][j].basis = false;
        }
    }

    cout << "\nОптимальный опорный план:" << endl;
    printPlan(cost);

    int totalCost = 0;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
            totalCost += cost[i][j].x * cost[i][j].val;
    cout << "Общая стоимость перевозок: " << totalCost << endl;

    int ai, bj, shm;
    cout << "\nВведите номер поставщика (от 1 до " << M << "): ";
    cin >> ai;
    cout << "Введите номер потребителя (от 1 до " << N << "): ";
    cin >> bj;
    cout << "Введите минимальное количество груза (shm): ";
    cin >> shm;

    int supplier = ai - 1;
    int consumer = bj - 1;

    if (a[supplier] < shm || b[consumer] < shm) {
        cout << "Ошибка: недостаточно запаса или спроса для выполнения фиксированной поставки." << endl;
        return 1;
    }

    a[supplier] -= shm;
    b[consumer] -= shm;

    cell newCost[M][N];
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            newCost[i][j].val = cost[i][j].val;
            newCost[i][j].x = 0;
            newCost[i][j].basis = false;
        }
    }

    int supply2[M], demand2[N];
    for (int i = 0; i < M; i++)
        supply2[i] = a[i];
    for (int j = 0; j < N; j++)
        demand2[j] = b[j];

    while (true) {
        int minCost = INT_MAX;
        int minI = -1, minJ = -1;
        for (int i = 0; i < M; i++) {
            if (supply2[i] > 0) {
                for (int j = 0; j < N; j++) {
                    if (demand2[j] > 0 && newCost[i][j].val <= minCost) {
                        minCost = newCost[i][j].val;
                        minI = i;
                        minJ = j;
                    }
                }
            }
        }
        if (minI == -1 || minJ == -1)
            break;
        int allocation = min(supply2[minI], demand2[minJ]);
        newCost[minI][minJ].x = allocation;
        newCost[minI][minJ].basis = true;
        supply2[minI] -= allocation;
        demand2[minJ] -= allocation;
    }

    while (true) {
        int u[M], v[N];
        computePotentials(newCost, u, v);
        int maxDelta = 0;
        int enterI = -1, enterJ = -1;
        for (int i = 0; i < M; i++) {
            for (int j = 0; j < N; j++) {
                if (!newCost[i][j].basis) {
                    int delta = u[i] + v[j] - newCost[i][j].val;
                    if (delta > maxDelta) {
                        maxDelta = delta;
                        enterI = i;
                        enterJ = j;
                    }
                }
            }
        }
        if (maxDelta <= 0)
            break;
        vector<pair<int, int>> cycle;
        cycle.push_back({ enterI, enterJ });
        bool cycleFound = findCycle(enterI, enterJ, enterI, enterJ, 0, cycle, newCost);
        if (!cycleFound) {
            cout << "Ошибка: цикл не найден при улучшении плана преобразованной задачи." << endl;
            break;
        }
        int theta = INT_MAX;
        for (size_t k = 1; k < cycle.size() - 1; k += 2) {
            int i = cycle[k].first;
            int j = cycle[k].second;
            theta = min(theta, newCost[i][j].x);
        }
        for (size_t k = 0; k < cycle.size() - 1; k++) {
            int i = cycle[k].first;
            int j = cycle[k].second;
            if (k % 2 == 0)
                newCost[i][j].x += theta;
            else
                newCost[i][j].x -= theta;
        }
        newCost[enterI][enterJ].basis = true;
        for (size_t k = 1; k < cycle.size() - 1; k += 2) {
            int i = cycle[k].first;
            int j = cycle[k].second;
            if (newCost[i][j].x == 0)
                newCost[i][j].basis = false;
        }
    }

    newCost[supplier][consumer].x += shm;
    newCost[supplier][consumer].basis = true;

    cout << "\nОптимальный опорный план для задачи с фиксированной поставкой:" << endl;
    printPlan(newCost);

    int totalCost2 = 0;
    for (int i = 0; i < M; i++)
        for (int j = 0; j < N; j++)
            totalCost2 += newCost[i][j].x * newCost[i][j].val;
    cout << "Общая стоимость перевозок (с учетом фиксированной поставки): " << totalCost2 << endl;

    return 0;
}
