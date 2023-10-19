// Laba3_base.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <cmath>
#include <array>
#include <type_traits>

/** класс для работы с трехдиагональной матрицей **/
template<typename Type, unsigned int N>
class ThreeDiagonalMatrix {
public:
    std::array<Type, N> A; // 1 диагональ
    std::array<Type, N> B; // 2 диагональ
    std::array<Type, N> C; // 3 диагональ
    //конструктор класса
    ThreeDiagonalMatrix(const std::array<Type,N>& A, const std::array<Type,N>& B, const std::array<Type, N>& C) {
        this->A = A;
        this->B = B;
        this->C = C;
    }
};

template<typename numeratorType, typename denominatorType>
using DivisType = decltype(std::declval<numeratorType>() / std::declval<denominatorType>());

template<typename Type>
using DiffType = decltype(std::declval<Type>() - std::declval<Type>());

//решение системы с трехдиагональной матрицей методом прогонки
template<typename mType, typename cType, unsigned int N>
std::array<mType, N> solve (const ThreeDiagonalMatrix<mType, N>& matrix, const std::array<cType, N>& column) 
{
    std::array<mType, N> P, Q, X;// коэффициенты прогонки и решение системы
    // считаем коээфициенты прогонки
    P[1] = -matrix.C[0] / matrix.B[0];
    Q[1] = column[0] / matrix.B[0];
    for (unsigned int k = 2; k < N; k++) {
        P[k] = -matrix.C[k - 1] / (matrix.A[k - 1] * P[k - 1] + matrix.B[k - 1]);
        Q[k] = (column[k - 1] - matrix.A[k - 1] * Q[k - 1]) / (matrix.A[k - 1] * P[k - 1] + matrix.B[k - 1]);
    }
    // находим решение системы
    X[N - 1] = (column[N - 1] - matrix.A[N - 1] * Q[N - 1]) / (matrix.A[N - 1] * P[N - 1] + matrix.B[N - 1]);
    for (unsigned int k = 2; k < N+1; k++) {
        X[N - k] = X[N - k + 1] * P[N - k + 1] + Q[ N - k + 1];
    }
    return X;
};

//класс для построения кубического сплайна
template<typename xType, typename yType, unsigned int N>
class CubicSpline {
private:
    using DeltaXType = DiffType<xType>;
    using AType = DivisType<DeltaXType, DiffType<DeltaXType>>;
    using DerivType = DivisType<DiffType<yType>, DeltaXType>;
    using Deriv2Type = DivisType<DiffType<DerivType>, DeltaXType>;
    
    Deriv2Type a; // значение второй производной функции на правом конце
    Deriv2Type b; // значение второй производной функции на левом конце
    std::array<xType, N + 1> points; // узлы
    std::array<yType, N + 1> f; // значения функции в узлах
    std::array<DeltaXType, N> h; // длины отрезков разбиения
    std::array<AType, N - 1> A; // 1 диагональ трехдиагональной матрицы
    std::array<AType, N - 1> B; // 2 диагональ трехдиагональной матрицы
    std::array<AType, N - 1> C; // 3 диагональ трехдиагональной матрицы
    std::array<Deriv2Type, N-1> D; // столбец свободных членов
public:
    std::array<std::array<AType, 4>, N> Coefs; // коэффициенты кубических сплайнов на отрезках разбиения
    //конструктор класса
    CubicSpline(const std::array<xType, N + 1>& points,  // Значения x
        const std::array<yType, N + 1>& values,  // значения y
        const Deriv2Type& first,  // значение для левой второй производной
        const Deriv2Type& second) {
        this->b = first;
        this->a = second;
        this->points = points;
        this->f = values;
        //считаем длины отрезков
        for (unsigned int k = 0; k < N; k++) {
            this->h[k] = points[k + 1] - points[k];
        }
        // задаем диагонали матрицы и столбец свободных членов
        DeltaXType sumh;
        sumh = h[0] + h[1];
        DerivType f1 = (f[2] - f[1]) / h[1];
        DerivType f0 = (f[1] - f[0]) / h[0];
        this->A[0] = 0;
        this->B[0] = 2;
        this->C[0] = h[1] / sumh;
        this->D[0] = 3 * (f1 - f0) / sumh - b * h[0] / 2 / sumh;
        for (unsigned int k = 1; k < N - 2; k++) {
            sumh = h[k] + h[k + 1];
            this->A[k]= h[k] / sumh;
            this->B[k] = 2;
            this->C[k] = h[k + 1] / sumh;
            f0 = f1;
            f1 = (f[k + 2] - f[k + 1]) / h[k + 1];
            this->D[k] = 3 * (f1 - f0) / sumh;
        }
        sumh = h[N - 2] + h[N - 1];
        this->A[N - 2] = h[N - 2] / sumh;
        this->B[N - 2] = 2;
        this->C[N - 2] = 0;
        f0 = f1;
        f1 = (f[N] - f[N - 1]) / h[N - 1];
        this->D[N - 2] = 3 * (f1 - f0) / sumh - 2 * a * h[N - 1] / sumh;
    };
    // подсчет коэффициентов кубических сплайнов
    void calculateCoefs() {
        ThreeDiagonalMatrix < AType, N - 1>  M(A, B, C);
        std::array<AType, N - 1> X;
        X = solve(M, D);
        Coefs[0][2] = b / 2;
        for (unsigned int k = 0; k < N - 1; k++) {
            Coefs[k + 1][2] = X[k];
        }
        for (unsigned int k = 0; k < N - 1; k++) {
            Coefs[k][1] = (f[k + 1] - f[k]) / h[k] - 2 * Coefs[k][2] * h[k] / 3 - Coefs[k + 1][2] * h[k] / 3;
            Coefs[k][3] = (Coefs[k + 1][2] - Coefs[k][2]) / h[k] / 3;
            Coefs[k][0] = f[k];
        }
        Coefs[N - 1][1] = (f[N] - f[N - 1]) / h[N - 1] - 2 * Coefs[N - 1][2] * h[N - 1] / 3 - a * h[N - 1] / 6;
        Coefs[N - 1][3] = (a - 2 * Coefs[N - 1][2]) / h[N - 1] / 6;
        Coefs[N - 1][0] = f[N - 1];
    };

    yType interpolate(const xType& x) const noexcept {
        yType y;
        int k = 0;
        //проверяем, на какой отрезок попала точка
        if (x == 0)
            k++;
        while (x > points[k])
        {
            k++;
        }
        y = Coefs[k-1][0] + Coefs[k-1][1] * (x - points[k-1]) + Coefs[k-1][2] * (x - points[k-1]) * (x - points[k-1]) + Coefs[k-1][3] * (x - points[k-1]) * (x - points[k-1]) * (x - points[k-1]);
        return y;
    };
    yType evaluate(const xType& l) const noexcept // l - правый конец отрезка интерполяции
    {
        yType y, dis;
        DeltaXType d;
        xType x;
        dis = 0;
        // считаем логарифм ошибки интерполяции
        for (unsigned int j = 0; j < 1000; j++)
        {
            d = l / 1000; //расстояние между точками
            x = d * j; // текущая точка
            y = abs(interpolate(x) - exp(x));
            if (y > dis)
                dis = y;
        };
        std::cout << "log_2 (N) = " << log(N) / log(2) << "    log_2 (discrepency) = " << log(dis)/ log(2) << ";";
        std::cout << "\n\n";
        return log(dis);
    };
};

// задание равномерного расположения N узлов на отрезке [0,l] и вычисление значений e^x в узлах
template<typename xType, typename yType, unsigned int N>
void UniformPoints(xType& endInterval, std::array<xType, N>& points, std::array<yType, N>& values)
{
    double l;
    l = endInterval / (N-1);
    values[0] = 1;
    for (unsigned int i = 1; i < N; i++)
    {
        points[i] = l * i;
        values[i] = exp(points[i]);
    };
}


int main()
{
    int N;
    double a, b, l, max, min, d;
    std::array<double, 6> delta = {}; //значения логарифма ошибки интерполяции в зависимости от числа узлов (естественный сплайн);
    std::array<double, 6> delta1 = {}; //значения логарифма ошибки интерполяции в зависимости от числа узлов (сплайн с реальными значениями второй производной на краях отрезка);
    l = 10;// длина отрезка
    
    //Естественный сплайн
    std::cout << "Natural spline" << '\n';
    //задаем значения второй производной на концах отрезка = 0
    a = 0;
    b = 0;
       
    N = 5;// число отрезков разбиения
    std::cout << "N=" << N << '\n';
    //задаем массивы точек разбиения (равномерно) и значений функции в них
    std::array<double, 6> values = {};
    std::array<double, 6> points = {};
    UniformPoints(l, points, values);
    //задаем кубический сплайн
    CubicSpline < double,double, 5> S(points, values, a,b);
    //считаем коэффициенты сплайна
    S.calculateCoefs();
    // считаем и печатаем логарифм ошибки интерполяции
    delta[0] = S.evaluate(l);

    N = 10;// число отрезков разбиения
    std::cout << "N=" << N << '\n';
    //задаем массивы точек разбиения (равномерно) и значений функции в них
    std::array<double, 11> values1 = {};
    std::array<double, 11> points1 = {};
    UniformPoints(l, points1, values1);
    //задаем кубический сплайн
    CubicSpline < double, double, 10> S1(points1, values1, a, b);
    //считаем коэффициенты сплайна
    S1.calculateCoefs();
    // считаем и печатаем логарифм ошибки интерполяции
    delta[1] = S1.evaluate(l);

    N = 20;// число отрезков разбиения
    std::cout << "N=" << N << '\n';
    //задаем массивы точек разбиения (равномерно) и значений функции в них
    std::array<double, 21> values2 = {};
    std::array<double, 21> points2 = {};
    UniformPoints(l, points2, values2);
    //задаем кубический сплайн
    CubicSpline < double, double, 20> S2(points2, values2, a, b);
    //считаем коэффициенты сплайна
    S2.calculateCoefs();
    // считаем и печатаем логарифм ошибки интерполяции
    delta[2] = S2.evaluate(l);

    N = 40;// число отрезков разбиения
    std::cout << "N=" << N << '\n';
    //задаем массивы точек разбиения (равномерно) и значений функции в них
    std::array<double, 41> values3 = {};
    std::array<double, 41> points3 = {};
    UniformPoints(l, points3, values3);
    //задаем кубический сплайн
    CubicSpline < double, double, 40> S3(points3, values3, a, b);
    //считаем коэффициенты сплайна
    S3.calculateCoefs();
    // считаем и печатаем логарифм ошибки интерполяции
    delta[3] = S3.evaluate(l);

    N = 80;// число отрезков разбиения
    std::cout << "N=" << N << '\n';
    //задаем массивы точек разбиения (равномерно) и значений функции в них
    std::array<double, 81> values4 = {};
    std::array<double, 81> points4 = {};
    UniformPoints(l, points4, values4);
    //задаем кубический сплайн
    CubicSpline < double, double, 80> S4(points4, values4, a, b);
    //считаем коэффициенты сплайна
    S4.calculateCoefs();
    // считаем и печатаем логарифм ошибки интерполяции
    delta[4] = S4.evaluate(l);

    N = 160;// число отрезков разбиения
    std::cout << "N=" << N << '\n';
    //задаем массивы точек разбиения (равномерно) и значений функции в них
    std::array<double, 161> values5 = {};
    std::array<double, 161> points5 = {};
    UniformPoints(l, points5, values5);
    //задаем кубический сплайн
    CubicSpline < double, double, 160> S5(points5, values5, a, b);
    //считаем коэффициенты сплайна
    S5.calculateCoefs();
    // считаем и печатаем логарифм ошибки интерполяции
    delta[5] = S5.evaluate(l);

    // оцениваем тангенсы углов наклона линейных участков зависимости логарифма ошибки от логарифма числа узлов
    std::cout << "\n";
    max = delta[1] - delta[0];
    min = delta[1] - delta[0];
    for (int i = 1; i < 5; i++) {
        d = (delta[i+1] - delta[i]);
        if (d > max) {
            max = d;
        }
        if (d < min) {
            min = d;
        }
    }
    // печаем оценку тангенсов
    std::cout << "max tg phi = " << max << ";\n";
    std::cout << "min tg phi = " << min << ";\n";

    //Сплайн со истинными значениями второй производной на концах отрезка
    std::cout << '\n' << "Real spline" << '\n';
    //задаем значения второй производной на концах отрезка = 0
    a = 1;
    b = exp(10);

    N = 5;// число отрезков разбиения
    std::cout << "N=" << N << '\n';
    //задаем массивы точек разбиения (равномерно) и значений функции в них
    std::array<double, 6> values6 = {};
    std::array<double, 6> points6 = {};
    UniformPoints(l, points6, values6);
    //задаем кубический сплайн
    CubicSpline < double, double, 5> S6(points6, values6, a, b);
    //считаем коэффициенты сплайна
    S6.calculateCoefs();
    // считаем и печатаем логарифм ошибки интерполяции
    delta1[0] = S6.evaluate(l);

    N = 10;// число отрезков разбиения
    std::cout << "N=" << N << '\n';
    //задаем массивы точек разбиения (равномерно) и значений функции в них
    std::array<double, 11> values7 = {};
    std::array<double, 11> points7 = {};
    UniformPoints(l, points7, values7);
    //задаем кубический сплайн
    CubicSpline < double, double, 10> S7(points7, values7, a, b);
    //считаем коэффициенты сплайна
    S7.calculateCoefs();
    // считаем и печатаем логарифм ошибки интерполяции
    delta1[1] = S7.evaluate(l);

    N = 20;// число отрезков разбиения
    std::cout << "N=" << N << '\n';
    //задаем массивы точек разбиения (равномерно) и значений функции в них
    std::array<double, 21> values8 = {};
    std::array<double, 21> points8 = {};
    UniformPoints(l, points8, values8);
    //задаем кубический сплайн
    CubicSpline < double, double, 20> S8(points8, values8, a, b);
    //считаем коэффициенты сплайна
    S8.calculateCoefs();
    // считаем и печатаем логарифм ошибки интерполяции
    delta1[2] = S8.evaluate(l);

    N = 40;// число отрезков разбиения
    std::cout << "N=" << N << '\n';
    //задаем массивы точек разбиения (равномерно) и значений функции в них
    std::array<double, 41> values9 = {};
    std::array<double, 41> points9 = {};
    UniformPoints(l, points9, values9);
    //задаем кубический сплайн
    CubicSpline < double, double, 40> S9(points9, values9, a, b);
    //считаем коэффициенты сплайна
    S9.calculateCoefs();
    // считаем и печатаем логарифм ошибки интерполяции
    delta1[3] = S9.evaluate(l);

    N = 80;// число отрезков разбиения
    std::cout << "N=" << N << '\n';
    //задаем массивы точек разбиения (равномерно) и значений функции в них
    std::array<double, 81> values10 = {};
    std::array<double, 81> points10 = {};
    UniformPoints(l, points10, values10);
    //задаем кубический сплайн
    CubicSpline < double, double, 80> S10(points10, values10, a, b);
    //считаем коэффициенты сплайна
    S10.calculateCoefs();
    // считаем и печатаем логарифм ошибки интерполяции
    delta1[4] = S10.evaluate(l);

    N = 160;// число отрезков разбиения
    std::cout << "N=" << N << '\n';
    //задаем массивы точек разбиения (равномерно) и значений функции в них
    std::array<double, 161> values11 = {};
    std::array<double, 161> points11 = {};
    UniformPoints(l, points11, values11);
    //задаем кубический сплайн
    CubicSpline < double, double, 160> S11(points11, values11, a, b);
    //считаем коэффициенты сплайна
    S11.calculateCoefs();
    // считаем и печатаем логарифм ошибки интерполяции
    delta1[5] = S11.evaluate(l);

    // оцениваем тангенсы углов наклона линейных участков зависимости логарифма ошибки от логарифма числа узлов
    std::cout << "\n";
    max = delta1[1] - delta1[0];
    min = delta1[1] - delta1[0];
    for (int i = 1; i < 5; i++) {
        d = (delta1[i + 1] - delta1[i]);
        if (d > max) {
            max = d;
        }
        if (d < min) {
            min = d;
        }
    }
    // печаем оценку тангенсов
    std::cout << "max tg phi = " << max << ";\n";
    std::cout << "min tg phi = " << min << ";\n";
}


