#include <iostream>
#include <cmath>
#include <array>

template<typename RealType, unsigned int N> //неопределенные коэффициенты
struct DerivativeCoef {
    RealType centralCoef;
    std::array<RealType, N> otherCoefs;
};

template<typename RealType, unsigned int N, typename L> // вычисление неопределенных коэффициентов для k-ой производной  
DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points, L k) noexcept
{
    //заполнение расширенной матрицы системы для поиска неопределенных коэффициентов 
    std::array<std::array<double, N + 2>, N + 1> A;
    A = {};// все зануляем
    A[0][0] = 1;
    A[k][N + 1] = 1;
    for (unsigned int l = 2; l < k + 1; l++)//вычисляем k!
    {
        A[k][N + 1] *= l;
    }

    for (unsigned int i = 1; i < N + 1; i++)
    {
        for (unsigned int j = 1; j < N + 1; j++)
        {
            A[i][j] = pow(points[j - 1], i);
            A[0][j] = 1;
        }
    }

    //метод Гаусса-Жордана, 0 столбец - ок уже
    double a;
    unsigned int s;
    for (unsigned int j = 1; j < N + 1; j++)
    {
        s = j;// находим ненулевой элемент в j-ом столбце
        while (A[s][j] == 0) {
            s++;
        }
        a = A[s][j];
        for (unsigned int l = j; l < N + 2; l++)//делаем 1 в A[s][j]
        {
            A[s][l] = A[s][l] / a;
        }
        for (unsigned int i1 = s + 1; i1 < N + 1; i1++)// делаем нули под A[s][j]
        {
            a = A[i1][j];
            for (unsigned int j1 = j; j1 < N + 2; j1++) {
                A[i1][j1] = A[i1][j1] - A[s][j1] * a;
            }
        }
        for (unsigned int l = j; l < N + 2; l++)//переставляем s строку и j строку
        {
            a = A[s][l];
            A[s][l] = A[j][l];
            A[j][l] = a;
        }
        for (unsigned int i1 = 0; i1 < j; i1++)// делаем нули выше A[j][j]
        {
            a = A[i1][j];
            for (unsigned int j1 = j; j1 < N + 2; j1++) {
                A[i1][j1] = A[i1][j1] - A[s][j1] * a;
            }
        }
    }

    DerivativeCoef<double, N> C; // внесем найденные неопределенные коэффициенты в структуру
    C.centralCoef = A[0][N + 1];

    for (unsigned int i = 1; i < N + 1; i++) {
        C.otherCoefs[i - 1] = A[i][N + 1];
    }
    return C;
}
template<typename RealType, unsigned int N, typename L> // вычисление всего, что требуется в лабе, для N узлов и k-ой производной
std::array<RealType, 16> calcDerivative(const std::array<RealType, N>& points, L k, std::array<double, 16>& steps) noexcept
{
    DerivativeCoef<double, N> A;
    std::array<double, 16> f; //значения первой производной в зависимости от шага
    std::array<double, 16> delta; //значения логарифма ошибки для первой производной в зависимости от шага
    std::array<double, 16> tg; //значения тангенса углов наклона линейных участков зависимости логарифма ошибки от логарифма шага
    double max, min;

    //поиск k-ой производной для e^x с числом узлов N=3
    std::cout << "\n\n\n" << k << "-order derivative" << "\n\n";
    std::cout << "\n" << "N = " << N << "\n\n";
    std::cout << "points = {" << points[0];
    for (int i = 1; i < N; i++) {
        std::cout << ", " << points[i];
    }
    std::cout << "}\n\n";
    //вычисляем неопределенные коэффициенты
    A = calcDerivativeCoef(points, k);
    //напечатаем неопределенные коэффициенты
    std::cout << "\n" << "centralCoef = " << A.centralCoef << '\n';
    std::cout << "otherCoefs = {";
    for (unsigned int i = 0; i < N - 1; i++) {
        std::cout << A.otherCoefs[i] << ',';
    }
    std::cout << A.otherCoefs[N - 1] << "}\n\n";
    // вычисляем k-ую производную для разных шагов и печатаем
    for (unsigned int i = 0; i < 16; i++) {
        std::cout << "h = " << steps[i] << ";   " << k << "-order derivative(1) = ";
        f[i] = exp(1) * A.centralCoef;
        for (unsigned int j = 0; j < N; j++) {
            f[i] += A.otherCoefs[j] * exp(1 + points[j] * steps[i]);
        }
        f[i] = f[i] / pow(steps[i], k);
        std::cout << f[i] << ";\n";
    }
    // вычисляем логарифм ошибки для разных шагов и печатаем
    std::cout << "\n\n";
    for (int i = 0; i < 16; i++) {
        std::cout << "ln h/ln10 = " << -i << ";   ln (delta " << k << "-order derivative(1)) = ";
        delta[i] = log(abs(f[i] - exp(1))) / log(10);
        std::cout << delta[i] << ";\n";
    }
    // вычисляем тангенсы углов наклона линейных участков зависимости логарифма ошибки от логарифма шага и печатаем
    std::cout << "\n\n";
    for (int i = 0; i < 15; i++) {
        tg[i] = (delta[i] - delta[i + 1]);
        std::cout << "tg phi = " << tg[i] << ";\n";
    }
    // оценим тангенсы
    std::cout << "\n\n";
    max = tg[0];
    min = tg[0];
    for (int i = 1; i < 15; i++) {
        if (tg[i] > max) {
            max = tg[i];
        }
        if (tg[i] < min) {
            min = tg[i];
        }
    }
    std::cout << "max tg phi = " << max << ";\n";
    std::cout << "min tg phi = " << min << ";\n";
    return delta;
}

int main() {

    //заполняем массив узлов симметрично по возможности
    std::array<double, 3> points3 = { -1,1,2 };
    std::array<double, 4> points4 = { -2,-1,1,2 };
    std::array<double, 5> points5 = { -2,-1,1,2,3 };
    std::array<double, 16> steps, delta1, delta2, delta3, delta4, delta5, delta6;

    //заполняем массив шагов
    for (int i = 0; i < 16; i++) {
        steps[i] = pow(10, -i);
    }
    
    // рассчитываем первую производную для N=3,4,5
    delta1 = calcDerivative(points3, 1, steps);
    delta2 = calcDerivative(points4, 1, steps);
    delta3 = calcDerivative(points5, 1, steps);
    
    // рассчитываем вторую производную для N=3,4,5 и рисуем графики
    delta4 = calcDerivative(points3, 2, steps);
    delta5 = calcDerivative(points4, 2, steps);
    delta6 = calcDerivative(points5, 2, steps);

    return 0;
}
