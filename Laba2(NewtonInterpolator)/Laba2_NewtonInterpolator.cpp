// Laba2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <cmath>
#include <array>

template<typename xType, typename yType, unsigned int N>
class NewtonInterpolator {
private:
    std::array<xType, N> points;//узлы
    std::array<yType, N> values;//значения функции в узлах
    std::array<double, N> coefs;//коэффициенты интерполяции
public:
    NewtonInterpolator(const std::array<xType, N>& points, const std::array<yType, N>& values) noexcept//конструктор класса
    {
        this->points = points;
        this->values = values;
        this->coefs = {};
    };
    int NewtonCoefs() noexcept//вычисление коэффициентов интерполяции
    {
        // создаем матрицу разделенных разностей
        std::array<std::array<double, N>, N> A;
        // создаем массив коэффициентов
        // заполняем 0 столбец - разделенные разности 0 порядка: значения в узлах
        for (unsigned int i = 0; i < N; i++)
        {
            A[i][0] = values[i];
        };
        coefs[0] = A[0][0];
        std::cout << "Interpolation Coefs: " << coefs[0];
        // последовательно заполняем остальные столбцы - разделенные разности (k+1) порядка
        for (unsigned int k = 0; k < N - 1; k++)
        {
            for (unsigned int i = 0; i < N - k - 1; i++)
            {
                A[i][k + 1] = (A[i + 1][k] - A[i][k]) / (points[i + k + 1] - points[i]);
            };
            coefs[k + 1] = A[0][k + 1];
            std::cout << " ;" << coefs[k + 1];
        };
        std::cout << "\n";
        return 0;
    }
    yType interpolate(const xType& x) noexcept// вычисление значения интерполянта Ньютона в точке x
    {
        double y;
        y = coefs[N - 1];
        for (unsigned int i = 0; i < N - 1; i++)
        {
            y = y * (x - points[N - i - 2]) + coefs[N - i - 2];
        };
        return y;
    };
    int UniformPoints(unsigned int N, xType& endInterval)  noexcept// задание равномерного расположения N узлов на отрезке [0,l] и вычисление значений e^x в узлах
    {
        points = {};
        double l;
        l = endInterval / (N - 1);
        values[0] = 1;
        std::cout << "point = " << " " << 0 << "; exp(point) = " << 1 << ";\n";
        for (unsigned int i = 1; i < N; i++)
        {
            points[i] = l * i;
            values[i] = exp(points[i]);
            std::cout << "point = " << " " << points[i] << "; exp(point) = " << values[i] << ";\n";
        };
        std::cout << "\n";
        return 0;
    }
    // задание Чебышевского расположения узлов на отрезке [0,l] и вычисление значения e^x в узлах
    int ChebishevPoints(unsigned int N, xType& endInterval)  noexcept
    {
        points = {};
        for (unsigned int i = 0; i < N; i++)
        {
            points[i] = endInterval / 2.0 * (1 + cos(acos(0.0) * (2 * (i + 1) - 1) / N));
            values[i] = exp(points[i]);
            std::cout << "point = " << " " << points[i] << "; exp(point) = " << values[i] << ";\n";
        };
        std::cout << "\n";
        return 0;
    }
};

int main()
{
    unsigned int N = 3;
    std::array<double, 3> points = {};
    std::array<double, 3> values = {};
    std::array<double, 6> l = {};//правые концы отрезков
    std::array<double, 6> dis_uniform = {};//ошибки интерполяции при равномерном распределении узлов
    std::array<double, 6> dis_Chebishev = {};//ошибки интерполяции при равномерном распределении узлов
    double x, d, f;
    //создаем класс интерполятора Ньютона для N=3
    NewtonInterpolator < double, double, 3> P(points, values);
    std::cout << "NEWTON INTERPOLATION"<< "\n\n";
    std::cout << "N = " << 3 << "\n\n";
    std::cout << "Uniform points" << "\n\n";
    // для каждого из отрезков
    for (unsigned int i = 0; i < 6; i++)
    {
        //определяем правый конец отрезка, располагаем узлы равномерно, вычисляем значение e^x в каждом узле
        l[i] = 2.0 / pow(2, i);
        std::cout << "length = " << l[i] << ";\n";
        P.UniformPoints(N, l[i]);
        // считаем коэффициенты интерполянта Ньютона
        P.NewtonCoefs();
        // считаем ошибку интерполяции, края отрезка не берем, там значения совпадают
        for (unsigned int j = 1; j < 1000; j++)
        {
            d = l[i] / 1000; //расстояние между точками
            x = d * j; // текущая точка
            f = abs(P.interpolate(x) - exp(x));
            if (f > dis_uniform[i])
                dis_uniform[i] = f;
        };
        std::cout << "discrepency_uniform = " << dis_uniform[i] << ";";
        std::cout << "\n\n";
    };
    std::cout << "Chebishev points" << "\n\n";
    for (unsigned int i = 0; i < 6; i++)
    {
        std::cout << "length = " << l[i] << ";\n";
        P.ChebishevPoints(N, l[i]);
        // считаем коэффициенты интерполянта Ньютона
        P.NewtonCoefs();
        // считаем ошибку интерполяции, края отрезка не берем, там значения совпадают
        for (unsigned int j = 1; j < 1000; j++)
        {
            d = l[i] / 1000; //расстояние между точками
            x = d * j; // текущая точка
            f = abs(P.interpolate(x) - exp(x));
            if (f > dis_Chebishev[i])
                dis_Chebishev[i] = f;
        };
        std::cout << "discrepency_Chebishev = " << dis_Chebishev[i] << ";";
        std::cout << "\n\n";
    };

    N = 4;
    std::array<double, 4> points4 = {};
    std::array<double, 4> values4 = {};
    dis_uniform = {};//ошибки интерполяции при равномерном распределении узлов
    dis_Chebishev = {};//ошибки интерполяции при равномерном распределении узлов

    //создаем класс интерполятора Ньютона для N=4
    NewtonInterpolator < double, double, 4> P4(points4, values4);
    std::cout << "N = " << 4 << "\n\n";
    std::cout << "Uniform points" << "\n\n";
    // для каждого из отрезков
    for (unsigned int i = 0; i < 6; i++)
    {
        std::cout << "length = " << l[i] << ";\n";
        P4.UniformPoints(N, l[i]);
        // считаем коэффициенты интерполянта Ньютона
        P4.NewtonCoefs();
        // считаем ошибку интерполяции, края отрезка не берем, там значения совпадают
        for (unsigned int j = 1; j < 1000; j++)
        {
            d = l[i] / 1000; //расстояние между точками
            x = d * j; // текущая точка
            f = abs(P4.interpolate(x) - exp(x));
            if (f > dis_uniform[i])
                dis_uniform[i] = f;
        };
        std::cout << "discrepency_uniform = " << dis_uniform[i] << ";";
        std::cout << "\n\n";
    };
    std::cout << "Chebishev points" << "\n\n";
    for (unsigned int i = 0; i < 6; i++)
    {
        std::cout << "length = " << l[i] << ";\n";
        P4.ChebishevPoints(N, l[i]);
        // считаем коэффициенты интерполянта Ньютона
        P4.NewtonCoefs();
        // считаем ошибку интерполяции, края отрезка не берем, там значения совпадают
        for (unsigned int j = 1; j < 1000; j++)
        {
            d = l[i] / 1000; //расстояние между точками
            x = d * j; // текущая точка
            f = abs(P4.interpolate(x) - exp(x));
            if (f > dis_Chebishev[i])
                dis_Chebishev[i] = f;
        };
        std::cout << "discrepency_Chebishev = " << dis_Chebishev[i] << ";";
        std::cout << "\n\n";
    };


    N = 5;
    std::array<double, 5> points5 = {};
    std::array<double, 5> values5 = {};
    dis_uniform = {};//ошибки интерполяции при равномерном распределении узлов
    dis_Chebishev = {};//ошибки интерполяции при равномерном распределении узлов

    //создаем класс интерполятора Ньютона для N=4
    NewtonInterpolator < double, double, 5> P5(points5, values5);
    std::cout << "N = " << 5 << "\n\n";
    std::cout << "Uniform points" << "\n\n";
    // для каждого из отрезков
    for (unsigned int i = 0; i < 6; i++)
    {
        std::cout << "length = " << l[i] << ";\n";
        P5.UniformPoints(N, l[i]);
        // считаем коэффициенты интерполянта Ньютона
        P5.NewtonCoefs();
        // считаем ошибку интерполяции, края отрезка не берем, там значения совпадают
        for (unsigned int j = 1; j < 1000; j++)
        {
            d = l[i] / 1000; //расстояние между точками
            x = d * j; // текущая точка
            f = abs(P5.interpolate(x) - exp(x));
            if (f > dis_uniform[i])
                dis_uniform[i] = f;
        };
        std::cout << "discrepency_uniform = " << dis_uniform[i] << ";";
        std::cout << "\n\n";
    };
    std::cout << "Chebishev points" << "\n\n";
    for (unsigned int i = 0; i < 6; i++)
    {
        std::cout << "length = " << l[i] << ";\n";
        P5.ChebishevPoints(N, l[i]);
        // считаем коэффициенты интерполянта Ньютона
        P5.NewtonCoefs();
        // считаем ошибку интерполяции, края отрезка не берем, там значения совпадают
        for (unsigned int j = 1; j < 1000; j++)
        {
            d = l[i] / 1000; //расстояние между точками
            x = d * j; // текущая точка
            f = abs(P5.interpolate(x) - exp(x));
            if (f > dis_Chebishev[i])
                dis_Chebishev[i] = f;
        };
        std::cout << "discrepency_Chebishev = " << dis_Chebishev[i] << ";";
        std::cout << "\n\n";
    };


    return 0;
}
