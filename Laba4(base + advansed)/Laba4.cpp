

#include <iostream>
#include <array>
#include <cmath>
#include <type_traits>

template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename T>
using Dif = decltype(std::declval<T>() - std::declval<T>());

/* Функция производит интегрирование на одном отрезке для квадратуры Гаусса на 3 узлах */
template<typename Callable, typename RealType>
decltype(auto) integrate3(
    const Callable& func,  // Интегрируемая функция
    const typename RealType& start,  // начало отрезка
    const typename RealType& end  // конец отрезка
) {
    RealType I,A,B;
    A = (end - start) / 2;
    B = (end + start) / 2;
    I = A * (5.0 / 9 * func(-A * sqrt(0.6) + B) + 8.0/9 * func(B) + 5.0 / 9 * func(A * sqrt(0.6) + B));
    return I;
};

/* Функция производит интегрирование на одном отрезке для квадратуры Гаусса на 4 узлах */
template<typename Callable, typename RealType>
decltype(auto) integrate4(
    const Callable& func,  // Интегрируемая функция
    const typename RealType& start,  // начало отрезка
    const typename RealType& end  // конец отрезка
) {
    RealType I, A, B, Z3, Z4,W3,W4;
    A = (end - start) / 2;
    B = (end + start) / 2;
    Z4 = sqrt(3.0 / 7 + 2.0 / 7 * sqrt(1.2));//четвертый узел на [-1,1]
    Z3 = sqrt(3.0 / 7 - 2.0 / 7 * sqrt(1.2));//третий узел на [-1,1]
    W4 = (18 - sqrt(30.0)) / 36;// веса первой и четвертой функций Лагранжа на [-1,1]
    W3 = (18 + sqrt(30.0)) / 36;// веса второй и третьей функций Лагранжа на [-1,1]
    I = A * (W4 * func(-A * Z4 + B) + W3 * func(-A * Z3 + B) + +W3 * func(A * Z3 + B)+ W4 * func(A * Z4 + B));
    return I;
};

/* Функция производит интегрирование, разбивая отрезок на подотрезки длиной не более dx */
template<typename Callable, typename RealType>
decltype(auto) integrate(
    const Callable& func,  // Интегрируемая функция
    const typename RealType& start,  // начало отрезка
    const typename RealType& end,  // конец отрезка
    const Dif<typename RealType>& dx,  // Длина подотрезка
    const int N  // число узлов для квадратуры Гаусса
) {
    RealType start_current, end_current, I;
    start_current = start;
    end_current = start + dx;
    I = 0;
    //интегрирование по каждому отрезку
    while (end_current <= end) {
        if (N == 3)
        {
            I += integrate3(func, start_current, end_current);
        }
        else {
            I += integrate4(func, start_current, end_current);
        }
        start_current += dx;
        end_current += dx;
        }
    //интегрирование на последнем отрезке разбиения
    if (N == 3)
    {
        I += integrate3(func, start_current, end);
    }
    else {
        I += integrate4(func, start_current, end);
    }
    return I;
};

// считает интеграл с заданной точностью (контроль точности по правилу Рунге) и делает экстраполяцию Ричардсона
template<typename Callable, typename RealType>
decltype(auto) integrate_with_Runge(
    const Callable& func,  // Интегрируемая функция
    const typename RealType& start,  // начало отрезка
    const typename RealType& end,  // конец отрезка
    const typename RealType& accuracy,  // заданная точность интегрирования
    const int N  // число узлов для квадратуры Гаусса
) {
    RealType accuracy_current, h, I2h, Ih, eh;
    //вычисляем 2^{2N}
    int A=64;
    if (N == 4) {
        A *= 4;
    }

    h = (end - start)/2.0;
    I2h = integrate(func, start, end, 2 * h, N);//интегрирование с шагом 2h
    Ih = integrate(func, start, end, h, N);//интегрирование с шагом h
    eh= (Ih - I2h) / (A - 1);
    accuracy_current = abs(eh);// оценка по правилу Рунге для Ih
    
    while (accuracy_current > accuracy) {
        h /= 2.0;
        I2h = Ih;
        Ih= integrate(func, start, end, h, N);
        eh = (Ih - I2h) / (A - 1);
        accuracy_current = abs(eh);
    }
    // экстраполяция Ричардсона
    Ih += eh;
    return Ih;
};
double g(double x) { return sin(x); }// подынтегральная функция: возвращает sin(x)

int main()
{
    double a, b, max, min, I, h,d, dd, accurancy;
    a = 0;//начало отрезка интегрирования
    b = 10;//конец отрезка интегрирования
    h = 4.0;
    I = 1 - cos(10);// точное значение интеграла
    std::array<double, 6> delta = {}; //значения логарифма ошибки интегрирования для квадратур Гаусса с 3 узлами в зависимости от шага;
    std::array<double, 6> delta1 = {}; //значения логарифма ошибки интегрирования для квадратур Гаусса с 4 узлами в зависимости от шага;

    //Для квадратуры Гаусса на 3 узлах
    //интегрируем с шагами 4, 2, 1, 0.5, 0.25 , 0.125 ; считаем ошибки интегрирования для каждой величины шага и печатаем 

    std::cout << "Gaussian quadrature with 3 nodes\n\n";
    for (unsigned i = 0; i < 6; i++)
    {
        std::cout << "h = " << h << ";   ln(h) / ln(2) = " << log(h)/log(2) << ";   ln(delta I) / ln(2) = ";
        delta[i]= log(abs(integrate(g, a, b, h, 3) - I)) / log(2);//логарифм ошибки
        std::cout << delta1[i] << "; \n";
        h /= 2.0;
    }

    // оцениваем тангенсы углов наклона линейных участков зависимости логарифма ошибки от логарифма числа узлов
    std::cout << "\n";
    max = delta[1] - delta[0];
    min = delta[1] - delta[0];
    for (int i = 1; i < 5; i++) {
        d = (delta[i + 1] - delta[i]);
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

    //Для квадратуры Гаусса на 4 узлах
    //интегрируем с шагами 4, 2, 1, 0.5, 0.25 , 0.125 ; считаем ошибки интегрирования для каждой величины шага и печатаем 

    std::cout << "\n\n Gaussian quadrature with 4 nodes\n\n";
    h = 4.0;
    for (unsigned i = 0; i < 6; i++)
    {
        std::cout << "h = " << h << ";   ln(h) / ln(2) = " << log(h) / log(2) << ";   ln(delta I) / ln(2) = ";
        delta1[i] = log(abs(integrate(g, a, b, h, 4) - I)) / log(2);//логарифм ошибки
        std::cout << delta1[i] << "; \n";
        h /= 2.0;
    }

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

    //Для квадратуры Гаусса на 3 узлах
    std::cout << "\n\nGaussian quadrature with 3 nodes\n\n";
    
    // выполняем интегрирования с заданной точностью 0.1; 0.01; 0.001; 0.0001; 0.00001; 0.000001; 0.0000001 и печатаем реальную ошибку
    accurancy = 0.1;
    for (unsigned i = 0; i < 7; i++)
    {
        std::cout << "accurancy = " << accurancy << ";   ln(accurancy) / ln(10) = " << log(accurancy) / log(10) << ";   ln(delta I) / ln(10) = ";
        dd = log(abs(integrate_with_Runge(g, a, b, accurancy, 3) - I)) / log(10);//логарифм ошибки
        std::cout << dd << "; \n";
        accurancy /= 10.0;
    }
   
    //Для квадратуры Гаусса на 4 узлах
    std::cout << "\n\nGaussian quadrature with 4 nodes\n\n";

    // выполняем интегрирования с заданной точностью 0.1; 0.01; 0.001; 0.0001; 0.00001; 0.000001; 0.0000001 и печатаем реальную ошибку
    accurancy = 0.1;
    for (unsigned i = 0; i < 7; i++)
    {
        std::cout << "accurancy = " << accurancy << ";   ln(accurancy) / ln(10) = " << log(accurancy) / log(10) << ";   ln(delta I) / ln(10) = ";
        dd = log(abs(integrate_with_Runge(g, a, b, accurancy, 4) - I)) / log(10);//логарифм ошибки
        std::cout << dd << "; \n";
        accurancy /= 10.0;
    } 
}

