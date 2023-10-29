// Laba5.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <iostream>
#include <array>
#include <cmath>
/**
    Решает уравнение Кеплера методом Ньютона
    * ecc - эксцентриситет, принадлежит (0, 1)
    * meanAnomaly - средняя аномалия, М (радианы)
    * maxIter - максимальное количество итераций
    * tol - точность, с которой нужно отыскать решение

    Рекомендуемое поведение. Если решение не нашлось за maxIter итераций - выбрасывать исключение.
    Если приближения к решению между итерациями меняются не более, чем на tol, то решение достигнуто.
**/
template<typename RealType, unsigned int N>
std:: array<RealType, N> keplerSolver(RealType& ecc, RealType& meanAnomaly, unsigned int& maxIter, RealType& tol) {//возвращает массив со значениями итераций 
    double E;// n приближение
    double E0 = meanAnomaly;// начальное приближение
    std::array<RealType, N> It = {};// создаем массив для итераций и зануляем его
    It[1] = E0;
    for (unsigned int i = 0; i < maxIter; i++) {//делаем итерации методом Ньютона
        E = E0 - (E0 - ecc * sin(E0) - meanAnomaly) / (1 - ecc * cos(E0));
        It[i + 2] = E;
        if (abs(E - E0) <= tol) {//Если приближения к решению между итерациями меняются не более, чем на tol, то решение достигнуто
            It[0] = E;
            return It;
        }
        else {// если нет, то продолжаем
            E0 = E;
        }
    }
    //Если решение не нашлось за maxIter итераций - выбрасываем исключение
    std::cout << "Тhe limit of the number of iterations has been exceeded" << "\n";
};

template<typename RealType>
struct SystemSolution {
    RealType x;
    RealType y;
    RealType Ac;//точность решения системы
};

template<typename A>
struct ArgumentGetter;

template<typename R, typename Arg>
struct ArgumentGetter<R(Arg)> {
    using Argument = Arg;
};

template<typename Callable, typename RealType>
decltype(auto) solve(
    const Callable& func,                                             // функция f
    const RealType& tau,                                              // шаг тау
    const typename ArgumentGetter<Callable>::Argument& initialGuess,  // начальное приближение
    const unsigned int nIteration                                     // количество итераций
) {
    RealType x;// n приближение
    RealType x0;// (n-1) приближение
    SystemSolution <double> S;
    
    x0 = initialGuess;
    std::cout << "\n Iterations:" << '\n';
    std::cout << x0 << '\n';
    for (unsigned int i = 0; i < nIteration; i++) {//делаем итерации методом простой итерации с релаксацией
        x = x0 - tau * func(x0);
        std::cout << x << '\n';
        S.Ac = sqrt((x - x0) * (x - x0) + (tan(x) - tan(x0)) * (tan(x) - tan(x0)));
        if (S.Ac <= 0.000001) {//Если приближения к решению системы между итерациями меняются не более, чем 10^{-6}, то решение достигнуто
            S.x = x;
            S.y = tan(x);
            return S;
        }
        else {// если нет, то продолжаем
            x0 = x;
        }
    }
    //Если решение не нашлось за maxIter итераций - выбрасываем исключение
    std::cout << "Тhe limit of the number of iterations has been exceeded" << "\n";
};
double g(double x) { return x*x+tan(x)*tan(x)-1; }// подынтегральная функция: возвращает f(x)

int main()
{
    double meanAnomaly, tol, E,x,x0,tau;
    E = 0;
    meanAnomaly = asin(1) / 2.0;
    tol = 0.000001;
    std::array<double, 9> It = {};//массив итераций для решения уравнения Кеплера
    unsigned int maxIter = 8;// максимальное число итераций
    unsigned int k = 1;
    std::cout << "Solution of the Kepler equation by Newton's method" << "\n";
    for (double ecc : {0.1, 0.2, 0.5, 0.8}) {//эксцентриситеты
        std::cout << "\n ecc = " << ecc << ";\n";
        It=keplerSolver <double, 9> (ecc, meanAnomaly, maxIter, tol);
        std::cout << "E = " << It[0] << ";\n\n";
        while (It[k] != It[0]) {
            std::cout << "log_10 |E_" << k-1 << " - E| = " << log(abs(It[k]-It[0]))/log(10) << "; \n";
            k++;
        }
        It = {};
        k = 1;
    }

    //Решение системы уравнений
    std::cout << "\n\nSolving the system of equations by simple iteration with relaxation" << "\n\n";
    //задаем параметр релаксации для решения системы 
    tau = 0.245;
    std::cout << "tau = " << tau << ";\n\n";
    //задаем начальное приближение для х
    x0 = asin(1) * 2 / 5;
    std::cout << "initialGuess = " << x0 << ";\n";
    SystemSolution <double> S;
    // решаем уравнение для х методом простой итерации с релаксацией
    S = solve(g, tau, x0, 8);
    //печатаем все решения системы
    std::cout << "\n System solutions:"  << '\n';
    std::cout << "x = " << S.x << ";  " << "y = " << S.y << ";\n";
    std::cout << "x = " << -S.x << ";  " << "y = " << -S.y << ";\n";
    std::cout << "Accuracy = " << S.Ac << ";\n";
}