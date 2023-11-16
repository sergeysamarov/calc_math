#include <iostream>
#include <array>
#include <cmath>
#include <type_traits>

/* Это таблица Бутчера для метода Рунге-Кутты 4 порядка */
struct RK4Table {
    static constexpr unsigned int stages = 4;
    static constexpr std::array<std::array<double, stages>, stages> table = { { {0,0,0,0}, {0.5,0,0,0},{0,0.5, 0,0},{0,0,1,0} } };
    static constexpr std::array<double, stages> cColumn = { 0,0.5,0.5,1 };
    static constexpr std::array<double, stages> bString = { 1.0 / 6,1.0 / 3,1.0 / 3,1.0 / 6 };
};

class Deq { // класс для решения задачи Коши для уравнения
public:

    static constexpr unsigned int dim = 1;  // размерность задачи
    using Argument = double;  // тип аргумента, тип t
    using State = double;  // состояние

    struct StateAndArg {
        State state;//значение y
        Argument arg;// значение t
    };

    /*** Вычисляет правую часть ДУ - функцию f***/
    double calc(const StateAndArg& stateAndArg) const {
        return stateAndArg.arg * stateAndArg.arg * stateAndArg.arg;
    };

    // решает задачу Коши для уравнения 1 порядка
    template<typename Table, unsigned int N>  // таблица бутчера, N - число точек, в которых получаем результат интегрирования, включая начальную
    std::array<StateAndArg, N> integrateDeq(
        const StateAndArg& initialState,// начальное условие
        const Argument& endTime,// конец отрезка интегрирования
        double step// шаг
        //const RHS& rhs,
    ) {
        StateAndArg B;//состояние для аргументов правой части при вычислении коэффициентов К
        std::array<StateAndArg, N> Res = {}; //массив со значениями решений во всех точках с шагом step 
        RK4Table A;// таблица Бутчера
        std::array<double, 4> K = {};//кэффициенты K

        Res[0] = initialState;//начальное условие
        for (unsigned int j = 0; j < N - 1; j++) {
            // считаем коэффициенты K_i
            K[0] = calc(Res[j]);
            for (unsigned int i = 1; i < 4; i++) {
                B.arg = Res[j].arg + step * A.cColumn[i];
                B.state = Res[j].state + step * A.table[i][i - 1] * K[i - 1];
                K[i] = calc(B);
            }
            // находим решение в следующей точке
            Res[j + 1].arg = Res[j].arg + step;
            Res[j + 1].state = Res[j].state + step * (K[0] + K[1] * 2 + K[2] * 2 + K[3]) / 6.0;
        }
        return Res;
    };

    //считает ошибку и печатает логарифм шага и логарифм ошибки, тангенс угла наклона и решение для шага h=1/128 
    template<typename RealType, unsigned int N>
    decltype(auto) findErr(
        const Deq::StateAndArg& initialState// начальное условие
    ) {
        RealType h, err, exactSolution;
        err = 0;//ошибка
        h = 5.0 / (N - 1);
        std::cout << "log_2 (h) = " << log(h) / log(2) << ";     log_2 (err) = ";
        std::array<Deq::StateAndArg, N> Res;// решение уравнения 
        Res = integrateDeq <RK4Table, N>(initialState, N, h);
        for (unsigned int j = 0; j < N; j++) {
            exactSolution = g(Res[j].arg);
            if (abs(Res[j].state - exactSolution) > err) {
                err = abs(Res[j].state - exactSolution);
            }
        }
        std::cout << log(err) / log(2) << ";\n";
        /*
        //распечатка тангенса угла наклона и решения для шага h=1/128
        if (N == 641) {
            std::cout << "\n\ntg(phi) = 0\n\nSolution equation for h=1/128\n\n";
            for (unsigned int j = 0; j < N; j++) {
                std::cout << "x = " << Res[j].arg << ";     y = " << Res[j].state << ";\n";
            }
        }
        */
        return err;
    };
    double g(double x) { return x * x * x * x / 4.0; }
};

class Oscillator {// класс для решения задачи Коши для осциллятора
public:

    static constexpr unsigned int dim = 2;  // размерность задачи
    using Argument = double;  // тип аргумента, тип t
    using State = std::array <double, dim>;  // состояние

    struct StateAndArg {
        State state;// значения y и z
        Argument arg;// значение х
    };

    /*** Вычисляет правые части системы ***/
    std::array <double, dim> calc(const StateAndArg& stateAndArg) const {
        return std::array<double, dim>{stateAndArg.state[1], -stateAndArg.state[0]};
    }

    // решает задачу Коши для системы 
    template<typename Table, unsigned int N>  // таблица бутчера, N - число точек, в которых получаем результат интегрирования, включая начальную
    std::array<StateAndArg, N> integrateSyst(
        const StateAndArg& initialState,// начальное условие
        const Argument& endTime,// конец отрезка интегрирования
        double step// шаг
        //const RHS& rhs,
    ) {
        StateAndArg B;//состояние для аргументов правых частей при вычислении коэффициентов К и М
        std::array<StateAndArg, N> Res = {};//массив со значениями решений во всех точках с шагом step
        RK4Table A;// таблица Бутчера
        std::array<double, 4> K = {};// коэффициенты K
        std::array<double, 4> M = {};// коэффициенты M
        Res[0] = initialState;//начальные условия
        for (unsigned int j = 0; j < N - 1; j++) {
            // считаем коэффициенты K_i, M_i
            K[0] = calc(Res[j])[0];
            M[0] = calc(Res[j])[1];
            for (unsigned int i = 1; i < 4; i++) {
                B.arg = Res[j].arg + step * A.cColumn[i];
                B.state[0] = Res[j].state[0] + step * A.table[i][i - 1] * K[i - 1];
                B.state[1] = Res[j].state[1] + step * A.table[i][i - 1] * M[i - 1];
                K[i] = calc(B)[0];
                M[i] = calc(B)[1];
            }
            // находим решение в следующей точке
            Res[j + 1].arg = Res[j].arg + step;
            Res[j + 1].state[0] = Res[j].state[0] + step * (K[0] + K[1] * 2 + K[2] * 2 + K[3]) / 6.0;
            Res[j + 1].state[1] = Res[j].state[1] + step * (M[0] + M[1] * 2 + M[2] * 2 + M[3]) / 6.0;
        }
        return Res;
    };
    //считает ошибку и печатает логарифм шага и логарифм ошибки, тангенс угла наклона и решение для шага h=1/128 
    template<typename RealType, unsigned int N>
    decltype(auto) findErrOs(
        const Oscillator::StateAndArg& initialState// начальное условие
    ) {
        RealType h, err, exactSolution;
        err = 0;//ошибка
        h = 5.0 / (N - 1);
        std::cout << "log_2 (h) = " << log(h) / log(2) << ";     log_2 (err) = ";
        std::array<Oscillator::StateAndArg, N> Res;// решение уравнения 
        Res = integrateSyst <RK4Table, N>(initialState, N, h);
        for (unsigned int j = 0; j < N; j++) {
            exactSolution = gOs(Res[j].arg);
            if (abs(Res[j].state[0] - exactSolution) > err) {
                err = abs(Res[j].state[0] - exactSolution);
            }
        }
        std::cout << log(err) / log(2) << ";\n";
        /*
        //распечатка решения для шага h=1/128
        if (N == 641) {
            std::cout << "\n\nSolution equation for Oscillator with h=1/128\n\n";
            for (unsigned int j = 0; j < N; j++) {
                std::cout << "x = " << Res[j].arg << ";     y = " << Res[j].state[0] << ";\n";
            }
        }*/
        return err;
    };
    double gOs(double x) { return sin(x); }
};



#include <iostream>

int main()
{
    Deq R;
    Deq::StateAndArg D;// начальное условие для уравнения
    D.arg = 0;
    D.state = 0;

    std::cout << "Runge-Kutta method for the equation" << "\n\n";

    //для шагов h = 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128 сосчитаем решение уравнения и найдем ошибку

    R.findErr <double, 6>(D);
    R.findErr <double, 11>(D);
    R.findErr <double, 21>(D);
    R.findErr <double, 41>(D);
    R.findErr <double, 81>(D);
    R.findErr <double, 161>(D);
    R.findErr <double, 321>(D);
    R.findErr <double, 641>(D);

    std::cout << "\n\nWe have found the exact solution\n\ntg(phi) = 0\n\n";

    std::cout << "\n\n";

    std::cout << "Runge-Kutta method for the Oscillator" << "\n\n";

    //для шагов h = 1, 1/2, 1/4, 1/8, 1/16, 1/32, 1/64, 1/128 сосчитаем решение уравнения осциллятора, найдем ошибку, 

    Oscillator RSyst;
    Oscillator::StateAndArg DSyst;// начальное условие для осциллятора
    DSyst.arg = 0;
    DSyst.state[0] = 0;
    DSyst.state[1] = 1;
    std::array<double, 8> delta = {}; //значения логарифма ошибки для разных шагов;

    delta[0] = log(RSyst.findErrOs <double, 6>(DSyst)) / log(2);
    delta[1] = log(RSyst.findErrOs <double, 11>(DSyst)) / log(2);
    delta[2] = log(RSyst.findErrOs <double, 21>(DSyst)) / log(2);
    delta[3] = log(RSyst.findErrOs <double, 41>(DSyst)) / log(2);
    delta[4] = log(RSyst.findErrOs <double, 81>(DSyst)) / log(2);
    delta[5] = log(RSyst.findErrOs <double, 161>(DSyst)) / log(2);
    delta[6] = log(RSyst.findErrOs <double, 321>(DSyst)) / log(2);
    delta[7] = log(RSyst.findErrOs <double, 641>(DSyst)) / log(2);

    // оценим тангенсы углов наклона линейных участков зависимости логарифма ошибки от логарифма шага
    double max, min, d;

    std::cout << "\n";
    max = -delta[1] + delta[0];
    min = -delta[1] + delta[0];
    for (int i = 1; i < 7; i++) {
        d = (-delta[i + 1] + delta[i]);
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