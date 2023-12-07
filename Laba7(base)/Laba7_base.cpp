
#include <iostream>
#include <array>
#include <cmath>
//#include <type_traits>

/* Это коэффициенты для метода на 4 шагах.Последний элемент - коэффициент перед h*f в правой части */
struct BDF4 {
    static constexpr unsigned int size = 4;
    static constexpr std::array<double, size + 1> alpha = { -48.0 / 25, 36.0 / 25, -16.0 / 25, 3.0 / 25, 12.0 / 25 };
};

/* Это таблица Бутчера для метода Рунге-Кутты 4 порядка, используем для разгона */
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

    //разгон, вычисляет стартовые значения y(h), y(2h), y(3h) методом Рунге-Кутты 4 порядка
    template<typename Table, unsigned int N>  // таблица бутчера, N - число узлов, включая начальный, N>=5
    std::array<StateAndArg, N> integrateDeq(
        const StateAndArg& initialState,// начальное условие
        double step// шаг
    ) {
        StateAndArg B;//состояние для аргументов правой части при вычислении коэффициентов К
        std::array<StateAndArg, N> Res = {}; //массив со значениями решений во всех точках с шагом step 
        RK4Table A;// таблица Бутчера
        std::array<double, 4> K = {};//кэффициенты K

        Res[0] = initialState;//начальное условие
        for (unsigned int j = 0; j < 4; j++) {
            // считаем коэффициенты K_i
            K[0] = calc(Res[j].arg);
            for (unsigned int i = 1; i < 4; i++) {
                B.arg = Res[j].arg + step * A.cColumn[i];
                B.state = Res[j].state + step * A.table[i][i - 1] * K[i - 1];
                K[i] = calc(B.arg);
            }
            // находим решение в следующей точке
            Res[j + 1].arg = Res[j].arg + step;
            Res[j + 1].state = Res[j].state + step * (K[0] + K[1] * 2 + K[2] * 2 + K[3]) / 6.0;
        }
        return Res;
    };

    //интегрирование уравнения
    template<typename BDF, unsigned int N>  // коэф. ФДН 4, N - число узлов, включая начальный, N>=5
    std::array<StateAndArg, N> integrateEq(
        const StateAndArg& initialState,// начальное условие
        double step// шаг
    ) {
        std::array<StateAndArg, N> Res = {}; //массив со значениями решений во всех точках с шагом step 
        BDF4 A;// коэффициенты для метода ФДН 4

        Res = integrateDeq <BDF4, N>(initialState, step);//разгон
        for (unsigned int j = 3; j < N-1; j++) {
            // находим решение в следующей точке
            Res[j + 1].arg = Res[j].arg + step;
            Res[j + 1].state = (48 * Res[j].state -36 * Res[j - 1].state +16 * Res[j - 2].state - 3 * Res[j - 3].state + 12 * step * calc(Res[j + 1].arg))/25.0;
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
        Res = integrateEq <RK4Table, N>(initialState, h);
        for (unsigned int j = 0; j < N; j++) {
            exactSolution = g(Res[j].arg);
            if (abs(Res[j].state - exactSolution) > err) {
                err = abs(Res[j].state - exactSolution);
            }
        }
        std::cout << log(err) / log(2) << ";\n";
        return err;
    };

    /*** Вычисляет правую часть ДУ - функцию f***/
    double calc(double& t) const {
        return t * t * t;
    };
    double g(double x) { return x * x * x * x / 4.0; }
};

//Реализация класса правой части дифференциального уравнения.То есть класс f(t, y) для y' = f(t,y). Здесь написан пример для уравнения осциллятора (x, v)' = (v, -x)
class Oscillator {// класс для решения задачи Коши для осциллятора

public:

    static constexpr unsigned int dim = 2;  // размерность задачи
    using Argument = double;  // тип аргумента, тип t
    using State = std::array<double, dim>;  // состояние

    struct StateAndArg {
        State state;// значения x и v
        Argument arg;// значение t
    };

    //разгон, вычисляет стартовые значения x(h),v(h); x(2h),v(2h); x(3h),v(3h) методом Рунге-Кутты 4 порядка
    template<typename Table, unsigned int N>  // таблица бутчера, N - число точек, в которых получаем результат интегрирования, включая начальную
    std::array<StateAndArg, N> integrateSyst(
        const StateAndArg& initialState,// начальное условие
        double step// шаг
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

    /*** Вычисляет правую часть ДУ - функцию f***/
    std::array<double, dim> calc(const StateAndArg& stateAndArg) const {
         return std::array<double, dim>{stateAndArg.state[1], -stateAndArg.state[0]};
    }

    //интегрирование осциллятора
    template<typename BDF, unsigned int N>  // коэф. ФДН 4, N - число узлов, включая начальный, N>=5
    std::array<StateAndArg, N> integrateOs(
        const StateAndArg& initialState,// начальное условие
        double step// шаг
    ) {
        std::array<StateAndArg, N> Res = {}; //массив со значениями решений во всех точках с шагом step 
        BDF4 A;// коэффициенты для метода ФДН 4
        double c1, c2, alpha;// константы, вычисляемые по 4 предыдущим шагам

        Res = integrateSyst <BDF4, N>(initialState, step);//разгон
        alpha = 12 * step * 0.04;
        for (unsigned int j = 3; j < N - 1; j++) {
            // находим решение в следующей точке
            Res[j + 1].arg = Res[j].arg + step;
            c1 = ( - 48 * Res[j].state[0] + 36 * Res[j - 1].state[0] - 16 * Res[j - 2].state[0] + 3 * Res[j - 3].state[0]) * 0.04;
            c2 =( - 48 * Res[j].state[1] + 36 * Res[j - 1].state[1] - 16 * Res[j - 2].state[1] + 3 * Res[j - 3].state[1]) * 0.04;
            Res[j + 1].state[0] = (-c1 - alpha * c2)/(1 + alpha * alpha);
            Res[j + 1].state[1] = (-c2 + alpha * c1) / (1 + alpha * alpha);
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
        std::cout << "log_10 (h) = " << log(h) / log(10) << ";     log_10 (err) = ";
        std::array<Oscillator::StateAndArg, N> Res;// решение уравнения 
        Res = integrateOs <BDF4, N>(initialState, h);
        for (unsigned int j = 0; j < N; j++) {
            exactSolution = gOs(Res[j].arg);
            if (abs(Res[j].state[0] - exactSolution) > err) {
                err = abs(Res[j].state[0] - exactSolution);
            }
        }
        std::cout << log(err) / log(10) << ";\n";
        return err;
    };
    double gOs(double x) { return sin(x); }
};

int main()
{
    Deq R;
    Deq::StateAndArg D;// начальное условие для уравнения
    D.arg = 0;
    D.state = 0;

    std::cout << "BDF 4 for the equation" << "\n\n";

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

    std::cout << "BDF 4 for the Oscillator" << "\n\n";

    //для шагов h = 0.1, 0.01, 0.001 сосчитаем решение уравнения осциллятора, найдем ошибку 

    Oscillator RSyst;
    Oscillator::StateAndArg DSyst;// начальное условие для осциллятора
    DSyst.arg = 0;
    DSyst.state[0] = 0;
    DSyst.state[1] = 1;

    std::array<double, 3> delta = {}; //значения логарифма ошибки для разных шагов;
   
    delta[0] = log(RSyst.findErrOs <double, 51>(DSyst)) / log(10);
    delta[1] = log(RSyst.findErrOs <double, 501>(DSyst)) / log(10);
    delta[2] = log(RSyst.findErrOs <double, 5001>(DSyst)) / log(10);

    // оценим тангенсы углов наклона линейных участков зависимости логарифма ошибки от логарифма шага
    double max, min, d;

    max = -delta[1] + delta[0];
    min = -delta[1] + delta[0];
    d = -delta[2] + delta[1];
    if (d > max) {
        max = d;
    }
    if (d < min) {
        min = d;
    }
        
    // печаем оценку тангенсов
    std::cout << "\n" << "max tg phi = " << max << ";\n";
    std::cout << "min tg phi = " << min << ";\n";
    
}

