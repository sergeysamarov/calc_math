// laba6(advanced).cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <array>
#include <cmath>
#include <type_traits>
#include <stdio.h>

/* Это таблица Бутчера для метода Дормана-Принса 5(4) */
struct DP45 {
    static constexpr unsigned int stages = 7;
    static constexpr std::array<std::array<double, stages>, stages> table = { { {0,0,0,0,0,0,0}, {0.2,0,0,0,0,0,0},{3.0 / 40,9.0 / 40, 0,0,0,0,0},{44.0 / 45,-56.0 / 15,32.0 / 9,0,0,0,0},
        {19372.0 / 6561,-25360.0 / 2187,64448.0 / 6561,-212.0 / 729,0,0,0}, {9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176, -5103.0 / 18656,0,0},
        {35.0 / 384, 0,500.0 / 1113, 125.0 / 192,-2187.0 / 6784,11.0 / 84,0}} };
    static constexpr std::array<double, stages> cColumn = { 0,0.2,0.3,0.8,8.0 / 9,1,1 };
    static constexpr std::array<double, stages> bString1 = { 35.0 / 384,0,500.0 / 1113,125.0 / 192,-2187.0 / 6784,11.0 / 84,0 };
    static constexpr std::array<double, stages> bString2 = { 5179.0 / 57600,0,7571.0 / 16695,393.0 / 640,-92097.0 / 339200,187.0 / 2100,1.0 / 40 };
    static constexpr unsigned int approximation = 5;
};

class Arenstorf {// класс для построения орбиты Аренсторфа
public:

    static constexpr unsigned int dim = 4;  // размерность задачи
    using Argument = double;  // тип аргумента, тип t
    using State = std::array <double, dim>;  // состояние

    struct StateAndArg {
        State state;// значения x,u,y,v
        Argument arg;// значение t
    };
    double step; //текущий шаг

    Arenstorf() {// конструктор класса
        step = 0.01;
    }
    /*** Вычисляет правые части системы ***/
    std::array <double, dim> calc(const StateAndArg& stateAndArg) const {
        double A, B, mu, nu;
        mu = 0.012277471;
        nu = 1 - mu;
        A = pow((stateAndArg.state[0] + mu) * (stateAndArg.state[0] + mu) + stateAndArg.state[2] * stateAndArg.state[2], 1.5);
        B = pow((stateAndArg.state[0] - nu) * (stateAndArg.state[0] - nu) + stateAndArg.state[2] * stateAndArg.state[2], 1.5);
        return std::array<double, dim>{stateAndArg.state[1], stateAndArg.state[0] + 2 * stateAndArg.state[3] - nu * (stateAndArg.state[0] + mu) / A - mu * (stateAndArg.state[0] - nu) / B,
            stateAndArg.state[3], stateAndArg.state[2] - 2 * stateAndArg.state[1] - nu * stateAndArg.state[2] / A - mu * stateAndArg.state[2] / B};
    }

    // считает значения в следующей точке интегрирования по двум методам DP4 и DP5
    std::array<StateAndArg, 2> integrateSyst(
        const StateAndArg& initialState,// текущее состояние
        double step// шаг
    ) {
        StateAndArg B;//состояние для аргументов правых частей при вычислении коэффициентов К
        StateAndArg Res1;//следующее состояние по методу DP4
        StateAndArg Res2;//следующее состояние по методу DP5
        DP45 A;// таблица Бутчера
        std::array<std::array<double, 7>, 4> K = { {{},{},{},{}} };// коэффициенты K
        // считаем коэффициенты К_1 для всех 4 переменных
        for (unsigned int l = 0; l < 4; l++) {
            K[l][0] = calc(initialState)[l];// l - номер компоненты К_1
        }
        //считаем следующие коэффициенты
        for (unsigned int i = 1; i < 7; i++) {
            // считаем состояние для аргументов правых частей при вычислении коэффициентов K
            B.arg = initialState.arg + step * A.cColumn[i];
            for (unsigned int l = 0; l < 4; l++) {
                B.state[l] = initialState.state[l];
                for (unsigned int j = 0; j < i; j++) {
                    B.state[l] += step * A.table[i][j] * K[l][j];
                }
            }
            // теперь сами коэф. K
            for (unsigned int l = 0; l < 4; l++) {
                K[l][i] = calc(B)[l];
            }
        }

        // находим решение в следующей точке для 2 методов DP4 и DP5
        Res1.arg = initialState.arg + step;
        Res2.arg = initialState.arg + step;
        for (unsigned int l = 0; l < 4; l++) {
            Res1.state[l] = initialState.state[l];
            Res2.state[l] = initialState.state[l];
            for (unsigned int j = 0; j < 7; j++) {
                Res1.state[l] += step * A.bString1[j] * K[l][j];
                Res2.state[l] += step * A.bString2[j] * K[l][j];
            }
            // печатаем оба решения в следующей точке
            std::cout << "\n\n";
            std::cout << "Res1.arg = " << Res1.arg << "  " << "Res1.state[" << l << "] = " << Res1.state[l] << "\n\n";
            std::cout << "Res2.arg = " << Res1.arg << "  " << "Res2.state[" << l << "] = " << Res2.state[l] << "\n\n";
        }
        return std::array<StateAndArg, 2>{Res1, Res2};
    };
    StateAndArg stepControl(//Контроль за шагом
        const StateAndArg& initialState// текущее состояние
    ) {
        double minStep = 0.001;//минимальный шаг
        double maxStep = 0.1;//максимальный шаг
        double tol = 0.000001;//точность
        double error = 0;
        double stepNew;// новый шаг
        std::array<StateAndArg, 2> CR;//два решения в следующей точке по DP4 и по DP5
        CR = integrateSyst(initialState, step);
        //вычисляем ошибку
        for (unsigned int j = 0; j < 4; j++) {
            if (abs(CR[0].state[j] - CR[1].state[j]) > error) {
                error = abs(CR[0].state[j] - CR[1].state[j]);
            }
        }
        //вычисляем и печатаем новый оптимальный шаг
        stepNew = pow(tol * step / 2.0 / error, 0.2) * step;
        std::cout << stepNew << "\n\n";
        if ((stepNew <= maxStep) && (stepNew >= minStep)) {
            step = stepNew;
        }
        else {
            if (stepNew > maxStep) {
                step = maxStep;
            }
            else {
                step = minStep;
            }
        }
        if (error > tol) {//если ошибка > требуемой точности, то пересчитываем результаты предыдущего шага
            CR = integrateSyst(initialState, step);
        }
        return CR[1];
    }
};

#include <iostream>

int main()
{
    Arenstorf R;
    Arenstorf::StateAndArg D, D1;// текущее состояние и следующее состояние
    //задаем начальное состояние
    D.arg = 0;
    D.state = { 0.994,0,0,-2.00158510637908252240537862224 };
    double T = 17.0652165601579625588917206249;//период

    std::cout << "\n\n";
    //считаем состояния системы и записываем данные в файл "orb.txt", чтобы потом нарисовать орбиту на питоне
    FILE* out;
    fopen_s(&out, "orb.txt", "a");//открываем файл для записи
    fprintf(out, "%.9f %.9f %.9f\n", D.arg, D.state[0], D.state[2]);//заносим начальное состояние
    while (D.arg < T) {//пока не пройдем весь период
        D1 = R.stepControl(D);//рассчитываем состояние
        D = D1;
        fprintf(out, "%.9f %.9f %.9f\n", D.arg, D.state[0], D.state[2]);//записываем данные по этому состоянию
    }
    fclose(out);//закрываем файл
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
