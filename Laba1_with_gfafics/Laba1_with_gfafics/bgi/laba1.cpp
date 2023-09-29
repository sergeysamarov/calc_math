#include <iostream>
#include <cmath>
#include <array>
#include "graphics.h"


//using namespace std;

template<typename RealType, unsigned int N> //�������������� ������������
struct DerivativeCoef {
    RealType centralCoef;
    std::array<RealType, N> otherCoefs;
};

template<typename RealType, unsigned int N, typename L> // ���������� �������������� ������������� ��� k-�� �����������  
DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points, L k) noexcept
{
    //���������� ����������� ������� ������� ��� ������ �������������� ������������� 
    std::array<std::array<double, N + 2>, N + 1> A;
    A = {};// ��� ��������
    A[0][0] = 1;
    A[k][N + 1] = 1;
    for (unsigned int l = 2; l < k + 1; l++)//��������� k!
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

    //����� ������-�������, 0 ������� - �� ���
    double a;
    unsigned int s;
    for (unsigned int j = 1; j < N + 1; j++)
    {
        s = j;// ������� ��������� ������� � j-�� �������
        while (A[s][j] == 0) {
            s++;
        }
        a = A[s][j];
        for (unsigned int l = j; l < N + 2; l++)//������ 1 � A[s][j]
        {
            A[s][l] = A[s][l] / a;
        }
        for (unsigned int i1 = s + 1; i1 < N + 1; i1++)// ������ ���� ��� A[s][j]
        {
            a = A[i1][j];
            for (unsigned int j1 = j; j1 < N + 2; j1++) {
                A[i1][j1] = A[i1][j1] - A[s][j1] * a;
            }
        }
        for (unsigned int l = j; l < N + 2; l++)//������������ s ������ � j ������
        {
            a = A[s][l];
            A[s][l] = A[j][l];
            A[j][l] = a;
        }
        for (unsigned int i1 = 0; i1 < j; i1++)// ������ ���� ���� A[j][j]
        {
            a = A[i1][j];
            for (unsigned int j1 = j; j1 < N + 2; j1++) {
                A[i1][j1] = A[i1][j1] - A[s][j1] * a;
            }
        }
    }

    DerivativeCoef<double, N> C; // ������ ��������� �������������� ������������ � ���������
    C.centralCoef = A[0][N + 1];

    for (unsigned int i = 1; i < N + 1; i++) {
        C.otherCoefs[i - 1] = A[i][N + 1];
    }
    return C;
}
template<typename RealType, unsigned int N, typename L> // ���������� �����, ��� ��������� � ����, ��� N ����� � k-�� �����������
std::array<RealType, 16> calcDerivative(const std::array<RealType, N>& points, L k, std::array<double, 16>& steps) noexcept
{
    DerivativeCoef<double, N> A;
    std::array<double, 16> f; //�������� ������ ����������� � ����������� �� ����
    std::array<double, 16> delta; //�������� ��������� ������ ��� ������ ����������� � ����������� �� ����
    std::array<double, 16> tg; //�������� �������� ����� ������� �������� �������� ����������� ��������� ������ �� ��������� ����
    double max, min;

    //����� k-�� ����������� ��� e^x � ������ ����� N=3
    std::cout << "\n\n\n" << k << "-order derivative" << "\n\n";
    std::cout << "\n" << "N = " << N << "\n\n";
    std::cout << "points = {" << points[0];
    for (int i = 1; i < N; i++) {
        std::cout << ", " << points[i];
    }
    std::cout << "}\n\n";
    //��������� �������������� ������������
    A = calcDerivativeCoef(points, k);
    //���������� �������������� ������������
    std::cout << "\n" << "centralCoef = " << A.centralCoef << '\n';
    std::cout << "otherCoefs = {";
    for (unsigned int i = 0; i < N - 1; i++) {
        std::cout << A.otherCoefs[i] << ',';
    }
    std::cout << A.otherCoefs[N - 1] << "}\n\n";
    // ��������� k-�� ����������� ��� ������ ����� � ��������
    for (unsigned int i = 0; i < 16; i++) {
        std::cout << "h = " << steps[i] << ";   " << k << "-order derivative(1) = ";
        f[i] = exp(1) * A.centralCoef;
        for (unsigned int j = 0; j < N; j++) {
            f[i] += A.otherCoefs[j] * exp(1 + points[j] * steps[i]);
        }
        f[i] = f[i] / pow(steps[i], k);
        std::cout << f[i] << ";\n";
    }
    // ��������� �������� ������ ��� ������ ����� � ��������
    std::cout << "\n\n";
    for (int i = 0; i < 16; i++) {
        std::cout << "ln h/ln10 = " << -i << ";   ln (delta " << k << "-order derivative(1)) = ";
        delta[i] = log(abs(f[i] - exp(1)))/log(10);
        std::cout << delta[i] << ";\n";
    }
    // ��������� �������� ����� ������� �������� �������� ����������� ��������� ������ �� ��������� ���� � ��������
    std::cout << "\n\n";
    for (int i = 0; i < 15; i++) {
        tg[i] = (delta[i] - delta[i + 1]);
        std::cout << "tg phi = " << tg[i] << ";\n";
    }
    // ������ ��������
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

void drawAxes(void)     //������ ���
{
    int xmax, ymax, i;
    char c[6];
    setcolor(15);           //���� - �����
    xmax = getmaxx();
    ymax = getmaxy();

    line(xmax / 2 - 250, ymax / 2, xmax / 2 + 300, ymax / 2);
    line(xmax - 140, ymax / 2 - 470, xmax -140, ymax / 2 + 370);

    for (i = -16; i < 2; i++)
    {
        line(xmax / 2 + i * 30+240, ymax / 2, xmax / 2 + i * 30+240, ymax / 2 - 5);
        
    }
    for (i = -15; i < 13; i++)
    {
        
        line(xmax - 140, ymax / 2 + i * 30, xmax - 140 + 5, ymax / 2 + i * 30);
    }
    setcolor(12);           //������� � ����� - �������
    for (i = -15; i < 0; i++)
    {
        itoa(i, c, 10);
        outtextxy(xmax + i * 30 - 148, ymax / 2 + 5, c);
    }
    for (i = 1; i < 3; i++)
    {
        itoa(i, c, 10);
        outtextxy(xmax + i * 30 - 143, ymax / 2 + 5, c);
    }
    itoa(0, c, 10);
    outtextxy(xmax  - 150, ymax / 2 + 5, c);
   for (i = -12; i < 0; i++)
    {
        itoa(i, c, 10);
        outtextxy(xmax - 140 + 15, ymax / 2 - i * 30-8, c);
    }
   for (i = 1; i < 16; i++)
   {
       itoa(i, c, 10);
       outtextxy(xmax - 140 + 10, ymax / 2 - i * 30 - 8, c);
   }
}

void drawFunction(std::array<double, 16>& delta)     //������ �������
{
    int xmax, ymax, i, j, i0, j0;
    double x, y;
    setcolor(11);               //���� - ���������
    xmax = getmaxx();
    ymax = getmaxy();

    j0 = delta[0] * 30 + ymax / 2; //y(x) ��� ��������� �����
    i0 = xmax / 2 +210;

    for (i = 1; i < 16; i++)
    {
        line(i0, j0, i0-30, delta[i] * 30 + ymax / 2);
        i0 = i0 - 30;
        j0 = delta[i] * 30 + ymax / 2;
    }
}


int main() {

    //��������� ������ ����� ����������� �� �����������
    std::array<double, 3> points3 = { -1,1,2 };
    std::array<double, 4> points4 = { -2,-1,1,2 };
    std::array<double, 5> points5 = { -2,-1,1,2,3 };
    std::array<double, 16> steps, delta1, delta2, delta3, delta4, delta5, delta6;
    
    int windows1 = initwindow(700, 1000, "����������� ��������� ������ ���������� f'(1) �� ��������� ���� ��� 3 �����");
    int windows2 = initwindow(700, 1000, "����������� ��������� ������ ���������� f'(1) �� ��������� ���� ��� 4 �����");
    int windows3 = initwindow(700, 1000, "����������� ��������� ������ ���������� f'(1) �� ��������� ���� ��� 5 �����");
    int windows4 = initwindow(700, 1000, "����������� ��������� ������ ���������� f''(1) �� ��������� ���� ��� 3 �����");
    int windows5 = initwindow(700, 1000, "����������� ��������� ������ ���������� f''(1) �� ��������� ���� ��� 4 �����");
    int windows6 = initwindow(700, 1000, "����������� ��������� ������ ���������� f''(1) �� ��������� ���� ��� 5 �����");
    //��������� ������ �����
    for (int i = 0; i < 16; i++) {
        steps[i] = pow(10, -i);
    }
    // ������������ ������ ����������� ��� N=3,4,5 � ������ �������
    
    delta1 = calcDerivative(points3, 1, steps);
    setcurrentwindow(windows1);
    drawAxes();
    drawFunction(delta1);
    delta2 = calcDerivative(points4, 1, steps);
    setcurrentwindow(windows2);
    drawAxes();
    drawFunction(delta2);
    delta3 = calcDerivative(points5, 1, steps);
    setcurrentwindow(windows3);
    drawAxes();
    drawFunction(delta3);
    // ������������ ������ ����������� ��� N=3,4,5 � ������ �������
    delta4 = calcDerivative(points3, 2, steps);
    setcurrentwindow(windows4);
    drawAxes();
    drawFunction(delta4);
    delta5 = calcDerivative(points4, 2, steps);
    setcurrentwindow(windows5);
    drawAxes();
    drawFunction(delta5);
    delta6 = calcDerivative(points5, 2, steps);
    setcurrentwindow(windows6);
    drawAxes();
    drawFunction(delta6);
    
    while (!kbhit())
    {
        delay(200);
    }
   
    getch();        //�������� �������
    closegraph();   //����� ������ � ����������� ������
   
    return 0;
}

