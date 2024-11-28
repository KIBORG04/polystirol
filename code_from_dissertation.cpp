#include <vcl.h>

#include <math.h>

#include <stdlib.h>

#pragma hdrstop
#define max(a, b)(((a) > (b)) ? (a) : (b)) // определяем максимальный изргументов
#include "main.h"

#include "USelect.h"
 // ------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm1 * Form1;
const N = 7,
    NN = 7000;
TColor cl;
bool change = false; // флаг, обозначающий наличие изменений на стр. 2
// ------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent * Owner): TForm(Owner) {
    cold = false;
    while (!Table1 -> Eof)
        if (Table1NAME -> Value == "Опыт1м")
            break;
        else
            Table1 -> Next();
    // Заполняем поля формы расчёта T0
    Edit9n -> Text = FloatToStr(Table1CAC1 -> Value);
    Edit10n -> Text = FloatToStr(Table1CAC1 -> Value);
    Edit9n -> Text = FloatToStr(Table1CAC1 -> Value);
    Edit13n -> Text = "373";
    Edit14n -> Text = DBMemo2 -> Lines -> operator[](0);
    Edit1n -> Text = DBMemo3 -> Lines -> operator[](0);
    cl = clBlack;
}
// ------------------------------------------------------------------------
// описание функции инициализации параметров модели и др. величин
void TForm1::Ini() {
    // Технологические параметры
    ps = StrToFloat(DBEdit7 -> Text); // плотность стирола
    pd = StrToFloat(DBEdit8 -> Text); // плотность дивинила
    pr = StrToFloat(DBEdit6 -> Text); // плотность растворителя
    ch = StrToFloat(Edit21 -> Text); // теплоемкость хладагента
    ph = StrToFloat(Edit22 -> Text); // плотность хладагента
    th1 = StrToFloat(DBEdit16 -> Text); // температура хладагента, стадия 1
    th2_1 = StrToFloat(DBEdit17 -> Text); // температура хладагента, стадия 2_1
    th2_2 = StrToFloat(DBEdit22 -> Text); // температура хладагента, стадия 2_2
    gh1 = StrToFloat(DBEdit21 -> Text); // расход хладагента, стадия 1
    gh2_1 = StrToFloat(DBEdit20 -> Text); // расход хладагента, стадия 2_1
    gh2_2 = StrToFloat(DBEdit19 -> Text); // расход хладагента, стадия 2_2
    ms = StrToFloat(DBEdit2 -> Text); // масса стирола
    md1 = StrToFloat(DBEdit3 -> Text); // масса дивинила, стадия 2_1
    md2 = StrToFloat(DBEdit4 -> Text); // масса дивинила, стадия 2_2
    mr = StrToFloat(DBEdit5 -> Text); // масса растворителя
    vr = mr / pr; // объем растворителя
    vs = ms / ps; // объем стирола
    vd1 = md1 / pd; // объем дивинила, стадия 2_1
    vd2 = md2 / pd; // объем дивинила, стадия 2_2
    vk = StrToFloat(Edit8 -> Text); // объем катализатора
    vt = vr + vs + vk; // объем реакционной смеси (первая ста-ия)
    Map = StrToFloat(Edit25 -> Text); // масса реактора
    cap = StrToFloat(Edit26 -> Text); // теплоемкость реактора
    F = StrToFloat(Edit27 -> Text); // площадь теплообмена
    pt = StrToFloat(Edit28 -> Text); // плотность реакционной смеси
    ct = StrToFloat(Edit29 -> Text); // теплоемкость реакционной смеси
    mms = StrToFloat(Edit4 -> Text); // молекулярная масса стирола
    mmd = StrToFloat(Edit5 -> Text); // молекулярная масса дивинила
    m0s = StrToFloat(DBEdit13 -> Text); // начальная концентрация стирола
    m0d1 = StrToFloat(DBEdit14 -> Text); // нач. концентрация дивинила (2_1)
    m0d2 = StrToFloat(DBEdit15 -> Text); // нач. концентрация дивинила (2_2)
    // коэффициенты модели
    A = StrToFloat(Edit30 -> Text);
    B = StrToFloat(Edit31 -> Text);
    C = StrToFloat(Edit32 -> Text);
    Km0 = StrToFloat(Edit33 -> Text); // константа скорости роста
    k2 = StrToFloat(Edit34 -> Text);
    b2 = StrToFloat(Edit36 -> Text);
    Ems = StrToFloat(Edit37 -> Text); //энергия активации роста цепи (стирол)
    Emd = StrToFloat(Edit38 -> Text); //энергия активации роста цепи (дивинил)
    Ktes = StrToFloat(Edit39 -> Text); // коэффициент тепловыделения (стирол)
    Kted = StrToFloat(Edit40 -> Text); // коэффициент тепловыделения (дивинил)
    // Массивы температуры
    // температура на стадии 1
    TMAX1 = 0;
    dTMAX1 = 0;
    Tex1 = new float[DBMemo1 -> Lines -> Count];
    for (int i = 0; i < DBMemo1 -> Lines -> Count; i++) {
        Tex1[i] = StrToFloat(DBMemo1 -> Lines -> Strings[i]);
        if (Tex1[i] > TMAX1)
            TMAX1 = Tex1[i];
        if (i > 0)
            if (dTMAX1 < Tex1[i] - Tex1[i - 1])
                dTMAX1 = Tex1[i] - Tex1[i - 1];
    }
    // температура на стадии 2_1
    TMAX2_1 = 0;
    dTMAX2_1 = 0;
    Tex2_1 = new float[DBMemo2 -> Lines -> Count];
    for (int i = 0; i < DBMemo2 -> Lines -> Count; i++) {
        Tex2_1[i] = StrToFloat(DBMemo2 -> Lines -> Strings[i]);
        if (Tex2_1[i] > TMAX2_1)
            TMAX2_1 = Tex2_1[i];
        if (i > 0)
            if (dTMAX2_1 < Tex2_1[i] - Tex2_1[i - 1])
                dTMAX2_1 = Tex2_1[i] - Tex2_1[i - 1];
    }
    // температура на стадии 2_2
    TMAX2_2 = 0;
    dTMAX2_2 = 0;
    Tex2_2 = new float[DBMemo3 -> Lines -> Count];
    for (int i = 0; i < DBMemo3 -> Lines -> Count; i++) {
        Tex2_2[i] = StrToFloat(DBMemo3 -> Lines -> Strings[i]);
        if (Tex2_2[i] > TMAX2_2)
            TMAX2_2 = Tex2_2[i];
        if (i > 0)
            if (dTMAX2_2 < Tex2_2[i] - Tex2_2[i - 1])
                dTMAX2_2 = Tex2_2[i] - Tex2_2[i - 1];
    }
    T10 = Tex1[0]; // начальная температура полимеризации (1)
    T20 = Tex2_1[0]; // начальная температура реакционной смеси (2_1)
    T30 = Tex2_2[0]; // начальная температура реакционной смеси (2_2)
    J0 = StrToFloat(DBEdit9 -> Text); // действительная КАЦ
    J0doz = StrToFloat(DBEdit10 -> Text); // КАЦ по данным дозировок
    Cac1 = StrToFloat(DBEdit11 -> Text); // КАЦ, стадия 2_1
    Cac2 = StrToFloat(DBEdit12 -> Text); // КАЦ, стадия 2_2
    if (PageControl1 -> ActivePage == TabSheet1) {
        e = StrToFloat(Edit6n -> Text); // точность оптимизации
        dt = StrToFloat(Edit7n -> Text); // шаг дискретизации
    }
    if (PageControl1 -> ActivePage == TabSheet6) {
        e = StrToFloat(Edit43 -> Text); // точность оптимизации
        dt = StrToFloat(Edit43 -> Text); // шаг дискретизации
    }
}
// ------------------------------------------------------------------------
//кнопка «Выход»
void __fastcall TForm1::Button2Click(TObject * Sender) {
    Form1 -> Close();
}
// ------------------------------------------------------------------------
// Описание функции расчёта по модели (без учета стадии инициирования)
int __fastcall TForm1::ModRas(float ml, // коэффициент идентификацииатематической модели
    float Km, // константа скорости роста цепи
    float J0, // КАЦ
    float St, // степень при КАЦ
    float Kte, // коэффициент тепловыделения
    float m0, // начальная концентрация мономера
    float vt, // текущий объем реакционной смеси
    float k, // коэффициент модели расчета Кm
    float b, // коэффициент модели расчета Кm
    float E1, // энергия активации
    float E2,
    float so_b, // количество сухого остатка
    float m, // масса мономера на текущей стадии
    float M, // масса реакционной смеси (все компоненты)
    float * Xm, // указатель на массив конверсии мономера
    float * T, // указатель на массив температуры реакции
    float gh, // расход хладагента
    float T_k, // температура начала следующей стадии
    float * Tmax, // максимальная температура
    float th, // температура хладагента
    bool f) {
    int i, n, c = 0;
    float L, // уровень загрузки реактора
    R = 8.32; // универсальная газовая постоянная
    float x1, y1, x2, y2, x3, y3, x4, y4; // коэффициенты Рунге-Кутта
    L = (vt - 3000) / 28500; // уровень заполнения реактора
    i = 0;
    * Tmax = T[0];
    while (1) {
        x1 = Km * pow(J0, St + ml * so) * (1 - Xm[i]);
        if (CBLim -> Checked) { // ограничитель
            if (x1 > 2 * (1 - Xm[i]) / dt)
                x1 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt; // 10% от общей "массы"
        }
        y1 = (Kte * x1 * m0 * vt - (Kf * F * L * (T[i] - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        x2 = Km * pow(J0, St + ml * so) * (1 - (Xm[i] + x1 * dt / 2));
        if (CBLim -> Checked) { // ограничитель
            if (x2 > 2 * (1 - Xm[i]) / dt)
                x2 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y2 = (Kte * x2 * m0 * vt - (Kf * F * L * ((T[i] + y1 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        x3 = Km * pow(J0, St + ml * so) * (1 - (Xm[i] + x2 * dt / 2));
        if (CBLim -> Checked) { // ограничитель
            if (x3 > 2 * (1 - Xm[i]) / dt)
                x3 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y3 = (Kte * x3 * m0 * vt - (Kf * F * L * ((T[i] + y2 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        x4 = Km * pow(J0, St + ml * so) * (1 - (Xm[i] + x3 * dt / 2));
        if (CBLim -> Checked) { // ограничитель
            if (x4 > 2 * (1 - Xm[i]) / dt)
                x4 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y4 = (Kte * x4 * m0 * vt - (Kf * F * L * ((T[i] + y3 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        Xm[i + 1] = Xm[i] + dt * (x1 + 2 * (x2 + x3) + x4) / 6;
        T[i + 1] = T[i] + dt * (y1 + 2 * (y2 + y3) + y4) / 6;
        if (T[i + 1] > * Tmax)
            *
            Tmax = T[i + 1];
        if (Xm[i + 1] > 1.001) {
            n = i + 1;
            break;
        }
        so = so_b + m * Xm[i + 1] / M;
        Km = k * (1 - b * so) * exp(-(E1 + E2 * so) / (R * T[i + 1])); // расчёт константыкорости роста цепи
        Kf = A - B * exp(C * so); // коэффициент теплопередачи
        // условие прекращения расчета для модели с текущей КАЦ
        if (Xm[i + 1] >= 0.98) {
            if (!f) {
                if (T[i + 1] <= T_k) {
                    n = i + 1;
                    break;
                }
            } else {
                c++;
                if (c > 2) {
                    n = i + 1;
                    break;
                }
            }
        }
        if (i > NN - 3) {
            n = i + 1;
            break;
        }
        i++;
    }
    return n;
}
// ------------------------------------------------------------------------
// Описание функции расчёта по модели (с учетом стадии инициирования)
int __fastcall TForm1::ModRas(float ml, // коэффициент идентификацииатематической модели
    float Ky, // константа скорости инициирования
    float Km, // константа скорости роста цепи
    float J0, // КАЦ
    float St, // степень при КАЦ
    float Kte, // коэффициент тепловыделения
    float m0, // начальная концентрация мономера
    float vt, // текущий объем реакционной смеси
    float Ky0, // постоянная константы скорости инииирования
    float Ei, // энергия активации реакции инициирования
    float k, // коэффициент модели расчета Кm
    float b, // коэффициент модели расчета Кm
    float Em, // енергия активации
    float so_b, // количество сухого остатка
    float m, // масса мономера на текущей стадии
    float M, // масса реакционной смеси (все компоненты)
    float * Xy, // указатель на массив конверсии инициатора
    float * Xm, // указатель на массив конверсии мономера
    float * T, // указатель на массив температуры реакции
    float gh, // расход хладагента
    float T_k, // температура начала следующей стадии
    float * Tmax,
    float th, // температура хладагента
    bool f,
    int i0) // начальное значение переменной i
{
    int i, n, c = 0;
    float L, // уровень загрузки реактора
    R = 8.32; // универсальная газовая постоянная
    float x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4; // коэффициенты Рунге-Кутта
    L = (vt - 3000) / 28500; // уровень заполнения реактора
    i = i0;
    * Tmax = 0;
    while (1) {
        z1 = Ky * m0 * (1 - Xy[i]) * (1 - Xm[i]);
        x1 = (Ky * J0 * (1 - Xy[i]) + Km * MyPow(J0 * Xy[i], St + ml * so)) * (1 - Xm[i]);
        if (CBLim -> Checked) {
            // ограничитель
            if (x1 > 2 * (1 - Xm[i]) / dt)
                x1 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt; // 10% от общей "массы"
        }
        y1 = (Kte * x1 * m0 * vt - (Kf * F * L * (T[i] - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);

        чш = Ky * m0 * (1 - (Xy[i] + z1 * dt / 2)) * (1 - Xm[i]);
        x2 = (Ky * J0 * (1 - (Xy[i] + z1 * dt / 2)) +
            Km * MyPow(J0 * (Xy[i] + z1 * dt / 2), St + ml * so)) * (1 - (Xm[i] + x1 * dt / 2));
        y2 = (Kte * x2 * m0 * vt - (Kf * F * L * ((T[i] + y1 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);

        z3 = Ky * m0 * (1 - (Xy[i] + z2 * dt / 2)) * (1 - Xm[i]);
        x3 = (Ky * J0 * (1 - (Xy[i] + z2 * dt / 2)) +
            Km * MyPow(J0 * (Xy[i] + z2 * dt / 2), St + ml * so)) * (1 - (Xm[i] + x2 * dt / 2));
        if (CBLim -> Checked) { // ограничитель
            if (x3 > 2 * (1 - Xm[i]) / dt)
                x3 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y3 = (Kte * x3 * m0 * vt - (Kf * F * L * ((T[i] + y2 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        z4 = Ky * m0 * (1 - (Xy[i] + z3 * dt / 2)) * (1 - Xm[i]);
        x4 = (Ky * J0 * (1 - (Xy[i] + z3 * dt / 2)) +
            Km * MyPow(J0 * (Xy[i] + z3 * dt / 2), St + ml * so)) * (1 - (Xm[i] + x3 * dt / 2));
        if (CBLim -> Checked) { // ограничитель
            if (x4 > 2 * (1 - Xm[i]) / dt)
                x4 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y4 = (Kte * x4 * m0 * vt - (Kf * F * L * ((T[i] + y3 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        Xy[i + 1] = Xy[i] + dt * (z1 + 2 * (z2 + z3) + z4) / 6;
        Xm[i + 1] = Xm[i] + dt * (x1 + 2 * (x2 + x3) + x4) / 6;
        T[i + 1] = T[i] + dt * (y1 + 2 * (y2 + y3) + y4) / 6;
        if (T[i + 1] > * Tmax)
            *
            Tmax = T[i + 1];
        if (Xm[i + 1] > 1.001 || Xy[i + 1] > 1.001) {
            n = i + 1;
            break;
        }
        so = so_b + m * Xm[i + 1] / M; // сухой остаток
        Km = k * (1 - b * so) * exp(-Em / (R * T[i + 1])); // расчёт константы скорости роста
        Kf = A - B * exp(C * so); // коэффициент теплопередачи
        Ky = Ky0 * exp(-Ei / (R * T[i + 1])); // константа скорости инициирования
        // условие прекращения расчета для модели с текущей КАЦ
        if (Xm[i + 1] >= 0.98) {
            if (!f) {
                if (T[i + 1] <= T_k) {
                    n = i + 1;
                    break;
                }
            } else {
                c++;
                if (c > 2) {
                    n = i + 1;
                    break;
                }
            }
        }
        if (i > NN - 3) {
            n = i + 1;
            break;
        }
        i++;
    }
    return n;
}
// ------------------------------------------------------------------------
// функция расчёта по модели с добавлением вероятности
// Описание функции расчёта по модели (с учетом стадии инициирования)
int __fastcall TForm1::ModRas_V(double ** pp,
    double * mm, // для определения концентрации мономера
    float ml, // коэффициент идентификации математической модели
    float Ky, // константа скорости инициирования
    float Km, // константа скорости роста цепи
    float J0, // КАЦ
    float St, // степень при КАЦ
    float Kte, // коэффициент тепловыделения
    float m0, // начальная концентрация мономера
    float vt, // текущий объем реакционной смеси
    float Ky0, // постоянная константы скорости инииирования
    float Ei, // энергия активации реакции инициирования
    float k, // коэффициент модели расчета Кm
    float b, // коэффициент модели расчета Кm
    float Em, // энергия активации
    float E2,
    float so_b, // количество сухого остатка с предыдущей стадии
    float m, // масса мономера на текущей стадии
    float M, // масса реакционной смеси (все компоненты)
    float * Xy, // указатель на массив конверсии инициатора
    float * Xm, // указатель на массив конверсии мономера
    float * T, // указатель на массив температуры реакции
    float gh, // расход хладагента
    float T_k, // температура начала следующей стадии
    float * Tmax,
    float th, // температура хладагента
    bool f, // булевская переменная
    int i0, // начальное значение переменной i
    int n1) // количество фракций для расчёта
{
    int i, n, c = 0, j = 0;
    float L, // уровень загрузки реактора
    R = 8.32; // универсальная газовая постоянная
    double x1 = 0, y1 = 0, z1 = 0, x2 = 0, y2 = 0, z2 = 0, x3 = 0, y3 = 0, z3 = 0, x4 = 0, y4 = 0, z4 = 0,
        p1 = 0, p2 = 0, p3 = 0, p4 = 0; // коэффициенты Рунге-Кутта
    L = (vt - 3000) / 28500; // уровень заполнения реактора
    i = i0;
    float Mux;
    if (Ky != 0)
        Mux = pow(J0, -0.5);
    else
        Mux = pow(J0, -0.75);
    while (1) {
        z1 = Ky * m0 * (1 - Xy[i]) * (1 - Xm[i]);
        x1 = (Ky * J0 * (1 - Xy[i]) + Km * MyPow(J0 * Xy[i], St + ml * so)) * (1 - Xm[i]);
        if (CBLim -> Checked) {
            // ограничитель
            if (x1 > 2 * (1 - Xm[i]) / dt)
                x1 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt; // 10% от общей "массы"
        }
        y1 = (Kte * x1 * m0 * vt - (Kf * F * L * (T[i] - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        z2 = Ky * m0 * (1 - (Xy[i] + z1 * dt / 2)) * (1 - Xm[i]);
        x2 = (Ky * J0 * (1 - (Xy[i] + z1 * dt / 2)) +
            Km * MyPow(J0 * (Xy[i] + z1 * dt / 2), St + ml * so)) * (1 - (Xm[i] + x1 * dt / 2));
        if (CBLim -> Checked) { // ограничитель
            if (x2 > 2 * (1 - Xm[i]) / dt)
                x2 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y2 = (Kte * x2 * m0 * vt - (Kf * F * L * ((T[i] + y1 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        z3 = Ky * m0 * (1 - (Xy[i] + z2 * dt / 2)) * (1 - Xm[i]);
        x3 = (Ky * J0 * (1 - (Xy[i] + z2 * dt / 2)) +
            Km * MyPow(J0 * (Xy[i] + z2 * dt / 2), St + ml * so)) * (1 - (Xm[i] + x2 * dt / 2));
        if (CBLim -> Checked) { // ограничитель
            if (x3 > 2 * (1 - Xm[i]) / dt)
                x3 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y3 = (Kte * x3 * m0 * vt - (Kf * F * L * ((T[i] + y2 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        z4 = Ky * m0 * (1 - (Xy[i] + z3 * dt / 2)) * (1 - Xm[i]);
        x4 = (Ky * J0 * (1 - (Xy[i] + z3 * dt / 2)) +
            Km * MyPow(J0 * (Xy[i] + z3 * dt / 2), St + ml * so)) * (1 - (Xm[i] + x3 * dt / 2));
        if (CBLim -> Checked) { // ограничитель
            if (x4 > 2 * (1 - Xm[i]) / dt)
                x4 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y4 = (Kte * x4 * m0 * vt - (Kf * F * L * ((T[i] + y3 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        Xy[i + 1] = Xy[i] + dt * (z1 + 2 * (z2 + z3) + z4) / 6;
        Xm[i + 1] = Xm[i] + dt * (x1 + 2 * (x2 + x3) + x4) / 6;
        T[i + 1] = T[i] + dt * (y1 + 2 * (y2 + y3) + y4) / 6;
        // вставляем блок с вероятностью
        mm[i + 1] = m0 * (1 - Xm[i + 1]); // расчёт концентрации мономера
        for (j = 0; j < n1; j++) {
            if (j == 0) {
                if (Ky != 0) {
                    p1 = -Km * (Mux) * mm[i] * pp[j][i] + Ky * J0 * (1 - Xy[i]) * mm[i];
                    p2 = -Km * (Mux) * mm[i] * (pp[j][i]) + Ky * J0 * (1 - Xy[i]) * mm[i];
                    p3 = -Km * (Mux) * mm[i] * (pp[j][i]) + Ky * J0 * (1 - Xy[i]) * mm[i];
                    p4 = -Km * (Mux) * mm[i] * (pp[j][i]) + Ky * J0 * (1 - Xy[i]) * mm[i];
                } else {
                    p1 = -Km * (Mux) * mm[i] * pp[j][i];
                    p2 = -Km * (Mux) * mm[i] * (pp[j][i]);
                    p3 = -Km * (Mux) * mm[i] * (pp[j][i]);
                    p4 = -Km * (Mux) * mm[i] * (pp[j][i]);
                }
            } else {
                p1 = Km * (Mux) * mm[i] * (pp[j - 1][i] - (pp[j][i]));
                p2 = Km * (Mux) * mm[i] * ((pp[j - 1][i]) - (pp[j][i]));
                p3 = Km * (Mux) * mm[i] * ((pp[j - 1][i]) - (pp[j][i]));
                p4 = Km * (Mux) * mm[i] * ((pp[j - 1][i]) - (pp[j][i]));
            }
            if (j == n1 - 1) {
                p1 = Km * (Mux) * mm[i] * pp[j - 1][i];
                p2 = Km * (Mux) * mm[i] * pp[j - 1][i];
                p3 = Km * (Mux) * mm[i] * pp[j - 1][i];
                p4 = Km * (Mux) * mm[i] * pp[j - 1][i];
            }
            pp[j][i + 1] = pp[j][i] + dt * (p1 + 2 * (p2 + p3) + p4) / 6;
            if (pp[j][i + 1] < 0)
                pp[j][i + 1] = 0;
        }
        if (T[i + 1] > * Tmax)
            *
            Tmax = T[i + 1];
        if (Xm[i + 1] > 1.001 || Xy[i + 1] > 1.001) {
            n = i + 1;
            break;
        }
        so = so_b + m * Xm[i + 1] / M; // сухой остаток
        Km = k * (1 - b * so) * exp(-(Em + E2 * so) / (R * T[i + 1])); // расчёт констанытыкорости роста
        Kf = A - B * exp(C * so); // коэффициент теплопередачи
        Ky = Ky0 * exp(-Ei / (R * T[i + 1])); // константа скорости инициирования
        // условие прекращения расчета для модели с текущей КАЦ
        if (Xm[i + 1] >= 0.99) {
            if (!f) {
                if (T[i + 1] <= T_k) {
                    n = i + 1;
                    break;
                }
            } else {
                c++;
                if (c > 200) {
                    n = i + 1;
                    break;
                }
            }
        }
        if (i > NN - 3) {
            n = i + 1;
            break;
        }
        i++;
    }
    return n;
}
// ------------------------------------------------------------------------
// функция расчёта по модели с целью идентификации параметров ml1 и ml2
(графики)
void TForm1::ModRas_Identif(float ml1, float ml2, bool graf) {
    int i, j, ii, n, n1;
    float t;
    float R = 8.32,
        tt1, tt2; // для запоминания времени окончания стадии
    if (!cold) {
        Chart1 -> SeriesList -> Clear();
        Chart2 -> SeriesList -> Clear();
    } else
        cold = !cold;
    for (i = 0; i < ListBox1 -> Count; i++) {
        t = 0;
        Table1 -> First();
        while (!Table1 -> Eof)
            if (Table1NAME -> Value == ListBox1 -> Items -> Strings[i])
                break;
            else
                Table1 -> Next();
        // описываем графики
        TPointSeries * Te = new TPointSeries(this);
        Te -> ParentChart = Chart2;
        Te -> Title = "Tэксп" + Table1NAME -> Value;
        Te -> Pointer -> Style = psCircle;
        TLineSeries * Tr = new TLineSeries(this);
        Tr -> ParentChart = Chart2;
        Tr -> Title = "Tрасч" + Table1NAME -> Value;
        Tr -> LinePen -> Width = 2;
        TLineSeries * X = new TLineSeries(this);
        X -> ParentChart = Chart1;
        X -> Title = "Xрасч" + Table1NAME -> Value;
        X -> LinePen -> Width = 2;
        // выделяем память под массивы данных
        Xms = new float[NN];
        Xmd_1 = new float[NN];
        Xmd_2 = new float[NN];
        T1 = new float[NN];
        T2_1 = new float[NN];
        T2_2 = new float[NN];
        // Инициализация параметров
        Ini();
        // Расчет первой стадии
        if (CB1 -> Checked) {
            Xms[0] = 0;
            T1[0] = Tex1[0];
            so = 0;
            Kms = Km0 * exp(-Ems / (R * T1[0]));
            Kf = A - B * exp(C * so);
            n = ModRas(ml1, Kms, J0, 0.5, Ktes, m0s, vt, Km0, 0, Ems, 0, so, ms,
                mr + ms, Xms, T1, gh1, Tex2_1[0], & Tmax, th1, false);
            if (graf)
                for (j = 0; j <= n; j++) {
                    X -> AddXY(t, Xms[j], "", clDefault);
                    Tr -> AddXY(t, T1[j], "", clDefault);
                    if (j < DBMemo1 -> Lines -> Count)
                        Te -> AddXY(j * Table1DT1 -> Value, Tex1[j], "", clDefault);
                    t += dt;
                }
        }
        // Расчет второй стадии 2_1
        if (CB2_1 -> Checked) {
            Xmd_1[0] = 0;
            T2_1[0] = Tex2_1[0];
            so = ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1);
            Kmd = k2 * (1 - b2 * so) * exp(-Emd / (R * T2_1[0]));
            Kf = A - B * exp(C * so);
            n = ModRas(ml2, Kmd, Cac1, 0.25, Kted, m0d1, vt + vd1, k2, b2, Emd, 0,
                so, md1, mr + ms + md1, Xmd_1, T2_1, gh2_1, Tex2_2[0], &
                Tmax, th2_1, false);
            tt1 = t; // запоминаем время окончания стадии
            if (graf)
                for (j = 0; j <= n; j++) {
                    X -> AddXY(t, Xmd_1[j], "", clDefault);
                    Tr -> AddXY(t, T2_1[j], "", clDefault);
                    if (j < DBMemo2 -> Lines -> Count)
                        Te -> AddXY(tt1 + j * Table1DT2_1 -> Value, Tex2_1[j], "",
                            clDefault);
                    t += dt;
                }
        }
        if (CB2_2 -> Checked) {
            // Расчет второй стадии 2_2
            Xmd_2[0] = 0;
            T2_2[0] = Tex2_2[0];
            so = (ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1) + md1 * 1 /
                (ms + mr + md1)) * (ms + mr + md1) / (ms + mr + md1 + md2);
            Kf = A - B * exp(C * so);
            Kmd = k2 * (1 - b2 * so) * exp(-Emd / (R * T2_2[0]));
            n = ModRas(ml2, Kmd, Cac2, 0.25, Kted, m0d2, vt + vd1 + vd2, k2, b2, Emd, 0,
                so, md2, mr + ms + md1 + md2, Xmd_2, T2_2, gh2_2, Tex2_2[0], &
                Tmax, th2_2, true);
            tt2 = t;
            if (graf)
                for (j = 0; j <= n; j++) {
                    X -> AddXY(t, Xmd_2[j], "", clDefault);
                    Tr -> AddXY(t, T2_2[j], "", clDefault);
                    if (j < DBMemo3 -> Lines -> Count)
                        Te -> AddXY(tt2 + j * Table1DT2_2 -> Value, Tex2_2[j], "",
                            clDefault);
                    t += dt;
                }
        }
        // Освобождаем память
        delete Xms;
        delete Xmd_1;
        delete Xmd_2;
        delete T1;
        delete T2_1;
        delete T2_2;
    }
}
// ------------------------------------------------------------------------
// функция расчёта по модели с целью идентификации параметров ml и Ky (граф.)
void TForm1::ModRas_Identif(float Ky0, float Ei, float ml1, float ml2, float E1, float E2, bool graf) {
    int i, j, ii, n;
    float t;
    float R = 8.32, Tmax,
        tt; // для запоминания времени окончания стадии
    if (!cold) {
        Chart1 -> SeriesList -> Clear();
        Chart2 -> SeriesList -> Clear();
    } else
        cold = !cold;
    for (i = 0; i < ListBox1 -> Count; i++) {
        t = 0;
        Table1 -> First();
        while (!Table1 -> Eof)
            if (Table1NAME -> Value == ListBox1 -> Items -> Strings[i])
                break;
            else
                Table1 -> Next();
        // описываем графики
        TPointSeries * Te = new TPointSeries(this);
        Te -> ParentChart = Chart2;
        Te -> Title = "Tэксп" + Table1NAME -> Value;
        Te -> Pointer -> Style = psCircle;
        TLineSeries * Tr = new TLineSeries(this);
        Tr -> ParentChart = Chart2;
        Tr -> Title = "Tрасч" + Table1NAME -> Value;
        Tr -> LinePen -> Width = 2;
        TLineSeries * X = new TLineSeries(this);
        X -> ParentChart = Chart1;
        X -> Title = "Xрасч" + Table1NAME -> Value;
        X -> LinePen -> Width = 2;
        TLineSeries * XX = new TLineSeries(this);
        XX -> ParentChart = Chart1;
        XX -> Title = "Xy_расч" + Table1NAME -> Value;
        XX -> LinePen -> Width = 2;
        // выделяем память под массивы данных
        Xy = new float[NN];
        Xyd = new float[NN];
        Xms = new float[NN];
        Xmd_1 = new float[NN];
        Xmd_2 = new float[NN];
        T1 = new float[NN];
        T2_1 = new float[NN];
        T2_2 = new float[NN];
        // Инициализация параметров
        Ini();
        // Расчет первой стадии
        if (CB1 -> Checked) {
            Xy[0] = 0;
            Xms[0] = 0;
            T1[0] = Tex1[0];
            so = 0;
            Kms = Km0 * exp(-Ems / (R * T1[0]));
            Ky = Ky0 * exp(-Ei / (R * T1[0]));
            Kf = A - B * exp(C * so);
            n = ModRas(ml1, Ky, Kms, J0, 0.5, Ktes, m0s, vt, Ky0, Ei, Km0, 0, Ems, so,
                ms, mr + ms, Xy, Xms, T1, gh1, Tex2_1[0], & Tmax, th1, false, 0);
            if (graf)
                for (j = 0; j <= n; j++) {
                    XX -> AddXY(t, Xy[j], "", clDefault);
                    X -> AddXY(t, Xms[j], "", clDefault);
                    Tr -> AddXY(t, T1[j], "", clDefault);
                    if (j < DBMemo1 -> Lines -> Count)
                        Te -> AddXY(j * Table1DT1 -> Value, Tex1[j], "", clDefault);
                    t += dt;
                }
        }
        // Расчет второй стадии 2_1
        if (CB2_1 -> Checked) {
            Xmd_1[0] = 0;
            T2_1[0] = Tex2_1[0];
            so = ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1);
            Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * T2_1[0]));
            Kf = A - B * exp(C * so);
            n = ModRas(ml2, Kmd, Cac1, 0.25, Kted, m0d1, vt + vd1, k2, b2, E1, E2, so, md1,
                mr + ms + md1, Xmd_1, T2_1, gh2_1, Tex2_2[0], & Tmax, th2_1, false);
            tt = t; // запоминаем время окончания стадии
            if (graf)
                for (j = 0; j <= n; j++) {
                    X -> AddXY(t, Xmd_1[j], "", clDefault);
                    Tr -> AddXY(t, T2_1[j], "", clDefault);
                    if (j < DBMemo2 -> Lines -> Count)
                        Te -> AddXY(tt + j * Table1DT2_1 -> Value, Tex2_1[j], "",
                            clDefault);
                    t += dt;
                }
        }
        if (CB2_2 -> Checked) {
            // Расчет второй стадии 2_2
            Xmd_2[0] = 0;
            T2_2[0] = Tex2_2[0];
            so = (ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1) +
                md1 * 1 / (ms + mr + md1)) * (ms + mr + md1) / (ms + mr + md1 + md2);
            Kf = A - B * exp(C * so);
            Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * T2_2[0]));
            n = ModRas(ml2, Kmd, Cac2, 0.25, Kted, m0d2, vt + vd1 + vd2, k2, b2,
                E1, E2, so, md2, mr + ms + md1 + md2, Xmd_2, T2_2, gh2_2, Tex2_2[0], &
                Tmax, th2_2, true);
            tt = t;
            if (graf)
                for (j = 0; j <= n; j++) {
                    X -> AddXY(t, Xmd_2[j], "", clDefault);
                    Tr -> AddXY(t, T2_2[j], "", clDefault);
                    if (j < DBMemo3 -> Lines -> Count)
                        Te -> AddXY(tt + j * Table1DT2_2 -> Value, Tex2_2[j], "",
                            clDefault);
                    t += dt;
                }
        }
        // Освобождаем память
        delete Xy;
        delete Xms;
        delete Xmd_1;
        delete Xmd_2;
        delete T1;
        delete T2_1;
        delete T2_2;
    }
}
// ------------------------------------------------------------------------
// функция расчёта по модели с целью идентификации параметров ml и Ky
float TForm1::ModRas_Identif(float Ky0, float Ei, float ml1, float ml2, float E1, float E2) {
    int i, j, ii, n, k;
    float S, a, b, aa = 0.5, bb = 0.5;
    float R = 8.32;
    float Tmax; // максимальное значение температуры
    S = 0;
    Tmax = 0;
    Sr = 0;
    Skv = 0;
    for (i = 0; i < ListBox1 -> Count; i++) {
        Table1 -> First();
        while (!Table1 -> Eof)
            if (Table1NAME -> Value == ListBox1 -> Items -> Strings[i])
                break;
            else
                Table1 -> Next();
        // выделяем память под массивы данных
        Xy = new float[NN];
        Xms = new float[NN];
        Xmd_1 = new float[NN];
        Xmd_2 = new float[NN];
        T1 = new float[NN];
        T2_1 = new float[NN];
        T2_2 = new float[NN];
        // Инициализация параметров
        Ini();
        // Расчет первой стадии
        if (CB1 -> Checked) {
            Xy[0] = 0;
            Xms[0] = 0;
            T1[0] = Tex1[0];
            so = 0;
            Kms = Km0 * exp(-Ems / (R * T1[0]));
            Ky = Ky0 * exp(-Ei / (R * T1[0]));
            Kf = A - B * exp(C * so);
            n = ModRas(ml1, Ky, Kms, J0, 0.5, Ktes, m0s, vt, Ky0, Ei, Km0, 0, Ems, so,
                ms, mr + ms, Xy, Xms, T1, gh1, Tex2_1[0], & Tmax, th1, false, 0);
            k = 0;
            for (j = 0; j < DBMemo1 -> Lines -> Count - 1; j++) {
                a = (Tex1[j + 1] - Tex1[j]) / Table1DT1 -> Value;
                b = Tex1[j] - a * (j * Table1DT1 -> Value);
                for (ii = 0; ii < (Table1DT1 -> Value / dt); ii++) {
                    if (RadioButton1 -> Checked) {
                        k++;
                        if (k < n) {
                            S = S + pow(T1[k] - (a * k * dt + b), 2) /
                                ((DBMemo1 -> Lines -> Count - 1) * Table1DT1 -> Value) * dt;
                            Skv = Skv + pow(T1[k] - (a * k * dt + b), 2) /
                                ((DBMemo1 -> Lines -> Count - 1) * Table1DT1 -> Value) * dt;
                            Sr = Sr + fabs(T1[k] - (a * k * dt + b)) / (a * k * dt + b - 273) /
                                ((DBMemo1 -> Lines -> Count - 1) * Table1DT1 -> Value) * dt;
                        } else {
                            S = S + pow(T1[n] - (a * n * dt + b), 2) /
                                ((DBMemo1 -> Lines -> Count - 1) * Table1DT1 -> Value) * dt;
                            Skv = Skv + pow(T1[n] - (a * k * dt + b), 2) /
                                ((DBMemo1 -> Lines -> Count - 1) * Table1DT1 -> Value) * dt;
                            Sr = Sr + fabs(T1[n] - (a * k * dt + b)) / (a * k * dt + b - 273) /
                                ((DBMemo1 -> Lines -> Count - 1) * Table1DT1 -> Value) * dt;
                        }
                    } else {
                        k++;
                        if (k < n)
                            S = S + (aa * pow((T1[k] - (a * k * dt + b)) /
                                (Tex1[0] - TMAX1), 2) + bb * pow(((T1[k] - T1[k - 1]) -
                                ((a * k * dt + b) - (a * (k - 1) * dt + b))) / dTMAX1, 2)) * dt;
                        else
                            S = S + (aa * pow((T1[n] - (a * n * dt + b)) /
                                (Tex1[0] - TMAX1), 2) + bb * pow(((T1[n] - T1[n - 1]) -
                                ((a * n * dt + b) - (a * (n - 1) * dt + b))) / dTMAX1, 2)) * dt;
                    }
                }
            }
        }
        // Расчет второй стадии 2_1
        if (CB2_1 -> Checked) {
            Xmd_1[0] = 0;
            T2_1[0] = Tex2_1[0];
            so = ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1);
            Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * T2_1[0]));
            Kf = A - B * exp(C * so);
            n = ModRas(ml2, Kmd, Cac1, 0.25, Kted, m0d1, vt + vd1, k2, b2, E1,
                E2, so, md1, mr + ms + md1, Xmd_1, T2_1, gh2_1, Tex2_2[0], & Tmax,
                th2_1, false);
            k = 0;
            for (j = 0; j < DBMemo2 -> Lines -> Count - 1; j++) {
                a = (Tex2_1[j + 1] - Tex2_1[j]) / Table1DT2_1 -> Value;
                b = Tex2_1[j] - a * (j * Table1DT2_1 -> Value);
                for (ii = 0; ii < (Table1DT2_1 -> Value / dt); ii++) {
                    if (RadioButton1 -> Checked) {
                        k++;
                        if (k < n) {
                            S = S + pow(T2_1[k] - (a * k * dt + b), 2) /
                                ((DBMemo2 -> Lines -> Count - 1) * Table1DT2_1 -> Value) * dt;
                            Skv = Skv + pow(T2_1[k] - (a * k * dt + b), 2) /
                                ((DBMemo2 -> Lines -> Count - 1) * Table1DT2_1 -> Value) * dt;
                            Sr = Sr + fabs(T2_1[k] - (a * k * dt + b)) / (a * k * dt + b - 273) /
                                ((DBMemo2 -> Lines -> Count - 1) * Table1DT2_1 -> Value) * dt;
                        } else {
                            S = S + pow(T2_1[n] - (a * n * dt + b), 2) /
                                ((DBMemo2 -> Lines -> Count - 1) * Table1DT2_1 -> Value) * dt;
                            Skv = Skv + pow(T2_1[n] - (a * k * dt + b), 2) /
                                ((DBMemo2 -> Lines -> Count - 1) * Table1DT2_1 -> Value) * dt;
                            Sr = Sr + fabs(T2_1[n] - (a * k * dt + b)) / (a * k * dt + b - 273) /
                                ((DBMemo2 -> Lines -> Count - 1) * Table1DT2_1 -> Value) * dt;
                        }
                    } else {
                        k++;
                        if (k < n)
                            S = S + (aa * pow((T2_1[k] - (a * k * dt + b)) /
                                (Tex2_1[0] - TMAX2_1), 2) + bb * pow(((T2_1[k] - T2_1[k - 1]) -
                                ((a * k * dt + b) - (a * (k - 1) * dt + b))) / dTMAX2_1, 2)) * dt;
                        else
                            S = S + (aa * pow((T2_1[n] - (a * n * dt + b)) /
                                (Tex2_1[0] - TMAX2_1), 2) + bb * pow(((T2_1[n] - T2_1[n - 1]) -
                                ((a * n * dt + b) - (a * (n - 1) * dt + b))) / dTMAX2_1, 2)) * dt;
                    }
                }
            }
        }
        if (CB2_2 -> Checked) {
            // Расчет второй стадии 2_2
            Xmd_2[0] = 0;
            T2_2[0] = Tex2_2[0];
            so = (ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1) +
                md1 * 1 / (ms + mr + md1)) * (ms + mr + md1) / (ms + mr + md1 + md2);
            Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * T2_2[0]));
            Kf = A - B * exp(C * so);
            n = ModRas(ml2, Kmd, Cac2, 0.25, Kted, m0d2, vt + vd1 + vd2, k2,
                b2, E1, E2, so, md2, mr + ms + md1 + md2, Xmd_2, T2_2, gh2_2, Tex2_2[0], &
                Tmax, th2_2, true);
            k = 0;
            for (j = 0; j < DBMemo3 -> Lines -> Count - 1; j++) {
                a = (Tex2_2[j + 1] - Tex2_2[j]) / Table1DT2_2 -> Value;
                b = Tex2_2[j] - a * (j * Table1DT2_2 -> Value);
                for (ii = 0; ii < (Table1DT2_2 -> Value / dt); ii++) {
                    if (RadioButton1 -> Checked) {
                        k++;
                        if (k < n) {
                            S = S + pow(T2_2[k] - (a * k * dt + b), 2) /
                                ((DBMemo3 -> Lines -> Count - 1) * Table1DT2_2 -> Value) * dt;
                            Skv = Skv + pow(T2_2[k] - (a * k * dt + b), 2) /
                                ((DBMemo3 -> Lines -> Count - 1) * Table1DT2_2 -> Value) * dt;
                            Sr = Sr + fabs(T2_2[k] - (a * k * dt + b)) / (a * k * dt + b - 273) /
                                ((DBMemo3 -> Lines -> Count - 1) * Table1DT2_2 -> Value) * dt;
                        } else {
                            S = S + pow(T2_2[n] - (a * n * dt + b), 2) /
                                ((DBMemo3 -> Lines -> Count - 1) * Table1DT2_2 -> Value) * dt;
                            Skv = Skv + pow(T2_2[k] - (a * k * dt + b), 2) /
                                ((DBMemo3 -> Lines -> Count - 1) * Table1DT2_2 -> Value) * dt;
                            Sr = Sr + fabs(T2_2[k] - (a * k * dt + b)) / (a * k * dt + b - 273) /
                                ((DBMemo3 -> Lines -> Count - 1) * Table1DT2_2 -> Value) * dt;
                        }
                    } else {
                        k++;
                        if (k < n)
                            S = S + (aa * pow((T2_2[k] - (a * k * dt + b)) /
                                (Tex2_2[0] - TMAX2_2), 2) + bb * pow(((T2_2[k] - T2_2[k - 1]) -
                                ((a * k * dt + b) - (a * (k - 1) * dt + b))) / dTMAX2_2, 2)) * dt;
                        else
                            S = S + (aa * pow((T2_2[n] - (a * n * dt + b)) /
                                (Tex2_2[0] - TMAX2_2), 2) + bb * pow(((T2_2[n] - T2_2[n - 1]) -
                                ((a * n * dt + b) - (a * (n - 1) * dt + b))) / dTMAX2_2, 2)) * dt;
                    }
                }
            }
        }
        // среднеквадратическое отклонение
        Edit41 -> Text = FloatToStr(pow(Skv, 0.5));
        // среднее относительное отклонение
        Edit3 -> Text = FloatToStr(Sr * 100);
        // Освобождаем память
        delete Xy;
        delete Xms;
        delete Xmd_1;
        delete Xmd_2;
        delete T1;
        delete T2_1;
        delete T2_2;
    }
    return S;
}
// ------------------------------------------------------------------------
// кнопка "Идентификация"
void __fastcall TForm1::Button4Click(TObject * Sender) {
    Button6Click(this);
}
// ------------------------------------------------------------------------
// Назад на одну запись
void __fastcall TForm1::BitBtn3Click(TObject * Sender) {
    if (Table1 -> RecordCount != 0) {
        Table1 -> Edit();
        Table1 -> Prior();
        if (Table1 -> Bof)
            Table1 -> Last();
    } else
        Table1 -> Insert();
    Edit45 -> Text = Table1J0_2 -> Value / 2;
    // Заполняем поля формы расчёта T0
    Edit9n -> Text = FloatToStr(Table1CAC1 -> Value);
    Edit10n -> Text = FloatToStr(Table1CAC1 -> Value);
    Edit9n -> Text = FloatToStr(Table1CAC1 -> Value);
    Edit14n -> Text = DBMemo2 -> Lines -> operator[](0);
    Edit1n -> Text = DBMemo3 -> Lines -> operator[](0);
}
// ------------------------------------------------------------------------
void __fastcall TForm1::BitBtn2Click(TObject * Sender) {
    Table1 -> Delete(); // удаление текущей записи из таблицы
}
// ------------------------------------------------------------------------
// Вперед на одну запись
void __fastcall TForm1::BitBtn4Click(TObject * Sender) {
    if (Table1 -> RecordCount != 0) {
        Table1 -> Edit();
        Table1 -> Next();
        if (Table1 -> Eof) // если последняя запись, то переходим на первую
            Table1 -> First();
    } else
        Table1 -> Insert();
    Edit45 -> Text = Table1J0_2 -> Value / 2;
    // Заполняем поля формы расчёта T0
    Edit9n -> Text = FloatToStr(Table1CAC1 -> Value);
    Edit10n -> Text = FloatToStr(Table1CAC1 -> Value);
    Edit9n -> Text = FloatToStr(Table1CAC1 -> Value);
    Edit14n -> Text = DBMemo2 -> Lines -> operator[](0);
    Edit1n -> Text = DBMemo3 -> Lines -> operator[](0);
}
// ------------------------------------------------------------------------
// кнопка "Добавить запись"
void __fastcall TForm1::BitBtn1Click(TObject * Sender) {
    Table1 -> Insert();
}
// ------------------------------------------------------------------------
void __fastcall TForm1::N1Click(TObject * Sender) {
    // Предварительная очистка ListBox-ов
    FSelect -> ListBox1 -> Clear();
    FSelect -> ListBox2 -> Clear();
    FSelect -> Top = Form1 -> Top;
    FSelect -> Left = Form1 -> Left;
    PageControl2 -> ActivePage = TabSheet3;
    // заполняем ListBox1 на форме выбора записями из базы данных
    Table1 -> First();
    while (!Table1 -> Eof) {
        FSelect -> ListBox1 -> Items -> Add(Table1NAME -> Value);
        Table1 -> Next();
    }
    Table1 -> First();
    FSelect -> ShowModal();
    // фокусируемся на первом элементе ListBox-а
    FSelect -> ListBox1 -> ItemIndex = 0;;
}
// ------------------------------------------------------------------------
void __fastcall TForm1::N2Click(TObject * Sender) {
    bool fRepeat = false;
    FSelect -> Top = Form1 -> Top;
    FSelect -> Left = Form1 -> Left;
    PageControl2 -> ActivePage = TabSheet3;
    FSelect -> ListBox2 -> Items = ListBox1 -> Items;
    // заполняем ListBox1 на форме выбора записями из базы данных
    FSelect -> ListBox1 -> Clear();
    Table1 -> First();
    while (!Table1 -> Eof) {
        for (int i = 0; i < FSelect -> ListBox2 -> Count; i++)
            if (FSelect -> ListBox2 -> Items -> Strings[i] == Table1NAME -> Value) {
                fRepeat = true;
                i = FSelect -> ListBox2 -> Count;
            }
        else
            fRepeat = false;
        if (!fRepeat)
            FSelect -> ListBox1 -> Items -> Add(Table1NAME -> Value);
        Table1 -> Next();
    }
    Table1 -> First();
    FSelect -> ShowModal();
}
// ------------------------------------------------------------------------
// контекстное меню
void __fastcall TForm1::PopupMenu1Popup(TObject * Sender) {
    if (ListBox1 -> Items -> Count != NULL) {
        N2 -> Enabled = true;
        N3 -> Enabled = true;
        N4 -> Enabled = true;
    } else {
        N2 -> Enabled = false;
        N3 -> Enabled = false;
        N4 -> Enabled = false;
    }
}
// ------------------------------------------------------------------------
void __fastcall TForm1::N4Click(TObject * Sender) {
    // заполняем ListBox1 на форме выбора записями из базы данных
    Table1 -> First();
    while (!Table1 -> Eof) {
        FSelect -> ListBox1 -> Items -> Add(Table1NAME -> Value);
        Table1 -> Next();
    }
    Table1 -> First();
    ListBox1 -> Clear();
}
// ------------------------------------------------------------------------
void __fastcall TForm1::N3Click(TObject * Sender) {
    if (ListBox1 -> ItemIndex >= 0) {
        FSelect -> ListBox1 -> Items -> Add(ListBox1 -> Items ->
            Strings[ListBox1 -> ItemIndex]);
        ListBox1 -> Items -> Delete(ListBox1 -> ItemIndex);
    }
}
// ------------------------------------------------------------------------
void __fastcall TForm1::N5Click(TObject * Sender) {
    bool fRepeat = false; // флаг нужен для проверки,есть ли в ListBox-е текущаяапись Edita
    for (int i = 0; i < ListBox1 -> Count; i++)
        if (ListBox1 -> Items -> Strings[i] == DBEdit25 -> Text) {
            fRepeat = true;
            break;
        }
    else
        fRepeat = false;
    if (!fRepeat)
        ListBox1 -> Items -> Add(DBEdit25 -> Text);
}
// ------------------------------------------------------------------------
void __fastcall TForm1::Button1Click(TObject * Sender) {
    if (ListBox1 -> Count != 0) {
        dt = StrToFloat(Edit7 -> Text); // шаг дискретизации
        Ky0 = StrToFloat(Edit13 -> Text); // коэффициенты идентификации
        Ei = StrToFloat(Edit52 -> Text);
        ml1 = StrToFloat(Edit9 -> Text);
        ml2 = StrToFloat(Edit11 -> Text);
        E1 = StrToFloat(Edit18 -> Text);
        E2 = StrToFloat(Edit20 -> Text);
        // вывод графиков
        if (Ky0 == 0) {
            ModRas_Identif(ml1, ml2, true);
            Edit2 -> Text = ModRas_Identif(ml1, ml2);
        } else {
            ModRas_Identif(Ky0, Ei, ml1, ml2, E1, E2, true);
            Edit2 -> Text = ModRas_Identif(Ky0, Ei, ml1, ml2, E1, E2);
        }
    } else
        ShowMessage("Выберете опыт");
}
// ------------------------------------------------------------------------
// двойной щелчок на графике
void __fastcall TForm1::Chart1DblClick(TObject * Sender) {
    Chart2 -> Visible = !Chart2 -> Visible;
    Chart1 -> Visible = !Chart1 -> Visible;
}
// ------------------------------------------------------------------------
void __fastcall TForm1::N7Click(TObject * Sender) {
    cold = false;
    for (int i = 0; i < Chart1 -> SeriesCount(); i++)
        Chart1 -> SeriesList -> operator[](i) -> Clear();
    for (int i = 0; i < Chart2 -> SeriesCount(); i++)
        Chart2 -> SeriesList -> operator[](i) -> Clear();
    for (int i = 0; i < Chart3 -> SeriesCount(); i++)
        Chart3 -> SeriesList -> operator[](i) -> Clear();
    Chart1 -> SeriesList -> Clear();
    Chart2 -> SeriesList -> Clear();
    Chart3 -> SeriesList -> Clear();
}
// ------------------------------------------------------------------------
// функция «заморозить»
void __fastcall TForm1::N6Click(TObject * Sender) {
    cold = true;
}
// ------------------------------------------------------------------------
void __fastcall TForm1::Button6Click(TObject * Sender) {
    float S = 0, S1 = 0;
    int p, q, c, i; // счётчики
    if (ListBox1 -> Count != 0) {
        Ky0 = StrToFloat(Edit13 -> Text); // коэффициент идентификации Ky0
        Ei = StrToFloat(Edit52 -> Text); // коэффициент идентификации Ei
        ml1 = StrToFloat(Edit9 -> Text); // коэффициент идентификации ml1
        ml2 = StrToFloat(Edit11 -> Text); // коэффициент идентификации ml2
        E1 = StrToFloat(Edit18 -> Text); // коэффициент идентификации Е1
        E2 = StrToFloat(Edit20 -> Text); // коэффициент идентификации Е2
        dml1 = StrToFloat(Edit10 -> Text); // шаг изменения коэффициента ml1
        dml2 = StrToFloat(Edit12 -> Text); // шаг изменения коэффициента ml2
        dKy0 = StrToFloat(Edit14 -> Text); // шаг изменения коэффициента Ky0
        dEi = StrToFloat(Edit53 -> Text); // шаг изменения коэффициента Ei
        dE1 = StrToFloat(Edit19 -> Text); // шаг изменения E1
        dE2 = StrToFloat(Edit23 -> Text); // шаг изменения E2
        dt = StrToFloat(Edit7 -> Text); // шаг дискретизации
        // первоначальный расчёт
        S = ModRas_Identif(Ky0, Ei, ml1, ml2, E1, E2);
        q = c = 0;
        while (1) {
            p = 0; // флаг попадания в цикл
            // ml1+
            if (CB_ml -> Checked) {
                if (CB1 -> Checked) {
                    while (1) {
                        S1 = ModRas_Identif(Ky0, Ei, ml1 + dml1, ml2, E1, E2);
                        if (S1 < S) {
                            S = S1;
                            ml1 += dml1;
                            q++;
                            p = 1;
                        } else
                            break;
                    }
                }
                if (CB2_1 -> Checked || CB2_2 -> Checked) {
                    // ml2+
                    while (1) {
                        S1 = ModRas_Identif(Ky0, Ei, ml1, ml2 + dml2, E1, E2);
                        if (S1 < S) {
                            S = S1;
                            ml2 += dml2;
                            q++;
                            p = 1;
                        } else
                            break;
                    }
                }
            }
            // Ky+
            if (CB_Ky -> Checked) {
                if (CB1 -> Checked) {
                    // Ei+
                    while (1) {
                        S1 = ModRas_Identif(Ky0, Ei + dEi, ml1, ml2, E1, E2);
                        if (S1 < S) {
                            S = S1;
                            Ei += dEi;
                            q++;
                            p = 1;
                        } else
                            break;
                    }
                    // Ky0+
                    while (1) {
                        S1 = ModRas_Identif(Ky0 + dKy0, Ei, ml1, ml2, E1, E2);
                        if (S1 < S) {
                            S = S1;
                            Ky0 += dKy0;
                            q++;
                            p = 1;
                        } else
                            break;
                    }
                }
            }
            // E1+
            if (CB_E -> Checked) {
                if (CB2_1 -> Checked || CB2_2 -> Checked) {
                    while (1) {
                        S1 = ModRas_Identif(Ky0, Ei, ml1, ml2, E1 + dE1, E2);
                        if (S1 < S) {
                            S = S1;
                            E1 += dE1;
                            q++;
                            p = 1;
                        } else
                            break;
                    }
                    // E2+
                    while (1) {
                        S1 = ModRas_Identif(Ky0, Ei, ml1, ml2, E1, E2 + dE2);
                        if (S1 < S) {
                            S = S1;
                            E2 += dE2;
                            q++;
                            p = 1;
                        } else
                            break;
                    }
                }
            }
            if (p == 0) {
                // ml1-
                if (CB_ml -> Checked) {
                    if (CB1 -> Checked) {
                        while (1) {
                            S1 = ModRas_Identif(Ky0, Ei, ml1 - dml1, ml2, E1, E2);
                            if (S1 < S) {
                                S = S1;
                                ml1 -= dml1;
                                q++;
                                p = 1;
                            } else
                                break;
                        }
                    }
                    if (CB2_1 -> Checked || CB2_2 -> Checked) {
                        // ml2-
                        while (1) {
                            S1 = ModRas_Identif(Ky0, Ei, ml1, ml2 - dml2, E1, E2);
                            if (S1 < S) {
                                S = S1;
                                ml2 -= dml2;
                                q++;
                                p = 1;
                            } else
                                break;
                        }
                    }
                }
                // Kyif(
                CB_Ky -> Checked) {
                if (CB1 -> Checked) {
                    // Eiwhile
                    (1) {
                        S1 = ModRas_Identif(Ky0, Ei - dEi, ml1, ml2, E1, E2);
                        if (S1 < S) {
                            S = S1;
                            Ei -= dEi;
                            q++;
                            p = 1;
                        } else
                            break;
                    }
                    // Ky0-
                    while (1) {
                        S1 = ModRas_Identif(Ky0 - dKy0, Ei, ml1, ml2, E1, E2);
                        if (S1 < S) {
                            S = S1;
                            Ky0 -= dKy0;
                            q++;
                            p = 1;
                        } else
                            break;
                    }
                }
            }
            if (CB_E -> Checked) {
                if (CB2_1 -> Checked || CB2_2 -> Checked) {
                    // E1-
                    while (1) {
                        S1 = ModRas_Identif(Ky0, Ei, ml1, ml2, E1 - dE1, E2);
                        if (S1 < S) {
                            S = S1;
                            E1 -= dE1;
                            q++;
                            p = 1;
                        } else
                            break;
                    }
                    // E2-
                    while (1) {
                        S1 = ModRas_Identif(Ky0, Ei, ml1, ml2, E1, E2 - dE2);
                        if (S1 < S) {
                            S = S1;
                            E2 -= dE2;
                            q++;
                            p = 1;
                        } else
                            break;
                    }
                }
            }
        }
        if (p == 0) { // Меняем шаг поиска настройки коэффициентов
            dml1 = dml1 / 2;
            dml2 = dml2 / 2;
            dKy0 = dKy0 / 2;
            dEi = dEi / 2;
            dE1 = dE1 / 2;
            dE2 = dE2 / 2;
            c++;
            if ((S <= e) || (c > 20))
                break;
        }
    }
    // количество итераций
    Edit3 -> Text = q;
    // величина критерия
    Edit2 -> Text = S;
    // коэффициенты идентификации
    if (CB_ml -> Checked) {
        if (CB1 -> Checked)
            Edit1 -> Text = ml1;
        if (CB2_1 -> Checked || CB2_2 -> Checked)
            Edit16 -> Text = ml2;
    }
    if (CB_Ky -> Checked) {
        if (CB1 -> Checked) {
            Edit15 -> Text = Ky0;
            Edit54 -> Text = Ei;
        }
    }
    if (CB_E -> Checked) {
        if (CB2_1 -> Checked || CB2_2 -> Checked) {
            Edit24 -> Text = E1;
            Edit35 -> Text = E2;
        }
    }
    // вывод графиков
    ModRas_Identif(Ky0, Ei, ml1, ml2, E1, E2, true);
    }
    else
        ShowMessage("Выберете опыт");
}
// ------------------------------------------------------------------------
// собственная функция возведения в степень
float TForm1::MyPow(float X, float Y) {
    float Answer;
    if (X < 1e-12 && X > -1e-12)
        Answer = 0;
    else
        Answer = pow(X, Y);
    return Answer;
}
// ------------------------------------------------------------------------
void __fastcall TForm1::CheckBox1Click(TObject * Sender) {
    if (CB_ml -> Checked) {
        Edit9DblClick(this);
        Edit11DblClick(this);
    }
    if (CB_Ky -> Checked) {
        Edit13DblClick(this);
        Edit52DblClick(this);
    }
    if (CB_E -> Checked) {
        Edit18DblClick(this);
        Edit20DblClick(this);
    }
}
// ------------------------------------------------------------------------
// функция расчёта по модели с целью идентификации J0
float TForm1::Identif_KAC(float J0) {
    int i, j, ii, n, k;
    float S, a, b;
    float R = 8.32;
    S = 0;
    // выделяем память под массивы данных
    Xy = new float[NN];
    Xms = new float[NN];
    T1 = new float[NN];
    // Инициализация параметров
    Ini();
    // Расчет первой стадии
    Xy[0] = 0;
    Xms[0] = 0;
    T1[0] = Tex1[0];
    so = 0;
    Kms = Km0 * exp(-Ems / (R * T1[0]));
    Ky = Ky0 * exp(-Ei / (R * T1[0]));
    Kf = A - B * exp(C * so);
    n = ModRas(ml1, Ky, Kms, J0, 0.5, Ktes, m0s, vt, Ky0, Ei, Km0, 0, Ems, so, ms, mr + ms,
        Xy, Xms, T1, gh1, Tex2_1[0], & Tmax, th1, false, 0);
    k = 0;
    for (j = 0; j < N_t; j++) {
        a = (Tex1[j + 1] - Tex1[j]) / Table1DT1 -> Value;
        b = Tex1[j] - a * (j * Table1DT1 -> Value);
        for (ii = 0; ii < (Table1DT1 -> Value / dt); ii++) {
            if (k == 0)
                S = S + pow(T1[k] - (a * k * dt + b), 2) * dt;
            k++;
            if (k < n)
                S = S + pow(T1[k] - (a * k * dt + b), 2) * dt;
            else
                S = S + pow(T1[n] - (a * n * dt + b), 2) * dt;
        }
    }
    // Освобождаем память
    delete Xy;
    delete Xms;
    delete T1;
    return S;
}
// -----------------------------------------------------------------------
void __fastcall TForm1::Button9Click(TObject * Sender) {
    int g = 0, c = 0, p;
    float J, J1; // первоначальный критерий
    float dK = 0.001;
    float dJ_J0_1, dJ_J0_2, H_J0, nabla;
    // ввод исходных данных
    N_t = StrToInt(Edit42 -> Text); // интервал идентификации
    dt = StrToFloat(Edit43 -> Text); // шаг дискретизации
    J0i = StrToFloat(Edit45 -> Text); // первоначальная КАЦ
    dJ0 = StrToFloat(Edit44 -> Text); // шаг изменения КАЦ
    e = StrToFloat(Edit48 -> Text); // точность оптимизации
    Ky0 = StrToFloat(Edit15 -> Text); // коэффициент идентификации Ky
    ml1 = StrToFloat(Edit1 -> Text); // коэффициент идентификации ml1
    Ei = StrToFloat(Edit54 -> Text);
    ml2 = StrToFloat(Edit16 -> Text); // коэффициент идентификации ml2
    E1 = StrToFloat(Edit24 -> Text); // коэффициент идентификации Е1
    E2 = StrToFloat(Edit35 -> Text); // коэффициент идентификации Е2
    if (RG -> ItemIndex == 0) {
        // первоначальный расчёт критерия
        J = Identif_KAC(J0i);
        // расчёт частных производных по искомым параметрам
        dJ_J0_1 = (Identif_KAC(J0i + dK) - J) / dK;
        // расчитываем величину "набла"
        nabla = sqrt(pow(dJ_J0_1, 2));
        // определяем величины приращения коэффициентов
        H_J0 = dJ0;
        // осуществляем приращения коэффициентов
        J0i = J0i - H_J0 * dJ_J0_1 / nabla;
        // увеличиваем счётчик
        g++;
        // начинаем цикл
        while (nabla > e) {
            // первоначальный расчёт критерия
            J = Identif_KAC(J0i);
            // расчёт частных производных по искомым параметрам
            dJ_J0_2 = (Identif_KAC(J0i + dK) - J) / dK;
            // расчитываем величину "набла"
            nabla = sqrt(pow(dJ_J0_2, 2));
            if (nabla < 0.001)
                break;
            // определяем величины приращения коэффициентов
            if (dJ_J0_2 * dJ_J0_1 > 0)
                H_J0 = 2 * H_J0;
            else
                H_J0 = 1 / 3.0 * H_J0;
            // осуществляем приращения коэффициентов
            J0i = J0i - H_J0 * dJ_J0_2 / nabla;
            // переприсваиваем значения частных производных
            dJ_J0_1 = dJ_J0_2;
            // увеличиваем счётчик
            g++;
            if (g > 300)
                break;
        }
    } else {
        // первоначальный расчёт
        J = Identif_KAC(J0i);
        g = c = 0;
        while (1) {
            p = 0; // флаг попадания в цикл
            // J0i+
            while (1) {
                J1 = Identif_KAC(J0i + dK);
                if (J1 < J) {
                    J = J1;
                    J0i += dK;
                    g++;
                    p = 1;
                } else
                    break;
            }
            if (p == 0) {
                // J0iwhile
                (1) {
                    J1 = Identif_KAC(J0i - dK);
                    if (J1 < J) {
                        J = J1;
                        J0i -= dK;
                        g++;
                        p = 1;
                    } else
                        break;
                }
            }
            if (p == 0) { // Меняем шаг поиска настройки коэффициентов
                dK = dK / 2;
                c++;
                if ((J <= e) || (c > 20))
                    break;
            }
        }
    }
    // количество итераций
    Edit49 -> Text = g;
    // точность
    Edit50 -> Text = nabla;
    // значение коэффициента
    Edit46 -> Text = J0i;
    // величина критерия
    Edit47 -> Text = Identif_KAC(J0i);
    // относительное отклонение
    Edit51 -> Text = fabs(J0 - J0i) / J0 * 100;
    // вывод графиков
    Identif_KAC(J0i, true);
    J0 = J0i;
}
// ------------------------------------------------------------------------
// функция расчёта по модели с целью идентификации J0 (графики)
void TForm1::Identif_KAC(float J0i, bool f) {
    int i, j, ii, n;
    float S, t;
    float R = 8.32, Tmax,
        tt; // для запоминания времени окончания стадии
    if (!cold) {
        Chart1 -> SeriesList -> Clear();
        Chart2 -> SeriesList -> Clear();
    } else
        cold = !cold;
    t = 0;
    // описываем графики
    TLineSeries * Tr1 = new TLineSeries(this);
    Tr1 -> ParentChart = Chart2;
    Tr1 -> Title = "Tрасч_дейст.КАЦ";
    Tr1 -> LinePen -> Width = 2;
    TLineSeries * Tr2 = new TLineSeries(this);
    Tr2 -> ParentChart = Chart2;
    Tr2 -> Title = "Tрасч_найд.КАЦ";
    Tr2 -> LinePen -> Width = 2;
    TPointSeries * Te = new TPointSeries(this);
    Te -> ParentChart = Chart2;
    Te -> Title = "Tэксп";
    Te -> Pointer -> Style = psCircle;
    TLineSeries * X1 = new TLineSeries(this);
    X1 -> ParentChart = Chart1;
    X1 -> Title = "Xрасч_дейст.КАЦ";
    X1 -> LinePen -> Width = 2;
    TLineSeries * X2 = new TLineSeries(this);
    X2 -> ParentChart = Chart1;
    X2 -> Title = "Xрасч_найд.КАЦ";
    X2 -> LinePen -> Width = 2;
    TLineSeries * XX1 = new TLineSeries(this);
    XX1 -> ParentChart = Chart1;
    XX1 -> Title = "Xy_расч_дейст.КАЦ";
    XX1 -> LinePen -> Width = 2;
    TLineSeries * XX2 = new TLineSeries(this);
    XX2 -> ParentChart = Chart1;
    XX2 -> Title = "Xy_расч_найд.КАЦ";
    XX2 -> LinePen -> Width = 2;
    // выделяем память под массивы данных
    Xy = new float[NN];
    Xms = new float[NN];
    T1 = new float[NN];
    // Инициализация параметров
    Ini();
    // для действительной КАЦ
    Xy[0] = 0;
    Xms[0] = 0;
    T1[0] = Tex1[0];
    so = 0;
    Kms = Km0 * exp(-Ems / (R * T1[0]));
    Ky = Ky0 * exp(-Ei / (R * T1[0]));
    Kf = A - B * exp(C * so);
    n = ModRas(ml1, Ky, Kms, J0, 0.5, Ktes, m0s, vt, Ky0, Ei, Km0, 0, Ems, so, ms, mr + ms, Xy,
        Xms, T1, gh1, Tex2_1[0], & Tmax, th1, false, 0);
    for (j = 0; j <= n; j++) {
        XX1 -> AddXY(t, Xy[j], "", clDefault);
        X1 -> AddXY(t, Xms[j], "", clDefault);
        Tr1 -> AddXY(t, T1[j], "", clDefault);
        if (j < DBMemo1 -> Lines -> Count)
            Te -> AddXY(j * Table1DT1 -> Value, Tex1[j], "", clDefault);
        t += dt;
    }
    // для найденного значения КАЦ
    Xy[0] = 0;
    Xms[0] = 0;
    T1[0] = Tex1[0];
    so = 0;
    t = 0;
    Kms = Km0 * exp(-Ems / (R * T1[0]));
    Ky = Ky0 * exp(-Ei / (R * T1[0]));
    Kf = A - B * exp(C * so);
    n = ModRas(ml1, Ky, Kms, J0i, 0.5, Ktes, m0s, vt, Ky0, Ei, Km0, 0, Ems, so, ms, mr + ms, Xy,
        Xms, T1, gh1, Tex2_1[0], & Tmax, th1, false, 0);
    for (j = 0; j <= n; j++) {
        XX2 -> AddXY(t, Xy[j], "", clDefault);
        X2 -> AddXY(t, Xms[j], "", clDefault);
        Tr2 -> AddXY(t, T1[j], "", clDefault);
        t += dt;
    }
    // Освобождаем память
    delete Xy;
    delete Xms;
    delete T1;
}
// кнопка «Расчёт» панели «Оценка КАЦ»
void __fastcall TForm1::Button10Click(TObject * Sender) {
    Ky0 = StrToFloat(Edit15 -> Text); // коэффициент идентификации Ky
    ml1 = StrToFloat(Edit1 -> Text); // коэффициент идентификации ml1
    Ei = StrToFloat(Edit54 -> Text);
    ml2 = StrToFloat(Edit16 -> Text); // коэффициент идентификации ml2
    E1 = StrToFloat(Edit24 -> Text); // коэффициент идентификации Е1
    E2 = StrToFloat(Edit35 -> Text); // коэффициент идентификации Е2
    if (StrToFloat(Edit46 -> Text) != NULL) {
        J0i = StrToFloat(Edit46 -> Text); // найденная КАЦ
        ModRas_Identif(J0i, Ky0, Ei, ml1, ml2, E1, E2, true); // расчёт сайденной КАЦ
    }
}
// ------------------------------------------------------------------------
void __fastcall TForm1::UpDown1Click(TObject * Sender, TUDBtnType Button) {
    UpDown1 -> Max = DBMemo1 -> Lines -> Count - 1;
}
// ------------------------------------------------------------------------
// функция расчёта по модели с целью идентификации параметра J0(графики)
void TForm1::ModRas_Identif(float J0, float Ky0, float Ei, float ml1, float ml2, float E1, float E2, bool graf) {
    int i, j, ii, n;
    float S, t;
    float R = 8.32, Tmax,
        tt; // для запоминания времени окончания стадии
    S = 0;
    if (!cold) {
        Chart1 -> SeriesList -> Clear();
        Chart2 -> SeriesList -> Clear();
    } else
        cold = !cold;
    // описываем графики
    TPointSeries * Te = new TPointSeries(this);
    Te -> ParentChart = Chart2;
    Te -> Title = "Tэксп" + Table1NAME -> Value;
    Te -> Pointer -> Style = psCircle;
    TLineSeries * Tr = new TLineSeries(this);
    Tr -> ParentChart = Chart2;
    Tr -> Title = "Tрасч" + Table1NAME -> Value;
    Tr -> LinePen -> Width = 2;
    TLineSeries * X = new TLineSeries(this);
    X -> ParentChart = Chart1;
    X -> Title = "Xрасч" + Table1NAME -> Value;
    X -> LinePen -> Width = 2;
    TLineSeries * XX = new TLineSeries(this);
    XX -> ParentChart = Chart1;
    XX -> Title = "Xy_расч" + Table1NAME -> Value;
    XX -> LinePen -> Width = 2;
    // выделяем память под массивы данных
    Xy = new float[NN];
    Xms = new float[NN];
    Xmd_1 = new float[NN];
    Xmd_2 = new float[NN];
    T1 = new float[NN];
    T2_1 = new float[NN];
    T2_2 = new float[NN];
    // Инициализация параметров
    Ini();
    // Расчет первой стадии
    if (CheckBox2 -> Checked) {
        Xy[0] = 0;
        Xms[0] = 0;
        T1[0] = Tex1[0];
        so = 0;
        Kms = Km0 * exp(-Ems / (R * T1[0]));
        Ky = Ky0 * exp(-Ei / (R * T1[0]));
        Kf = A - B * exp(C * so);
        n = ModRas(ml1, Ky, Kms, J0, 0.5, Ktes, m0s, vt, Ky0, Ei, Km0, 0, Ems, so, ms,
            mr + ms, Xy, Xms, T1, gh1, Tex2_1[0], & Tmax, th1, false, 0);
        if (graf)
            for (j = 0; j <= n; j++) {
                XX -> AddXY(t, Xy[j], "", clDefault);
                X -> AddXY(t, Xms[j], "", clDefault);
                Tr -> AddXY(t, T1[j], "", clDefault);
                if (j < DBMemo1 -> Lines -> Count)
                    Te -> AddXY(j * Table1DT1 -> Value, Tex1[j], "", clDefault);
                t += dt;
            }
    }
    // Расчет второй стадии 2_1
    if (CheckBox3 -> Checked) {
        Xmd_1[0] = 0;
        T2_1[0] = Tex2_1[0];
        so = ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1);
        Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * T2_1[0]));
        Kf = A - B * exp(C * so);
        n = ModRas(ml2, Kmd, Cac1, 0.25, Kted, m0d1, vt + vd1, k2, b2, E1, E2, so,
            md1, mr + ms + md1, Xmd_1, T2_1, gh2_1, Tex2_2[0], & Tmax, th2_1,
            false);
        tt = t; // запоминаем время окончания стадии
        if (graf)
            for (j = 0; j <= n; j++) {
                X -> AddXY(t, Xmd_1[j], "", clDefault);
                Tr -> AddXY(t, T2_1[j], "", clDefault);
                if (j < DBMemo2 -> Lines -> Count)
                    Te -> AddXY(tt + j * Table1DT2_1 -> Value, Tex2_1[j], "",
                        clDefault);
                t += dt;
            }
    }
    if (CheckBox4 -> Checked) {
        // Расчет второй стадии 2_2
        Xmd_2[0] = 0;
        T2_2[0] = Tex2_2[0];
        so = (ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1) +
            md1 * 1 / (ms + mr + md1)) * (ms + mr + md1) / (ms + mr + md1 + md2);
        Kf = A - B * exp(C * so);
        Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * T2_2[0]));
        n = ModRas(ml2, Kmd, Cac2, 0.25, Kted, m0d2, vt + vd1 + vd2, k2, b2,
            E1, E2, so, md2, mr + ms + md1 + md2, Xmd_2, T2_2, gh2_2, Tex2_2[0], &
            Tmax, th2_2, true);
        tt = t;
        if (graf)
            for (j = 0; j <= n; j++) {
                X -> AddXY(t, Xmd_2[j], "", clDefault);
                Tr -> AddXY(t, T2_2[j], "", clDefault);
                if (j < DBMemo3 -> Lines -> Count)
                    Te -> AddXY(tt + j * Table1DT2_2 -> Value, Tex2_2[j], "",
                        clDefault);
                t += dt;
            }
    }
    // Освобождаем память
    delete Xy;
    delete Xms;
    delete Xmd_1;
    delete Xmd_2;
    delete T1;
    delete T2_1;
    delete T2_2;
}
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// при выходе из программы
void __fastcall TForm1::FormClose(TObject * Sender, TCloseAction & Action) {
    delete Tex1;
    delete Tex2_1;
    delete Tex2_2;
    delete MMRg;
    delete MMRgNM;
    delete XX;
    delete TT;
    delete XXi1;
    delete XXi2;
    delete Tx;
    delete Tr;
    delete Tet;
    delete Xg;
    delete Xet;
}
// ------------------------------------------------------------------------
// кнопка "Определение управляющего воздействия"
void __fastcall TForm1::ButClick(TObject * Sender) {
    float S = 0;
    int p, q, c, i, j; // счётчики
    float t;
    // присваиваем начальное значение переменным для расчёта дисперсии
    N11 = N2_1 = N2_2 = N33 = 0;
    int n1, n2, n3; // количество итераций расчёт по модели
    // присваиваем начальное значение переменной one
    one = false;
    // Исходные данные
    Ky0 = StrToFloat(Edit15 -> Text); // коэффициент идентификации Ky0
    Ei = StrToFloat(Edit54 -> Text); // коэффициент идентификации Ei
    ml1 = StrToFloat(Edit1 -> Text); // коэффициент идентификации ml1
    dt = StrToFloat(Edit17 -> Text); // шаг дискретизации
    ml2 = StrToFloat(Edit16 -> Text); // коэффициент идентификации
    E1 = StrToFloat(Edit24 -> Text); // коэффициент идентификации E 1
    E2 = StrToFloat(Edit35 -> Text); // коэффициент идентификации E2
    /*работа с графиками*/
    if (!cold) {
        Chart1 -> SeriesList -> Clear();
        Chart2 -> SeriesList -> Clear();
        Chart3 -> SeriesList -> Clear();
    } else
        cold = !cold;
    float k, Mn, Mw, disp;
    // ЭТАЛОННЫЙ ПРОФИЛЬ
    n1 = RasUpr(StrToFloat(Edit63 -> Text), StrToFloat(DBEdit34 -> Text), &
        k, & Mn, & Mw, & disp, 1);
    Edit68 -> Text = FloatToStrF(Mn, ffGeneral, 7, 0);
    Edit69 -> Text = FloatToStrF(Mw, ffGeneral, 7, 0);
    Edit70 -> Text = FloatToStrF(k, ffGeneral, 5, 3);
    Edit77 -> Text = FloatToStrF(disp, ffGeneral, 7, 0);
    // ПРОФИЛЬ БЕЗ УПРАВЛЕНИЯ
    n2 = RasUpr(StrToFloat(Edit63 -> Text), StrToFloat(DBEdit34 -> Text), &
        k, & Mn, & Mw, & disp, 2);
    Edit71 -> Text = FloatToStrF(Mn, ffGeneral, 7, 0);
    Edit72 -> Text = FloatToStrF(Mw, ffGeneral, 7, 0);
    Edit73 -> Text = FloatToStrF(k, ffGeneral, 5, 3);
    Edit78 -> Text = FloatToStrF(disp, ffGeneral, 7, 0);
    // ПРОФИЛЬ С УПРАВЛЕНИЕМ
    if (StrToFloat(Edit63 -> Text) < StrToFloat(DBEdit34 -> Text))
        n3 = RasUpr(StrToFloat(Edit63 -> Text), StrToFloat(DBEdit34 -> Text), &
            k, & Mn, & Mw, & disp, 3);
    else
        n3 = RasUpr(StrToFloat(Edit63 -> Text), StrToFloat(DBEdit34 -> Text), &
            k, & Mn, & Mw, & disp, 4);
    Edit74 -> Text = FloatToStrF(Mn, ffGeneral, 7, 0);
    Edit75 -> Text = FloatToStrF(Mw, ffGeneral, 7, 0);
    Edit76 -> Text = FloatToStrF(k, ffGeneral, 5, 3);
    Edit82 -> Text = FloatToStrF(disp, ffGeneral, 7, 0);
}
// ------------------------------------------------------------------------
// функция расчёта управляющего воздействия при изменении КАЦ
int __fastcall TForm1::RasUpr(float J0, float J0d, float * kk, float * Mn,
    float * Mw, float * disp, int sw) {
    float vts = 0;
    bool m1bool, m2bool, m11bool, m22bool; /*булевские переменные, необходимыеля организации вычислений минимального и максимального значимых значенийаспределения */
    bool one;
    bool enm; // булевская переменная, отвечающая за учёт низкомолекулярныхракций при расчёте ММР по второму способу управления
    int i, j, ii, s, sk = 1, n, iz = 1, n1, kz = 0, q;
    int m1 = 0, m2 = 0, m11 = 0, m22 = 0, mSr1 = 0, mSr2 = 0; // переменные для масштабирова-ия распределений
    float dif, mdif = 0, t1 = 0;
    float min = 0, max = 0, dmas;
    float kp = 1; // количество прогонов
    float R = 8.32,
        Nm = 0; // количество деактивированных молекул АЦ (моль)
    // объявляем переменную для хранения интервала разбиения
    float dp = 0;
    // объявляем переменную для хранения количества фракций в сшитом продукте
    int nf = 0;
    // считываем значения
    dp = StrToFloat(Edit64 -> Text);
    nf = StrToInt(Edit65 -> Text);
    kp = StrToInt(Edit90 -> Text);
    // объявляем переменную для хранения шага разбиения ММР сшитого продуктао массе
    float t = 0;
    // описываем цвета графиков
    TColor cl[4];
    cl[0] = clRed;
    cl[1] = clBlack;
    cl[2] = cl[3] = clGreen;
    // выделяем память под массивы данных
    Xy = new float[NN];
    Xyd = new float[NN];
    Xms = new float[NN];
    Xmd_1 = new float[NN];
    Xmd_2 = new float[NN];
    T1 = new float[NN];
    T2_1 = new float[NN];
    T2_2 = new float[NN];
    // описываем графики
    MMRg = new TLineSeries(this);
    MMRg -> ParentChart = Chart3;
    MMRg -> Title = "MMP " + IntToStr(sw) + "_" + Table1NAME -> Value;
    MMRg -> SeriesColor = cl[sw - 1];
    MMRg -> LinePen -> Width = 2;
    MMRgNM = new TLineSeries(this);
    MMRgNM -> ParentChart = Chart3;
    MMRgNM -> SeriesColor = cl[sw - 1];
    MMRgNM -> LinePen -> Width = 2;
    MMRan = new TLineSeries(this);
    MMRan -> ParentChart = Chart3;
    MMRan -> Title = "MMPсгл" + IntToStr(sw) + "_" + Table1NAME -> Value;
    MMRan -> LinePen -> Width = 3;
    MMRan -> SeriesColor = cl[sw - 1];
    XX = new TLineSeries(this);
    XX -> ParentChart = Chart1;
    XX -> Title = "X_" + IntToStr(sw) + "_" + Table1NAME -> Value;
    XX -> LinePen -> Width = 2;
    XX -> SeriesColor = cl[sw - 1];
    TT = new TLineSeries(this);
    TT -> ParentChart = Chart2;
    TT -> Title = "T_" + IntToStr(sw) + "_" + Table1NAME -> Value;
    TT -> LinePen -> Width = 2;
    TT -> SeriesColor = cl[sw - 1];
    // конверсия инициатора
    XXi1 = new TLineSeries(this);
    XXi1 -> ParentChart = Chart1;
    XXi1 -> Title = "Х иниц.1 " + IntToStr(sw) + "_" + Table1NAME -> Value;
    XXi1 -> LinePen -> Width = 2;
    XXi1 -> SeriesColor = cl[sw - 1];
    XXi2 = new TLineSeries(this);
    XXi2 -> ParentChart = Chart1;
    XXi2 -> Title = "Х иниц.2 " + IntToStr(sw) + "_" + Table1NAME -> Value;
    XXi2 -> LinePen -> Width = 2;
    XXi2 -> SeriesColor = cl[sw - 1];
    i = 1; // параметр, требуемый для формирования шага расчёта ММР
    // Инициализация параметров
    Ini();
    n1 = StrToInt(Edit87 -> Text);
    int i0 = StrToInt(Edit89 -> Text);
    enm = CheckBox10 -> Checked;
    double ** pp = new double * [n1 + 1];
    double * mm = new double[NN];
    double * nmk = new double[NN];
    for (i = 0; i < n1 + 1; i++)
        pp[i] = new double[NN];
    // задаём начальные условия
    pp[0][0] = 0;
    mm[0] = m0s;
    for (i = 1; i < n1 + 1; i++) {
        mm[i] = 0;
        pp[i][0] = 0;
    }
    if (CB1u -> Checked) {
        // Расчет первой стадии
        Xy[0] = 0;
        // Обнуляем добавочный массив конверсии инициатора
        for (ii = 0; ii < 25 / dt; ii++) {
            Xyd[ii] = 0;
            nmk[ii] = 0;
        }
        Xms[0] = 0;
        T1[0] = Tex1[0];
        so = 0;
        Kms = Km0 * exp(-Ems / (R * T1[0]));
        Ky = Ky0 * exp(-Ei / (R * T1[0]));
        Kf = A - B * exp(C * so);
        switch (sw) {
        case 1: {
            n = ModRas_V(pp, mm, ml1, Ky, Kms, J0d, 0.5, Ktes, m0s, vt, Ky0, Ei,
                Km0, 0, Ems, 0, so, ms, mr + ms, Xy, Xms, T1, gh1, Tex2_1[0], & Tmax, th1,
                true, 0, n1);
            break;
        }
        case 2: {
            n = ModRas_V(pp, mm, ml1, Ky, Kms, J0, 0.5, Ktes, m0s, vt, Ky0, Ei,
                Km0, 0, Ems, 0, so, ms, mr + ms, Xy, Xms, T1, gh1, Tex2_1[0], & Tmax, th1,
                true, 0, n1);
            break;
        }
        case 3: {
            Edit66 -> Text = FloatToStr(J0d - J0);
            n = ModRasV1(pp, mm, ml1, Ky, Kms, J0, J0d, 0.5, Ktes, m0s, vt, Ky0,
                Ei, Km0, 0, Ems, 0, so, ms, mr + ms, Xy, Xyd, Xms, T1, gh1, Tex2_1[0], &
                Tmax, th1, true, i0, n1);
            break;
        }
        case 4: {
            msd = RasMSD(J0d, J0);
            mrd = mr / ms * msd;
            Edit67 -> Text = FloatToStr(msd);
            Edit61 -> Text = FloatToStr(mrd);
            n = ModRasV2(pp, nmk, mm, ml1, Ky, Kms, J0, J0d, 0.5, Ktes, m0s, vt,
                Ky0, Ei, Km0, 0, Ems, 0, so, ms, mr + ms, msd, mrd, Xy, Xms, T1, gh1,
                Tex2_1[0], & Tmax, th1, true, i0, n1, enm);
            break;
        }
        };
        // вывод графиков первой стадии
        for (j = 0; j <= n; j++) {
            XXi1 -> AddXY(t, Xy[j], "", cl[sw - 1]);
            if (sw == 3)
                XXi2 -> AddXY(t, Xyd[j], "", cl[sw - 1]);
            XX -> AddXY(t, Xms[j], "", cl[sw - 1]);
            TT -> AddXY(t, T1[j], "", cl[sw - 1]);
            t += dt;
        }
        /*вычисляем минимальное и максимальное значимые величины распределения,аралельно сдвигая весь массив распределений по времени*/
        m1bool = m2bool = false;
        for (int k = 0; k <= n1; k++) {
            // ищем минимальное значение
            if (pp[k][n] > dp && !m1bool) {
                m1 = k;
                m1bool = true;
            }
            // ищем максимальное значение
            if (pp[k][n] < dp && m1bool && !m2bool) {
                m2 = k;
                m2bool = true;
            }
            // переприсвоение
            pp[k][0] = pp[k][n];
        }
    }
    // РАСЧЁТ СТАДИИ 2_1
    if (CB2_1u -> Checked) {
        // начальные условия
        mm[0] = m0d1;
        Xy[0] = 1;
        for (i = 1; i < NN; i++) {
            mm[i] = 0;
        }
        Xmd_1[0] = 0;
        T2_1[0] = Tex2_1[0];
        so = ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1);
        Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * T2_1[0]));
        Kf = A - B * exp(C * so);
        switch (sw) {
        case 1: {
            n = ModRas_V(pp, mm, ml2, 0, Kmd, Cac1, 0.25, Kted, m0d1, vt + vd1,
                0, 0, k2, b2, E1, E2, so, md1, mr + ms + md1, Xy, Xmd_1, T2_1,
                gh2_1, Tex2_2[0], & Tmax, th2_1, true, 0, n1);
            break;
        }
        case 2: {
            J0 = J0 * vt / (vt + md1 / pd);
            n = ModRas_V(pp, mm, ml2, 0, Kmd, J0, 0.25, Kted, m0d1, vt + vd1,
                0, 0, k2, b2, E1, E2, so, md1, mr + ms + md1, Xy, Xmd_1, T2_1, gh2_1,
                Tex2_2[0], & Tmax, th2_1, true, 0, n1);
            break;
        }
        case 3: {
            n = ModRasV1(pp, mm, ml2, 0, Kmd, Cac1, Cac1, 0.25, Kted,
                m0d1, vt + vd1, 0, 0, k2, b2, E1, E2, so, md1, mr + ms + md1, Xy, Xyd, Xmd_1,
                T2_1, gh2_1, Tex2_2[0], & Tmax, th2_1, true, 0, n1);
            break;
        }
        case 4: {
            /*расчёт добавочных значений бутадиена*/
            md1d = 1.0 * md1 * (msd / ms);
            md2d = 1.0 * md2 * (msd / ms);
            Edit55 -> Text = FloatToStr(md1d);
            Edit56 -> Text = FloatToStr(md2d);
            so = (ms + msd) * 1 / (mr + ms + msd + mrd) * (ms + mr + msd + mrd) /
                (ms + msd + mrd + mr + md1 + md1d);
            T2_1[0] = T2_1[0] - 1;
            Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * T2_1[0]));
            Kf = A - B * exp(C * so);
            m0d1 = (md1 + md1d) / mmd / (vt + msd / ps + mrd / pr + (md1 + md1d) / pd);
            mm[0] = m0d1;
            J0 = J0d * (vt + msd / ps + mrd / pr) / (vt + msd / ps + mrd / pr + (md1 + md1d) / pd);
            n = ModRasV2(pp, nmk, mm, ml2, 0, Kmd, J0, Cac1, 0.25, Kted, m0d1,
                vt + msd / ps + mrd / pd + vd1 + md1d / pd, 0, 0, k2, b2, E1, E2,
                so, md1 + md1d, mr + ms + mrd + msd + md1 + md1d, msd, mrd, Xy, Xmd_1,
                T2_1, gh2_1, Tex2_2[0], & Tmax, th2_1, true, 0, n1, enm);
            break;
        }
        };
        // вывод графиков второй стадии (первая порция дивинила)
        for (j = 0; j <= n; j++) {
            XX -> AddXY(t, Xmd_1[j], "", cl[sw - 1]);
            TT -> AddXY(t, T2_1[j], "", cl[sw - 1]);
            t += dt;
        }
    }
    // РАСЧЁТ ВТОРОЙ СТАДИИ 2_2
    if (CB2_2u -> Checked) {
        // начальные условия
        for (int k = 0; k <= n1; k++)
            pp[k][0] = pp[k][n];
        mm[0] = m0d2;
        for (i = 1; i < NN; i++)
            mm[i] = 0;
        Xmd_2[0] = 0;
        T2_2[0] = Tex2_2[0];
        so = (ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1) +
            md1 * 1 / (ms + mr + md1)) * (ms + mr + md1) / (ms + mr + md1 + md2);
        Kf = A - B * exp(C * so);
        Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * T2_2[0]));
        switch (sw) {
        case 1: {
            n = ModRas_V(pp, mm, ml2, 0, Kmd, Cac2, 0.25, Kted, m0d2,
                vt + vd1 + vd2, 0, 0, k2, b2, E1, E2, so, md2, mr + ms + md1 + md2,
                Xy, Xmd_2, T2_2, gh2_2, Tex2_2[0], & Tmax, th2_2, true, 0, n1);
            break;
        }
        case 2: {
            J0 = J0 * (vt + md1 / pd) / (vt + md1 / pd + md2 / pd);
            n = ModRas_V(pp, mm, ml2, 0, Kmd, J0, 0.25, Kted, m0d2,
                vt + vd1 + vd2, 0, 0, k2, b2, E1, E2, so, md2, mr + ms + md1 + md2,
                Xy, Xmd_2, T2_2, gh2_2, Tex2_2[0], & Tmax, th2_2, true, 0, n1);
            break;
        }
        case 3: {
            n = ModRasV1(pp, mm, ml2, 0, Kmd, Cac2, Cac2, 0.25, Kted, m0d2, vt, 0, 0,
                k2, b2, E1, E2, so, md2, mr + ms + md1 + md2, Xy, Xyd, Xmd_2, T2_2,
                gh2_2, Tex2_2[0], & Tmax, th2_2, true, i0, n1);
            break;
        }
        case 4: {
            // изменяем процентное содержание сухого остатка за счётобавления второй порции бутадиена
            so = ((ms + msd) * 1 / (mr + ms + mrd + msd) * (ms + mr + mrd + msd) /
                    (ms + mr + mrd + msd + md1 + md1d) + (md1 + md1d) * 1 /
                    (ms + mr + md1 + md1d + mrd + msd)) * (ms + mr + md1 + mrd + msd) /
                (ms + mr + msd + mrd + md1 + md2 + md1d + md2d);
            Kf = A - B * exp(C * so);
            T2_2[0] = T2_2[0] - 1;
            Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * T2_2[0]));
            m0d2 = (md2 + md2d) / mmd / (vt + msd / ps + mrd / pr +
                (md1 + md1d + md2 + md2d) / pd);
            mm[0] = m0d2;
            J0 = J0 * (vt + msd / ps + mrd / pr + (md1 + md1d) / pd) /
                (vt + msd / ps + mrd / pr + (md1 + md1d + md2 + md2d) / pd);
            n = ModRasV2(pp, nmk, mm, ml2, 0, Kmd, J0, Cac2, 0.25, Kted, m0d2,
                vt + msd / ps + mrd / pr + vd1 + vd2 + (md1d + md2d) / pd, 0, 0, k2, b2, E1, E2,
                so, md2 + md2d, mr + ms + msd + mrd + md1 + md2 + md1d + md2d, msd,
                mrd, Xy, Xmd_2, T2_2, gh2_2, Tex2_2[0], & Tmax, th2_2, true, i0, n1, enm);
            break;
        }
        };
        // вывод графиков для стадии 2_2
        for (j = 0; j <= n; j++) {
            XX -> AddXY(t, Xmd_2[j], "", cl[sw - 1]);
            TT -> AddXY(t, T2_2[j], "", cl[sw - 1]);
            t += dt;
        }
    }
    /*вычисляем минимальное и максимальное значимые величины распределения*/
    m11bool = m22bool = false;
    for (int k = 0; k <= n1; k++) {
        // ищем минимальное значение
        if (pp[k][n] > dp && !m11bool) {
            m11 = k;
            m11bool = true;
        }
        // ищем максимальное значение
        if (pp[k][n] < dp && m11bool && !m22bool) {
            m22 = k;
            m22bool = true;
        }
    }
    float P = 0, Mns = 0, Mws1 = 0, Mws2 = 0, ds = 0;
    float minM = 0;
    int km = 0; // индекс массива pp[][], которому соответствует среднечис-енная молекулярная масса
    if (CB2_1u -> Checked) {
        MMRg -> AddXY(0, 0, "", cl[sw - 1]);
        for (int k = m11; k <= m22; k++) {
            kz = ((k - m11) * (m2 - m1) + m1 * (m22 - m11)) / (m22 - m11);
            MMRg -> AddXY(kz * mms * 1000 +
                (k - kz) * mmd * 1000, pp[k][n] * 1000, "", cl[sw - 1]);
            P = P + pp[k][n];
            Mns = Mns + (kz * mms * 1000 + (k - kz) * mmd * 1000) * pp[k][n];
            Mws1 = Mws1 + pow(kz * mms * 1000 + (k - kz) * mmd * 1000, 2) * pp[k][n];
            Mws2 = Mws2 + (kz * mms * 1000 + (k - kz) * mmd * 1000) * pp[k][n];
        }
        MMRg -> AddXY(m2 * mms * 1000 + m22 * mmd * 1000, 0, "", cl[sw - 1]);
    } else {
        for (int k = 0; k <= n1; k++) {
            MMRg -> AddXY(k * mms * 1000, pp[k][n] * 1000, "", cl[sw - 1]);
            MMRgNM -> AddXY(k * mms * 1000, nmk[k] * 1000, "", cl[sw - 1]);
            P = P + pp[k][n];
            Mns = Mns + (k * mms * 1000) * pp[k][n];
            Mws1 = Mws1 + pow(k * mms * 1000, 2) * pp[k][n];
            Mws2 = Mws2 + (k * mms * 1000) * pp[k][n];
        }
    }
    * Mn = Mns / P;
    * Mw = Mws1 / Mws2;
    * kk = Mws1 / Mws2 / (Mns / P);
    // Расчёт соеднеквадратического отклонения
    ds = 0;
    if (CB2_1u -> Checked) {
        // определение концентрации (моль/л) цепочек со среднечисленной моле-улярной массой
        minM = * Mn;
        for (int k = m11; k <= m22; k++) {
            kz = ((k - m11) * (m2 - m1) + m1 * (m22 - m11)) / (m22 - m11);
            if (fabs((kz * mms * 1000 + (k - kz) * mmd * 1000) - * Mn) < minM) {
                km = k;
                minM = fabs((kz * mms * 1000 + (k - kz) * mmd * 1000) - * Mn);
            }
        }
        for (int k = m11; k <= m22; k++) {
            kz = ((k - m11) * (m2 - m1) + m1 * (m22 - m11)) / (m22 - m11);
            ds = ds + pow((kz * mms * 1000 + (k - kz) * mmd * 1000) * pp[k][n] /
                pp[km][n] - * Mn, 2);
        }
        if (N2_1 == 0)
            N2_1 = m22 - m11;
        ds = pow(ds / (N2_1 - 1), 2);
    } else {
        // определение количества цепочек со среднечисленной молекулярнойассой
        minM = * Mn;
        for (int k = 0; k <= n1; k++) {
            kz = ((k - m11) * (m2 - m1) + m1 * (m22 - m11)) / (m22 - m11);
            if (fabs((k * mms * 1000) - * Mn) < minM) {
                km = k;
                minM = fabs((k * mms * 1000) - * Mn);
            }
        }
        for (int k = 0; k <= n1; k++)
            ds = ds + pow((kz * mms * 1000) * pp[k][n] / pp[km][n] - * Mn, 2);
        if (N11 == 0)
            N11 = n1;
        ds = pow(ds / (N11 - 1), 0.5);
    }
    // СТАДИЯ СШИВКИ!!!!!
    // объявляем и выделяем память под массив целочисленных элементов
    int * KolCh = new int[m22 - m11 + 1];
    float * MolMas, // указатель на массив, хранящий данные о всехепочках (последовательно)
        * NewMolMas, // указатель на массив, хранящий данные о всехшитых цепочках
        * cF, // указатель на массив для хранения мольныхонцентраций каждой фракции
        * mF; // указатель на массив для хранения молекулярныхесов фракций (оба массива предназначеныля ММР сшитого полимера)
    int KolEl = 0; // счётчик количества элементов в массиве MolMas
    int ss = 0; // счётчик элементов массива MolMas
    float S = 0; // для накопления суммы
    int nm = 0; // индекс количества элементов в массиве сшитых молекул
    float mu1 = 0, mu2 = 0; // моменты нормального логарифмическогоаспределения
    if (CB3 -> Checked) {
        // Блок формирования массива
        for (int k = m11; k <= m22; k++) {
            KolCh[k - m11] = pp[k][n] / dp; // определяем количество элементарныхракций в каждой фракции
            KolEl = KolEl + KolCh[k - m11]; // определяем общее количестволементарных фракций
        }
        // динамически формируем массив для хранения всех цепочек
        MolMas = new float[KolEl];
        // динамически формируем массив для хранения всех сшитых цепочек
        NewMolMas = new float[KolEl / 4 + 4];
        // а также массив для хранениия концентраций фракций
        cF = new float[KolEl / 4 + 4];
        // и массив для хранения молекулярных весов фракций
        mF = new float[KolEl / 4 + 4];
        for (int k = m11; k <= m22; k++) {
            kz = ((k - m11) * (m2 - m1) + m1 * (m22 - m11)) / (m22 - m11); // пересчётасштаба одного распределения в другое
            for (int kk = 1; kk <= KolCh[k - m11]; kk++) {
                MolMas[ss] = (kz * mms + (k - kz) * mmd) * 1000; // определяемоличество элементарных фракций в каждой фракции
                ss++;
            }
        }
        // Блок реализации случайного механизма взаимодействия молекул
        nm = 0; // начальное значение индекса массива сшитых молекул
        randomize();
        while (KolEl >= 4) {
            S = 0;
            for (int four = 1; four <= 4; four++) {
                i = random(KolEl - 1);
                S = S + MolMas[i];
                for (int k = i; k < (KolEl - 1); k++)
                    MolMas[k] = MolMas[k + 1];
                KolEl--;
            }
            NewMolMas[nm] = S;
            nm++;
        }
        // Блок формирования ММР сшитого продукта
        ss = 0;
        min = 4 * (m1 * mms + (m11 - m1) * mmd) * 1000;
        max = 6 * (m2 * mms + (m22 - m2) * mmd) * 1000;
        dmas = (max - min) / nf;
        for (float mas = min; mas < max; mas = mas + dmas) {
            S = 0;
            q = 0;
            for (int k = 0; k < nm; k++)
                if (NewMolMas[k] > mas && NewMolMas[k] < mas + dmas) {
                    S = S + NewMolMas[k];
                    q++;
                }
            if (S != 0) {
                S = S / q;
                cF[ss] = dp * q;
                mF[ss] = S;
            } else {
                cF[ss] = 0;
                mF[ss] = (2 * mas + dmas) / 2;
            }
            ss++;
        }
        // Блок вывода ММР сшитого продукта
        // вывод нефильтрованного графика
        for (int k = 0; k < ss; k++)
            if (CheckBox11 -> Checked)
                MMRg -> AddXY(mF[k], cF[k] * 1000, "", cl[sw - 1]);
        // добавляем усреднённое решение
        for (int or = 1; or <= kp; or++)
            for (int k = 0; k < ss - 1; k++) {
                cF[k] = (cF[k] + cF[k + 1]) / 2;
                mF[k] = (mF[k] + mF[k + 1]) / 2;
            }
        // расчёт характеристик
        P = 0;
        Mns = 0;
        Mws1 = 0;
        Mws2 = 0;
        for (int k = 0; k < ss; k++) {
            P = P + cF[k];
            Mns = Mns + mF[k] * cF[k];
            Mws1 = Mws1 + pow(mF[k], 2) * cF[k];
            Mws2 = Mws2 + (mF[k]) * cF[k];
        }
        * Mn = Mns / P;
        * Mw = Mws1 / Mws2;
        * kk = Mws1 / Mws2 / (Mns / P);
        // расчёт дисперсии
        ds = 0;
        // определение концентрации цепочек со среднечисленной молекулярнойассой
        minM = * Mn;
        for (int k = 0; k < ss; k++) {
            if (fabs(mF[k] - * Mn) < minM) {
                km = k;
                minM = fabs(mF[k] - * Mn);
            }
        }
        j = 0;
        for (int k = 0; k < ss; k++) {
            if (cF[k] > 0) {
                ds = ds + pow((mF[k] - * Mn) * cF[k] / cF[km], 2);
                j++;
            }
        }
        if (N33 == 0)
            N33 = j;
        ds = pow(ds / (N33 - 1), 0.5);
        for (int k = 0; k < ss; k++)
            MMRan -> AddXY(mF[k], cF[k] * 1000, "", cl[sw - 1]);
        // очищаем память
        delete[] MolMas;
        delete[] NewMolMas;
        delete[] cF;
        delete[] mF;
    }
    // сохранение значения дисперсии
    * disp = ds;
    // Освобождаем память
    delete[] Xy;
    delete[] Xyd;
    delete[] Xms;
    delete[] T1;
    delete[] Xmd_1;
    delete[] Xmd_2;
    delete[] T2_1;
    delete[] T2_2;
    for (int k = 0; k < n1 + 1; k++)
        delete[] pp[k];
    delete[] pp;
    delete[] nmk;
    delete[] mm;
    delete[] KolCh;
    return n;
}
//----------------------------------------------------------------------
// Функция для расчёта модели и ММР при добавлении катализатор
int __fastcall TForm1::ModRasV1(
    double ** pp, // для расчёта с вероятностью
    double * mm, // для определения концентрации мономера
    float ml, // коэффициент идентификации матем. модели
    float Ky, // константа скорости инициирования
    float Km, // константа скорости роста цепи
    float J0, // КАЦ
    float J0d, // действительное значение КАЦ
    float St, // степень при КАЦ
    float Kte, // коэффициент тепловыделения
    float m0, // начальная концентрация мономера
    float vt, // текущий объем реакционной смеси
    float Ky0, // постоянная константы скорости инииирования
    float Ei, // энергия активации реакции инициирования
    float k, // коэффициент модели расчета Кm
    float b, // коэффициент модели расчета Кm
    float Em, // енергия активации
    float E2, // второй коэффициент энергии активации
    float so_b, // количество сухого остатка с предыдущей стадии
    float m, // масса мономера на текущей стадии
    float M, // масса реакционной смеси (все компоненты)
    float * Xy, // указатель на массив конверсии инициатора
    float * Xyd, // указатель на массив конверсии инициатора (2)
    float * Xm, // указатель на массив конверсии мономера
    float * T, // указатель на массив температуры реакции
    float gh, // расход хладагента
    float T_k, // температура начала следующей стадии
    float * Tmax,
    float th, // температура хладагента
    bool f, // булевская переменная
    int i0, // время начала расчётов при добавлении катализ.
    int n1) {
    int i = 0, n, c = 0, g = 0, j = 0;
    float L, // уровень загрузки реактора
    R = 8.32; // универсальная газовая постоянная
    float x1 = 0, y1 = 0, z1 = 0, x2 = 0, y2 = 0, z2 = 0, x3 = 0, y3 = 0, z3 = 0, x4 = 0, y4 = 0, z4 = 0; // коэф-ициенты Рунге-Кутта
    float z1d = 0, z2d = 0, z3d = 0, z4d = 0, p1 = 0, p2 = 0, p3 = 0, p4 = 0; // коэффициенты Рунге-утта для системы с учётом добавления катализатора
    L = (vt - 3000) / 28500; // уровень заполнения реактора
    float J0dob = 0;
    float Mux = 0;
    if (Ky != 0)
        Mux = pow(J0, -0.5);
    else
        Mux = pow(J0, -0.75);
    while (1) {
        if (i >= i0 / dt && !one) {
            J0dob = (J0d - J0);
            Mux = pow(J0d, -0.5);
            // меняем значение булевской переменной для закрытия входа в блокересчёта
            one = true;
        }
        z1 = Ky * m0 * (1 - Xy[i]) * (1 - Xm[i]);
        if (i >= i0 / dt) // компонент, отвечающий за начало рас-ёта следующего уравнения после прохождения i0 мин
            z1d = Ky * m0 * (1 - Xyd[i]) * (1 - Xm[i]);
        x1 = (Ky * (J0 * (1 - Xy[i]) + J0dob * (1 - Xyd[i])) +
            Km * MyPow(J0 * Xy[i] + J0dob * Xyd[i], St + ml * so)) * (1 - Xm[i]);
        // ограничитель
        if (CBLim -> Checked) {
            if (x1 > 2 * (1 - Xm[i]) / dt)
                x1 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt; // 10% от общей "массы"
        }
        y1 = (Kte * x1 * m0 * vt - (Kf * F * L * (T[i] - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        z2 = Ky * m0 * (1 - (Xy[i] + z1 * dt / 2)) * (1 - Xm[i]);
        if (i >= i0 / dt)
            z2d = Ky * m0 * (1 - (Xyd[i] + z1d * dt / 2)) * (1 - Xm[i]);
        x2 = (Ky * (J0 * (1 - (Xy[i] + z1 * dt / 2)) + J0dob * (1 - (Xyd[i] + z1d * dt / 2))) +
                Km * MyPow(J0 * (Xy[i] + z1 * dt / 2) + J0dob * (Xyd[i] + z1d * dt / 2), St + ml * so)) *
            (1 - (Xm[i] + x1 * dt / 2));
        // ограничитель
        if (CBLim -> Checked) {
            if (x2 > 2 * (1 - Xm[i]) / dt)
                x2 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y2 = (Kte * x2 * m0 * vt - (Kf * F * L * ((T[i] + y1 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        z3 = Ky * m0 * (1 - (Xy[i] + z2 * dt / 2)) * (1 - Xm[i]);
        if (i >= i0 / dt)
            z3d = Ky * m0 * (1 - (Xyd[i] + z2d * dt / 2)) * (1 - Xm[i]);
        x3 = (Ky * (J0 * (1 - (Xy[i] + z2 * dt / 2)) + J0dob * (1 - (Xyd[i] + z2d * dt / 2))) +
                Km * MyPow(J0 * (Xy[i] + z2 * dt / 2) + J0dob * (Xyd[i] + z2d * dt / 2), St + ml * so)) *
            (1 - (Xm[i] + x2 * dt / 2));
        // ограничитель
        if (CBLim -> Checked) {
            if (x3 > 2 * (1 - Xm[i]) / dt)
                x3 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y3 = (Kte * x3 * m0 * vt - (Kf * F * L * ((T[i] + y2 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        z4 = Ky * m0 * (1 - (Xy[i] + z3 * dt / 2)) * (1 - Xm[i]);
        if (i >= i0 / dt)
            z4d = Ky * m0 * (1 - (Xyd[i] + z3d * dt / 2)) * (1 - Xm[i]);
        x4 = (Ky * (J0 * (1 - (Xy[i] + z3 * dt / 2)) + J0dob * (1 - (Xyd[i] + z3d * dt / 2))) +
                Km * MyPow(J0 * (Xy[i] + z3 * dt / 2) + J0dob * (Xyd[i] + z3d * dt / 2), St + ml * so)) *
            (1 - (Xm[i] + x3 * dt / 2));
        // ограничитель
        if (CBLim -> Checked) {
            if (x4 > 2 * (1 - Xm[i]) / dt)
                x4 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y4 = (Kte * x4 * m0 * vt - (Kf * F * L * ((T[i] + y3 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        Xy[i + 1] = Xy[i] + dt * (z1 + 2 * (z2 + z3) + z4) / 6;
        if (i >= i0 / dt)
            Xyd[i + 1] = Xyd[i] + dt * (z1d + 2 * (z2d + z3d) + z4d) / 6;
        Xm[i + 1] = Xm[i] + dt * (x1 + 2 * (x2 + x3) + x4) / 6;
        T[i + 1] = T[i] + dt * (y1 + 2 * (y2 + y3) + y4) / 6;
        // вставляем блок с вероятностью
        mm[i + 1] = m0 * (1 - Xm[i + 1]); // расчёт концентрации мономера
        for (j = 0; j < n1; j++) {
            if (j == 0) {
                p1 = -Km * (Mux) * mm[i] * pp[j][i] + Ky * J0 * (1 - Xy[i]) * mm[i] +
                    Ky * J0dob * (1 - Xyd[i]) * mm[i];
                p2 = -Km * (Mux) * mm[i] * (pp[j][i]) + Ky * J0 * (1 - Xy[i]) * mm[i] +
                    Ky * J0dob * (1 - Xyd[i]) * mm[i];
                p3 = -Km * (Mux) * mm[i] * (pp[j][i]) + Ky * J0 * (1 - Xy[i]) * mm[i] +
                    Ky * J0dob * (1 - Xyd[i]) * mm[i];
                p4 = -Km * (Mux) * mm[i] * (pp[j][i]) + Ky * J0 * (1 - Xy[i]) * mm[i] +
                    Ky * J0dob * (1 - Xyd[i]) * mm[i];
            } else {
                p1 = Km * (Mux) * mm[i] * (pp[j - 1][i] - (pp[j][i]));
                p2 = Km * (Mux) * mm[i] * ((pp[j - 1][i]) - (pp[j][i]));
                p3 = Km * (Mux) * mm[i] * ((pp[j - 1][i]) - (pp[j][i]));
                p4 = Km * (Mux) * mm[i] * ((pp[j - 1][i]) - (pp[j][i]));
            }
            if (j == n1 - 1) {
                p1 = Km * (Mux) * mm[i] * pp[j - 1][i];
                p2 = Km * (Mux) * mm[i] * pp[j - 1][i];
                p3 = Km * (Mux) * mm[i] * pp[j - 1][i];
                p4 = Km * (Mux) * mm[i] * pp[j - 1][i];
            }
            pp[j][i + 1] = pp[j][i] + dt * (p1 + 2 * (p2 + p3) + p4) / 6;
            if (pp[j][i + 1] < 0)
                pp[j][i + 1] = 0;
        }
        if (T[i + 1] > * Tmax)
            *
            Tmax = T[i + 1];
        // необходимо проверить необходимость нижележащего блока
        // начало блока
        if (Xm[i + 1] > 1.001 || Xy[i + 1] > 1.001) {
            n = i + 1;
            break;
        }
        // конец блока
        so = so_b + m * Xm[i + 1] / M; // сухой остаток
        Km = k * (1 - b * so) * exp(-(Em + E2 * so) / (R * T[i + 1])); // расчёт констанытыкорости роста цепи
        Kf = A - B * exp(C * so); // коэффициент теплопередачи
        Ky = Ky0 * exp(-Ei / (R * T[i + 1])); // константа скорости инициирования
        // условие прекращения расчета для модели с текущей КАЦ
        if (Xm[i + 1] >= 0.99) {
            if (!f) {
                if (T[i + 1] <= T_k) {
                    n = i + 1;
                    break;
                }
            } else {
                c++;
                if (c > 200) {
                    n = i + 1;
                    break;
                }
            }
        }
        if (i > NN - 3) {
            n = i + 1;
            break;
        }
        i++;
    }
    return n;
}
//------------------------------------------------------------------------
// Функция для расчёта модели и ММР при добавлении шихты
int __fastcall TForm1::ModRasV2(
    double ** pp, // для расчёта с вероятностью
    double * nmk, // массив для хранения низкомолекуляр. каучука
    double * mm, // для определения концентрации мономера
    float ml, // коэффициент идентификации математ. модели
    float Ky, // константа скорости инициирования
    float Km, // константа скорости роста цепи
    float J0, // КАЦ
    float J0d, // действительное значение КАЦ
    float St, // степень при КАЦ
    float Kte, // коэффициент тепловыделения
    float m0, // начальная концентрация мономера
    float vt, // текущий объем реакционной смеси
    float Ky0, // постоянная константы скорости инииирования
    float Ei, // энергия активации реакции инициирования
    float k, // коэффициент модели расчета Кm
    float b, // коэффициент модели расчета Кm
    float Em, // енергия активации
    float E2,
    float so_b, // количество сухого остатка с предыдущей стадии
    float m, // масса мономера на текущей стадии
    float M, // масса реакционной смеси (все компоненты)
    float msd, // добавочное количество стирола
    float mrd, // добавочное количество растворителя
    float * Xy, // указатель на массив конверсии инициатора
    float * Xm, // указатель на массив конверсии мономера
    float * T, // указатель на массив температуры реакции
    float gh, // расход хладагента
    float T_k, // температура начала следующей стадии
    float * Tmax,
    float th, // температура хладагента
    bool f, // булевская переменная
    int i0, // время начала расчётов при добавлении катализат.
    int n1,
    bool enm) // булевская переменная, отвечающая за учёт низкомоле-улярных фракций при расчёте ММР
{
    int i = 0, n, c = 0, g = 0, j = 0, iz = 0;
    float L, // уровень загрузки реактора
    R = 8.32; // универсальная газовая постоянная
    float x1 = 0, y1 = 0, z1 = 0, x2 = 0, y2 = 0, z2 = 0, x3 = 0, y3 = 0, z3 = 0, x4 = 0, y4 = 0, z4 = 0; // коэф-ициенты Рунге-Кутта
    float p1 = 0, p2 = 0, p3 = 0, p4 = 0;
    L = (vt - 3000) / 28500; // уровень заполнения реактора
    float J0dob = 0, cs = 0, Nm = 0;
    float Mux = 0;
    float * buf;
    buf = new float[n1 + 1];
    if (Ky != 0)
        Mux = pow(J0, -0.5);
    else
        Mux = pow(J0, -0.75);
    while (1) {
        if (i >= i0 / dt && !one) {
            // текущее значение концентрации стирола, моль/л
            cs = m0 * (1 - Xm[i]);
            // пересчитываем значение конверсии
            Xm[i] = (m - cs * vt * mms) / (m + msd);
            // рассчёт количества деактивированных молекул
            Nm = (J0doz - J0) * msd / m; // рассчитываем процент количествоеактивированных молекул при добавлении новой порции стирола
            for (int kk = 0; kk < n1; kk++) {
                buf[kk] = pp[kk][i] * Nm / J0;
                pp[kk][i] = pp[kk][i] - buf[kk];
            }
            // пересчёт КАЦ
            J0 = J0d;
            // пересчёт объёма
            vt = vt + msd / ps + (mrd / pr);
            // пересчёт массы стирола
            m = m + msd;
            // пересчёт общей массы в реакторе
            M = M + msd + mrd;
            // пересчёт уровня заполнения реактора
            L = (vt - 3000) / 28500;
            // пересчёт температуры полимеризационной массы
            T[i] = ((ms + mr) * T[i] + (msd + mrd) * T[0]) / (ms + mr + msd + mrd);
            Mux = pow(J0 - 0.07 * J0, -0.5);
            // меняем значение булевской переменной для закрытияхода в блок пересчёта
            one = true;
        }
        z1 = Ky * m0 * (1 - Xy[i]) * (1 - Xm[i]);
        x1 = (Ky * J0 * (1 - Xy[i]) + Km * MyPow(J0 * Xy[i], St + ml * so)) * (1 - Xm[i]);
        // ограничитель
        if (CBLim -> Checked) {
            if (x1 > 2 * (1 - Xm[i]) / dt)
                x1 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt; // 10% от общей "массы"
        }
        y1 = (Kte * x1 * m0 * vt - (Kf * F * L * (T[i] - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        z2 = Ky * m0 * (1 - (Xy[i] + z1 * dt / 2)) * (1 - Xm[i]);
        x2 = (Ky * J0 * (1 - (Xy[i] + z1 * dt / 2)) +
            Km * MyPow(J0 * (Xy[i] + z1 * dt / 2), St + ml * so)) * (1 - (Xm[i] + x1 * dt / 2));
        // ограничитель
        if (CBLim -> Checked) {
            if (x2 > 2 * (1 - Xm[i]) / dt)
                x2 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y2 = (Kte * x2 * m0 * vt - (Kf * F * L * ((T[i] + y1 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        z3 = Ky * m0 * (1 - (Xy[i] + z2 * dt / 2)) * (1 - Xm[i]);
        x3 = (Ky * J0 * (1 - (Xy[i] + z2 * dt / 2)) +
            Km * MyPow(J0 * (Xy[i] + z2 * dt / 2), St + ml * so)) * (1 - (Xm[i] + x2 * dt / 2));
        // ограничитель
        if (CBLim -> Checked) {
            if (x3 > 2 * (1 - Xm[i]) / dt)
                x3 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y3 = (Kte * x3 * m0 * vt - (Kf * F * L * ((T[i] + y2 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        z4 = Ky * m0 * (1 - (Xy[i] + z3 * dt / 2)) * (1 - Xm[i]);
        x4 = (Ky * J0 * (1 - (Xy[i] + z3 * dt / 2)) +
            Km * MyPow(J0 * (Xy[i] + z3 * dt / 2), St + ml * so)) * (1 - (Xm[i] + x3 * dt / 2));
        // ограничитель
        if (CBLim -> Checked) {
            if (x4 > 2 * (1 - Xm[i]) / dt)
                x4 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y4 = (Kte * x4 * m0 * vt - (Kf * F * L * ((T[i] + y3 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        Xy[i + 1] = Xy[i] + dt * (z1 + 2 * (z2 + z3) + z4) / 6;
        Xm[i + 1] = Xm[i] + dt * (x1 + 2 * (x2 + x3) + x4) / 6;
        T[i + 1] = T[i] + dt * (y1 + 2 * (y2 + y3) + y4) / 6;
        // вставляем блок с вероятностью
        mm[i + 1] = m0 * (1 - Xm[i + 1]); // расчёт концентрации мономера
        for (j = 0; j < n1; j++) {
            if (j == 0) {
                p1 = -Km * (Mux) * mm[i] * pp[j][i] + Ky * J0 * (1 - Xy[i]) * mm[i];
                p2 = -Km * (Mux) * mm[i] * (pp[j][i]) + Ky * J0 * (1 - Xy[i]) * mm[i];
                p3 = -Km * (Mux) * mm[i] * (pp[j][i]) + Ky * J0 * (1 - Xy[i]) * mm[i];
                p4 = -Km * (Mux) * mm[i] * (pp[j][i]) + Ky * J0 * (1 - Xy[i]) * mm[i];
            } else {
                p1 = Km * (Mux) * mm[i] * (pp[j - 1][i] - (pp[j][i]));
                p2 = Km * (Mux) * mm[i] * ((pp[j - 1][i]) - (pp[j][i]));
                p3 = Km * (Mux) * mm[i] * ((pp[j - 1][i]) - (pp[j][i]));
                p4 = Km * (Mux) * mm[i] * ((pp[j - 1][i]) - (pp[j][i]));
            }
            if (j == n1 - 1) {
                p1 = Km * (Mux) * mm[i] * pp[j - 1][i];
                p2 = Km * (Mux) * mm[i] * pp[j - 1][i];
                p3 = Km * (Mux) * mm[i] * pp[j - 1][i];
                p4 = Km * (Mux) * mm[i] * pp[j - 1][i];
            }
            pp[j][i + 1] = pp[j][i] + dt * (p1 + 2 * (p2 + p3) + p4) / 6;
            if (pp[j][i + 1] < 0)
                pp[j][i + 1] = 0;
        }
        if (T[i + 1] > * Tmax)
            *
            Tmax = T[i + 1];
        // начало блока
        if (Xm[i + 1] > 1.001 || Xy[i + 1] > 1.001) {
            n = i + 1;
            break;
        }
        // конец блока
        so = so_b + m * Xm[i + 1] / M; // сухой остаток
        Km = k * (1 - b * so) * exp(-(Em + E2 * so) / (R * T[i + 1])); // расчёт констанытыкорости роста
        Kf = A - B * exp(C * so); // коэффициент теплопередачи
        Ky = Ky0 * exp(-Ei / (R * T[i + 1])); // константа скорости инициирования
        // условие прекращения расчета для модели с текущей КАЦ
        if (Xm[i + 1] >= 0.99) {
            if (!f) {
                if (T[i + 1] <= T_k) {
                    n = i + 1;
                    break;
                }
            } else {
                c++;
                if (c > 200) {
                    n = i + 1;
                    break;
                }
            }
        }
        if (i > NN - 3) {
            n = i + 1;
            break;
        }
        i++;
    }
    // вывод низкомолекулярной фракции
    if (enm)
        for (int kk = 0; kk < n1; kk++)
            nmk[kk] = buf[kk];
    delete buf;
    return n;
}
// ------------------------------------------------------------------------
// функция расчёта добавочного значения стирола
float __fastcall TForm1::RasMSD(float J0d, float J0) {
    float alfa = 0;
    float msd = 0;
    alfa = (1 / ps + mr / ms / pr);
    msd = vt * (J0 - J0d) / (alfa * J0d + vt / ms * (J0doz - J0));
    return msd;
}
// ------------------------------------------------------------------------
// контекстное меню графика
void __fastcall TForm1::N8Click(TObject * Sender) {
    Chart1 -> Visible = true;
    Chart2 -> Visible = false;
    Chart3 -> Visible = false;
}
// ------------------------------------------------------------------------
void __fastcall TForm1::N9Click(TObject * Sender) {
    Chart1 -> Visible = false;
    Chart2 -> Visible = true;
    Chart3 -> Visible = false;
}
// ------------------------------------------------------------------------
void __fastcall TForm1::N10Click(TObject * Sender) {
    Chart1 -> Visible = false;
    Chart2 -> Visible = false;
    Chart3 -> Visible = true;
}
// ------------------------------------------------------------------------
// кнопка "Расчёт"
void __fastcall TForm1::Button1nClick(TObject * Sender) {
    float S, J0;
    T_lim = StrToFloat(Edit13n -> Text);
    J0 = StrToFloat(Edit10n -> Text);
    J0et = StrToFloat(Edit9n -> Text);
    dT = StrToFloat(Edit18n -> Text);
    dt = StrToFloat(Edit7n -> Text);
    a = StrToFloat(Edit11n -> Text);
    b = StrToFloat(Edit12n -> Text);
    ml2 = StrToFloat(Edit16 -> Text); // коэффициент идентификации
    E1 = StrToFloat(Edit24 -> Text); // коэффициент идентификации E1
    E2 = StrToFloat(Edit35 -> Text); // коэффициент идентификации E2
    if (!cold) {
        Chart1 -> SeriesList -> Clear();
        Chart2 -> SeriesList -> Clear();
    } else
        cold = !cold;
    if (CBIni -> Checked) {
        T01 = StrToFloat(Edit14n -> Text);
        T02 = StrToFloat(Edit1n -> Text);
        cl = clBlack;
        S = Optima(T01, T02, J0, J0et, true, true, Edit35n, Edit41n);
    }
    if (CBOpt -> Checked) {
        T01 = StrToFloat(Edit14n -> Text);
        T02 = StrToFloat(Edit1n -> Text);
        if (Edit2n -> Text != "")
            T01 = StrToFloat(Edit2n -> Text);
        if (Edit3n -> Text != "")
            T02 = StrToFloat(Edit3n -> Text);
        cl = clBlue;
        S = Optima(T01, T02, J0, J0et, true, true, Edit15n, Edit16n);
    }
    if (CBObr -> Checked) {
        T01 = StrToFloat(Edit14n -> Text);
        T02 = StrToFloat(Edit1n -> Text);
        if (Edit19n -> Text != "")
            T01 = StrToFloat(Edit19n -> Text);
        if (Edit20n -> Text != "")
            T02 = StrToFloat(Edit20n -> Text);
        cl = clGreen;
        S = Optima(T01, T02, J0, J0et, true, true, Edit23n, Edit24n);
    }
}
// кнопка "Идентификация"
void __fastcall TForm1::Button4nClick(TObject * Sender) {
    if (!cold) {
        Chart1 -> SeriesList -> Clear();
        Chart2 -> SeriesList -> Clear();
    } else
        cold = !cold;
    Opt_Ras();
}
// кнопка "Обратный расчёт"
void __fastcall TForm1::Button6nClick(TObject * Sender) {
    if (!cold) {
        Chart1 -> SeriesList -> Clear();
        Chart2 -> SeriesList -> Clear();
    } else
        cold = !cold;
    Ob_Ras();
}
// ------------------------------------------------------------------------
// функция для расчёта оптимальной начальной температуры
float __fastcall TForm1::Optima(float T01, // начальная температура (2_1)
        float T02, // начальная температура (стадия 2_2)
        float J0, // КАЦ текущая на стадии 2_1
        float J0et, // КАЦ эталонная на стадии 2_1
        bool graf, // переменная включения/отключения графиков
        bool error, // переменная включения расчёта ошибки
        TEdit * MyEd1, // окно для сохранения ошибки
        TEdit * MyEd2) {
        int i, j, ii, n2, n1;
        float S, t1, t2, er;
        float R = 8.32,
            tt1 = 0, tt2 = 0; // для запоминания времени окончания стадии
        S = 0;
        // Присваиваем начальное значение времени полимеризации
        t1 = 0;
        t2 = 0;
        // описываем графики
        if (graf) {
            Tx = new TPointSeries(this);
            Tx -> ParentChart = Chart2;
            Tx -> Title = "Tэксп" + Table1NAME -> Value;
            Tx -> Pointer -> Style = psCircle;
            // массив T для хранения температуры при расчёте для действительной J0
            Tr = new TLineSeries(this);
            Tr -> ParentChart = Chart2;
            Tr -> Title = "Tрасч" + Table1NAME -> Value;
            Tr -> LinePen -> Width = 2;
            // массив T для хранения температуры при расчёте для эталонной J0
            Tet = new TLineSeries(this);
            Tet -> ParentChart = Chart2;
            Tet -> Title = "Tрасч" + Table1NAME -> Value;
            Tet -> LinePen -> Width = 2;
            // массив Х для хранения конверсии при расчёте для действительной J0
            Xg = new TLineSeries(this);
            Xg -> ParentChart = Chart1;
            Xg -> Title = "Xрасч" + Table1NAME -> Value;
            Xg -> LinePen -> Width = 2;
            // массив Х для хранения конверсии при расчёте для эталонной J0
            Xet = new TLineSeries(this);
            Xet -> ParentChart = Chart1;
            Xet -> Title = "Xet_расч" + Table1NAME -> Value;
            Xet -> LinePen -> Width = 2;
        }
        // выделяем память под массивы данных
        X = new float[NN];
        T2 = new float[NN];
        Xe = new float[NN];
        Te = new float[NN];
        // Инициализация параметров
        Ini();
        // Расчет стадии 2_1
        if (CB2_1n -> Checked) {
            // расчёт для текущих значений T0 и J0
            X[0] = 0;
            T2[0] = T01;
            so = ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1);
            Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * T2[0]));
            Kf = A - B * exp(C * so);
            n1 = ModRasT0(J0, Kmd, m0d1, vt + vd1, so, md1, mr + ms + md1,
                X, T2, gh2_1, T02, & Tmax, th2_1, false, 2);
            // расчёт для эталонных значений T0 и J0
            Xe[0] = 0;
            Te[0] = Tex2_1[0];
            so = ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1);
            Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * Te[0]));
            Kf = A - B * exp(C * so);
            n2 = ModRasT0(J0et, Kmd, m0d1, vt + vd1, so, md1, mr + ms + md1,
                Xe, Te, gh2_1, Tex2_2[0], & Tmax_et, th2_1, false, 2);
            // Расчёт величины критерия S и ошибки
            S = 0;
            er = 0;
            if (n1 > n2) {
                for (j = 0; j <= n2; j++) {
                    if (j > n1)
                        Xe[j] = Xe[n1];
                    S = S + pow(X[j] - Xe[j], 2);
                    if (j != 0)
                        er = er + fabs(X[j] - Xe[j]) * 100 / Xe[j];
                }
                er = er / n2;
            } else {
                for (j = 0; j <= n1; j++) {
                    if (j > n2)
                        Xe[j] = Xe[n2];
                    S = S + pow(X[j] - Xe[j], 2);
                    if (j != 0)
                        er = er + fabs(X[j] - Xe[j]) * 100 / Xe[j];
                }
                er = er / n1;
            }
            S = a * S + b * pow(Tmax - T_lim, 2); // добавляем к критерию ограничениео максимальной температуре.
            // вывод значения ошибки
            if (error)
                MyEd1 -> Text = FloatToStr(er);
            // вывод графиков
            // оптимальный профиль
            if (graf) {
                for (j = 0; j <= n1; j++) {
                    Xg -> AddXY(t1, X[j], "", cl);
                    Tr -> AddXY(t1, T2[j], "", cl);
                    t1 += dt;
                }
                // эталонный профиль
                for (j = 0; j <= n2; j++) {
                    Xet -> AddXY(t2, Xe[j], "", clRed);
                    Tet -> AddXY(t2, Te[j], "", clRed);
                    t2 += dt;
                }
            }
        }
        if (CB2_2n -> Checked) {
            // Расчет стадии 2_2
            // расчёт для текущих значений T0 и J0
            X[0] = 0;
            T2[0] = T02;
            so = (ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1) +
                md1 * 1 / (ms + mr + md1)) * (ms + mr + md1) / (ms + mr + md1 + md2);
            Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * T2[0]));
            Kf = A - B * exp(C * so);
            J0 = J0 * (vt + vd1) / (vt + vd1 + vd2); // пересчёт КАЦ для второй порцииивинила
            n1 = ModRasT0(J0, Kmd, m0d2, vt + vd1 + vd2, so, md2, mr + ms + md1 + md2,
                X, T2, gh2_2, Tex2_2[0], & Tmax, th2_2, true, 2);
            // расчёт для эталонных значений T0 и J0
            Xe[0] = 0;
            Te[0] = Tex2_2[0];
            so = (ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1) +
                md1 * 1 / (ms + mr + md1)) * (ms + mr + md1) / (ms + mr + md1 + md2);
            Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * Te[0]));
            Kf = A - B * exp(C * so);
            J0et = J0et * (vt + vd1) / (vt + vd1 + vd2); // пересчёт КАЦ для второй порцииивинила
            n2 = ModRasT0(J0et, Kmd, m0d2, vt + vd1 + vd2, so, md2, mr + ms + md1 + md2,
                Xe, Te, gh2_2, Tex2_2[0], & Tmax_et, th2_2, true, 2);
            // Расчёт величины критерия S и ошибки
            S = 0;
            er = 0;
            if (n1 > n2) {
                for (j = 0; j <= n2; j++) {
                    if (j > n1)
                        Xe[j] = Xe[n1];
                    S = S + pow(X[j] - Xe[j], 2);
                    if (j != 0)
                        er = er + fabs(X[j] - Xe[j]) * 100 / Xe[j];
                }
                er = er / n2;
            } else {
                for (j = 0; j <= n1; j++) {
                    if (j > n2)
                        Xe[j] = Xe[n2];
                    S = S + pow(X[j] - Xe[j], 2);
                    if (j != 0)
                        er = er + fabs(X[j] - Xe[j]) * 100 / Xe[j];
                }
                er = er / n1;
            }
            S = a * S + b * pow(Tmax - T_lim, 2); // добавляем к критерию ограничениео максимальной температуре.
            // вывод значения ошибки
            if (error)
                MyEd2 -> Text = FloatToStr(er);
            // вывод графиков
            // оптимальный профиль
            if (graf) {
                for (j = 0; j <= n1; j++) {
                    Xg -> AddXY(t1, X[j], "", cl);
                    Tr -> AddXY(t1, T2[j], "", cl);
                    t1 += dt;
                }
                // эталонный профиль
                for (j = 0; j <= n2; j++) {
                    Xet -> AddXY(t2, Xe[j], "", clRed);
                    Tet -> AddXY(t2, Te[j], "", clRed);
                    t2 += dt;
                }
            }
        }
        // Освобождаем память
        delete X;
        delete T2;
        delete Xe;
        delete Te;
        return S;
    }
    // ------------------------------------------------------------------------
    -- --
// функция метода покоординатного спуска для нахождения Т0
void __fastcall TForm1::Opt_Ras() {
    float S = 0, S1 = 0, e, J0;
    int p, q, c; // счётчики
    // получение исходных параметров
    T01 = StrToFloat(Edit14n -> Text);
    T02 = StrToFloat(Edit1n -> Text);
    T_lim = StrToFloat(Edit13n -> Text);
    J0 = StrToFloat(Edit10n -> Text);
    J0et = StrToFloat(Edit9n -> Text);
    dT = StrToFloat(Edit18n -> Text);
    dt = StrToFloat(Edit7n -> Text);
    a = StrToFloat(Edit11n -> Text);
    b = StrToFloat(Edit12n -> Text);
    e = StrToFloat(Edit6n -> Text);
    ml2 = StrToFloat(Edit16 -> Text); // коэффициент идентификации
    E1 = StrToFloat(Edit24 -> Text); // коэффициент идентификации E1
    E2 = StrToFloat(Edit35 -> Text); // коэффициент идентификации E2
    // идентификация T0 на стадии 2_1
    if (CB2_1n -> Checked) {
        // первоначальный расчёт
        S = Optima(T01, T02, J0, J0et, false, false, Edit35n, Edit41n);
        c = 0;
        while (1) {
            p = 0; // флаг попадания в цикл
            // T01+
            while (1) {
                S1 = Optima(T01 + dT, T02, J0, J0et, false, false, Edit35n, Edit41n);
                if (S1 < S) {
                    S = S1;
                    T01 += dT;
                    p = 1;
                } else
                    break;
            }
            // T01-
            if (p == 0) {
                while (1) {
                    S1 = Optima(T01 - dT, T02, J0, J0et, false, false, Edit35n, Edit41n);
                    if (S1 < S) {
                        S = S1;
                        T01 -= dT;
                        p = 1;
                    } else
                        break;
                }
            }
            if (p == 0) { // Меняем шаг поиска настроек
                dT = dT / 2;
                c++;
                if ((S <= e) || (c > 5))
                    break;
            }
        }
    }
    // идентификация T0 на стадии 2_2
    if (CB2_2n -> Checked) {
        // первоначальный расчёт
        S = Optima(T01, T02, J0, J0et, false, false, Edit35n, Edit41n);
        c = 0;
        while (1) {
            p = 0; // флаг попадания в цикл
            // T02+
            while (1) {
                S1 = Optima(T01, T02 + dT, J0, J0et, false, false, Edit35n, Edit41n);
                if (S1 < S) {
                    S = S1;
                    T02 += dT;
                    p = 1;
                } else
                    break;
            }
            // T01-
            if (p == 0) {
                while (1) {
                    S1 = Optima(T01, T02 - dT, J0, J0et, false, false, Edit35n, Edit41n);
                    if (S1 < S) {
                        S = S1;
                        T02 -= dT;
                        p = 1;
                    } else
                        break;
                }
            }
            if (p == 0) { // Меняем шаг поиска настроек
                dT = dT / 2;
                c++;
                if ((S <= e) || (c > 5))
                    break;
            }
        }
    }
    // вывод оптимальных значений начальной температуры для 1 и 2 порций диви-ила
    if (CB2_1n -> Checked)
        Edit2n -> Text = FloatToStr(T01);
    if (CB2_2n -> Checked)
        Edit3n -> Text = FloatToStr(T02);
    // вывод графиков
    CBOpt -> Checked = true;
    Optima(T01, T02, J0, J0et, true, true, Edit15n, Edit16n);
}
// описание функции осуществляющей обратный расчёт и возврат величины най-енной температуры
float __fastcall TForm1::Ob_Ras(float J0et, // эталонное значение КАЦ
    float J0, // действительное значение КАЦ
    float T0et, // начальное значение температуры полимеризации (берёмначение эталонной Т0)
    float vt, // текущий объём реакционной смеси (vt+vd1[+vd2])
    float Kf, // коэффициент теплопередачи (расчитывается перед вызо-ом текущей функции)
    float Kmd, // константа скорости роста дивинильных цепей
    float so_d, // значение сухого остатка на текущей стадии (расчиты-ается до вызова данной функции)
    float m0d, // начальная концентрация мономера
    float m, // масса мономера добавляемая на текущей стадии
    float M, // масса реакционной смеси (mr+ms+md1[+md2])
    float gh, // расход хладагента
    float th, // температура хладагента
    bool graf) // флаг необходимости построения графиков
{
    float so, Xmd, T;
    float L = 0, // уровень загрузки реактора
        R = 8.32; // универсальная газовая постоянная
    // коэффициенты Рунге-Кутта
    float x1, y1, x2, y2, x3, y3, x4, y4;
    // время
    float t = 0;
    int I, j; // счётчик для элементов эталоного массива
    // начальные условия для расчёта системы ДУ (2 стадия)
    Xr = new float[NN];
    // описываем графики
    if (graf) {
        Tet = new TLineSeries(this);
        Tet -> ParentChart = Chart2;
        Tet -> Title = "Tet" + Table1NAME -> Value;
        Tet -> LinePen -> Width = 2;
        Tobr = new TLineSeries(this);
        Tobr -> ParentChart = Chart2;
        Tobr -> Title = "Tobr" + Table1NAME -> Value;
        Tobr -> LinePen -> Width = 2;
    }
    // обнуляем начальный элемент массива
    Xr[0] = 0;
    Xmd = 0;
    T = T0et;
    so = so_d;
    I = 1;
    L = (vt - 3000) / 28500;
    j = 0;
    t = 0;
    if (graf)
        Tet -> AddXY(t, T, "", clBlack);
    t = t + dt;
    // Решаем систему описывающую полимеризацию второй стадии методом Рунге-утта 4 порядка
    while (1) {
        x1 = Kmd * pow(J0et, 0.25 + ml2 * so_d) * (1 - Xmd);
        y1 = (Kted * x1 * m0d * (vt) - (Kf * F * L * (T - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + (vt) * pt * ct);
        x2 = Kmd * pow(J0et, 0.25 + ml2 * so_d) * (1 - (Xmd + x1 * dt / 2));
        y2 = (Kted * x2 * m0d * (vt) - (Kf * F * L * ((T + y1 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + (vt) * pt * ct);
        x3 = Kmd * pow(J0et, 0.25 + ml2 * so_d) * (1 - (Xmd + x2 * dt / 2));
        y3 = (Kted * x3 * m0d * (vt) - (Kf * F * L * ((T + y2 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + (vt) * pt * ct);
        x4 = Kmd * pow(J0et, 0.25 + ml2 * so_d) * (1 - (Xmd + x3 * dt / 2));
        y4 = (Kted * x4 * m0d * (vt) - (Kf * F * L * ((T + y3 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + (vt) * pt * ct);
        // расчёт переменных дифференциального уравнения
        Xmd = Xmd + dt * (x1 + 2 * (x2 + x3) + x4) / 6;
        T = T + dt * (y1 + 2 * (y2 + y3) + y4) / 6;
        // расчёт вспомогательных переменных ММ
        so_d = so + m * Xmd / M;
        Kmd = k2 * (1 - b2 * so_d) * exp(-(E1 + E2 * so_d) / (R * T));
        Kf = A - B * exp(C * so_d);
        // построение графиков (при необходимости)
        if (graf)
            Tet -> AddXY(t, T, "", clBlack);
        // приращение времени
        t = t + dt;
        // заполнение массива конверсий
        Xr[I] = Xmd;
        // приращение индекса массива
        I++;
        // условие прекращения расчета для модели с текущей КАЦ
        if (Xmd > 0.99)
            break;
    }
    // осуществляем обратный расчёт
    t = t - dt;
    I = I - 1;
    // построение графиков
    if (graf)
        Tobr -> AddXY(t, T, "", clRed);
    t = t - dt;
    // Решаем систему описывающую полимеризацию второй стадии методом Рунге-утта 4 порядка
    while (1) {
        so_d = so + m * Xr[I] / M;
        Kmd = (k2 * (1 - b2 * so_d)) * exp(-(E1 + E2 * so_d) / (R * T));
        Kf = A - B * exp(C * so_d);
        y1 = (Kted * Kmd * pow(J0, 0.25 + ml2 * so_d) * (1 - Xr[I]) * m0d * (vt) -
            (Kf * F * L * (T - th) * gh * ch * ph / (Kf * F * L + gh * ch * ph))) / (Map * cap + (vt) * pt * ct);
        y2 = (Kted * Kmd * pow(J0, 0.25 + ml2 * so_d) * (1 - Xr[I]) * m0d * (vt) - (Kf * F * L * ((Ty1 *
            dt / 2) - th) * gh * ch * ph / (Kf * F * L + gh * ch * ph))) / (Map * cap + (vt) * pt * ct);
        y3 = (Kted * Kmd * pow(J0, 0.25 + ml2 * so_d) * (1 - Xr[I]) * m0d * (vt) - (Kf * F * L * ((Ty2 *
            dt / 2) - th) * gh * ch * ph / (Kf * F * L + gh * ch * ph))) / (Map * cap + (vt) * pt * ct);
        y4 = (Kted * Kmd * pow(J0, 0.25 + ml2 * so_d) * (1 - Xr[I]) * m0d * (vt) - (Kf * F * L * ((Ty3 *
            dt / 2) - th) * gh * ch * ph / (Kf * F * L + gh * ch * ph))) / (Map * cap + (vt) * pt * ct);
        T = T - dt * (y1 + 2 * (y2 + y3) + y4) / 6;
        // построение графиков
        if (graf)
            Tobr -> AddXY(t, T, "", clRed);
        // приращение времени
        t = t - dt;
        // инкремент индекса массива для хранения конверсии инициатора
        I--;
        // условие прекращения расчета для модели с текущей КАЦ
        if (Xmd <= 0 || I < 1)
            break;
    }
    // очистка памяти
    delete Xr;
    // возвращение оптимальной начальной температуры
    return T;
}
// ------------------------------------------------------------------------
-- -- -
// перегружаем функцию Ob_ras()
void __fastcall TForm1::Ob_Ras() {
        float Tobr1 = 0, Tobr2 = 0, J0;
        bool graf;
        float R = 8.32; // универсальная газовая постоянная
        if (CBObr -> Checked)
            graf = true;
        else
            graf = false;
        // получение исходных параметров
        J0 = StrToFloat(Edit10n -> Text);
        J0et = StrToFloat(Edit9n -> Text);
        dt = StrToFloat(Edit7n -> Text);
        ml2 = StrToFloat(Edit16 -> Text); // коэффициент идентификации
        E1 = StrToFloat(Edit24 -> Text); // коэффициент идентификации E1
        E2 = StrToFloat(Edit35 -> Text); // коэффициент идентификации E2
        // Инициализация параметров
        Ini();
        // Расчет стадии 2_1
        if (CB2_1n -> Checked) {
            so = ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1);
            Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * Tex2_1[0]));
            Kf = A - B * exp(C * so);
            Tobr1 = Ob_Ras(J0et, J0, Tex2_1[0], vt + vd1, Kf, Kmd, so, m0d1, md1,
                mr + ms + md1, gh2_1, th2_1, graf);
            Edit19n -> Text = Tobr1;
        }
        if (CB2_2n -> Checked) {
            so = (ms * 1 / (mr + ms) * (ms + mr) / (ms + mr + md1) +
                md1 * 1 / (ms + mr + md1)) * (ms + mr + md1) / (ms + mr + md1 + md2);
            Kmd = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * Tex2_2[0]));
            Kf = A - B * exp(C * so);
            J0 = J0 * (vt + vd1) / (vt + vd1 + vd2); // пересчёт КАЦ для второй порцииивинила
            J0et = J0et * (vt + vd1) / (vt + vd1 + vd2); // пересчёт КАЦ для второй порцииивинила
            Tobr2 = Ob_Ras(J0et, J0, Tex2_2[0], vt + vd1 + vd2, Kf, Kmd, so, m0d2, md2,
                mr + ms + md1 + md2, gh2_2, th2_2, graf);
            Edit20n -> Text = Tobr2;
        }
    }
    // ------------------------------------------------------------------------
    -- -
    // Описание функции расчёта по модели при расчёте Т0
    int __fastcall TForm1::ModRasT0(float J0, // КАЦ
        float Km, // константа скорости роста цепи
        float m0, // начальная концентрация мономера
        float vt, // текущий объем реакционной смеси
        float so_b, // количество сухого остатка с предыдущей стадии
        float m, // масса мономера на текущей стадии (md1 или md2)
        float M, // масса реакционной смеси (все компоненты)
        float * Xm, // указатель на массив конверсии мономера
        float * T, // указатель на массив температуры реакции
        float gh, // расход хладагента
        float T_k, // температура начала следующей стадии
        float * Tmax, // максимальная температура
        float th, // температура хладагента
        bool f, // булевская переменная
        int lim) // количество итераций после выполнения условия за-ершения цикла
{
    int i, n, c = 0, g = 0;
    float L, // уровень загрузки реактора
    R = 8.32; // универсальная газовая постоянная
    float x1, y1, x2, y2, x3, y3, x4, y4; // коэффициенты Рунге-Кутта
    L = (vt - 3000) / 28500; // уровень заполнения реактора
    i = 0;
    * Tmax = T[0];
    while (1) {
        x1 = Km * pow(J0, 0.25 + ml2 * so) * (1 - Xm[i]);
        if (CBLim -> Checked) { // ограничитель
            if (x1 > 2 * (1 - Xm[i]) / dt)
                x1 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt; // 10% от общей "массы"
        }
        y1 = (Kted * x1 * m0 * vt - (Kf * F * L * (T[i] - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        x2 = Km * pow(J0, 0.25 + ml2 * so) * (1 - (Xm[i] + x1 * dt / 2));
        if (CBLim -> Checked) { // ограничитель
            if (x2 > 2 * (1 - Xm[i]) / dt)
                x2 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y2 = (Kted * x2 * m0 * vt - (Kf * F * L * ((T[i] + y1 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        x3 = Km * pow(J0, 0.25 + ml2 * so) * (1 - (Xm[i] + x2 * dt / 2));
        if (CBLim -> Checked) { // ограничитель
            if (x3 > 2 * (1 - Xm[i]) / dt)
                x3 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y3 = (Kted * x3 * m0 * vt - (Kf * F * L * ((T[i] + y2 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        x4 = Km * pow(J0, 0.25 + ml2 * so) * (1 - (Xm[i] + x3 * dt / 2));
        if (CBLim -> Checked) { // ограничитель
            if (x4 > 2 * (1 - Xm[i]) / dt)
                x4 = 2 * (1 - Xm[i]) / dt - 2 * 0.1 * (1 - Xm[i]) / dt;
        }
        y4 = (Kted * x4 * m0 * vt - (Kf * F * L * ((T[i] + y3 * dt / 2) - th) * gh * ch * ph /
            (Kf * F * L + gh * ch * ph))) / (Map * cap + vt * pt * ct);
        Xm[i + 1] = Xm[i] + dt * (x1 + 2 * (x2 + x3) + x4) / 6;
        T[i + 1] = T[i] + dt * (y1 + 2 * (y2 + y3) + y4) / 6;
        if (T[i + 1] > * Tmax)
            *
            Tmax = T[i + 1];
        if (Xm[i + 1] > 1.001) {
            n = i + 1;
            break;
        }
        so = so_b + m * Xm[i + 1] / M;
        Km = k2 * (1 - b2 * so) * exp(-(E1 + E2 * so) / (R * T[i + 1])); // расчёт константыкорости роста цепи
        Kf = A - B * exp(C * so); // коэффициент теплопередачи
        // условие прекращения расчета для модели с текущей КАЦ
        if (Xm[i + 1] >= 0.98) {
            if (!f) {
                if (T[i + 1] <= T_k) {
                    n = i + 1;
                    break;
                }
            } else {
                c++;
                if (c > lim) {
                    n = i + 1;
                    break;
                }
            }
        }
        if (i > NN - 3) {
            n = i + 1;
            break;
        }
        i++;
    }
    return n;
}