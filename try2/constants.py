# Количество фракций
FACTIONS = 6

# Начальная температура для каждого опыта
T0 = [310.4, 309.6, 312.7, 309.3, 313.0, 311.2]
# Начальная температура для каждого опыта
T1 = [307.4, 310.4, 309.5, 313.3, 314.5, 309.4]
# Начальная концентрация стирола, моль/л.
ms0 = [0.745, 0.7304, 0.7129, 0.7106, 0.6683, 0.7168]
# Концентрация катализатора , моль/л.
Jk = [0.00478, 0.0045, 0.00416, 0.00361, 0.0041, 0.00402]
# Температура хладагента Тh, °С
Th = [290, 289.3, 288.5, 290, 296, 288.5]
# Расход хладагента в реакции стирола Gh, л/мин
Ghs = [1500, 2000, 2000, 1666, 1500, 1500]
# Расход хладагента в реакции дивинила 1 Gh, л/мин
Ghd1 = [1666, 1666, 1366, 1666, 1866, 1166]
# Расход хладагента в реакции дивинила 2 Gh, л/мин
Ghd2 = [1666, 1666, 1366, 1666, 1866, 1166]
# Молекулярная масса стирола, кг/моль
mms = 0.092
MMs = [mms * i for i in range(FACTIONS)]
# Молекулярная масса дивинила, кг/моль
MMd = 0.054

# Динамическая вязкость
ng = [0.125, 0.129, 0.13, 0.147, 0.129, 0.137]

# Масса дивинила на стадии 1
Md1 =[1765, 1687, 1701, 1664, 1558, 1702]
# Масса дивинила на стадии 2
Md2 =[837, 715, 715, 722, 663, 710]
# Плотность дивинила
Pd = [0.638, 0.637, 0.637, 0.635, 0.637, 0.637]
# Объем дивинила на стадии 1
Vd1 =[Md1[i] / Pd[i] for i in range(len(Md1))]
# Объем дивинила на стадии 2
Vd2 = [Md2[i] / Pd[i] for i in range(len(Md2))]

# Масса растворителя
Mr = [11400, 11140, 11370, 11500, 11270, 11290]
# Плотность растворителя
Pr = [0.748, 0.748, 0.739, 0.75, 0.744, 0.739]
# Объем растворителя
Vr = [Mr[i] / Pr[i] for i in range(len(Mr))]

# Масса стирола
Ms = [1140, 1080, 1085, 1081, 998, 1085]
# Плотность стирола
Ps = [0.91, 0.916, 0.92, 0.919, 0.921, 0.92]
# Объем стирола
Vs = [Ms[i] / Ps[i] for i in range(len(Mr))]

# Объем катализатора
Vk = 139.2
# Объём реакционной смеси
Vt = [(Vk + Vs[i] + Vr[i]) for i in range(len(Md1))]
# Расчет степени заполнения
L = [((Vt[i] - 3000) / 28500) for i in range(len(Vt))]

# Концентрация активных центровполимеризации на первой и второй стадии синтеза
# Jac = np.array([(Jk[i] * ((Vt[i]) / (Vt[i] + Vd2[i]))) for i in range(len(Vt))])
#Jacd1 = [0.00409, 0.00385, 0.00359, 0.00311, 0.00355, 0.00345]
Jacd1 = [(Jk[i] * (Vt[i] / (Vt[i] + Vd1[i]))) for i in range(len(Vt))]
#Jacd2 = [0.00383, 0.00364, 0.00338, 0.00294, 0.00336, 0.00327]
Jacd2 = [(Jk[i] * ((Vt[i]) / (Vt[i] + Vd2[i]))) for i in range(len(Vt))]
# Начальная концентрация дивинила на первой и второй стаддии
# md0 = np.array([((pd[i] * Vd[i]) / (MMd * (Vt[i] + Vd[i]))) for i in range(len(pd))])
#md10 = [1.682, 1.666, 1.6402, 1.6107, 1.542, 1.648]
md10 = [((Pd[i] * Vd1[i]) / (MMd * (Vt[i] + Vd1[i]))) for i in range(len(Vd1))]
#md20 = [0.758, 0.6664, 0.6511, 0.659, 0.6234, 0.649]
md20 = [((Pd[i] * Vd2[i]) / (MMd * (Vt[i] + Vd2[i]))) for i in range(len(Vd2))]

# Постоянные констант скорости инициирования и роста цепи [л/(мин*моль)]
Ki0, Ks0, Kd0 = 0.835 * 10**10, 5.76 * 10**11, 1.997 * 10**12
# Энергия активации инициирования и роста полистирольных цепей [Дж/моль]
Ei, Es, Ed = 59962, 71184, 80850
# Универсальная газовая постоянная
#R = 8.31446
R = 8.32
#
b = 3.77

# Параметр зависимости, учитывающей влияние концентрации полимера на порядок реакции по катализатору
ws = -0.4142
#
wd = -0.184
# Масса аппарата (реактора), кг
Map = 16820
# Теплоёмкость материала аппарата, Дж/(кг*К)
Cap = 0.46
# Теплоёмкость реакционной массы, кДж/(кг*К)
Cpm = 1.8436
# Плотность реакционной массы, кг/л
dpm = 0.763
# Теплоёмкость (кДж/(кг*К)) и плотность хладагента (кг/л)
Ch, ph = 4.19, 1

# Коэффициенты идентификации [Дж/(м2*К*мин)]
A, B, C = 14.17, 1.007, 10.3823

# Площадь поверхности теплообмена в реакторе, m^2
Fst = 42.0
# Коэффициент тепловедение стирола и дивинила
Kted, Ktes = 80.9, 74.87

# Экспериментальные значения температуры
Texp1 = [
        [
        310.4,310.9,311.5,312.4,313.6,315.2,316.7,318.5,320.7,322.3,324.4,326.6,328.2,329.8,331.0,332.0,332.7,333.3,333.7,333.7,333.7,333.5,333.4,333.2,332.8,332.6,332.3,332.0,331.7,331.4,
        ],
        [309.6,310.0,310.3,310.7,311.4,312.0,313.0,314.0,315.0,316.2,317.5,319.0,320.4,322.3,323.7,325.4,326.8,327.9,329.0,330.0,330.7,331.2,331.6,331.8,332.0,332.0,331.9,331.9,331.6,331.4,
        ],
        [312.7,312.9,313.4,314.0,315.0,316.0,317.3,319.2,320.8,322.4,324.8,326.7,328.6,330.1,331.4,332.7,333.6,334.3,334.6,334.3,334.2,333.9,333.5,333.3,333.0,332.7,332.4,332.2,331.8,331.5,
        ],
        [309.3,309.4,309.6,310.0,310.3,310.9,311.6,312.2,313.3,314.2,315.7,317.0,318.3,319.8,321.6,323.1,324.5,326.3,327.4,328.2,329.2,329.9,330.4,330.8,331.1,331.1,331.1,331.0,330.9,330.7,
        ],
        [313.0,313.2,313.6,314.2,315.6,316.7,318.2,320.3,322.9,325.4,328.2,330.0,331.0,331.7,333.9,334.5,334.9,335.1,335.2,335.1,334.9,334.7,334.5,334.3,334.1,333.9,333.7,333.4,333.1,332.8,
        ],
        [311.2,311.4,311.8,312.6,313.4,314.4,315.5,316.6,317.7,319.4,321.0,322.4,324.2,325.7,327.4,328.8,330.1,331.0,332.0,332.5,333.1,333.5,333.7,333.8,333.9,333.8,333.7,333.5,333.5,333.1,
        ],
    ]

Texp2 = [
        [307.4,307.5,307.6,307.7,308.2,308.8,309.6,310.5,311.5,312.5,313.7,315.0,316.5,318.2,320.2,322.6,325.6,329.2,334.6,341.6,353.0,363.2,368.8,370.8,371.4,370.8,
        ],
        [310.4,310.6,311.3,312.2,313.3,314.4,315.9,317.4,319.4,321.4,323.9,327.3,332.1,338.6,349.6,361.4,368.9,372.7,373.9,373.9,373.5,372.6,371.8,
        ],
        [309.5,309.6,309.6,309.6,311.6,312.6,313.6,314.9,316.3,317.8,319.5,321.4,323.7,326.8,330.4,336.1,344.1,354.9,365.3,370.0,372.6,373.6,373.9,
        ],
        [313.3,314.3,315.4,316.7,318.5,320.5,323.2,325.8,329.5,335.9,343.8,353.9,361.3,367.4,371.0,373.1,373.5,373.8,373.6,373.1,372.6,
        ],
        [314.5,315.7,317.2,318.7,320.6,323.2,326.1,329.7,334.6,342.1,353.4,364.3,370.2,372.6,373.2,372.9,372.0,371.1,370.2,369.4,368.6,
        ],
        [309.4,309.6,309.8,310.2,310.8,312.1,313.1,314.3,315.7,317.1,318.9,321.2,323.7,326.9,331.0,337.9,346.8,356.4,363.0,368.1,371.2,372.7,373.3,
        ],
    ]