import numpy as np
import json, codecs
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import interp1d

data = {"gamma_w": 1.06, "h":1312, "d_in":95, "q_liq":588,
        "eps":9.0E-5, "p_wh": 24.8, "MD": [0, 656, 1312],
        "TVD": [0, 608.5, 1027], "T0": 27.5}

# with open('initial_data.json', 'w', encoding='utf-8') as f:
#     json.dump(data, f, ensure_ascii=False, indent=4)


#Входные данные
gamma_w = data["gamma_w"]         #удельная плотность воды
h = data["h"]                     # глубина скважины (м)
d_in = data["d_in"] * 10**(-3)    #диаметр (м)
q_liq = data["q_liq"]             #дебит (м3/сут)
eps = data["eps"]                 #шероховатость (м)
p_wh = data["p_wh"]               #буферное давление (атм)
MD = data["MD"]                   #(м)
TVD = data["TVD"]                 #(м)
T0 = data["T0"]                   #устьевая температура (градусы Цельсия)
g = 9.8                           #ускорение свободного падения (м/с2)


# Плотность, кг/м3
def density(T_K, gamma=gamma_w):
    rho_0 = 1000 * gamma  # плотность воды в ст.у. (кг/м3)
    return rho_0 / (1 + (T_K - 273) / 10000 * (0.269 * (T_K - 273) ** 0.637 - 0.8))


# Вязкость мПа*с
def viscosity(P_MPa, T_K):
    A = 109.574
    B = 1.1217
    mu = A * (1.8 * T_K - 460) ** (-B) * (0.9994 + 0.0058 * P_MPa + 0.6534 * 10**(-4) * P_MPa ** 2)
    return mu


# Скорость потока (м/с)
def velocity(q_liq, d=d_in):
    q = q_liq / 86400  # м3/с
    s = np.pi * d ** 2 / 4
    return q / s


# Расчет числа Рейнольдса
def Reynolds(q_liq, mu, rho, d=d_in):  # кг/(м*с2*Па)
    # q_liq (м3/сут)
    # mu  (мПа*с)
    # rho (кг/м3)
    # d  (м)
    v = velocity(q_liq, d)
    return rho * v * d / mu * 1000


# Коэффициент трения
def friction(q_liq, mu, rho, d=d_in, roughness=eps):
    Re = Reynolds(q_liq, mu, rho, d)
    if Re < 3000:
        return 64 / Re
    else:
        return (1.14 - 2 * np.log10(roughness / d + 21.25 / (Re ** 0.9))) ** (-2)


#Вычисление зенитного угла
#Косинус угла определяется как производная от TVD по МD
#Определяем зависимость TVD(MD)
H_vertical = interp1d(MD, TVD, kind='linear', fill_value="extrapolate")

def cos_alpha(h, H_vert = H_vertical,dh = 1): # h - координата вдоль скважины (MD)
    return (H_vert(h+dh)-H_vert(h))/dh


def grad_P(q_liq, P_MPa, T_C, cos_alpha, d=d_in, roughness=eps):
    rho = density(T_C + 273)
    mu = viscosity(P_MPa, T_C + 273)
    f = friction(q_liq, mu, rho, d, roughness)
    v = velocity(q_liq, d)
    return (rho * g * cos_alpha - f * rho * v**2 / (2*d)) / 10**6  #Давление будем считать в МПа


grad_T = 50/max(TVD)


def grad_PT(PT, x):
    dPdx = grad_P(q_liq, PT[0], PT[1], cos_alpha(x))
    dTdx = grad_T * cos_alpha(x)
    return [dPdx, dTdx]

# Граничные условия
p_wh = data["p_wh"] * 0.101325 # атмосферы в МПа
PTwh = [p_wh, T0]

x = np.linspace(0, h, 200)
solution = odeint(grad_PT, PTwh, x)

P = solution[:,0]/0.101325  #МПа->атм
T = solution[:,1]  #Градусы Цельсия

plt.plot(P,x)
plt.xlabel("P (атм)")
plt.ylabel("x, длина скважины (м)")
ax = plt.gca()
ax.invert_yaxis()
plt.title("Распределение давления")
plt.show()


# Зависимость забойного давления от дебита (от 10 до 500 м3/сут)
Q = np.linspace(10, 500, 100)
P_wf = []

for q in Q:
    def grad_PT(PT, x):
        dPdx = grad_P(q, PT[0], PT[1], cos_alpha(x))
        dTdx = grad_T * cos_alpha(x)
        return [dPdx, dTdx]
    solution = odeint(grad_PT, PTwh, x)
    p_ = solution[:,0]
    P_wf.append(p_[len(x)-1]/0.101325) #МПа->атм

plt.title("Зависимость забойного давления от дебита")
plt.plot(Q, P_wf)
plt.xlabel("Q, м3/сут")
plt.ylabel("P_wf, атм")
plt.show()

# Зависимость забойного давления от температуры T0(от 10 до 50 Цельсия)
t0 = np.linspace(10, 50, 100)
P_wf_T = []

for twh in t0:
    def grad_PT(PT, x):
        dPdx = grad_P(q_liq, PT[0], PT[1], cos_alpha(x))
        dTdx = grad_T * cos_alpha(x)
        return [dPdx, dTdx]
    Ptwh = [p_wh, twh]
    solution = odeint(grad_PT, Ptwh, x)
    p_ = solution[:, 0]
    P_wf_T.append(p_[len(x) - 1] / 0.101325)  # МПа->атм

plt.title("Зависимость забойного давления от температуры T0")
plt.plot(t0, P_wf_T)
plt.xlabel("t0, градусы Цельсия")
plt.ylabel("P_wf (атм)")
plt.show()


#Запись результата в файл
result = {"ex1": {"h": list(x[:]), "p_wf": list(P[:])},
          "ex2": {"q_liq": list(Q[:]), "p_wf": list(P_wf[:])},
          "ex3": {"t": list(t0[:]),  "p_wf": list(P_wf_T[:])}
          }
with open('result.json', 'wb') as f:
    json.dump(result, codecs.getwriter('utf-8')(f), ensure_ascii=False)
