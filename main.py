import json
import matplotlib.pyplot as plt
from math import log as ln

with open('Data.txt',encoding='utf-8') as json_file:
    data = json.load(json_file)
    print('Source data:')
    print('Gamma_w: ' + str(data['gamma_w']))
    gamma_w = data['gamma_w']
    print('Height: ' + str(data['h']))
    h = data['h']
    print('Inner Diameter: ' + str(data['d_in']))
    d_in = data['d_in']
    print('Qliq: ' + str(data['q_liq']))
    q_liq = data['q_liq']
    print('Roughness: ' + str(data['eps']))
    eps=data['eps']
    print('Buffer pressure: ' + str(data['p_wh']))
    p_wh=data['p_wh']
    print('MD: ' + str(data['MD']))
    MD = data['MD']
    print('TVD: ' + str(data['TVD']))
    TVD = data['TVD']
    print('T0: ' + str(data['T0']))
    T0 = data['T0']
    print('')

def dens(gamma_W,T):
    #Плотность воды
    #gamma_W - плотность в ст.у. г/см3
    #T - температура в градусах Цельсия
    T = T+273.15
    ds = gamma_W/(1+(T-273.15)/10000*(0.269*((T-273.15)**0.637)-0.8))*1000
    return ds

def visc(P,T):
    # Вязкость воды
    # P - давление в МПа
    # T - температура в градусах Цельсия
    T = T + 273.15
    A = 109.574
    B =  1.1217
    vsc = A*((1.8*T-460)**-B)*(0.9994+0.0058*P+(0.63534*10**-4)*P**2)/1000
    return vsc

def Rein(q_liq,d_m,visc0,dens0):
    #Число Рейнольдса
    # q_liq - дебит жидкости в м3/сут
    # d_m - диаметр скважины в мм
    # visc0 - вязкость жидкости на поверхности  в Па*с
    # dens0 - плотность жидкости на поверхности в кг/м3
    d_m=d_m/1000
    r=d_m/2
    V = q_liq/(24*60*60)/(r**2)/3.1415927
    Re = dens0*V*d_m/visc0
    return Re

def friction(Re,eps, d_m):
    #Определение коэффициента трения
    # Re - число Рейнольдса
    # esp - шероховатость абс.м.
    # d_m - диаметр в мм
    d_m = d_m/1000
    if Re <=3000:
        f=64/Re
    else:
        f=1/(1.14-2*ln(eps/d_m+21.25/(Re**0.9)))**2
    return f

def cosA(x,MD, TVD):
    # Определение косинуса угла скважины
    # x - текущая высота
    # MD - инклинометрия по глубине скважины вдоль скважины
    # TVD - инклинометрия по глубине скважины по вертикали
    if x<MD[1]:
        cA=(TVD[1]-TVD[0])/(MD[1]-MD[0])
    elif x<=MD[2]:
        cA=(TVD[2]-TVD[1])/(MD[2]-MD[1])
    return cA

def pr_grad(gamma_w,P,T,eps,d_m,x,q_liq,MD,TVD):
    # Определение градиента
    # gamma_w - плотность в нормальных условиях
    # MD - инклинометрия по глубине скважины вдоль скважины
    # TVD - инклинометрия по глубине скважины по вертикали
    ds = dens(gamma_w,T)
    vs = visc(P,T)
    Re = Rein(q_liq, d_in, vs, ds)
    f = friction(Re, eps, d_in)
    csA = cosA(x, MD, TVD)
    d_m=d_m/1000
    r=d_m/2
    v = q_liq/(24*60*60)/(r**2)/3.1415927
    Pgr = (ds*9.81*csA-(f*ds*v**2)/(d_m*2))/1000000
    return Pgr

# 0m data
P0 = p_wh*101325/1000000
visc0 = visc(P0,T0)
dens0 = dens(gamma_w,T0)
cA = cosA(0,MD,TVD)
p_gr = pr_grad(gamma_w, P0, T0, eps, d_in, 0, q_liq, MD, TVD)
t_gr = 50/(TVD[2]-TVD[0])
P = P0+p_gr*1
T = T0+1*cA*t_gr


print('Calculated 0m data:')
print('Viscosity0: ' + str(visc0) + ' Pa*s')
print('Density0: ' + str(dens0) + ' kg/m3')
print('Cosine alpha: ' + str(cA))
print('Pressure gradient: ' + str(p_gr) + ' MPa/m')
print('')

Res = []

if int(MD[1]) % 10 == 0:
    cap = int(MD[1])+10
else:
    cap = int(MD[1])

for x in range(1,int(MD[1])+1):
    p_gr = pr_grad(gamma_w,P,T,eps,d_in,x,q_liq,MD,TVD)
    cA = cosA(x, MD, TVD)
    P = P + p_gr*1
    T = T+1*cA*t_gr
    Res.append([x,P,T])
if x == MD[1]:
    x = x
else:
    p_gr = pr_grad(gamma_w,P,T,eps,d_in,x,q_liq,MD,TVD)
    cA = cosA(x, MD, TVD)
    P = P + (p_gr*(MD[1]-x))
    T = T = T + ((MD[1]-x)*cA*t_gr)
    Res.append([MD[1],P,T])
    x = x + 1

if int(MD[2]) % 10 == 0:
    cap = int(MD[2])+10
else:
    cap = int(MD[2])

for x1 in range(x,int(MD[2])+1):
    p_gr = pr_grad(gamma_w, P, T, eps, d_in, x1, q_liq, MD, TVD)
    cA = cosA(x1, MD, TVD)
    P = P + p_gr * 1
    T = T + 1 * cA * t_gr
    Res.append([x1, P, T])

if x1 == MD[2]:
    x1 = x1
else:
    p_gr = pr_grad(gamma_w,P,T,eps,d_in,x1,q_liq,MD,TVD)
    cA = cosA(x1, MD, TVD)
    P = P + (p_gr*(MD[2]-x1))
    T = T = T + ((MD[2]-x1)*cA*t_gr)
    Res.append([MD[2],P,T])

h = []
P = []
T = []
for i in range(0,len(Res)):
    h.append(Res[i][0])
    P.append(Res[i][1])
    T.append(Res[i][2])

plt.figure(1)
plt.plot(P, h, label="Pressure")
plt.xlabel("P, pressure, MPa")
plt.ylabel("h, well deepness, m")
ax = plt.gca()
ax.invert_yaxis()
plt.legend()
plt.title("Pressure distribution");

plt.figure(2)
plt.plot(T, h, label="Temperature", color='orange')
plt.xlabel("T, temperature, Celsius")
plt.ylabel("h, well deepness, m")
ax = plt.gca()
ax.invert_yaxis()
plt.legend()
plt.title("Temperature distribution")
print(Res)

Res_h_1 = h
Res_P_wf_1 = P
# Fluid quantity variation (2-nd task)

P_q_res=[]
P_T_res=[]
Q_var=[]
T_var=[]
for q_liq in range(10,510,10):
    # 0m data
    P0 = p_wh*101325/1000000
    visc0 = visc(P0,T0)
    dens0 = dens(gamma_w,T0)
    cA = cosA(0,MD,TVD)
    p_gr = pr_grad(gamma_w, P0, T0, eps, d_in, 0, q_liq, MD, TVD)
    t_gr = 50/(TVD[2]-TVD[0])
    P = P0+p_gr*1
    T = T0+1*cA*t_gr

    Res = []

    # Distribution calculation

    if int(MD[1]) % 10 == 0:
        cap = int(MD[1])+10
    else:
        cap = int(MD[1])

    for x in range(1,int(MD[1])+1):
        p_gr = pr_grad(gamma_w,P,T,eps,d_in,x,q_liq,MD,TVD)
        cA = cosA(x, MD, TVD)
        P = P + p_gr*1
        T = T+1*cA*t_gr
        Res.append([x,P,T])
    if x == MD[1]:
        x = x
    else:
        p_gr = pr_grad(gamma_w,P,T,eps,d_in,x,q_liq,MD,TVD)
        cA = cosA(x, MD, TVD)
        P = P + (p_gr*(MD[1]-x))
        T = T = T + ((MD[1]-x)*cA*t_gr)
        Res.append([MD[1],P,T])
        x = x + 1

    if int(MD[2]) % 10 == 0:
        cap = int(MD[2])+10
    else:
        cap = int(MD[2])

    for x1 in range(x,int(MD[2])+1):
        p_gr = pr_grad(gamma_w, P, T, eps, d_in, x1, q_liq, MD, TVD)
        cA = cosA(x1, MD, TVD)
        P = P + p_gr * 1
        T = T + 1 * cA * t_gr
        Res.append([x1, P, T])

    if x1 == MD[2]:
        x1 = x1
    else:
        p_gr = pr_grad(gamma_w,P,T,eps,d_in,x1,q_liq,MD,TVD)
        cA = cosA(x1, MD, TVD)
        P = P + (p_gr*(MD[2]-x1))
        T = T = T + ((MD[2]-x1)*cA*t_gr)
        Res.append([MD[2],P,T])
    Q_var.append(q_liq)
    P_q_res.append(Res[len(Res)-1][1])


print(P_q_res)

plt.figure(3)
plt.plot(Q_var, P_q_res, label="Well pressure",color='blue')
plt.xlabel("Q, fluid quantity, m3/day")
plt.ylabel("P, well pressure, MPa")
plt.legend()
plt.title("Pressure well distribution (Fluid quantity)")

Res_Q_2 = Q_var
Res_P_wf_2 = P_q_res

# Temperature variation (3-rd task)

for T0 in range(10,50):
    # 0m data
    P0 = p_wh*101325/1000000
    visc0 = visc(P0,T0)
    dens0 = dens(gamma_w,T0)
    cA = cosA(0,MD,TVD)
    p_gr = pr_grad(gamma_w, P0, T0, eps, d_in, 0, q_liq, MD, TVD)
    t_gr = 50/(TVD[2]-TVD[0])
    P = P0+p_gr*1
    T = T0+1*cA*t_gr

    Res = []

    # Distribution calculation

    if int(MD[1]) % 10 == 0:
        cap = int(MD[1])+10
    else:
        cap = int(MD[1])

    for x in range(1,int(MD[1])+1):
        p_gr = pr_grad(gamma_w,P,T,eps,d_in,x,q_liq,MD,TVD)
        cA = cosA(x, MD, TVD)
        P = P + p_gr*1
        T = T+1*cA*t_gr
        Res.append([x,P,T])
    if x == MD[1]:
        x = x
    else:
        p_gr = pr_grad(gamma_w,P,T,eps,d_in,x,q_liq,MD,TVD)
        cA = cosA(x, MD, TVD)
        P = P + (p_gr*(MD[1]-x))
        T = T = T + ((MD[1]-x)*cA*t_gr)
        Res.append([MD[1],P,T])
        x = x + 1

    if int(MD[2]) % 10 == 0:
        cap = int(MD[2])+10
    else:
        cap = int(MD[2])

    for x1 in range(x,int(MD[2])+1):
        p_gr = pr_grad(gamma_w, P, T, eps, d_in, x1, q_liq, MD, TVD)
        cA = cosA(x1, MD, TVD)
        P = P + p_gr * 1
        T = T + 1 * cA * t_gr
        Res.append([x1, P, T])

    if x1 == MD[2]:
        x1 = x1
    else:
        p_gr = pr_grad(gamma_w,P,T,eps,d_in,x1,q_liq,MD,TVD)
        cA = cosA(x1, MD, TVD)
        P = P + (p_gr*(MD[2]-x1))
        T = T = T + ((MD[2]-x1)*cA*t_gr)
        Res.append([MD[2],P,T])
    T_var.append(T0)
    P_T_res.append(Res[len(Res)-1][1])

print(P_q_res)

plt.figure(4)
plt.plot(T_var, P_T_res, label="Well pressure",color='red')
plt.xlabel("T, temperature, Celsius")
plt.ylabel("P, well pressure, MPa")
plt.legend()
plt.title("Pressure well distribution (Temperature)")

plt.show()


Res_T_3 = T_var
Res_P_wf_3 = P_T_res

Result = {
        "1 task": {
            "h": str(Res_h_1),
            "P": str(Res_P_wf_1)
        },
        "2 task": {
            "q_liq": str(Res_Q_2),
            "P_wf": str(Res_P_wf_2)
        },
        "3 task": {
            "q_liq": str(Res_T_3),
            "P_wf": str(Res_P_wf_3)
        },
        }
file = open("result.json","w")
json.dump(Result,file,indent=6)
file.close()