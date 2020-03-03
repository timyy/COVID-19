# -*- coding: utf-8 -*-

# Drawing animated GIFs with matplotlib
# https://eli.thegreenplace.net/2016/drawing-animated-gifs-with-matplotlib/
# https://github.com/maxberggren/blog-notebooks/blob/master/SweEbola.md
# https://zhuanlan.zhihu.com/p/104091330?utm_source=qq
# SEIR 模型

"""
Created on Sat Jun 4 22:22:22 2020
@author: timyy
"""

import scipy.integrate as spi
import numpy as np
import matplotlib.pyplot as plt

def SI():
    # N为人群总数
    N = 10000
    # β为传染率系数
    beta = 0.25
    # gamma为恢复率系数
    gamma = 0
    # I_0为感染者的初始人数
    I_0 = 1
    # S_0为易感者的初始人数
    S_0 = N - I_0
    # T为传播时间
    T = 150

    # INI为初始状态下的数组
    INI = (S_0,I_0)


    def funcSI(inivalue,_):
        Y = np.zeros(2)
        X = inivalue
        # 易感个体变化
        Y[0] = - (beta * X[0] * X[1]) / N + gamma * X[1]
        # 感染个体变化
        Y[1] = (beta * X[0] * X[1]) / N - gamma * X[1]
        return Y

    T_range = np.arange(0,T + 1)

    RES = spi.odeint(funcSI,INI,T_range)


    plt.plot(RES[:,0],color = 'darkblue',label = 'Susceptible',marker = '.')
    plt.plot(RES[:,1],color = 'red',label = 'Infection',marker = '.')
    plt.title('SI Model')
    plt.legend()
    plt.xlabel('Day')
    plt.ylabel('Number')
    plt.show()

def SIS():
    # N为人群总数
    N = 10000
    # β为传染率系数
    beta = 0.25
    # gamma为恢复率系数
    gamma = 0.05
    # I_0为感染者的初始人数
    I_0 = 1
    # S_0为易感者的初始人数
    S_0 = N - I_0
    # T为传播时间
    T = 150

    # INI为初始状态下的数组
    INI = (S_0,I_0)


    def funcSIS(inivalue,_):
        Y = np.zeros(2)
        X = inivalue
        # 易感个体变化
        Y[0] = - (beta * X[0]) / N * X[1] + gamma * X[1]
        # 感染个体变化
        Y[1] = (beta * X[0] * X[1]) / N - gamma * X[1]
        return Y

    T_range = np.arange(0,T + 1)

    RES = spi.odeint(funcSIS,INI,T_range)


    plt.plot(RES[:,0],color = 'darkblue',label = 'Susceptible',marker = '.')
    plt.plot(RES[:,1],color = 'red',label = 'Infection',marker = '.')
    plt.title('SIS Model')
    plt.legend()
    plt.xlabel('Day')
    plt.ylabel('Number')
    plt.show()

def SIR():
    # N为人群总数
    N = 10000
    # β为传染率系数
    beta = 0.25
    # gamma为恢复率系数
    gamma = 0.05
    # I_0为感染者的初始人数
    I_0 = 1
    # R_0为治愈者的初始人数
    R_0 = 0
    # S_0为易感者的初始人数
    S_0 = N - I_0 - R_0
    # T为传播时间
    T = 150

    # INI为初始状态下的数组
    INI = (S_0,I_0,R_0)


    def funcSIR(inivalue,_):
        Y = np.zeros(3)
        X = inivalue
        # 易感个体变化
        Y[0] = - (beta * X[0] * X[1]) / N
        # 感染个体变化
        Y[1] = (beta * X[0] * X[1]) / N - gamma * X[1]
        # 治愈个体变化
        Y[2] = gamma * X[1]
        return Y

    T_range = np.arange(0,T + 1)

    RES = spi.odeint(funcSIR,INI,T_range)


    plt.plot(RES[:,0],color = 'darkblue',label = 'Susceptible',marker = '.')
    plt.plot(RES[:,1],color = 'red',label = 'Infection',marker = '.')
    plt.plot(RES[:,2],color = 'green',label = 'Recovery',marker = '.')
    plt.title('SIR Model')
    plt.legend()
    plt.xlabel('Day')
    plt.ylabel('Number')
    plt.show()

def SIRS():
    # N为人群总数
    N = 10000
    # β为传染率系数
    beta = 0.25
    # gamma为恢复率系数
    gamma = 0.05
    # Ts为抗体持续时间
    Ts = 7
    # I_0为感染者的初始人数
    I_0 = 1
    # R_0为治愈者的初始人数
    R_0 = 0
    # S_0为易感者的初始人数
    S_0 = N - I_0 - R_0
    # T为传播时间
    T = 150

    # INI为初始状态下的数组
    INI = (S_0,I_0,R_0)


    def funcSIRS(inivalue,_):
        Y = np.zeros(3)
        X = inivalue
        # 易感个体变化
        Y[0] = - (beta * X[0] * X[1]) / N + X[2] / Ts
        # 感染个体变化
        Y[1] = (beta * X[0] * X[1]) / N - gamma * X[1]
        # 治愈个体变化
        Y[2] = gamma * X[1] - X[2] / Ts
        return Y

    T_range = np.arange(0,T + 1)

    RES = spi.odeint(funcSIRS,INI,T_range)


    plt.plot(RES[:,0],color = 'darkblue',label = 'Susceptible',marker = '.')
    plt.plot(RES[:,1],color = 'red',label = 'Infection',marker = '.')
    plt.plot(RES[:,2],color = 'green',label = 'Recovery',marker = '.')
    plt.title('SIRS Model')
    plt.legend()
    plt.xlabel('Day')
    plt.ylabel('Number')
    plt.show()

def SEIR():
    # N为人群总数
    N = 10000
    # β为传染率系数
    beta = 0.6
    # gamma为恢复率系数
    gamma = 0.1
    # Te为疾病潜伏期
    Te = 14
    # I_0为感染者的初始人数
    I_0 = 1
    # E_0为潜伏者的初始人数
    E_0 = 0
    # R_0为治愈者的初始人数
    R_0 = 0
    # S_0为易感者的初始人数
    S_0 = N - I_0 - E_0 - R_0
    # T为传播时间
    T = 150

    # INI为初始状态下的数组
    INI = (S_0,E_0,I_0,R_0)


    def funcSEIR(inivalue,_):
        Y = np.zeros(4)
        X = inivalue
        # 易感个体变化
        Y[0] = - (beta * X[0] * X[2]) / N
        # 潜伏个体变化
        Y[1] = (beta * X[0] * X[2]) / N - X[1] / Te
        # 感染个体变化
        Y[2] = X[1] / Te - gamma * X[2]
        # 治愈个体变化
        Y[3] = gamma * X[2]
        return Y

    T_range = np.arange(0,T + 1)

    RES = spi.odeint(funcSEIR,INI,T_range)


    plt.plot(RES[:,0],color = 'darkblue',label = 'Susceptible',marker = '.')
    plt.plot(RES[:,1],color = 'orange',label = 'Exposed',marker = '.')
    plt.plot(RES[:,2],color = 'red',label = 'Infection',marker = '.')
    plt.plot(RES[:,3],color = 'green',label = 'Recovery',marker = '.')

    plt.title('SEIR Model')
    plt.legend()
    plt.xlabel('Day')
    plt.ylabel('Number')
    plt.show()
    
def SEIRS():
    # N为人群总数
    N = 10000
    # β为传染率系数
    beta = 0.6
    # gamma为恢复率系数
    gamma = 0.1
    # Ts为抗体持续时间
    Ts = 7
    # Te为疾病潜伏期
    Te = 14
    # I_0为感染者的初始人数
    I_0 = 1
    # E_0为潜伏者的初始人数
    E_0 = 0
    # R_0为治愈者的初始人数
    R_0 = 0
    # S_0为易感者的初始人数
    S_0 = N - I_0 - E_0 - R_0
    # T为传播时间
    T = 150

    # INI为初始状态下的数组
    INI = (S_0,E_0,I_0,R_0)


    def funcSEIRS(inivalue,_):
        Y = np.zeros(4)
        X = inivalue
        # 易感个体变化
        Y[0] = - (beta * X[0] * X[2]) / N + X[3] / Ts
        # 潜伏个体变化
        Y[1] = (beta * X[0] * X[2]) / N - X[1] / Te
        # 感染个体变化
        Y[2] = X[1] / Te - gamma * X[2]
        # 治愈个体变化
        Y[3] = gamma * X[2] - X[3] / Ts
        return Y

    T_range = np.arange(0,T + 1)

    RES = spi.odeint(funcSEIRS,INI,T_range)


    plt.plot(RES[:,0],color = 'darkblue',label = 'Susceptible',marker = '.')
    plt.plot(RES[:,1],color = 'orange',label = 'Exposed',marker = '.')
    plt.plot(RES[:,2],color = 'red',label = 'Infection',marker = '.')
    plt.plot(RES[:,3],color = 'green',label = 'Recovery',marker = '.')

    plt.title('SEIRS Model')
    plt.legend()
    plt.xlabel('Day')
    plt.ylabel('Number')
    plt.show()    

if __name__ == '__main__':

    # gif_name = r'D:\Timyy\project\python\AniDraw\zombie.gif'
    SI()
    SIS()
    SIR()
    SIRS()
    SEIR()
    SEIRS()



#    anim = FuncAnimation(virus.imgplot, update, frames=np.arange(0, 10), interval=200)
#    if len(sys.argv) > 1 and sys.argv[1] == 'save':
        # anim.save('line.gif', dpi=80, writer='imagemagick')
        # anim.save('line.gif', dpi=80, writer='ffmpeg')
#        anim.save(gif_name, dpi=300, writer='pillow')        
#    else:
        # plt.show() will just loop the animation forever.
#        plt.show()