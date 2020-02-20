# -*- coding: utf-8 -*-

# 
# Drawing animated GIFs with matplotlib
# https://eli.thegreenplace.net/2016/drawing-animated-gifs-with-matplotlib/
# https://github.com/maxberggren/blog-notebooks/blob/master/SweEbola.md
# 传染病 SEIR 模型
# pip install matplotlib
# pip install PIL
# pip install aabbtree
# pip install pp 
# ( in python3 ,should download zip file form https://www.parallelpython.com/downloads.php)
# ( and use python setup.py install)
# 2020/02/13 add pp to parallel cal.

"""
Created on Tue Feb 4 22:22:22 2020
@author: timyy@126.com
timyy.github.com/COVID-19
©2020 Timyy, All rights reserved.
"""

import numpy as np
import sys
import math
import time
import random
import pp

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.image as mpimg
import matplotlib.cm as cm

from PIL import Image
from aabbtree import AABB
from aabbtree import AABBTree

rcParams['font.family'] = 'serif'
rcParams['font.family'] = 'SimHei'
rcParams['font.size'] = 10
rcParams['figure.figsize'] = 12, 9

gMaxInfectOther = 10  #每个人最大传播数
gInfectionRate = 0.1
gInfectDis = 80  #感染距离
gExposedDay = 10  # 潜伏期
gSelfRecoveryDay = 50
gSelfRecoverRate = 0.
gDeathDay = 50
gDeathRate = 0.

N = 2000  #总数
gN0 = 10  # 初始感染人数
MaxSize = 10000  #地图尺寸
MoveStep = 160  # 每个人移动距离

personList = []
MaxTime = 300  # 执行时间
keyFrames = []
frames = 30.0  #录像频数
tree = AABBTree()
gif_name = r'D:\Timyy\project\python\covid-19\ShowVirus.1.3.3.gif'
gtmpFolder = r'D:\Timyy\project\python\covid-19\tmp\\'
gMethod = 1
gUsePP = False # 使用并行算法。python 是单线程的。这个方法用了多线程，使用多个CPU计算

class Hosptal:
    def __init__(self):
        self.HosptalBedMax = 4  # 医院床位数
        self.HosptalBedRemain = 4  # 床位剩余
        self.HosptalRecoveryTime = 10  #康复时间
        self.HosptalResponse = 10  #响应时间

    def heal(self):
        if self.HosptalBedRemain > 0:
            self.HosptalBedRemain -= 1
            return True
        else:
            return False

    def healed(self):
        if self.HosptalBedRemain < self.HosptalBedMax:
            self.HosptalBedRemain += 1


gHosptal = Hosptal()


class Person:
    def __init__(self, id=0):
        self.id = id
        self.status = 0
        self.startInfect = 0
        self.startExposed = 0
        self.startHeal = 0
        self.InfOther = 0
        # 初始化直接放数，不放0了。
        self.posionX = random.randint(0, MaxSize - 1)
        self.posionY = random.randint(0, MaxSize - 1)

    def move(self, moveStep):
        # 随机移动一下，出边了要拉回来。
        # 取模的方法有点问题，出右边了，会到左边去。相对而言会造成数据向0点集中。应该给反弹回去。
        self.posionX = math.fabs(self.posionX +
                                 random.randint(-moveStep, moveStep))
        if (self.posionX > MaxSize):
            self.posionX -= moveStep
        self.posionY = math.fabs(self.posionY +
                                 random.randint(-moveStep, moveStep))
        if (self.posionY > MaxSize):
            self.posionY -= moveStep

    def dis(self, a, b):
        return math.fabs(a.posionX - b.posionX) + math.fabs(a.posionY -
                                                            b.posionY)

    def infectOther(self, day, persionlist):
        global tree
        if ((self.status == 1 or self.status == 2)
                and self.InfOther < gMaxInfectOther):
            if (gMethod == 1):
                aabb1 = AABB([(self.posionX - gInfectDis // 2,
                               self.posionX + gInfectDis // 2),
                              (self.posionY - gInfectDis // 2,
                               self.posionY + gInfectDis // 2)])
                result = tree.overlap_values(aabb1)
                if (len(result) > 0):
                    for apersonid in result:
                        aperson = persionlist[apersonid]
                        if (aperson.status == 0):
                            # 易感人群才会被感染
                            if (aperson.infect(day)):                                
                                self.InfOther += 1
            else:
                #这是个N*2的做法，换为AABBTree方法，nlog(n),二叉树
                for aperson in persionlist:
                    if (aperson.status == 0):
                        # 易感人群才会被感染
                        if (self.dis(self, aperson) < gInfectDis):
                            #todo 这里计算距离可以做优化，排个序用快查，用box。整个系统的速度瓶颈
                            if (aperson.infect(day)):
                                self.InfOther += 1

    def infect(self, day):
        # 感染， 有机率的感染
        if (self.status == 0):
            if (random.random() <= gInfectionRate):
                # 在感染机率内，则感染
                self.status = 2
                self.startExposed = day
                return True
        return False

    def update(self, day):
        # 潜伏期的
        if (self.status == 2):
            if (day - self.startExposed >= gExposedDay):
                self.status = 1
                self.startInfect = day
        # 已经感染发病的
        if (self.status == 1):
            #被医院冶好
            if (day - self.startInfect >= gHosptal.HosptalResponse):
                if (gHosptal.heal()):
                    self.status = 3
                    self.startHeal = day
            #auto recovery
            if (self.status == 1):
                if (day - self.startInfect >= gSelfRecoveryDay):
                    if (random.random() <= gSelfRecoverRate):
                        self.status = 4
                        self.startHeal = 0
            #death
            if (self.status == 1):
                if (day - self.startInfect >= gDeathDay):
                    if (random.random() <= gDeathRate):
                        self.status = 5

        # 已经治好的,release hosptal after recoverytime.
        if (self.status == 3):
            if ((self.startHeal > 0) and
                (day - self.startHeal >= gHosptal.HosptalRecoveryTime)):
                gHosptal.healed()
                self.startHeal = 0
def infectOtherPP(day ,personlist, start, end):
    # 并行算法
    for i in range(start, end):
        apersion = personlist[i]
        apersion.infectOther(day, personlist)

def create_gif(image_list, gif_name, gif_duration=0.1):
    frames = []
    for image_name in image_list:
        frames.append(Image.open(image_name))
    # Save them as frames into a gif
    frames[0].save(
        gif_name, save_all=True, append_images=frames, duration=gif_duration)


def virus(bSaveGif=False):
    def scatterLegend(personlist, x, y):
        type0 = []
        type1 = []
        type2 = []
        type3 = []
        type4 = []
        type5 = []
        for aperson in personlist:
            if aperson.status == 0:
                type0.append(np.array((aperson.posionX, aperson.posionY)))
            elif aperson.status == 1:
                type1.append(np.array((aperson.posionX, aperson.posionY)))
            elif aperson.status == 2:
                type2.append(np.array((aperson.posionX, aperson.posionY)))
            elif aperson.status == 3:
                type3.append(np.array((aperson.posionX, aperson.posionY)))
            elif aperson.status == 4:
                type4.append(np.array((aperson.posionX, aperson.posionY)))
            elif aperson.status == 5:
                type5.append(np.array((aperson.posionX, aperson.posionY)))

        type0 = np.array(type0)
        type1 = np.array(type1)
        type2 = np.array(type2)
        type3 = np.array(type3)
        type4 = np.array(type4)
        type5 = np.array(type5)

        # 会有空的情况，就要做一下
        handles = []
        labels = []
        if (len(type0) != 0):
            g0 = plt.scatter(
                type0[:, x], type0[:, y], color='darkblue', marker='.')
            handles.append(g0)
            labels.append('易感者')
        if (len(type1) != 0):
            g1 = plt.scatter(type1[:, x], type1[:, y], c='red', marker='.')
            handles.append(g1)
            labels.append('感染者')
        if (len(type2) != 0):
            g2 = plt.scatter(type2[:, x], type2[:, y], c='orange', marker='.')
            handles.append(g2)
            labels.append('潜伏者')
        if (len(type3) != 0):
            g3 = plt.scatter(type3[:, x], type3[:, y], c='green', marker='.')
            handles.append(g3)
            labels.append('康复者')
        if (len(type4) != 0):
            g4 = plt.scatter(
                type4[:, x], type4[:, y], c='lightgreen', marker='.')
            handles.append(g4)
            labels.append('自愈者')
        if (len(type5) != 0):
            g5 = plt.scatter(type5[:, x], type5[:, y], c='gray', marker='x')
            handles.append(g5)
            labels.append('死亡者')
        #plt.legend(handles=[g0, g1, g2,  g3], labels=['Susceptible', 'Infection', 'Exposed', 'Recovery'])
        plt.legend(handles=handles, labels=labels, loc=1)

    def drawFig():
        plt.cla()
        plt.plot(nSuscept, color='darkblue', label='Susceptible', marker='.')
        plt.plot(nInfect, color='red', label='Infection', marker='.')
        plt.plot(nExposed, color='orange', label='Exposed', marker='.')
        plt.plot(nRecovery, color='green', label='Recovery', marker='.')
        plt.plot(nSelfRecovery, color='lightgreen', label='SelfRecovery', marker='.')
        plt.plot(nDeath, color='grey', label='Death', marker='.')
        #plt.title('SEIR Model')
        plt.legend(loc=1)
        plt.xlabel('Day')
        plt.ylabel('Number')
    
    def DoStep(day, personlist):
        nDaySuscept = 0
        nDayExposed = 0
        nDayInfect = 0
        nDayRecovery = 0
        nDaySelfRecovery = 0
        nDayDeath = 0

        global tree
        for aperson in personlist:
            if aperson.status == 5: # has death
                pass
            else:
                aperson.move(MoveStep)
        
        if gUsePP:
            parts = 4
            start = 0
            end = len(personlist)
            step = int((end-start) / parts + 1) 

            jobs = []

            for index in range(parts):
                starti = start + index*step
                endi= min(start+(index+1)*step, end)
                # Submit a job which will calculate partial sum
                # part_sum - the function
                # (starti, endi) - tuple with arguments for part_sum
                # () - tuple with functions on which function part_sum depends
                # () - tuple with module names which must be
                #      imported before part_sum execution
                jobs.append(job_server.submit(infectOtherPP, (day, personlist,starti, endi)))
            #job_server.print_stats()
            
        else:
            for aperson in personlist:
                aperson.infectOther(day, personlist)
        for aperson in personlist:
            aperson.update(day)
            if aperson.status == 0:
                nDaySuscept += 1
            if aperson.status == 1:
                nDayInfect += 1
            if aperson.status == 2:
                nDayExposed += 1
            if aperson.status == 3:
                nDayRecovery += 1
            if aperson.status == 4:
                nDaySelfRecovery += 1
            if aperson.status == 5:
                nDayDeath += 1

        nSuscept.append(nDaySuscept)
        nInfect.append(nDayInfect)
        nExposed.append(nDayExposed)
        nRecovery.append(nDayRecovery)
        nSelfRecovery.append(nDaySelfRecovery)
        nDeath.append(nDayDeath)
#
    nSuscept = []
    nExposed = []
    nInfect = []
    nRecovery = []
    nSelfRecovery = []
    nDeath = []

    # 初始化每个人的位置，随机放置
    for i in range(0, N - 1):
        aperson = Person(i)
        personList.append(aperson)
        aabb1 = AABB([(aperson.posionX - gInfectDis // 2,
                       aperson.posionX + gInfectDis // 2),
                      (aperson.posionY - gInfectDis // 2,
                       aperson.posionY + gInfectDis // 2)])
        tree.add(aabb1, aperson.id)

    # 初始感染
    for i in range(0, gN0 - 1):
        personList[i].status = 1
    for i in range(0, MaxTime):
        DoStep(i, personList)
        if (i % int(MaxTime / frames) == 0):
            plt.subplot(2, 1, 1)  #两行一列的子图，目前在第一个位置画。
            plt.cla()
            label = 'Day {0}'.format(i)
            plt.title(label)
            scatterLegend(personList, 0, 1)
            plt.subplot(2, 1, 2)  #两行一列的子图，目前在第二个位置画。
            drawFig()
            if (bSaveGif):
                filename =gtmpFolder +  r"outbreak" + str(i) + ".png"
                plt.savefig(filename, dpi=150, bbox_inches='tight')
                keyFrames.append(filename)
            else:
                plt.pause(0.01)

if __name__ == '__main__':

    # Create jobserver

    if gUsePP:
        job_server = pp.Server()
        job_server.set_ncpus() 
    #    anim = FuncAnimation(virus.imgplot, update, frames=np.arange(0, 10), interval=200)
    t0 = time.time()
    if len(sys.argv) > 1 and sys.argv[1] == 'save':
        virus(True)
        image_list = keyFrames
        create_gif(image_list, gif_name, gif_duration=0.3)
    else:
        virus()
    print((time.time()-t0))

        # anim.save('line.gif', dpi=80, writer='imagemagick')
        # anim.save('line.gif', dpi=80, writer='ffmpeg')
#        anim.save(gif_name, dpi=300, writer='pillow')
#    else:
# plt.show() will just loop the animation forever.
#        plt.show()