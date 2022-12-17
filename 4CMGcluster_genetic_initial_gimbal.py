import numpy as np
import random as rnd
import math
from pyquaternion import Quaternion
from texttable import Texttable
import matplotlib.pyplot as plt


def skew(vector):
    vector = np.array(vector)

    return np.array([[0, -vector.item(2), vector.item(1)],
                     [vector.item(2), 0, -vector.item(0)],
                     [-vector.item(1), vector.item(0), 0]])


def cross(x1, x2):
    return np.cross(x1, x2)


def cos(x):
    return math.cos(x)


def sin(x):
    return math.sin(x)


def quaternion_to_euler(q):
    x = q[1]
    y = q[2]
    z = q[3]
    w = q[0]
    t0 = +2.0 * (w * x + y * z)
    t1 = +1.0 - 2.0 * (x * x + y * y)
    roll = math.atan2(t0, t1)
    t2 = +2.0 * (w * y - z * x)
    t2 = +1.0 if t2 > +1.0 else t2
    t2 = -1.0 if t2 < -1.0 else t2
    pitch = math.asin(t2)
    t3 = +2.0 * (w * z + x * y)
    t4 = +1.0 - 2.0 * (y * y + z * z)
    yaw = math.atan2(t3, t4)
    RPY = np.array([[roll], [pitch], [yaw]])
    return RPY


def mainCMG_code(gene1,gene2,gene3,gene4):


    w = np.array([[0], [0], [0]])
    J = 1 * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])  # % Mass=24kg , Side of Cube=0.5m
    Td = 0 * np.array([[0], [10], [0]])

    Tc = np.array([[0], [0], [0]])
    w_prev = w;
    dt = 0.01
    u = np.array([[0], [0], [0]])
    ini = np.array([[1], [0], [0]])

    qref = Quaternion(1, 0, 0, 0)
    q = qref

    q_prev = q.unit

    b = math.radians(54.73)
    # d1 = math.radians(-70)
    # d2 = 0
    # d3 = math.radians(75)
    # d4 = 0

    d1 = gene1
    d2 = gene2
    d3 = gene3
    d4 = gene4


    dprev = np.array([[d1], [d2], [d3], [d4]])

    Integral = np.array([[0], [0], [0]])
    axiserrorprev = np.array([[0], [0], [0]])

    gfactor = 1
    Kp = 20 / gfactor
    Ki = 0.00001 / gfactor
    Kwmega = 100 / gfactor

    endtime = 10 / dt
    timeinsec = np.arange(0, endtime * dt - dt, dt)
    stop = 3 / dt
    plot_motion = 0

    h0 = 1
    # Favoid = np.array([[0], [0], [0]])
    manip = np.array([])
    # dd = np.empty((0, 4), float)
    attitude_angle = np.empty((3, 0), float)
    # wplot = np.empty((3, 0), float)
    # hplot = np.empty((3, 0), float)
    # errorplot = np.empty((3, 0), float)
    # d_dotplot = np.empty((4, 0), float)

    for dur in range(1, int(endtime)):
        if dur > stop:
            Td = np.array([[0], [0], [0]])
            qref = Quaternion(0.70711, -0.70711, 0, 0)

        wdot = np.dot(np.linalg.inv(J), (np.dot(np.dot(-skew(w), J), w) + Td + Tc))
        w = wdot * dt + w_prev
        wquaternion = Quaternion(np.insert(w, 0, 0))
        qd = 0.5 * q * wquaternion
        # if ( w[1]!=0 and w[2]!=0 and w[3]!=0 ):
        if np.all(w):
            temp1 = (w / np.linalg.norm(w)) * np.sin(np.linalg.norm(w) * dt / 2)
            interpquat = Quaternion(np.cos(np.linalg.norm(w) * dt / 2), temp1[0], temp1[1], temp1[2])
            q = q_prev * interpquat

        RPY = quaternion_to_euler(q)
        attitude_angle = np.append(attitude_angle, RPY, axis=1)
        qerror = (qref.unit).conjugate * q.unit
        th = 2 * np.arccos(qerror[0])
        axiserror = np.array([[qerror[1]], [qerror[2]], [qerror[3]]])

        if qerror[0] < 0:
            axiserror = -1.0 * axiserror

        Integral = Integral + axiserror
        if dur == 1:
            axiserrorprev = axiserror

        u = -Kp * axiserror - Ki * Integral - Kwmega * w

        A = np.array([[-cos(b) * cos(d1), sin(d2), cos(b) * cos(d3), -sin(d4)],
                      [-sin(d1), -cos(b) * cos(d2), sin(d3), cos(b) * cos(d4)],
                      [sin(b) * cos(d1), sin(b) * cos(d2), sin(b) * cos(d3), sin(b) * cos(d4)]])

        h = h0 * np.array([[-cos(b) * sin(d1), -cos(d2), cos(b) * sin(d3), cos(d4)],
                           [cos(d1), -cos(b) * sin(d2), -cos(d3), cos(b) * sin(d4)],
                           [sin(b) * sin(d1), sin(b) * sin(d2), sin(b) * sin(d3), sin(b) * sin(d4)]])

        hx = np.sum(h[0, :])
        hy = np.sum(h[1, :])
        hz = np.sum(h[2, :])
        Hvector = np.array([[hx], [hy], [hz]])
        Hdot = -u - cross(w.transpose(), Hvector.transpose()).transpose()

        # mmm = np.linalg.det(np.dot(A,A.transpose()))
        mmm = np.linalg.det(np.dot(A, A.transpose()))
        manip = np.append(manip, np.sqrt(mmm))

        # pi=math.pi
        # l_gain = 1
        # mi = 0.010
        # e0 = 0.01
        # freq = 100
        # phi1 = 0
        # phi2 = pi/2
        # phi3 = pi
        # lamda = l_gain * math.exp(mi * mmm)
        # e1 = math.fabs(e0 * sin(freq * dur + phi1))
        # e2 = math.fabs(e0 * sin(freq * dur + phi2))
        # e3 = math.fabs( e0 * sin(freq * dur + phi3))
        # E =skew([e1, e2, e3])
        # A_hash = np.dot (A.T , np.linalg.inv(  np.dot(A, A.T) + lamda * E))
        A_hash = np.linalg.pinv(A)

        d_dot = np.dot(A_hash, Hdot)  # or np.dot(A.T,np.linalg.inv(np.dot(A,A.T)))

        # %%%%%%%%%%%% SATURATION %%%%%%%%%%%%%%
        ddmax = np.amax(np.absolute(d_dot))
        ddlim = np.deg2rad(50)
        ddgamma = ddmax / ddlim
        if ddgamma > 1:
            d_dot = d_dot / ddgamma

        d = dt * d_dot + dprev

        d1 = d[0]
        d2 = d[1]
        d3 = d[2]
        d4 = d[3]
        dprev = d

        g1 = np.array([[sin(b)], [0], [cos(b)]])
        g2 = np.array([[0], [sin(b)], [cos(b)]])
        g3 = np.array([[-sin(b)], [0], [cos(b)]])
        g4 = np.array([[0], [-sin(b)], [cos(b)]])

        d1_vect = d_dot[0] * g1
        d2_vect = d_dot[1] * g2
        d3_vect = d_dot[2] * g3
        d4_vect = d_dot[3] * g4

        Tc = np.array(
            cross(h[:, 0], d1_vect.transpose()).transpose() + cross(h[:, 1], d2_vect.transpose()).transpose() + cross(
                h[:, 2], d3_vect.transpose()).transpose() + cross(h[:, 3], d4_vect.transpose()).transpose())

        q_prev = q
        w_prev = w

        # dd = np.append(dd, np.rad2deg(d.transpose()), axis=0)
        # wplot = np.append(wplot, w, axis=1)
        # hplot = np.append(hplot, Hvector, axis=1)
        # d_dotplot = np.append(d_dotplot, d_dot, axis=1)
        # errorplot = np.append(errorplot, quaternion_to_euler(qerror), axis=1)

    return np.min(manip)



def cf(pop):
    f=[0]*len(pop)
    for i in range(0, len(pop)):
        chromosome = pop[i, :]
        a = chromosome[0]
        b = chromosome[1]
        c = chromosome[2]
        d = chromosome[3]

        # f[i] = a - b
        # f[i]=abs(f[i]) #minimize this function
        minmanip=mainCMG_code(a,b,c,d)
        f[i]=1/minmanip

    return f

def pie(f):
    p = [0] * len(f)
    r = [0] * len(p)
    ss=np.sum(f)
    for i in range(0, len(f)):
        p[i]=f[i]/ss
    for i in range(0, len(p)):
        if i==0:
            r[i] = p[i]
        else:
            r[i] = r[i-1] + p[i]
    return r

def selection(r,pop):
    prob=np.random.uniform(0,1, size=(len(pop), 1))
    indices=[0] * len(prob)
    for j in range(0,len(prob)):
        i=1
        while i<len(r):
            if prob[j]>r[len(r)-1]:
                indices[j]=len(r)-1
                break
            elif prob[j]<r[i]:
                indices[j]=i-1
                break
            i=i+1

    selectedpop=pop[indices]
    return selectedpop

def crossmut(crossprob,mutationprob,pop,LB,UB,numofgenes):
    #crossover
    for i in range(0, round(len(pop) / 2) ,2):
        if (np.random.uniform(0, 1) < crossprob ):
            rindex1=rnd.randint(0, numofgenes-1)
            rindex2 = rnd.randint(0, numofgenes-1)
            item1=pop[i,rindex1]
            item2=pop[i+1,rindex2]
            pop[i, rindex1]=item2
            pop[i+1, rindex2] = item1
    #mutation
    for i in range(0, len(pop)):
        if (np.random.uniform(0, 1) < mutationprob):
            pop[i, rnd.randint(0, numofgenes-1)]=np.round(np.random.uniform(LB,UB)).astype(int)

    return pop

########### PARAMETERS ##########
population_size=40
generations_limit=200
numofgenes=4
stall_limit=20
UB=math.pi
LB=-math.pi
crossprob = 0.1
mutationprob = 0.3
################################
stall=0
rep=0
errortol=0.02
printtable = Texttable()

pop=np.random.uniform(LB,UB,size=(population_size,numofgenes))
f=cf(pop)
while rep<generations_limit and stall<stall_limit:
    r=pie(f)
    pop=selection(r,pop)
    pop=crossmut(crossprob,mutationprob,pop,LB,UB,numofgenes)
    f = cf(pop)

    if rep != 0:
        if abs(np.min(f) - argument_prev) < errortol:
            stall = stall + 1
        else:
            stall = 0

    # %%%%% PRINT STUFF %%%%%%%
    printtable.add_rows([['Generation', 'Stall'], [rep,stall]])
    print(printtable.draw())
    # %%%%%%%%%%%%%%%%%%%%%%%%%%

    argument_prev = np.min(f)
    rep=rep+1



print(pop)
print('best solution is given by child')
print(pop[np.argmin(f)])
print('best solution is equal to')
print( f[np.argmin(f)] )