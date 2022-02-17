import numpy as np
import matplotlib.pyplot as plt
from bin.AA import AA,IA
import math
#
def plotGraph(targetFun, useAA, title):
    sx = -2.0
    ex = +2.0
    numStep = 20
    stepSizeHalf = (ex - sx) * 0.5 / (numStep-1)
    x = np.linspace(sx, ex, numStep)
    yl = np.zeros(numStep)
    yh = np.zeros(numStep)
    plt.style.use('dark_background')
    fig, ax = plt.subplots(figsize=(10.0, 6.0))
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(0.0, 1.6)
    for (i,_x) in enumerate(np.nditer(x)):
        ia = IA(_x-stepSizeHalf, _x+stepSizeHalf)
        if useAA:
            ia = AA(ia)
            ia = targetFun(ia)
            ia = ia.toIA()
        else:
            ia = targetFun(ia)
        yl[i] = ia.low()
        yh[i] = ia.high()
        x0 = _x-stepSizeHalf
        x1 = _x+stepSizeHalf
        y0 = ia.low()
        y1 = ia.high()
        ax.add_patch(plt.Rectangle((x0,y0), x1-x0, y1-y0, fc='#A0A0A0', ec='black'))
    x = np.linspace(sx, ex, 200)
    y = np.zeros(200)
    for (i,_x) in enumerate(np.nditer(x)):
        y[i] = targetFun(_x)
    plt.title(title)
    plt.plot(x, y, color='red')
    plt.savefig(f"{title}.png",transparent=True)
    plt.show()

def g(x):
    if  isinstance(x, AA) or isinstance(x, IA):
        return (x * x - x + 0.5).sqrt() / (x * x + 0.5).sqrt()
    else:
        return np.sqrt(x * x - x + 0.5) / np.sqrt(x * x + 0.5)

def h(x):
    return g(g(x))

plotGraph(g, False, 'g(x) with IA')
plotGraph(h, False, 'h(x) with IA')
plotGraph(g, True,  'g(x) with AA')
plotGraph(h, True,  'h(x) with AA')