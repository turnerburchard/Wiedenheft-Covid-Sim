##### Here is a model for the spread of a pandemic. The spread of a disease
##### in a cell will go as dI/dt = Co*exp(-r(t-to)^2). It might not,
##### though. Pick a cell size between 30 and 60. Let MATLAB decide how to
##### size the cells. Cells should overlap. If someone in your cell is a 2,
##### you should be much more likely to become a 2 than if no one in your
##### cell is a 2. How likely you are to become a 2 will have an Arrhenius
##### dependence. Each cell will have a temperature. The spread condition
##### will be an individual event. A value of '0' will imply that the
##### individual does not have covid. A value of '1' will indicate that
##### the individual has covid and is contagious. A value of '2' will
##### indicate that the individual is removed from the stream. Base this
##### spread on the Arrhenius rate law :-)

import numpy as np
import random

def PanSim(temp):
    GoodFit = 0
    R2D2Best = 1e9
    
    Threshold = 0.0675
    CloseVal = 9
    FarVal = 13
    CellOdds = 5

    while GoodFit == 0:
        nIndividuals = 0
        k1Max = 76
        k2Max = 76

        America = []
        
        for i in range(k1Max * k2Max):
            nSize = random.randint(5, 10)
            Row = np.zeros([nSize,nSize,2])
            America.append(Row)
            
            nIndividuals = nIndividuals + nSize*nSize

        print('Number of People in Simulation: ' + str(nIndividuals))

        CellOne = America[1]

        CellOne[1,1,0] = 1
        CellOne[1,1,1] = 1
        America[1] = CellOne

        CellTwo = America[2]
        CellTwo[2,1,0] = 1
        America[2] = CellTwo

        CellThree = America[4]
        CellThree[1,2,0] = 1
        America[4] = CellThree
        
        rngMax = 1e6
        Contagious = 1.5*2.314497e-16
        gbContagious = 1.5*2.314497e-16
        ConTime, Do, C, R, Tmin, Tmax = 14, 6.8668, 1e-8, 8.314, 1380, 1570
        nDays = 253
        time = list(range(1,251)) #all integers 1 to 250 inclusive. Should this be different for indexing?
        T = np.zeros((nDays, 1))

        QVol = 1.7555*475e3*np.ones((nDays, 1))
        Qgb = 1.25*532e3*np.ones((nDays, 1))

        DayVec = np.zeros((nDays,1))
        dIdtVec = np.zeros((nDays,1))
        InterdIdtVec = np.zeros((nDays,1))
        IntradIdtVec = np.zeros((nDays,1))

        InterEx, IntraEx, InterTrans, IntraTrans, InterExpAtt, IntraExpAtt = 0, 0, 0, 0,0,0

        for day in range(0, nDays):
            DailyDiff = 100
            ECOunt = 0
            if day > 1:
                QVol[day] = QVol[day - 1]
                Qgb[day] = Qgb[day-1]0 

        GoodFit = 1


def main():
    PanSim("init")

main()