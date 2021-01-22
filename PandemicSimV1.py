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

def PanSim():
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

        America = np.empty([k1Max, k2Max])

        America = list(America)
        
        for i in range(len(America)):
            nSize = random.randint(5, 10)
            Row = np.zeros([nSize,nSize,2])
            #America.append(Row)
            #np.append(America, Row)
            America[i] = Row
            
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
            ECount = 0
            if day > 1:
                QVol[day] = QVol[day - 1]
                Qgb[day] = Qgb[day-1] 
            while DailyDiff > 10:
                 
                 #Maybe there should be more here? lots of comments
                print('Day ' + str(day))

                InterdIdt = np.zeros(1)
                IntradIdt = np.zeros(1)

                for k1 in (0, k1Max):
                    for k2 in (0, k2Max):
                        InterStack = America[k1-1][k2-1] # Real possible off-by-one error here
                        InterMatInf = InterStack[:,:,0] 
                        InterMatRec = InterStack[:,:,1]
                        #My America doesn't have enough dimensions to index this? 
                        #I think it would require a 5-dimensional array

                        GateOne = [i for i in InterMatInf if i == 1]

                        GateTwo = [i for i in InterMatRec if i>=2]

                        SumInf = 0
                        for i in range(len(GateOne)):
                            if GateOne[i] == 1 and GateTwo[i] == 1:
                                SumInf += 1
                        
                        FracInf = SumInf/len(InterMatInf)


                        if SumInf > 0:
                            for i in range(-1, 1): #is this right?
                                for j in range(-1, 1):
                                    if day <651:
                                        T[day] = 2850 
                                    
                                    InterExpAtt += 1
                                    if (FracInf > Threshold) and (1/(1 - FracInf)*np.exp(-Qgb(day)/(R*T(day))) > gbContagious):
                                        InterEx += 1
                                        dx = i
                                        dy = j
                                        if (k1 + dx) > k1Max or (k1 + dx) == 0:
                                            dx = 0
                                        if (k2 + dy) > k2Max or (k2 + dy) == 0:
                                            dy = 0

                                        if (abs(dx) + abs(dy)) == 1:
                                            Enya = CloseVal
                                        elif (abs(dx) + abs(dy)) == 2:
                                            Enya = FarVal

                                            if random.randint(1, Enya) == random.randint(1, Enya):
                                                InterTrans += 1
                                                NeighborStack = America[k1+dx][k2+dy] # this part will be messy,  need indexing context
                                                NeighborMat = NeighborStack[:,:,0]
                                                NeighborRec = NeighborStack[:,:,1]

                                                NeiShape = np.shape(NeighborMat)

                                                NeiX = random.randint(1, NeiShape[0])
                                                NeiY = random.randint(1, NeiShape[1])

                                                if NeighborMat[NeiX][NeiY] == 0:
                                                    NeighborMat[NeiX][NeiY] = 1
                                                    NeighborRec[NeiX][NeiY] = 1

                                                    InterdIdt += 1

                                                    NeighborStack[:,:,0] = NeighborMat
                                                    NeighborStack[:,:,1] = NeighborRec
                                                    America[k1+dx][k2+dy] = NeighborStack

                for k in range(len(America)):

                    TheMatrix = America[k]
                    CelldIdt = 0
                    DisMat = TheMatrix[:,:,0]

                    if DisMat == 1:
                        SumOne = DisMat  

                    if SumOne > 0:
                        Col, Row = np.shape(TheMatrix[:,:,0])

                        for i in range(Col):
                            for j in range(Row):

                                if TheMatrix[i][j][0] == 1:
                                    TheMatrix[i][j][1] += 1

                                    if TheMatrix[i][j][1] > ConTime:
                                        TheMatrix[i][j][0] = 2
                                    
                        CellStart = sum([ii for ii in TheMatrix[:,:,0] if ii == 1])

                        Rho = SumOne/len(TheMatrix[:,:,0])

                        if CellStart >= round(np.size(TheMatrix[:,:,0])/2):
                            nEnd = np.size(TheMatrix[:,:,0])
                        else:
                            nEnd = 2*CellStart
                        
                        for n in range(nEnd):
                            IntraExpAtt += 1
                            xInf = random.randint(0,Row) #should this be 0 or 1?
                            yInf = random.randint(0,Col) #again

                            if DisMat[xInf][yInf] == 0:
                                DisMat[xInf][yInf] = 1
                                IntradIdt += 1
                                CelldIdt += 1

                if CelldIdt > 2*CellStart:
                    print('Too many infections')

            dIdt = InterdIdt + IntradIdt

            if dIdt > 0:
                TheMatrix[:][:][0] = DisMat      
                America[k] = TheMatrix

        dIdt = InterdIdt + IntradIdt 

        DayVec[day] = day
        InterdIdtVec[day] = InterdIdt
        IntradIdtVec[day] = IntradIdt
        dIdtVec[day] = dIdt


        if day > 7:
            dIdt7 = np.mean(dIdtVec[day-7:day])
            DailyDiff = abs(dIdt7 - didt7Data[day]) #didt7Data is not defined before here?

            if DailyDiff > 10:
                ECount += 1

                if ECount < 16:

                    if dIdt7 > didt7Data(day):

                        if np.sum(dIdtVec) < 215:
                            QVol[day] = 1.0125*QVol[day]
                        else:
                            QVol[day] = 1.0125*QVol[day]
                            Qgb[day] = 1.0125*Qgb[day]
                        
                        if ECount > 3:
                            print("Threshold too low: Increasing Threshold\n")
                            Threshold = 1.01*Threshold
                            if ECount > 7:
                                CellOdds +- 1
                                print('Cell Contagiousness Decreased\n')
                        
                    elif dIdt7 < didt7Data[day]:
                        if sum(dIdtVec) < 215:
                            QVol[day] = 0.9875*QVol[day]
                        else:
                            QVol[day] = 0.9875*QVol[day]
                            Qgb[day] = 0.9875*Qgb[day]
                        
                        if ECount > 3:
                            print('Threshold too high: Decreasing Threshold\n')

                            Threshold = 0.99*Threshold

                            if ECount > 7:
                                CellOdds -= 1
                                print('Cell Contagiousness Increased\n')
                    
                    if CellOdds < 1:
                        CellOdds = 1
                    
                else:
                    QVol[day] = QVol[day-1] #Need to be careful with indexing errors here?
                    Qgb[day] = Qgb[day -1]
                    print('Energy change too extreme\n')
                    DailyDiff = 5
        else:
            DailyDiff = 5

        GoodFit = 1 #temporary way to exit the initial while loop

    dIdt7Sim = np.zeros([len(dIdtVec)-7,1])
    for n in range(7, len(dIdtVec)): #off by one error?
        dIdt7Sim[n-7] = np.mean(dIdtVec[n-7:n])
    
    SimNum = len(dIdt7Sim)

    if SimNum >= len(didt7Data):
        SimComp = dIdt7Sim[1:len(didt7Data)]

    #Weird end statements here, not sure that I'm in the correct indentations

                        








                        


        


def main():
    PanSim()

main()