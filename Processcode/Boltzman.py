import numpy as np
import math
def LoadMatrix():
    aaDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
            "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
            "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
            "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}
    List = aaDict.keys()
    List.sort()
    for amino1 in List:
        for amino2 in List:
            aaDict[amino1][amino2] = np.zeros((20,20))
    return aaDict,List
def loadResidueTypeNum():
    ResidueType = '/Users/xg666/Desktop/NeighborChain/ResidueTypeNum.npy'
    ResidueTypeNumber = np.load(ResidueType).item()
    allResidue = sum(ResidueTypeNumber.values())
    return ResidueTypeNumber,allResidue
# calculate pi,which represent Amino Acid i's probability occur in Neighbors
def CalculatePa(pAnum,allAnum):
    pa = float(pAnum)/float(allAnum)
    return pa
# calculate pij,which represent the probability about i is neighbor of j
def CalculatePij(Nij,allNij):
    pij = float(Nij)/float(allNij)
    return pij
# calculate Pij(theta,phi,r)
def CalcuPijcoordinate(Nijcoordinate,Nij):
    pijcoordinate = float(Nijcoordinate)/float(Nij)
    return pijcoordinate
# calculate pi*pj
def CalculatePiPj(pi,pj):
    pipj = pi*pj
    return pipj
#load the coodninatd dict and number dict
def LoadCoordinateNumber(filec):
    coordict = np.load(filec).item()
    #numdict = np.load(filen).item()
    return coordict
#sum Ni,i = range(0,20)
def sumNi(numdict,List):
    sumN1_20 = 0
    for amino in List:
        sumN1_20 += sum(numdict[amino].values())
    return sumN1_20
def Calculatepij(Nij,allAnum):
    pij = float(Nij) /float(allAnum)
    return pij
# have foud smallest value is 2
def FindSmallNijCoordinate(coordict,List):
    small_value = 100
    for amino1 in List:
        for amino2 in List:
            for m in range(0,20):
                for n in range(0,20):
                    if coordict[amino1][amino2][m][n] < small_value:
                        small_value = coordict[amino1][amino2][m][n]
                        #print small_value
    return small_value

#calculate the number about a amino acid the frequency of occurence
def SumMatrixs(matrixs):
    allNumberAmino = 0
    for i in range(0,20):
        #print matrixs[i]
        #print sum(sum(matrixs[i]))
        allNumberAmino+= sum(sum(matrixs[i]))
    return allNumberAmino

# calculate the occurence about A1 and A2
def sumMatrixaiaj(matrix):
    numAiAj = sum(sum(matrix))
    return numAiAj
# calculate all the amino acid's frequency of accurence
def calculateallnumber(coordict,List):
    allAnum = 0
    for amino1 in List:
        matrixsamnio1 = coordict[amino1].values()
        allAnumamino1 = SumMatrixs(matrixsamnio1)
        allAnum += allAnumamino1
    return allAnum

def zerosMatrix():
    Dict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
            "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
            "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
            "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}
    List = Dict.keys()
    List.sort()
    for amino1 in List:
        for amino2 in List:
            Dict[amino1][amino2] = np.zeros((20,20))
    return Dict


#calculateenergy
def CalculateEnergy(filec):
    coordict = LoadCoordinateNumber(filec)
    aaDict,List = LoadMatrix()
    allnumNeighbor = calculateallnumber(coordict,List)
    expDict = zerosMatrix()
    objDict = zerosMatrix()
    objexpDict = zerosMatrix()
    ResidueTypeNumbermatrix,allResidue = loadResidueTypeNum()
    #molecule = 0
    sumpi = 0
    sumpij = 0

    sumexp = 0
    sumobj = 0
    #amino1 is the Amino Acid i
    for amino1 in List:
        pinum = ResidueTypeNumbermatrix[amino1]
        pi = CalculatePa(pinum,allResidue)
        sumpi = sumpi + pi
        #amino2 is Amino Acid j
        for amino2 in List:
            sumpijcoordinate = 0
            # calculate the data about amino acid i
            # calculate the data about amino acid is
            #pjmatrixs = coordict[amino2].values()
            pjnum = ResidueTypeNumbermatrix[amino2]
            pj = CalculatePa(pjnum,allResidue)
            #calculate pipj
            pipj = CalculatePiPj(pi,pj)
            matrixNij = coordict[amino1][amino2]
            Nij = sumMatrixaiaj(matrixNij)
            pij = Calculatepij(Nij,allnumNeighbor)
            sumpij = sumpij + pij
            sumNIJ = 0
            #print pij,pipj,pi,pj,allAnum,Nij
            for m in range(0,20):
                #molecule += np.cos(m*(np.pi/20)) - np.cos((m+1)*(np.pi/20))
                for n in range(0,20):
                    theta = (m *np.pi)/20
                    phi = (n * np.pi )/10
                    #theta = min(int(math.floor(theta*20/np.pi)),19)
                    #phi = min(int(math.floor(phi*10/np.pi) + 10),19)
                    Nijcoordinate = coordict[amino1][amino2][m][n]
                    if Nijcoordinate == 0:
                        Nijcoordinate = 0.0001
                    sumNIJ += Nijcoordinate
                    pijcoordinate = CalcuPijcoordinate(Nijcoordinate,Nij)
                    sumpijcoordinate += pijcoordinate
                    objected = pij * pijcoordinate
                    expection = (pipj * (np.cos(m*(np.pi/20)) - np.cos((m+1)*(np.pi/20))))/40
                    expDict[amino1][amino2][m][n] = expection
                    objDict[amino1][amino2][m][n] = objected
                    objexpDict[amino1][amino2][m][n] = objected/expection
                    sumexp+=expection
                    sumobj+=objected
                    energy = -math.log(objected/expection)

                    #print objected,expection,objected/expection
                    #print pij,pipj,pijcoordinates
                    #print energy
                    #energy = - energy
                    print energy,sumexp,sumobj
                    print 'hhhhhhhhhh'
                    aaDict[amino1][amino2][m][n] = energy
        #print Nij,sumNIJ,float(sumNIJ)/float(Nij)
        print sumpijcoordinate,sumpij
    return aaDict,expDict,objDict,objexpDict
def savefile(file,result):
    savefile = file
    np.save(savefile,result)


if __name__ == "__main__":
    file = '/Users/xg666/Desktop/NeighborChain/resultnpy/coordinatesum.npy'
    fileenergymatrix = '/Users/xg666/Desktop/NeighborChain/resultnpy/energymatrixsum.npy'
    fileexp = '/Users/xg666/Desktop/NeighborChain/resultnpy/expsum.npy'
    fileobj = '/Users/xg666/Desktop/NeighborChain/resultnpy/objsum.npy'
    fileobjexp = '/Users/xg666/Desktop/NeighborChain/resultnpy/objexpsum.npy'
    aaDict,expDict,objDict,objexpDict = CalculateEnergy(file)
    savefile(fileenergymatrix,aaDict)
    savefile(fileexp,expDict)
    savefile(fileobj,objDict)
    savefile(fileobjexp,objexpDict)

