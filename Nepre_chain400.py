import os
import math
import numpy as np
import AminoAcid as AA
import gc
import sys
import math
import csv

def listdirInMac(path):
    os_list = os.listdir(path)
    for item in os_list:
        if item.startswith('.') and os.path.isfile(os.path.join(path, item)):
            os_list.remove(item)
    return os_list

def LoadRadius():
    radiusDict = {"ALA":0,"VAL":0,"LEU":0,"ILE":0,"PHE":0,\
                  "TRP":0,"MET":0,"PRO":0,"GLY":0,"SER":0,\
                  "THR":0,"CYS":0,"TYR":0,"ASN":0,"GLN":0,\
                  "HIS":0,"LYS":0,"ARG":0,"ASP":0,"GLU":0,}

    f = open("./mean_radius.txt")
    for line in f.readlines():
        temp = line.strip().split()
        if(temp[0] != "Name"):
            radiusDict[temp[1]] = float(temp[0])
    return radiusDict

def pearson(rmsd,energy):
    size = np.shape(rmsd)[0]
    x = np.empty(shape=[2,size])
    for i in range(size):
        x[0][i] = rmsd[i]
    for j in range(size):
        x[1][j] = energy[j]
    y = np.corrcoef(x)
    return y[0][1]

def load_EnergyMatrixsum400():
    filesum400 = './Energymatrixs/chainenergymatrixsum3_400.npy'
    aaDictsum400 = np.load(filesum400).item()
    return aaDictsum400

def load_EnergyMatrixsum():
    filesum = './Energymatrixs/chainenergymatrixsum3.npy'
    aaDictsum = np.load(filesum).item()
    return aaDictsum


def extract_Data(line):
    """
    This part will extracted data from line according to the standard
    PDB file format(Version 3.3.0, Nov.21, 2012)
    """
    res = []

    line = line.strip()
    #record_name
    res.append(line[0:4].strip(' '))

    #atom_serial
    res.append(line[6:11].strip(' '))

    #atom_name
    res.append(line[12:16].strip(' '))

    #alternate_indicator
    res.append(line[16])

    #residue_name
    res.append(line[17:20].strip(' '))

    #chain_id
    res.append(line[21].strip(' '))

    #residue_num
    res.append(line[22:26].strip(' '))

    #xcor
    res.append(line[30:38].strip(' '))

    #ycor
    res.append(line[38:46].strip(' '))

    #zcor
    res.append(line[46:54].strip(' '))

    return res

def getlines_for_eachchain(f):
    chains = []
    chainA = []
    with open(f,'r') as file:
        for line in file.readlines():
            if (line[0:4] == 'ATOM'):
                chainA.append(line)
            elif (line[:3].strip() == 'TER'):
                chainA.append(line)
                chains.append(chainA)
                chainA = []
    return chains

def ProcessPDB(chainlines,matrix):
    #df = open(file,'r')
    radiusDict = LoadRadius()
    CurrentAANitrogen = None
    CurrentAACA = None
    Currentresidue_num = None
    EachAA = []
    CurrentAA = None

    for line in chainlines:
        if(line[0:4] != "ATOM"):
            continue
        element_list = extract_Data(line)
        record_name = element_list[0]
        atom_name = element_list[2]
        residue_name = element_list[4]
        alternate_indicator = element_list[3]
        residue_num = element_list[-4]
        xcor = float(element_list[-3])
        ycor = float(element_list[-2])
        zcor = float(element_list[-1])

        if(atom_name == "H"):
            continue
        if(residue_name not in matrix):
            continue

        if(CurrentAA == None):
            CurrentAA = AA.AminoAcid(residue_name)
            Currentresidue_num = residue_num
            if(atom_name == "N" or atom_name == "CA"):
                if(alternate_indicator == "B"):
                    continue
                if(atom_name == "N"):
                    CurrentAANitrogen = np.array([xcor,ycor,zcor])
                else:
                    CurrentAACA = np.array([xcor,ycor,zcor])
            if(residue_name == "GLY" or atom_name not in {"N","CA","C","O","O1","02"}):
                if(alternate_indicator != " "):
                    #If cases like "AASN or BASN" appears, we only add A
                    if(alternate_indicator == "A" and line[15] == "1"):
                        CurrentAA.SumCenters(xcor,ycor,zcor)
                    else:
                        continue
                else:
                    CurrentAA.SumCenters(xcor,ycor,zcor)
        else:
            #If another amino acid begins
            if(residue_num != Currentresidue_num):
                state = CurrentAA.CalculateCenter()
                if(state == False):
                    CurrentAA = AA.AminoAcid(residue_name)
                    Currentresidue_num = residue_num
                    continue

                CurrentAA.InputCAN(CurrentAANitrogen,CurrentAACA)
                EachAA.append(CurrentAA)
                del CurrentAA
                CurrentAA = AA.AminoAcid(residue_name)

                Currentresidue_num = residue_num
                if(atom_name == "N" or atom_name == "CA"):
                    if(alternate_indicator == "B"):
                        continue
                    if(atom_name == "N"):
                        CurrentAANitrogen = np.array([xcor,ycor,zcor])
                    else:
                        CurrentAACA = np.array([xcor,ycor,zcor])
                if(residue_name == "GLY" or atom_name not in {"N","CA","C","O","O1","02"}):
                    if(alternate_indicator != " "):
                    #If cases like "AASN or BASN" appears, we only add A
                        if(alternate_indicator == "A" and line[15] == "1"):
                            CurrentAA.SumCenters(xcor,ycor,zcor)
                        else:
                            continue
                    else:
                        CurrentAA.SumCenters(xcor,ycor,zcor)
            #If still the same amino acid
            else:
                if(atom_name == "N" or atom_name == "CA"):
                    if(alternate_indicator == "B"):
                        continue
                    if(atom_name == "N"):
                        CurrentAANitrogen = np.array([xcor,ycor,zcor])
                    else:
                        CurrentAACA = np.array([xcor,ycor,zcor])
                if(residue_name == "GLY" or atom_name not in {"N","CA","C","O","O1","02"}):
                    if(alternate_indicator != " "):
                    #If cases like "AASN or BASN" appears, we only add A
                        if(alternate_indicator == "A" and line[15] == "1"):
                            CurrentAA.SumCenters(xcor,ycor,zcor)
                        else:
                            continue
                    else:
                        CurrentAA.SumCenters(xcor,ycor,zcor)

    state = CurrentAA.CalculateCenter()
    if(state != False):
        CurrentAA.CalculateCenter()
        CurrentAA.InputCAN(CurrentAANitrogen,CurrentAACA)
        EachAA.append(CurrentAA)
    return EachAA

def change_chain(chains,radiusDict):
    changedchains = []
    for g in range(len(chains)):
        chain_processed = ProcessPDB(chains[g],radiusDict)
        changedchains.append(chain_processed)
    return changedchains

def calculate_Energy(changedchains,radiusDict,matrix,matrix400):
    E = 0
    E400 = 0
    length = len(changedchains)
    #print length
    for i in range(length):
        #print i
        for j in range(length):
            #print j
            if i != j:
                chainA = changedchains[i]
                chainB = changedchains[j]
                for m in range(len(chainA)):
                    chainA[m].EstablishCoordinate()
                    for n in range(len(chainB)):
                        dis = chainA[m].DistanceBetweenAA(chainB[n].center)
                        radiusSum = radiusDict[chainA[m].name] + radiusDict[chainB[n].name] + 3
                        if(dis <= radiusSum):
                            #print 'distance <=sum'
                            rho,theta,phi = chainA[m].ChangeCoordinate(chainB[n].center)
                            theta = min(int(math.floor(theta*20/np.pi)),19)
                            phi = min(int(math.floor(phi*10/np.pi) + 10),19)
                            #print matrixsum[chainA[m].name][chainB[n].name]
                            E += matrix[chainA[m].name][chainB[n].name][theta][phi] / rho
                            E400 += matrix400[chainA[m].name][chainB[n].name]
    return E,E400

def energymain(pdb):
    matrix400 = load_EnergyMatrixsum400()
    matrix = load_EnergyMatrixsum()
    #print matrixsum
    radiusDict = LoadRadius()
    chains = getlines_for_eachchain(pdb)
    length =len(chains)
    print length
    changedchains = change_chain(chains,radiusDict)
    E = calculate_Energy(changedchains,radiusDict,matrix,matrix400)
    return E



if __name__ == "__main__":

    args = sys.argv[1:]
    pdb = args[0]
    print pdb
    #f = open(pdb)
    Esum = energymain(pdb)
    print "Nepre Potential Energy(distance is radiusSum)",Esum
