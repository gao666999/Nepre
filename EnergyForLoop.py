import os
import math
import numpy as np
import AminoAcidv2 as AA
import gc
import sys
import math
#import csv
#import matplotlib.pyplot as plt
#from matplotlib.pyplot import savefig


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


def load_EnergyMatrix():
    aaDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
            "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
            "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
            "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}

    List = aaDict.keys()
    List.sort()

    f1 = open("./Energymatrixs/radius.npy")
    for amino1 in List:
        for amino2 in List:
            aaDict[amino1][amino2] = np.load(f1)
    f1.close()
    return aaDict


def load_EnergyMatrix400():
    filesum = './Energymatrixs/energymatrixsum3file7_400.npy'
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


def process_pdb_file(pdbfile,matrix):

    CurrentAANitrogen = None
    CurrentAACA = None
    Currentresidue_num = None
    EachAA = []
    CurrentAA = None
    f = open(pdbfile)
    lines = f.readlines()
    f.close()
    #pdbfile = lines
    for line in lines:
        if (line[0:4] != "ATOM"):
            continue
        element_list = extract_Data(line)
        record_name = element_list[0]
        atom_name = element_list[2]
        residue_name = element_list[4]
        alternate_indicator = element_list[3]
        residue_num = element_list[-4]
        chain_id = element_list[-5]
        xcor = float(element_list[-3])
        ycor = float(element_list[-2])
        zcor = float(element_list[-1])

        if (atom_name == "H"):
            continue
        if (residue_name not in matrix):
            continue

        if (CurrentAA == None):
            CurrentAA = AA.AminoAcid(residue_name,residue_num,chain_id)
            Currentresidue_num = residue_num
            if (atom_name == "N" or atom_name == "CA"):
                if (alternate_indicator == "B"):
                    continue
                if (atom_name == "N"):
                    CurrentAANitrogen = np.array([xcor, ycor, zcor])
                else:
                    CurrentAACA = np.array([xcor, ycor, zcor])
            if (residue_name == "GLY" or atom_name not in {"N", "CA", "C", "O", "O1", "02"}):
                if (alternate_indicator != " "):
                    # If cases like "AASN or BASN" appears, we only add A
                    if (alternate_indicator == "A"):
                        CurrentAA.SumCenters(xcor, ycor, zcor)
                    else:
                        continue
                else:
                    CurrentAA.SumCenters(xcor, ycor, zcor)
        else:
            # If another amino acid begins
            if (residue_num != Currentresidue_num):
                state = CurrentAA.CalculateCenter()
                if (state == False):
                    CurrentAA =AA.AminoAcid(residue_name,residue_num,chain_id)
                    Currentresidue_num = residue_num
                    continue

                CurrentAA.InputCAN(CurrentAANitrogen, CurrentAACA)
                CurrentAA.EstablishCoordinate()
                # Amino Acid check
                EachAA.append(CurrentAA)
                del CurrentAA
                CurrentAA = AA.AminoAcid(residue_name,residue_num,chain_id)

                Currentresidue_num = residue_num
                if (atom_name == "N" or atom_name == "CA"):
                    if (alternate_indicator == "B"):
                        continue
                    if (atom_name == "N"):
                        CurrentAANitrogen = np.array([xcor, ycor, zcor])
                    else:
                        CurrentAACA = np.array([xcor, ycor, zcor])
                if (residue_name == "GLY" or atom_name not in {"N", "CA", "C", "O", "O1", "02"}):
                    if (alternate_indicator != " "):
                        # If cases like "AASN or BASN" appears, we only add A
                        if (alternate_indicator == "A"):
                            CurrentAA.SumCenters(xcor, ycor, zcor)
                        else:
                            continue
                    else:
                        CurrentAA.SumCenters(xcor, ycor, zcor)
            # If still the same amino acid
            else:
                if (atom_name == "N" or atom_name == "CA"):
                    if (alternate_indicator == "B"):
                        continue
                    if (atom_name == "N"):
                        CurrentAANitrogen = np.array([xcor, ycor, zcor])
                    else:
                        CurrentAACA = np.array([xcor, ycor, zcor])
                if (residue_name == "GLY" or atom_name not in {"N", "CA", "C", "O", "O1", "02"}):
                    if (alternate_indicator != " "):
                        # If cases like "AASN or BASN" appears, we only add A
                        if (alternate_indicator == "A"):
                            CurrentAA.SumCenters(xcor, ycor, zcor)
                        else:
                            continue
                    else:
                        CurrentAA.SumCenters(xcor, ycor, zcor)

    state = CurrentAA.CalculateCenter()
    if (state != False):
        #CurrentAA.CalculateCenter()
        CurrentAA.InputCAN(CurrentAANitrogen, CurrentAACA)
        CurrentAA.EstablishCoordinate()
        EachAA.append(CurrentAA)
    return EachAA



def calculate_Energy(wholefile, loopfile, matrix,matrix400):
    #print wholefile,loopfile
    whole_residues = process_pdb_file(wholefile,matrix)
    loop_residues = process_pdb_file(loopfile,matrix)
    radiusDict = LoadRadius()
    CurrentAAlist = []
    Energy = 0
    Energy400 = 0
    #avoid to calculate a residue more than once
    for currentResidue in loop_residues:
        CurrentAAlist.append([currentResidue.name,currentResidue.idnumber,currentResidue.chainID])
        for Residue_nn in whole_residues:
            if ([Residue_nn.name,Residue_nn.idnumber,Residue_nn.chainID] in CurrentAAlist):
                #print [Residue_nn.name,Residue_nn.idnumber,Residue_nn.chainID]
            #if(currentResidue.idnumber == Residue_nn.idnumber and currentResidue.chainID == Residue_nn.chainID and currentResidue.name == Residue_nn.name):
                continue
            else:
                distance = currentResidue.DistanceBetweenAA(Residue_nn.center)
                radiusSum = radiusDict[currentResidue.name] + radiusDict[Residue_nn.name] + 3
                if(distance <= radiusSum):
                    #print currentResidue.idnumber,currentResidue.chainID
                    #print Residue_nn.idnumber,Residue_nn.chainID
                    rho_w,theta_w,phi_w = currentResidue.ChangeCoordinate(Residue_nn.center)
                    #print rho_w,theta_w,phi_w
                    theta_w = min(int(math.floor(theta_w*20/np.pi)),19)
                    phi_w = min(int(math.floor(phi_w*10/np.pi) + 10),19)
                    rho_w2,theta_w2,phi_w2 = Residue_nn.ChangeCoordinate(currentResidue.center)
                    theta_w2 = min(int(math.floor(theta_w2*20/np.pi)),19)
                    phi_w2 = min(int(math.floor(phi_w2*10/np.pi) + 10),19)
                    Energy += matrix[currentResidue.name][Residue_nn.name][theta_w][phi_w] / rho_w
                    Energy += matrix[Residue_nn.name][currentResidue.name][theta_w2][phi_w2] / rho_w2

                    Energy400 += matrix400[currentResidue.name][Residue_nn.name]
                    Energy400 += matrix400[Residue_nn.name][currentResidue.name]

    return Energy,Energy400

def GetEnergy(wholepdbfile,loopfile):
    matrix = load_EnergyMatrix()
    matrix400 = load_EnergyMatrix400()
    wholepdbfile = ''
    loopfile = ''
    E,E400 = calculate_Energy(wholepdbfile,loopfile,matrix,matrix400)


if __name__ == "__main__":

    #args = sys.argv[1:]
    #pdb = args[0]
    matrix = load_EnergyMatrix()
    matrix400 = load_EnergyMatrix400()
    wholepdbfile = './testdata/6dek.pdb'
    loopfile = './testdata/6dekloop.pdb'
    E = calculate_Energy(wholepdbfile,loopfile,matrix,matrix400)
    print "Nepre Potential Energy(Radius)"
    print E
