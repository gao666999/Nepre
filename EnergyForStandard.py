import os
import math
import numpy as np
import AminoAcidv2 as AA
import gc
import sys
import math
import csv
import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig


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
def calculate_Energy(wholefile, matrix,matrix400):
    #Scan over. Each amino acid is stored as an object in EachAA. Next step is to calculate the energy, results will be saved in EnergyList.
    EachAA = process_pdb_file(wholefile,matrix)
    radiusDict = LoadRadius()
    #Store the energy
    E = 0
    E400 = 0
    for m in range(len(EachAA)):
        EachAA[m].EstablishCoordinate()
        for n in range(len(EachAA)):
            if ([EachAA[m].name,EachAA[m].idnumber,EachAA[m].chainID]==[EachAA[n].name,EachAA[n].idnumber,EachAA[n].chainID]):
            #if(m == n):
                continue
            else:
                dis = EachAA[m].DistanceBetweenAA(EachAA[n].center)
                radiusSum = radiusDict[EachAA[m].name] + radiusDict[EachAA[n].name]
                if(dis <= radiusSum):#If the distance between two amino acid less than 10, we believe the two amino acid have interaction
                    print [EachAA[m].name,EachAA[m].idnumber,EachAA[m].chainID]
                    print [EachAA[n].name,EachAA[n].idnumber,EachAA[n].chainID]
                    rho,theta,phi = EachAA[m].ChangeCoordinate(EachAA[n].center)
                    theta = min(int(math.floor(theta*20/np.pi)),19)
                    phi = min(int(math.floor(phi*10/np.pi) + 10),19)

                    E += matrix[EachAA[m].name][EachAA[n].name][theta][phi] / rho
                    E400 += matrix400[EachAA[m].name][EachAA[n].name]


    return E,E400


if __name__ == "__main__":

    #args = sys.argv[1:]
    #pdb = args[0]
    pdb = '/Users/xg666/Desktop/loop/getTER/file7/pdb1pft.ent'
    matrix = load_EnergyMatrix()
    matrix400 = load_EnergyMatrix400()
    E = calculate_Energy(pdb,matrix,matrix400)
    print "Nepre Potential Energy(Radius)"
    print pdb,E
