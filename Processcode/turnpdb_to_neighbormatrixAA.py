import os
import math
import numpy as np
import AminoAcid as AA
#import gc
import sys
import math
#import csv
import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig
from multiprocessing import Pool
import multiprocessing

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


def ProcessPDB(file,matrix):
    df = open(file,'r')
    radiusDict = LoadRadius()
    CurrentAANitrogen = None
    CurrentAACA = None
    Currentresidue_num = None
    EachAA = []
    CurrentAA = None


    for line in df.readlines():
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
    df.close()
    return EachAA


def load_coordinate_number_matrix():
    aaDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
            "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
            "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
            "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}
    List = aaDict.keys()
    #print(type(List))
    #print(List)
    List=list(List)
    List.sort()
    #f1 = open("/Users/xg666/Desktop/gitgao/looppredict/Radius/radius.npy")
    #dualArray = np.zeros((20,20))
    for amino1 in List:
        for amino2 in List:
            aaDict[amino1][amino2] = 0
            #print aaDict[amino1]
    return aaDict


def judge_Neighbor(EachAAlist,resultaaDict,radiusDict):
    #radiusDict = LoadRadius()
    for m in range(len(EachAAlist)):
        EachAAlist[m].EstablishCoordinate()
        for n in range(len(EachAAlist)):
            if (m == n):
                continue
            dis = EachAAlist[m].DistanceBetweenAA(EachAAlist[n].center)
            radiusSum = radiusDict[EachAAlist[m].name] + radiusDict[EachAAlist[n].name] + 3
            if(dis <= radiusSum):
                #print (dis)
                #print(aaDict[chainA[m].name][chainB[n].name][theta][phi])
                #rho,theta,phi = EachAAlist[m].ChangeCoordinate(EachAAlist[n].center)
                #theta = min(int(math.floor(theta*20/np.pi)),19)
                #phi = min(int(math.floor(phi*10/np.pi) + 10),19)
                resultaaDict[EachAAlist[m].name][EachAAlist[n].name] += 1
    return resultaaDict


def main(resultpath,filepath,filename):
    file = os.path.join(filepath,filename)
    radiusDict = LoadRadius()
    aaDict = load_coordinate_number_matrix()
    EachAAlist = ProcessPDB(file,aaDict)
    #print EachAAlist
    resultDict = judge_Neighbor(EachAAlist,aaDict,radiusDict)
    #resultDict = static_Neighbor(changedchains,aaDict,radiusDict)
    savefile = os.path.join(resultpath,str(filename[:-4])+str('.npy'))
    #print savefile
    np.save(savefile,resultDict)
    #return resultDict

def multiprocess_run(filepath,resultpath):
    allnames = listdirInMac(filepath)
    #allnames= ['pdb1pft.ent']
    filelength = len(allnames)
    #for filename in filenames:
    processnum = 3
    #perSize = int(math.ceil(float(filelength)/float(processnum)))
    pool = multiprocessing.Pool(processes = processnum)
    for process_n in range(filelength):
        filename = allnames[process_n]
        print (filename)
        '''
        if perSize*(process_n + 1) > len(allnames):
            rightside = len(allnames)
            leftside = perSize*process_n
        else:
            leftside = perSize * process_n
            rightside = perSize * (process_n + 1)
        print leftside,rightside
        '''
        #filenames = allnames[leftside,rightside]
        #main(resultpath,filepath,filename)
        #get_startid_endid_inifile(filepath,reportpath,pdbpath)
        pool.apply_async(main,(resultpath,filepath,filename))
    pool.close()
    pool.join()

if __name__ == "__main__":
    #args = sys.argv[1:]
    #PDBname = args[0]
    #print (PDBname)
    #filepath = os.path.join('/public/home/xgao/Desktop/zdocktest/newTERzdockdata/',PDBname)
    #resultpath =os.path.join('/public/home/xgao/Desktop/zdocktest/zdocknpydata/',PDBname)
    filepath = '/Users/xg666/Desktop/loop/getTER/file7'
    resultpath = '/Users/xg666/Desktop/loop/getTER/file7npy'
    if not os.path.exists(resultpath):
        os.makedirs(resultpath)
    multiprocess_run(filepath,resultpath)

