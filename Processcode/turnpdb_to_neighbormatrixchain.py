# -*- coding: utf-8 -*-
import os
import math
import numpy as np
import AminoAcid as AA
import gc
import sys
import string
from multiprocessing import Pool
import multiprocessing
import numpy as np 
np.seterr(divide='ignore', invalid='ignore')

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

    f = open("/public/home/xgao/Desktop/gitgao/looppredict/Radius/mean_radius.txt")
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



def processAAforchian(chain,aaDict):
    CurrentAANitrogen = None
    CurrentAACA = None
    Currentresidue_num = None
    EachAA = []
    CurrentAA = None
    for line in chain:
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
        if (residue_name not in aaDict):
            continue
        if (CurrentAA == None):
            CurrentAA = AA.AminoAcid(residue_name, residue_num, chain_id)
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
                    CurrentAA = AA.AminoAcid(residue_name, residue_num, chain_id)
                    Currentresidue_num = residue_num
                    continue

                CurrentAA.InputCAN(CurrentAANitrogen, CurrentAACA)
                CurrentAA.EstablishCoordinate()
                # Amino Acid check
                EachAA.append(CurrentAA)
                del CurrentAA
                CurrentAA = AA.AminoAcid(residue_name, residue_num, chain_id)

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
# process every chain,to make erery AA has a
def change_chain(chains,aaDict):
    changedchains = []
    #words = string.uppercase
    for g in range(len(chains)):
        #chainname = 'chain' + str(words[g])
        #print name
        chain_processed = processAAforchian(chains[g],aaDict)
        changedchains.append(chain_processed)
    return changedchains

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
            #print aaDict[amino1]#[amino2]
            #aaDict[amino1][amino2] = np.load(f1)
            aaDict[amino1][amino2] = np.zeros((20,20))
            #print aaDict[amino1]
    return aaDict


def judge_Neighbor(chainA,chainB,aaDict,radiusDict):
    #radiusDict = LoadRadius()
    for m in range(len(chainA)):
        chainA[m].EstablishCoordinate()
        for n in range(len(chainB)):
            dis = chainA[m].DistanceBetweenAA(chainB[n].center)
            radiusSum = radiusDict[chainA[m].name] + radiusDict[chainB[n].name] + 3
            if(dis <= radiusSum):
                print (dis)
                #print(aaDict[chainA[m].name][chainB[n].name][theta][phi]) 
                rho,theta,phi = chainA[m].ChangeCoordinate(chainB[n].center)
                theta = min(int(math.floor(theta*20/np.pi)),19)
                phi = min(int(math.floor(phi*10/np.pi) + 10),19)
                aaDict[chainA[m].name][chainB[n].name][theta][phi] += 1
    return aaDict


def static_Neighbor(changedchains,aaDict,radiusDict):
    length = len(changedchains)
    for i in range(length):
        for j in range(length):
            if i != j:
                chainA = changedchains[i]
                chainB = changedchains[j]
                aaDict = judge_Neighbor(chainA,chainB,aaDict,radiusDict)

                #judge_Neighbor(changedchains[i],changedchains[j],aaDict)
    return aaDict

def getlines_for_eachchain(file):
    chains = []
    chainA = []
    with open(file,'r') as file:
        for line in file.readlines():
            if (line[0:4] == 'ATOM'):
                chainA.append(line)
            elif (line[:3].strip() == 'TER'):
                print (line)
                chainA.append(line)
                chains.append(chainA)
                chainA = []
    return chains

#def separate_files(path):


def main(resultpath,filepath,filename):
    file = os.path.join(filepath,filename)
    radiusDict = LoadRadius()
    aaDict = load_coordinate_number_matrix()
    chains = getlines_for_eachchain(file)
    print (len(chains))
    changedchains = change_chain(chains,aaDict)
    resultDict = static_Neighbor(changedchains,aaDict,radiusDict)
    savefile = os.path.join(resultpath,str(filename[:-4])+str('.npy'))
    np.save(savefile,resultDict)
    #return resultDict

def multiprocess_run(filepath,resultpath):
    allnames = listdirInMac(filepath)[80000:]
    filelength = len(allnames)
    #for filename in filenames:
    processnum = 24
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
    filepath = '/public/home/xgao/Desktop/CreateZdock/NegativeSample/'
    resultpath = '/public/home/xgao/Desktop/CreateZdock/Nnpy/'
    if not os.path.exists(resultpath):
        os.makedirs(resultpath)
    multiprocess_run(filepath,resultpath)
