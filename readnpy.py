import numpy as np
import math
def LoadMatrix():
    aaDict={"ALA":{},"VAL":{},"LEU":{},"ILE":{},"PHE":{},\
            "TRP":{},"MET":{},"PRO":{},"GLY":{},"SER":{},\
            "THR":{},"CYS":{},"TYR":{},"ASN":{},"GLN":{},\
            "HIS":{},"LYS":{},"ARG":{},"ASP":{},"GLU":{},}
    List = aaDict.keys()
    List = list(List)
    List.sort()
    for amino1 in List:
        for amino2 in List:
            aaDict[amino1][amino2] = 0
    return aaDict,List

def LoadCoordinateNumber(filec):
    coordict = np.load(filec).item()
    #numdict = np.load(filen).item()
    return coordict
def main(filec):
    aaDict,List = LoadMatrix()
    coordict = LoadCoordinateNumber(filec)
    for amino1 in List:
        for amino2 in List:
            print coordict[amino1][amino2]


if __name__ == "__main__":
    filec = '/Users/xg666/Desktop/loop/getTER/file7npy/pdb1pft.npy'
    main(filec)