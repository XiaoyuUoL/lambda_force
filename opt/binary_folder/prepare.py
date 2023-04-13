import numpy as np
from openbabel import openbabel
from openbabel import pybel

import input
import g16
import orca

## get Cartesian coordinates from smiles
def FromSmiles(smiles):
    mol = pybel.readstring('smi', smiles)
    mol.make3D(forcefield=input.FF, steps=input.StepsInit)

    ff = openbabel.OBForceField.FindForceField(input.FF)
    ff.Setup(mol.OBMol)
    ff.FastRotorSearch()
    ff.SteepestDescent(input.StepsOpt)
    ff.GetCoordinates(mol.OBMol)
    
    name = []
    coord = []
    for atom in mol.atoms:
        name.append(atom.atomicnum)
        coord.append(list(atom.coords))
    
    return name,coord,mol.charge,mol.spin

## prepare inputs of QC calculation
# name is for S0sp calculation of single atom
def QCinput(task, name=None):
    if (input.QCFlag.lower() == 'g16'):
        g16.InitPara(task)
        if (task == 'S0opt'):
            g16.g16para['name'],g16.g16para['coord'],g16.g16para['charge'],g16.g16para['multi'] = FromSmiles(input.SystemSmiles)
        elif (task in ['S0freq', 'S1force', 'S1opt']):
            g16.g16para['name'], g16.g16para['coord'] = g16.FchkRead('S0opt', 'coord')
            g16.g16para['coord'] *= input.au2Angstrom
        elif (task == 'S0sp' and name != None):
            g16.g16para['name'] = [name]
            g16.g16para['coord'] = [[0.0, 0.0, 0.0]]
            if (name % 2):
                g16.g16para['multi'] = 2
            else:
                g16.g16para['multi'] = 1
        else:
            print('error: task is not valid for prepare.QCinput')
            exit()
        g16.GjfGen(g16.g16para)
    
    elif (input.QCFlag.lower() == 'orca'):
        orca.InitPara(task)
        if (task == 'S0opt'):
            orca.orcapara['name'],orca.orcapara['coord'],orca.orcapara['charge'],orca.orcapara['multi'] = FromSmiles(input.SystemSmiles)
            if ('BOD' in input.Properties):
                orcapara['keywords'].append('printMOs')
        elif (task in ['S0freq', 'S1force', 'S1opt']):
            orca.orcapara['name'],orca.orcapara['coord'] = orca.LogRead('S0opt', 'coord')
            orca.orcapara['coord'] *= input.au2Angstrom
        elif (task == 'S0sp' and name != None):
            orca.orcapara['name'] = [name]
            orca.orcapara['coord'] = [[0.0, 0.0, 0.0]]
            print(name, name % 2)
            if (name % 2):
                orca.orcapara['multi'] = 2
            else:
                orca.orcapara['multi'] = 1
        else:
            print('error: task is not valid for QCinput')
            exit()
        orca.InpGen(orca.orcapara)
    
    else:
        print('"QCFlag" now has only two options: g16 or orca')
        exit()
    
# check if there is no imaginary frequency (ensure the correction of S0opt task)
def CheckFreq():
    if (input.QCFlag.lower() == 'g16'):
        AtomMass, ModeFreq, ModeQ, ModeVect = g16.FchkRead('S0freq', 'vibration')
    elif (input.QCFlag.lower() == 'orca'):
        AtomMass, ModeFreq, ModeQ, ModeVect = orca.LogRead('S0freq', 'vibration')
    else:
        print('"QCFlag" now has only two options: g16 or orca')
        exit()
    
    # normal termination of S0opt
    if (ModeFreq[0] > input.ImagFreq):
        return 0
    # rerun S0 opt with distorted initial geometry
    else:
        # make a distortion based on imaginary modes
        disp = np.zeros_like(np.reshape(ModeVect[0, :], (-1, 3)))
        for freq,vect in zip(ModeFreq,ModeVect):
            if (freq < 0.):
                vect /= np.sqrt(AtomMass)
                vect /= np.sqrt(np.sum(vect * vect))
                disp += np.reshape(vect, (-1, 3)) * np.random.rand()
        disp /= np.sqrt(np.sum(disp * disp))    # |disp| = 0.5

        task = 'S0opt'
        if (input.QCFlag.lower() == 'g16'):
            g16.InitPara(task)
            g16.g16para['name'],g16.g16para['coord'] = g16.FchkRead('S0freq', 'coord')
            g16.g16para['coord'] *= input.au2Angstrom
            g16.g16para['coord'] += disp
            g16.GjfGen(g16.g16para)

        else:
            orca.InitPara(task)
            orca.orcapara['name'],orca.orcapara['coord'] = orca.LogRead('S0freq', 'coord')
            orca.orcapara['coord'] *= input.au2Angstrom
            orca.orcapara['coord'] += disp
            orca.InpGen(orca.orcapara)

        return -1