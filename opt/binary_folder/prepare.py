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
    
    return name,coord

## prepare inputs of QC calculation and barkla submission (S0opt)
def Prepare1():
    task = 'S0opt'

    if (input.QCFlag.lower() == 'g16'):
        g16.InitPara(task)
        g16.g16para['name'],g16.g16para['coord'] = FromSmiles(input.SystemSmiles)
        g16.GjfGen(g16.g16para)
    
    elif (input.QCFlag.lower() == 'orca'):
        orca.InitPara(task)
        orca.orcapara['name'],orca.orcapara['coord'] = FromSmiles(input.SystemSmiles)
        orca.InpGen(orca.orcapara)
    
    else:
        print('"QCFlag" now has only two options: g16 or orca')
        exit()

## prepare inputs of QC calculation and barkla submission (S0freq)
def Prepare2():
    task = 'S0freq'
    if (input.QCFlag.lower() == 'g16'):

        # S0 freq
        g16.InitPara(task)
        g16.g16para['name'], g16.g16para['coord'] = g16.FchkRead('{}/S0opt'.format(input.FilePath), 'coord')
        g16.g16para['coord'] *= input.au2Angstrom
        g16.GjfGen(g16.g16para)

    elif (input.QCFlag.lower() == 'orca'):

        # S0 freq
        orca.InitPara(task)
        orca.orcapara['name'],orca.orcapara['coord'] = orca.LogRead('{}/S0opt'.format(input.FilePath), 'coord')
        orca.orcapara['coord'] *= input.au2Angstrom
        orca.InpGen(orca.orcapara)

    else:
        print('"QCFlag" now has only two options: g16 or orca')
        exit()

# check if there is no imaginary frequency (ensure the correction of S0 opt task)
# return 0
def CheckFreq():
    if (input.QCFlag.lower() == 'g16'):
        AtomMass, ModeFreq, ModeQ, ModeVect = g16.FchkRead('{}/S0freq'.format(input.FilePath), 'vibration')
    elif (input.QCFlag.lower() == 'orca'):
        AtomMass, ModeFreq, ModeQ, ModeVect = orca.LogRead('{}/S0freq'.format(input.FilePath), 'vibration')
    else:
        print('"QCFlag" now has only two options: g16 or orca')
        exit()
    
    # normal termination of S0 opt
    if (ModeFreq[0] > 0.):
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
            g16.g16para['name'],g16.g16para['coord'] = g16.FchkRead('{}/S0freq'.format(input.FilePath), 'coord')
            g16.g16para['coord'] *= input.au2Angstrom
            g16.g16para['coord'] += disp
            g16.GjfGen(g16.g16para)

        else:
            orca.InitPara(task)
            orca.orcapara['name'],orca.orcapara['coord'] = orca.LogRead('{}/S0freq'.format(input.FilePath), 'coord')
            orca.orcapara['coord'] *= input.au2Angstrom
            orca.orcapara['coord'] += disp
            orca.InpGen(orca.orcapara)

        return -1

## prepare inputs of QC calculation and barkla submission (S1force)
def Prepare3():
    if (input.QCFlag.lower() == 'g16'):
        
        g16.g16para['name'], g16.g16para['coord'] = g16.FchkRead('{}/S0opt'.format(input.FilePath), 'coord')
        g16.g16para['coord'] *= input.au2Angstrom

        # S1 force
        task = 'S1force'
        g16.InitPara(task)
        g16.GjfGen(g16.g16para)

        
    elif (input.QCFlag.lower() == 'orca'):
                
        orca.orcapara['name'],orca.orcapara['coord'] = orca.LogRead('{}/S0opt'.format(input.FilePath), 'coord')
        orca.orcapara['coord'] *= input.au2Angstrom

        # S1 force
        task = 'S1force'
        orca.InitPara(task)
        orca.InpGen(orca.orcapara)

    else:
        print('"QCFlag" now has only two options: g16 or orca')
        exit()

## prepare inputs of QC calculation and barkla submission (S1opt)
def Prepare4():
    if (input.QCFlag.lower() == 'g16'):
        
        g16.g16para['name'], g16.g16para['coord'] = g16.FchkRead('{}/S0opt'.format(input.FilePath), 'coord')
        g16.g16para['coord'] *= input.au2Angstrom

        # S1 opt
        task = 'S1opt'
        g16.InitPara(task)
        g16.GjfGen(g16.g16para)
        
    elif (input.QCFlag.lower() == 'orca'):
                
        orca.orcapara['name'],orca.orcapara['coord'] = orca.LogRead('{}/S0opt'.format(input.FilePath), 'coord')
        orca.orcapara['coord'] *= input.au2Angstrom

        # S1 opt
        task = 'S1opt'
        orca.InitPara(task)
        orca.InpGen(orca.orcapara)

    else:
        print('"QCFlag" now has only two options: g16 or orca')
        exit()