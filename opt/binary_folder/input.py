import yaml

# constants
c = 2.99792458e-5                                                               # speed of light (unit: cm/fs)

# atomic unit
hbar = 1.0
au2Angstrom = 0.52917721090380                                                  # length
au2eV = 27.21138624598853                                                       # energy
au2fs = 2.418884326585747e-2                                                    # time
au2aum = 5.48579909065e-4                                                       # mass

# read information from calculator.yml
CalInput = open('calculator.yml', 'rt')
CalDict = yaml.safe_load(CalInput)
CalInput.close()

# parameters for conversion from smiles to xyz
FF = CalDict['specification']['MMFF']                                           # forcefield for mol.localopt()
StepsInit =  20                                                                 # steps for initial conformation
StepsOpt =  1000                                                                # steps for mol.localopt()
#print('FF:', FF)

# QC package ('g16' or 'orca')
QCFlag = CalDict['specification']['Software']
if (QCFlag.lower() != 'g16' and QCFlag.lower() != 'orca'):
    print('"QCFlag" now has only two options: g16 or orca')
    exit()
#print('QCFlag:', QCFlag)

# computation/calculation details
Memory = CalDict['specification']['Memory']                                     # unit: MB
ProcNumber = CalDict['specification']['ProcNumber']
Functional = CalDict['specification']['Functional']
BasisSet = CalDict['specification']['BasisSet']
NRoots = CalDict['specification']['NRoots']
IRoot = CalDict['specification']['IRoot']
if (IRoot > NRoots):
    print('error in calculator.yml: IRoot is larger than NRoots')
    exit()
#print('Memory:', Memory)
#print('ProcNumber:', ProcNumber)
#print('Functional:', Functional)
#print('BasisSet:', BasisSet)
#print('NRoots:', NRoots)
#print('IRoot:', IRoot)

# basis function of elements for different basis set
BasisFunc = {}

# dushin parameters
dq = 1e-5                                                                       # dq for numerical Wilson Matrix
de = 1e-5                                                                       # threshold for discarding eigenvectors of GMatrix

# Freq check for S0opt
ImagFreq = CalDict['specification']['ImagFreq']
#print('ImagFreq:', ImagFreq)

# Calculation properties
Properties = CalDict['provides']
#print('Properties', Properties)

# read information from molecule.yml
MolInput = open('molecule.yml', 'rt')
MolDict = yaml.safe_load(MolInput)
MolInput.close()

# system & input/output
SystemName = MolDict['id']
SystemSmiles = MolDict['smiles']