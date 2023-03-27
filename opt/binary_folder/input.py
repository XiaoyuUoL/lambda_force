import yaml

# constants
c = 2.99792458e-5                                        # speed of light (unit: cm/fs)

# atomic unit
hbar = 1.0
au2Angstrom = 0.52917721090380                           # length
au2eV = 27.21138624598853                                # energy
au2fs = 2.418884326585747e-2                             # time
au2aum = 5.48579909065e-4                                # mass

# read information from calculator.yml
CalInput = open('calculator.yml', 'rt')
CalDict = yaml.safe_load(CalInput)
CalInput.close()

# parameters for conversion from smiles to xyz
FF = CalDict['specifications']['MMFF']                   # forcefield for mol.localopt()
StepsInit =  20                                          # steps for initial conformation
StepsOpt =  1000                                         # steps for mol.localopt()
#print('FF:', FF)

# QC package ('g16' or 'orca')
QCFlag = CalDict['specifications']['Software']
if (QCFlag.lower() != 'g16' and QCFlag.lower() != 'orca'):
    print('"QCFlag" now has only two options: g16 or orca')
    exit()
#print('QCFlag:', QCFlag)

# computation/calculation details
Memory = CalDict['specifications']['Memory']             # MB
ProcNumber = CalDict['specifications']['ProcNumber']
Functional = CalDict['specifications']['Functional']
BasisSet = CalDict['specifications']['BasisSet']
#print('Memory:', Memory)
#print('ProcNumber:', ProcNumber)
#print('Functional:', Functional)
#print('BasisSet:', BasisSet)

# dushin parameters
dq = 1e-5                                                # dq for numerical Wilson Matrix
de = 1e-5                                                # threshold for discarding eigenvectors of GMatrix

# Freq check for S0opt
ImagFreq = CalDict['specifications']['ImagFreq']
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
