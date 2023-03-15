# constants
c = 2.99792458e-5                  # speed of light (unit: cm/fs)

# atomic unit
hbar = 1.0
au2Angstrom = 0.52917721090380     # length
au2eV = 27.21138624598853          # energy
au2fs = 2.418884326585747e-2       # time
au2aum = 5.48579909065e-4          # mass

# QC package
#QCFlag = 'g16' # or 'orca'
QCFlag = 'orca'

# computation/calculation details
memeory = 10000 # MB
ProcNumber = 4
if (QCFlag.lower() == 'g16'):
    Functional = 'm062x'
    BasisSet = '3-21g*'
elif (QCFlag.lower() == 'orca'):
    Functional = 'm06-2x'
    BasisSet = '3-21g*'
else:
    print('"QCFlag" now has only two options: g16 or orca')
    exit()

# basis function
#BasisFunc = [[1, 2, 2], [3, 10, 9], [11, 18, 19], [19, 36, 29], [37, 54, 39]]  # 3-21g* ([start, end, number])
#BasisFunc = [[1, 2, 2], [3, 10, 15], [11, 18, 19]] # 6-31g(d) 
#BasisFunc = [[1, 2, 5], [3, 10, 14], [11, 18, 18], [19, 36, 32]] # def2-SVP
#BasisFunc = [[1, 2, 6], [3, 10, 31], [11, 18, 37]]  # def2-TZVP ([start, end, number])

# system & input/output
SystemName = ''
SystemSmiles = ''
FilePath = ''

# parameters for conversion from smiles to xyz
FF = 'mmff94'                      # forcefield for mol.localopt()
StepsInit =  20                    # steps for initial conformation
StepsOpt =  1000                   # steps for mol.localopt()

# dushin parameters
dq = 1e-5                          # dq for numerical Wilson Matrix
de = 1e-5                          # threshold for discarding eigenvectors of GMatrix