import numpy as np

import input

para = {
    'file': '',
    'procnumber': input.ProcNumber,
    'keywords': [],
    'desc': '',
    'charge': 0,
    'multi': 1,
    'name': [],
    'coord': [],
    'module': '',
    'option': {}
}

# common keywords/module/options for different calculation task
# can be changed
keywords = {
    'S0sp': [input.Functional, input.BasisSet, 'tightSCF'],
    'S0opt': [input.Functional, input.BasisSet, 'opt', 'tightSCF'],
    'S0freq': [input.Functional, input.BasisSet, 'AnFreq', 'tightSCF'],
    'S1force': [input.Functional, input.BasisSet, 'engrad', 'tightSCF'],
    'S1opt': [input.Functional, input.BasisSet, 'opt', 'tightSCF'],
    'S1nac': [input.Functional, input.BasisSet, 'tightSCF'],
    'soc': [input.Functional, input.BasisSet, 'tightSCF']
}

modules = {
    'S0sp': '',
    'S0opt': '',
    'S0freq': '',
    'S1force': 'tddft',
    'S1opt': 'tddft',
    'S1nac': 'tddft',
    'soc': 'tddft'
}

options = {
    'S0sp': {},
    'S0opt': {},
    'S0freq': {},
    'S1force': {
        'tda': 'false',
        'nroots': '5',
        'iroot': '1',
        'printlevel': '4'
    },
    'S1opt': {
        'tda': 'false',
        'nroots': '5',
        'iroot': '1',
        'printlevel': '4'
    },
    'S1nac': {
        'tda': 'false',
        'nroots': '5',
        'iroot': '1',
        'nacme': 'true',
        'printlevel': '4'
    },
    'soc': {
        'tda': 'false',
        'nroots': '5',
        'dosoc': 'true',
        'printlevel': '4'
    }
}

def InitPara(task):
    para['desc'] = task
    para['keywords'] = keywords[task]
    para['module'] = modules[task]
    para['option'] = options[task]
    para['file'] = task

# generate input file for orca calculation
def InpGen(para):
    fout = open('{}.inp'.format(para['file']), 'w')
    fout.writelines('%pal\n')
    fout.writelines('  nprocs {:d}\n'.format(para['procnumber']))
    fout.writelines('end\n')
    keywords = '!'
    for keyword in para['keywords']:
        keywords += ' ' + keyword
    fout.writelines('{}\n'.format(keywords))
    if (para['module'] != ''):
        fout.writelines('%tddft\n')
        for key,value in para['option'].items():
            fout.writelines('  {} {}\n'.format(key, value))
        fout.writelines('end\n')
    fout.writelines('*xyz {:d} {:d}\n'.format(para['charge'], para['multi']))
    for name, coord in zip(para['name'], para['coord']):
        fout.writelines('{}{:14.7f}{:14.7f}{:14.7f}\n'.format(name, coord[0], coord[1], coord[2]))
    fout.writelines('*\n')
    fout.close()

## read log file. keyword options:
# 'energy-gd': return energy of electronic ground state
# 'energy-ex': return energies of electronic excited state
# 'coord': return atom name and atom coordinate
# 'gradient': return energy gradient
# 'vibration': return normal mode information
# 'zindices': return indice of zmatrix
# 'basis': return number of basis functions
# 'orbital': return AO/MO information
# 'nac': return nonadiabatic coupling between S0 and S1
# 'soc': return spin-orbit couplings among singlet and triplet
def LogRead(logfile, keyword):
    fin = open('{}.log'.format(logfile), 'r')

    # 'energy': return energy of electronic ground state
    if (keyword.lower() == 'energy-gd'):
        line = fin.readline()
        e = None
        while (len(line) != 0):
            words = line.rstrip().split()
            # read ground state energy
            if (len(words) > 2 and words[0] == 'Total' and words[1] == 'Energy'):
                e = float(words[3])
            line = fin.readline()
        fin.close()

        if (e != None):
            return e  # unit: Hartree

        print('error of orca.LogRead(): No ground state energy information in {}.log'.format(logfile))
        exit()

    # 'energy': return energies of electronic excited state
    if (keyword.lower() == 'energy-ex'):
        e = []
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            # read excited state energies
            if (len(words) > 2 and words[0] == 'TD-DFT' and words[1] == 'EXCITED'):
                line = fin.readline()
                while (len(line) != 0):
                    words = line.rstrip().split()
                    if (len(words) > 1 and words[0] == 'STATE'):
                        e.append(float(words[3]))
                    elif (len(words) > 1 and words[0] == 'TD-DFT-EXCITATION'):
                        break
                    line = fin.readline()
            line = fin.readline()
        fin.close()

        if (len(e) != 0):
            return e   # unit: Hartree

        print('error of orca.LogRead(): No excited state energy information in {}.log'.format(logfile))
        exit()

    # 'coord': return atom name and atom coordinate
    elif (keyword.lower() == 'coord'):
        line = fin.readline()
        Name = []
        Coord = []
        while (len(line) != 0):
            words = line.rstrip().split()
            # read atomic number and Cartesian coordinates
            if (len(words) == 3 and words[0] == 'CARTESIAN' and words[1] == 'COORDINATES' and words[2] == '(A.U.)'):
                Name = []
                Coord = []
                fin.readline()
                fin.readline()
                line = fin.readline()
                while (line != '\n'):
                    data = line.rstrip().split()
                    Name.append(data[1])
                    Coord.append([data[5], data[6], data[7]])  # unit: Bohr
                    line = fin.readline()
            line = fin.readline()
        fin.close()

        if (len(Name) * len(Coord) != 0):
            return Name, np.array(Coord, dtype=float)

        print('error of orca.LogRead(): No coordinate information in {}.log'.format(logfile))
        exit()

    # 'gradient': return energy gradient
    elif (keyword.lower() == 'gradient'):
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            # read gradient
            if (len(words) == 2 and words[0] == 'CARTESIAN' and words[1] == 'GRADIENT'):
                Gradient = []
                fin.readline()
                fin.readline()
                line = fin.readline()
                while (line != '\n'):
                    data = line.rstrip().split()
                    Gradient += [data[3], data[4], data[5]]
                    line = fin.readline()
                fin.close()
                return np.array(Gradient, dtype=float)  # unit: Hartree/Bohr (3N 1D-array)
            line = fin.readline()
        fin.close()

        print('error of orca.LogRead(): No gradient information in {}'.format(logfile))
        exit()

    # 'vibration': return normal mode information
    elif (keyword.lower() == 'vibration'):
        IfGetAtom = False
        AtomNumber = 0
        AtomMass = []
        IfGetVib = False
        ModeNumber = 0
        ModeFreq = []
        ZeroIndex = []
        IfGetMode = False
        ModeVectTmp = np.array([], dtype=float)
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            # read atom number and atomic mass
            if (len(words) == 3 and words[0] == 'CARTESIAN' and words[1] == 'COORDINATES' and words[2] == '(A.U.)'):
                IfGetAtom = True
                fin.readline()
                fin.readline()
                line = fin.readline()
                while (line != '\n'):
                    AtomNumber += 1
                    data = line.rstrip().split()
                    AtomMass+= [data[4], data[4], data[4]]
                    line = fin.readline()
            # read frequency
            elif (len(words) == 2 and words[0] == 'VIBRATIONAL' and words[1] == 'FREQUENCIES' and IfGetAtom):
                IfGetVib = True
                fin.readline()
                fin.readline()
                fin.readline()
                fin.readline()
                for i in np.arange(AtomNumber * 3):
                    data = fin.readline().rstrip().split()
                    freq = float(data[1])
                    if (freq != 0):
                        ModeFreq.append(freq)
                        ModeNumber += 1
                    else:
                        ZeroIndex.append(i)
            # read vibration vector
            elif (len(words) == 2 and words[0] == 'NORMAL' and words[1] == 'MODES' and IfGetVib):
                IfGetMode = True
                ModeVectTmp = np.zeros((AtomNumber * 3, AtomNumber * 3), dtype=float)
                fin.readline()
                fin.readline()
                fin.readline()
                fin.readline()
                fin.readline()
                fin.readline()
                for i in np.arange((AtomNumber * 3 - 1)//6 + 1):
                    index = np.array(fin.readline().rstrip().split(), dtype=int)
                    for j in np.arange(AtomNumber * 3):
                        ModeVectTmp[index, j] = np.array(fin.readline().rstrip().split()[1:], dtype=float)
            line = fin.readline()
        fin.close()

        if (not IfGetMode):
            print('error of orca.LogRead(): No normal mode information in {}.log'.format(logfile))
            exit()

        AtomMass = np.array(AtomMass, dtype=float) # unit: amu
        ModeFreq = np.array(ModeFreq)
        ModeFreq *= 2 * np.pi * input.c * input.au2fs  # unit: Hartree

        ModeVect = np.zeros((ModeNumber, AtomNumber * 3), dtype=float)
        index = 0
        for i in np.arange(AtomNumber * 3):
            if i not in ZeroIndex:
                ModeVect[index, :] = ModeVectTmp[i, :]
                index += 1
        for i in np.arange(AtomNumber * 3):
            ModeVect[:, i] *= np.sqrt(AtomMass[i])
        for i in np.arange(ModeNumber):
                ModeVect[i, :] /= np.sqrt(np.sum(ModeVect[i, :] * ModeVect[i, :]))

        ModeQ = np.zeros_like(ModeFreq)
        for i,freq in enumerate(ModeFreq):
            if (ModeFreq[0] > 0.):
                ModeQ[i] = np.sqrt(freq * (1. / input.au2aum) / input.hbar)

        return AtomMass, ModeFreq, ModeQ, ModeVect

    # return redundant internal coordinates
    elif (keyword.lower() == 'zindices'):
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            if (len(words) == 3 and words[0] == 'Redundant' and words[1] == 'Internal' and words[2] == 'Coordinates'):
                ZIndices = []
                fin.readline()
                fin.readline()
                fin.readline()
                fin.readline()
                fin.readline()
                line = fin.readline()
                while ('---' not in line):
                    data = line.rstrip().split()
                    if (data[1][0] == 'B'):
                        ZIndices.append([int(data[2][:-2]), int(data[3][:-1]), -1, -1])
                    elif (data[1][0] == 'A'):
                        ZIndices.append([int(data[2][:-2]), int(data[3][:-2]), int(data[4][:-1]), -1])
                    else:
                        ZIndices.append([int(data[2][:-2]), int(data[3][:-2]), int(data[4][:-2]), int(data[5][:-1])])
                    line = fin.readline()
                fin.close()
                return ZIndices
            line = fin.readline()
        fin.close()

        print('error of orca.LogRead(): No redundant internal coordinates information in {}.log'.format(logfile))
        exit()

    # return number of basis functions
    elif (keyword.lower() == 'basis'):
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            if(len(words) > 4 and words[0] == 'Number' and words[1] == 'of' and words[2] == 'basis' and words[3] == 'functions'):
                return int(words[-1])
            line = fin.readline()
        fin.close()

        print('error of orca.LogRead(): No basis functions number information in {}.log'.format(logfile))
        exit()

    # return orbital information (AO, MO) ##### only close shell now
    elif (keyword.lower() == 'orbital'):
        def GetBasisFunc(elements):
            BasisFunc = []
            for element in elements:
                BasisFunc.append(input.BasisFunc[element])
            
            return np.array(BasisFunc)

        # calculate the power of symmetry matrix
        def MatrixPower(M, x):
            e,v = np.linalg.eigh(M)
            return np.matmul(v, np.matmul(np.diag(np.power(e, x)), v.T))

        IfGetON = False
        OccNumber = None
        IfGetEle = False
        Elements = []
        BasisFunc = []
        BasisNumber = None
        IfGetMO = False
        MOCoeff = []
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            if(len(words) > 2 and words[0] == 'Number' and words[1] == 'of' and words[2] == 'Electrons' and not IfGetON):
                IfGetON = True
                OccNumber = int(words[-1]) // 2
            elif(len(words) > 2 and words[0] == 'CARTESIAN' and words[1] == 'COORDINATES' and words[2] == '(A.U.)' and not IfGetEle):
                IfGetEle = True
                fin.readline()
                fin.readline()
                line = fin.readline()
                while (line != '\n'):
                    data = line.rstrip().split()
                    Elements += str(int(float(data[2])))
                    line = fin.readline()
                BasisFunc = GetBasisFunc(Elements)
                BasisNumber = np.sum(BasisFunc)
            elif(len(words) == 2 and words[0] == 'MOLECULAR' and words[1] == 'ORBITALS'):
                IfGetMO = True
                MOCoeff = np.zeros((BasisNumber, BasisNumber), dtype=float)
                fin.readline()
                for i in np.arange((BasisNumber - 1)//6 + 1):
                    index = np.array(fin.readline().rstrip().split(), dtype=int)
                    fin.readline()
                    fin.readline()
                    fin.readline()
                    for j in np.arange(BasisNumber):
                        MOCoeff[j, index] = np.array(fin.readline().rstrip().split()[2:], dtype=float)
            line = fin.readline()
        fin.close()

        if (not IfGetMO):
            print('error of orca.LogRead(): No orbital information in {}.log'.format(logfile))
            exit()

        AOoverlap = MatrixPower(np.matmul(MOCoeff, MOCoeff.T), -1.0)

        # indices of basis function for each atom
        indices = np.array([[0] + list(np.cumsum(BasisFunc)[:-1]), list(np.cumsum(BasisFunc))]).T

        return indices, AOoverlap, MOCoeff[:, :OccNumber], MOCoeff[:, OccNumber:]

    # return nonadiabatic coupling between S0 and S1
    elif (keyword.lower() == 'nac'):
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            # read nonadiabatic coupling
            if (len(words) == 2 and words[0] == 'CARTESIAN' and words[1] == 'NON-ADIABATIC' and words[2] == 'COUPLINGS'):
                NAC = []
                fin.readline()
                fin.readline()
                fin.readline()
                line = fin.readline()
                while (line != '\n'):
                    data = line.rstrip().split()
                    NAC += [data[3], data[4], data[5]]
                    line = fin.readline()
                fin.close()
                return np.array(NAC, dtype=float)  # unit: 1/Bohr (3N 1D-array)
            line = fin.readline()
        fin.close()

        print('error of orca.LogRead(): No nonadiabatic coupling information in {}'.format(logfile))
        exit()

    # return spin-orbit couplings among singlet and triplet
    elif (keyword.lower() == 'soc'):
        RootNumber = None
        IfGetRN = False
        HSOC = []
        IfGetHSR = False
        IfGetHSI = False
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            # read number of roots
            if(len(words) > 6 and words[0] == 'Number' and words[1] == 'of' and words[2] == 'roots'):
                IfGetRN = True
                RootNumber = int(words[-1])
                HSOC = np.zeros((RootNumber * 4 + 1, RootNumber * 4 + 1), dtype = complex)
            # read real part of HSOC
            elif(len(words) == 2 and words[0] == 'Real' and words[1] == 'part:' and IfGetRN):
                IfGetHSR = True
                for i in np.arange((RootNumber * 4) // 6 + 1):
                    fin.readline()
                    for j in np.arange(RootNumber * 4 + 1):
                        line = fin.readline()
                        data = np.array(line.rstrip().split(), dtype = float)
                        HSOC[j, i * 6 : min(i * 6 + 6, RootNumber * 4 + 1)] = data[1:]
            # read image part of HSOC
            elif(len(words) == 2 and words[0] == 'Image' and words[1] == 'part:' and IfGetRN):
                IfGetHSI = True
                for i in np.arange((RootNumber * 4) // 6 + 1):
                    fin.readline()
                    for j in np.arange(RootNumber * 4 + 1):
                        line = fin.readline()
                        data = np.array(line.rstrip().split(), dtype = float)
                        HSOC[j, i * 6 : min(i * 6 + 6, RootNumber * 4 + 1)] += data[1:] * complex(0.0, 1.0)
            line = fin.readline()
        fin.close()

        if (not IfGetHSR or not IfGetHSI):
            print('error of orca.LogRead(): No nonadiabatic coupling information in {}.log'.format(logfile))
            exit()

        SOC = []
        for i in np.arange(RootNumber):
            for j in np.arange(RootNumber):
                # Sindex, Tindex, SOC(MS=0), SOC(MS=-1), SOC(MS=+1) (unit: Hartree)
                SOC.append([i, j + 1, HSOC[j + RootNumber + 1, i], HSOC[j + RootNumber * 2 + 1, i], HSOC[j + RootNumber * 3 + 1, i]])
        
        return SOC

    else:
        print('error of orca.LogRead(): option is needed')
        exit()