# input/output of orca
import input
import numpy as np

orcapara = {
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
    'S0opt': [input.Functional, input.BasisSet, 'opt', 'tightSCF'],
    'S0freq': [input.Functional, input.BasisSet, 'AnFreq', 'tightSCF'],
    'S1force': [input.Functional, input.BasisSet, 'engrad', 'tightSCF'],
    'S1opt': [input.Functional, input.BasisSet, 'opt', 'tightSCF']
}

modules = {
    'S0opt': '',
    'S0freq': '',
    'S1force': 'tddft',
    'S1opt': 'tddft'
}

options = {
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
    }
}

def InitPara(task):
    orcapara['desc'] = task
    orcapara['keywords'] = keywords[task]
    orcapara['module'] = modules[task]
    orcapara['option'] = options[task]
    orcapara['file'] = task

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

# read information from orca log files
# keyword options:
# 'energy-gd': return energy of ground state
# 'energy-ex': return energy of excited states
# 'coord': return atom name and atom coordinate
# 'gradient': return energy gradient
# 'vibration': return normal mode information
# 'zindices': return indice of zmatrix
def LogRead(logfile, keyword):
    fin = open('{}.log'.format(logfile), 'r')
    
    # 'energy': return energy of electronic ground state
    if (keyword.lower() == 'energy-gd'):
        line = fin.readline()
        e = 0.0
        while (len(line) != 0):
            words = line.rstrip().split()
            # read ground state energy
            if (len(words) > 2 and words[0] == 'Total' and words[1] == 'Energy'):
                e = float(words[3])
            line = fin.readline()
        fin.close()

        if (e != 0.0):
            return e  # unit: Hartree
        
        print('orca: No ground state energy information in {}.log'.format(logfile))
        exit()

    # 'energy': return energies of electronic excited state
    if (keyword.lower() == 'energy-ex'):
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            # read excited state energy
            if (len(words) > 2 and words[0] == 'TD-DFT' and words[1] == 'EXCITED'):
                e = []
                line = fin.readline()
                while (len(line) != 0):
                    words = line.rstrip().split()
                    # read total energy
                    if (len(words) > 1 and words[0] == 'STATE'):
                        e.append(float(words[3]))
                    elif (len(words) > 1 and words[0] == 'TD-DFT-EXCITATION'):
                        break
                    line = fin.readline()
            line = fin.readline()
        fin.close()

        if (len(e) != 0):
            return e   # unit: Hartree

        print('orca: No ground state energy information in {}.log'.format(logfile))
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

        if (len(Name) != 0):
            return Name, np.array(Coord, dtype=float)
        
        print('orca: No coordinate information in {}.log'.format(logfile))
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
                return np.array(Gradient, dtype=float)  # unit: Hartree/Bohr
            line = fin.readline()
        fin.close()
        
        print('orca: No gradient information in {}'.format(logfile))
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
            print('orca: No vibration information in {}.log'.format(logfile))
            exit()

        AtomMass = np.array(AtomMass, dtype=float) # unit: amu
        ModeFreq = np.array(ModeFreq)
        ModeFreq *= 2 * np.pi * input.c * input.au2fs  # unit: Hartree

        ModeVect = np.zeros((ModeNumber, AtomNumber * 3), dtype=float)
        index = 0
        for i in np.arange(AtomNumber * 3):
            if i not in ZeroIndex:
                ModeVect[index, :] = ModeVectTmp[i, :]
                ModeVect[index, :] /= np.sqrt(np.sum(ModeVectTmp[i, :] * ModeVectTmp[i, :]))
                index += 1
        for i in np.arange(AtomNumber * 3):
            ModeVect[:, i] *= np.sqrt(AtomMass[i])
        
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

        print('orca: No internal coordinate information in {}.log'.format(logfile))
        exit()

    else:
        print('options are needed')
        exit()