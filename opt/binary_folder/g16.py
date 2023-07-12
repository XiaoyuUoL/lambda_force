import numpy as np

import input

para = {
    'file': '',
    'procnumber': input.ProcNumber,
    'memory': '{:>d}MW'.format(input.Memory),
    'chkfile': '',
    'keywords': [],
    'desc': '',
    'charge': 0,
    'multi': 1,
    'name': [],
    'coord': []
}

# common keywords for different calculation task
# can be changed (e.g., 'opt=loose' for Loose optimization)
keywords = {
    'S0sp': [input.Functional, input.BasisSet],
    'S0opt': [input.Functional, input.BasisSet, 'opt=loose'],
    'S0freq': [input.Functional, input.BasisSet, 'freq'],
    'S1force': [input.Functional, input.BasisSet, 'force',
                'td(nstates={:d},root={:d})'.format(input.NRoots, input.IRoot)],
    'S1opt': [input.Functional, input.BasisSet, 'opt',
              'td(nstates={:d},root={:d})'.format(input.NRoots, input.IRoot)],
    'S1nac': [input.Functional, input.BasisSet,
              'td(nstates={:d},root={:d},NonAdiabaticCoupling)'.format(
              input.NRoots, input.IRoot)],
}

def InitPara(task):
    para['desc'] = task
    para['keywords'] = keywords[task]
    para['file'] = '{}'.format(task)
    para['chkfile'] = '{}.chk'.format(task)

# generate input file for g16 calculation
def GjfGen(para):
    fout = open('{}.gjf'.format(para['file']), 'w')
    fout.writelines('%nprocshared={:d}\n'.format(para['procnumber']))
    fout.writelines('%mem={}\n'.format(para['memory']))
    fout.writelines('%chk={}\n'.format(para['chkfile']))
    keywords = '#p'
    for keyword in para['keywords']:
        keywords += ' ' + keyword
    fout.writelines('{}\n'.format(keywords))
    fout.writelines('\n')
    fout.writelines('{}\n'.format(para['desc']))
    fout.writelines('\n')
    fout.writelines('{:d} {:d}\n'.format(para['charge'], para['multi']))
    for name, coord in zip(para['name'], para['coord']):
        fout.writelines('{}{:14.7f}{:14.7f}{:14.7f}\n'.format(
            name, coord[0], coord[1], coord[2]))
    fout.writelines('\n')
    fout.writelines('\n')
    fout.close()

## read fchk file. keyword options:
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
def FchkRead(fchkfile, keyword):
    fin = open('{}.fchk'.format(fchkfile), 'r')

    # 'energy-gd': return energy of electronic ground state
    if (keyword.lower() == 'energy-gd'):
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            # read SCF energy
            if (len(words) > 2 and words[:2] == ['SCF', 'Energy']):
                fin.close()
                return float(words[3])                                          # unit: Hartree
            line = fin.readline()
        fin.close()

        print('error of g16.ReadFchk(): No ground energy information in {}.fchk'.format(fchkfile))
        exit()

    # 'energy': return energies of electronic excited state
    elif (keyword.lower() == 'energy-ex'):
        StateNumber = None
        IfGetEI = False
        ExcitedInfo = []
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            # read excited state information
            if (len(words) > 3 and words[:3] == ['ETran', 'state', 'values']):
                IfGetEI = True
                StateNumber = int(words[5]) // 16
                for i in range(int((StateNumber * 16 - 1) / 5) + 1):
                    ExcitedInfo += fin.readline().rstrip().split()
            line = fin.readline()
        fin.close()

        if (not IfGetEI):
            print('error of g16.ReadFchk(): No excited energy information in {}.fchk'.format(fchkfile))
            exit()

        ExcitedInfo = np.reshape(np.array(ExcitedInfo, dtype=float),
            (StateNumber, 16))

        # ExcitedInfo[:, 0]: energies of excited states                         (unit: Hartree)
        # ExcitedInfo[:, 1:4]: transition electric dipole moments               (unit: a.u.)
        # ExcitedInfo[:, 4:7]: transition velocity dipole moments               (unit: a.u.)
        # ExcitedInfo[:, 7:10]: transition magnetic dipole moments              (unit: a.u.)
        # ExcitedInfo[:, 10:]: transition velocity quadrupole moments           (unit: a.u.)
        return ExcitedInfo[:, 0]                                                # unit: Hartree

    # 'coord': return atom name and atom coordinate
    elif (keyword.lower() == 'coord'):
        IfGetAN = False
        AtomNumber = None
        IfGetName = False
        Name = []
        IfGetCoord = False
        Coord = []
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            # read number of atoms
            if (len(words) > 3 and words[:3] == ['Number', 'of', 'atoms']):
                IfGetAN = True
                AtomNumber = int(words[4])
            # read atomic number (name)
            elif (len(words) > 2 and words[:2] == ['Atomic', 'numbers']
                and IfGetAN):
                IfGetName = True
                for i in range(int((AtomNumber - 1) / 6) + 1):
                    Name += fin.readline().rstrip().split()
            # read Cartesian coordinates
            elif (len(words) > 3 and words[:3] == ['Current', 'cartesian',
                'coordinates'] and IfGetAN):
                IfGetCoord = True
                for i in np.arange(int((AtomNumber * 3 - 1) / 5) + 1):
                    Coord += fin.readline().rstrip().split()
            line = fin.readline()
        fin.close()

        if (not IfGetName or not IfGetCoord):
            print('error of g16.ReadFchk(): No coordinate information in {}.fchk'.format(fchkfile))
            exit()

        Coord = np.reshape(np.array(Coord, dtype=float), (AtomNumber, 3))       # unit: Bohr
        return Name, Coord

    # 'gradient': return energy gradient
    elif (keyword.lower() == 'gradient'):
        IfGetAN = False
        AtomNumber = None
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            # read number of atoms
            if (len(words) > 3 and words[:3] == ['Number', 'of', 'atoms']):
                IfGetAN = True
                AtomNumber = int(words[4])
            # read gradient
            elif (len(words) > 2 and words[:2] == ['Cartesian', 'Gradient']
                and IfGetAN):
                Gradient = []
                for i in np.arange(int((AtomNumber * 3 - 1) / 5) + 1):
                    Gradient += fin.readline().rstrip().split()
                fin.close()
                return np.array(Gradient, dtype=float)                          # unit: Hartree/Bohr, 3N 1D-array
            line = fin.readline()
        fin.close()

        print('error of g16.ReadFchk(): No gradient information in {}.fchk'.format(fchkfile))
        exit()

    # 'vibration': return normal mode information
    elif (keyword.lower() == 'vibration'):
        IfGetAN = False
        AtomNumber = None
        IfGetMN = False
        ModeNumber = None
        IfGetMass = False
        AtomMass = []
        IfGetModeInfo = False
        ModeInfo = []
        IfGetModeVect = False
        ModeVect = []
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            # read number of atoms
            if (len(words) > 3 and words[:3] == ['Number', 'of', 'atoms']):
                IfGetAN = True
                AtomNumber = int(words[4])
            # read number of normal modes
            elif (len(words) > 4 and words[:4] == ['Number', 'of', 'Normal',
                'Modes']):
                IfGetMN = True
                ModeNumber = int(words[5])
            # read information of atom masses (3N 1D-array)
            elif (len(words) > 1 and words[0] == 'Vib-AtMass' and IfGetAN):
                IfGetMass = True
                for i in np.arange(int((AtomNumber - 1) / 5) + 1):
                    masses = fin.readline().rstrip().split()
                    for mass in masses:
                        AtomMass += [mass, mass, mass]
            # read information of normal modes: 
            # frequency, reduce mass, force constant and IR intensity
            elif (len(words) > 1 and words[0] == 'Vib-E2' and IfGetMN):
                IfGetModeInfo = True
                for i in np.arange(int((ModeNumber * 4 - 1) / 5) + 1):
                    ModeInfo += fin.readline().rstrip().split()
            # read information of normal modes: vibration vector
            elif (len(words) > 1 and words[0] == 'Vib-Modes' and IfGetAN
                and IfGetMN):
                IfGetModeVect = True
                for i in np.arange(int((ModeNumber * AtomNumber * 3 - 1) / 5) + 1):
                    ModeVect += fin.readline().rstrip().split()
            line = fin.readline()
        fin.close()

        if (not IfGetMass or not IfGetModeInfo or not IfGetModeVect):
            print('error of g16.ReadFchk(): No normal mode information in {}.fchk'.format(fchkfile))
            exit()

        AtomMass = np.array(AtomMass, dtype=float)                              # unit: amu

        ModeInfo = np.reshape(np.array(ModeInfo[:ModeNumber * 4], dtype=float),
            (4, ModeNumber))
        ModeFreq = ModeInfo[0] * 2 * np.pi * input.c * input.au2fs              # unit: Hartree
        ReduceMass = ModeInfo[1]                                                # unit: amu
        #ForceConst = ModeInfo[2]
        #IRIntensity = ModeInfo[3]

        ModeVect = np.reshape(np.array(ModeVect, dtype=float),
            (ModeNumber, AtomNumber * 3))
        for i in np.arange(ModeNumber):
            ModeVect[i, :] /= np.sqrt(ReduceMass[i])                            # displacement (not mass-weigthed) but |Q|=1
        for i in np.arange(AtomNumber * 3):
            ModeVect[:, i] *= np.sqrt(AtomMass[i])                              # mass-weigthed displacement

        ModeQ = np.zeros_like(ModeFreq)
        for i,freq in enumerate(ModeFreq):
            if (ModeFreq[0] > 0.):
                ModeQ[i] = np.sqrt(freq * (1. / input.au2aum) / input.hbar)

        return AtomMass, ModeFreq, ModeQ, ModeVect

    # return redundant internal coordinate indices
    elif (keyword.lower() == 'zindices'):
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            # read redundant internal coordinate indices
            if (len(words) > 4 and words[:4] == ['Redundant', 'internal',
                'coordinate', 'indices']):
                ZCoordNumber = int(words[6]) // 4
                ZIndices = []
                for i in np.arange(int((ZCoordNumber * 4 - 1) / 6) + 1):
                    ZIndices += fin.readline().rstrip().split()
                fin.close()
                ZIndices = np.reshape(np.array(ZIndices, dtype=int),
                    (ZCoordNumber, 4)) - 1                                      # number of '-1': 0 for dihedral, 1 for angle and 2 for length
                # remove invalid indices
                for i in np.arange(ZCoordNumber - 1, -1, -1):
                    if (ZIndices[i, -1] < -1):
                        ZIndices = np.delete(ZIndices, i, 0)
                return ZIndices
            line = fin.readline()
        fin.close()

        print('error of g16.ReadFchk(): No redundant internal coordinate indices information in {}.fchk'.format(fchkfile))
        exit()

    # return number of basis functions
    elif (keyword.lower() == 'basis'):
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            if(len(words) > 4 and words[:4] == ['Number', 'of', 'basis',
                'functions']):
                fin.close()
                return int(words[-1])
            line = fin.readline()
        fin.close()

        print('error of g16.ReadFchk(): No basis functions number information in {}.fchk'.format(fchkfile))
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
        IfGetAN = False
        Elements = []
        BasisFunc = []
        BasisNumber = None
        IfGetMO = False
        MOCoeff = []
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            # read number of alpha/beta electrons (index for HOMO)
            if(len(words) > 4 and words[:4] == ['Number', 'of', 'alpha',
                'electrons']):
                IfGetON = True
                OccNumber = int(words[-1])
            # read atomic number of atoms
            elif(len(words) > 2 and words[:2] == ['Atomic', 'numbers']):
                IfGetAN = True
                AtomNumber = int(words[-1])
                for i in range(int((AtomNumber - 1) / 6) + 1):
                    Elements += fin.readline().rstrip().split()
                BasisFunc = GetBasisFunc(Elements)
                BasisNumber = np.sum(BasisFunc)
            # read MO coefficients
            elif(len(words) > 3 and words[:3] == ['Alpha', 'MO', 'coefficients']
                and IfGetON and IfGetAN):
                IfGetMO = True
                for i in np.arange(int((BasisNumber * BasisNumber - 1) / 5) + 1):
                    MOCoeff += fin.readline().rstrip().split()
            line = fin.readline()
        fin.close()

        if (not IfGetMO):
            print('error of g16.ReadFchk(): No orbital information in {}.fchk'.format(fchkfile))
            exit()

        MOCoeff = np.reshape(np.array(MOCoeff, dtype=float),
            (BasisNumber, BasisNumber)).T
        AOoverlap = MatrixPower(np.matmul(MOCoeff, MOCoeff.T), -1.0)
        # indices of basis function for each atom
        indices = np.array([[0] + list(np.cumsum(BasisFunc)[:-1]),
            list(np.cumsum(BasisFunc))]).T 

        return indices, AOoverlap, MOCoeff[:, :OccNumber], MOCoeff[:, OccNumber:]

    # return nonadiabatic coupling
    elif (keyword.lower() == 'nac'):
        IfGetAN = False
        AtomNumber = None
        IFGetNAC = False
        NAC = []
        line = fin.readline()
        while (len(line) != 0):
            words = line.rstrip().split()
            # read number of atoms
            if (len(words) > 3 and words[:3] == ['Number', 'of', 'atoms']):
                IfGetAN = True
                AtomNumber = int(words[4])
            # read nonadiabatic coupling
            elif(len(words) > 2 and words[:2] == ['Nonadiabatic', 'coupling']
                and IfGetAN):
                IFGetNAC = True
                for i in np.arange(int((AtomNumber * 3 - 1) / 5) + 1):
                    NAC += fin.readline().rstrip().split()
            line = fin.readline()
        fin.close()

        if (not IFGetNAC):
            print('error of g16.ReadFchk(): No nonadiabatic coupling information in {}.fchk'.format(fchkfile))
            exit()

        NAC = np.array(NAC, dtype=float)                                        # unit: 1/Bohr (3N 1D-array)
        return NAC

    else:
        print('error of g16.ReadFchk(): option is needed')
        exit()