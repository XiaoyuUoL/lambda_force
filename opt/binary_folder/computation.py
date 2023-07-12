import numpy as np
import os
import yaml

import input
import g16
import orca
import prepare

def QCRun(task):
    if (input.QCFlag.lower() == 'g16'):
        os.system('g16 {}.gjf'.format(task))
        os.system('formchk {}.chk > /dev/null'.format(task))
    elif (input.QCFlag.lower() == 'orca'):
        os.system('$ORCA/orca {}.inp > {}.log'.format(task, task))
    else:
        print('"QCFlag" now has only two options, "g16" or "orca"')
        exit()

# prepare QC input and run QC calculation
def QCCalculate(task):
    if (task in ['S0opt', 'S0freq', 'S1force', 'S1opt', 'S1nac', 'soc']):
        prepare.QCinput(task)
        QCRun(task)
    elif (task == 'CheckFreq'):
        result = prepare.CheckFreq()
        if (result == -1):
            QCRun(task)
        return result
    else:
        print('error in computation.QCcalculate(): task is not valid')
        exit()

# Diff of ZCoord (angle/dihedral with unit of length)
def ZCoordDiff(qCoord0, qCoord, zIndices):
    def Length(Coord):
        v = Coord[0] - Coord[1]
        return np.sqrt(np.sum(v * v))

    def Angle(Coord):
        v1 = Coord[0] - Coord[1]
        l1 = np.sqrt(np.sum(v1 * v1))
        v2 = Coord[2] - Coord[1]
        l2 = np.sqrt(np.sum(v2 * v2))
        cos = np.sum(v1 * v2) / (l1 * l2)
        if (cos < -1.):
            cos = -1.
        elif (cos > 1.):
            cos = 1.
        return np.arccos(cos)

    def Dihedral(Coord):
        v1 = Coord[0] - Coord[1]
        v2 = Coord[2] - Coord[1]
        pv1 = np.sum(v1 * v2) / np.sum(v2 * v2) * v2 - v1
        v1 = Coord[3] - Coord[2]
        v2 = Coord[1] - Coord[2]
        pv4 = np.sum(v1 * v2) / np.sum(v2 * v2) * v2 - v1
        cos = np.sum(pv1 * pv4) / np.sqrt(np.sum(pv1 * pv1) * np.sum(pv4 * pv4))
        if (cos < -1.):
            cos = -1.
        elif (cos > 1.):
            cos = 1.
        return np.arccos(cos)

    zCoordNum = len(zIndices)
    zCoordDiff = np.zeros(zCoordNum, dtype=float)

    for i,zIndex in enumerate(zIndices):
        qcoord0 = []
        qcoord = []
        for Index in zIndex:
            if (Index == -1):
                break
            qcoord0.append(qCoord0[Index])
            qcoord.append(qCoord[Index])
        # bond
        if (len(qcoord) == 2):
            zCoordDiff[i] = Length(qcoord) - Length(qcoord0)
        # angle
        elif (len(qcoord) == 3):
            zCoordDiff[i] = Angle(qcoord) - Angle(qcoord0)
        # dihedral
        else:
            zCoordDiff[i] = Dihedral(qcoord) - Dihedral(qcoord0)

    return zCoordDiff

# calculate lambda/HR via 4p/displacement approaches
# using notation in [JCP, 2001, 115, 9103.]
# return R = cpp^T * Gp^-1/2 * a^T in equation (20)
def Lambda4p():
    def WilsonMatrix(qCoord, zIndices, Mass, Mode):
        AtomNum = len(qCoord)
        zNum = len(zIndices)

        def WilsonBMatrix(qCoord, zIndices):
            def Vector(coord1, coord2):
                dq = coord2 - coord1
                r = np.sqrt(np.sum(dq * dq))
                return dq / r, r

            def Angle(Coord):
                v1 = Coord[0] - Coord[1]
                l1 = np.sqrt(np.sum(v1 * v1))
                v2 = Coord[2] - Coord[1]
                l2 = np.sqrt(np.sum(v2 * v2))
                cos = np.sum(v1 * v2) / (l1 * l2)
                if (cos < -1.):
                    cos = -1.
                elif (cos > 1.):
                    cos = 1.
                return np.arccos(cos)

            BMatrix = np.zeros((zNum, 3 * AtomNum), dtype=float)
            # analytical results (Molecular Vibrations, E. Bright Wilson, etc.)
            for i,zindex in enumerate(zIndices):
                # length
                if (zindex[2] == -1):
                    i1,i2 = zindex[0],zindex[1]
                    e12,r12 = Vector(qCoord[i2], qCoord[i1])
                    BMatrix[i, i1*3:i1*3+3] = -e12
                    BMatrix[i, i2*3:i2*3+3] = e12
                # angle
                elif (zindex[3] == -1):
                    i1,i2,i3 = zindex[0],zindex[1],zindex[2]
                    e21,r21 = Vector(qCoord[i1], qCoord[i2])
                    e23,r23 = Vector(qCoord[i3], qCoord[i2])
                    theta = Angle([qCoord[i1], qCoord[i2], qCoord[i3]])
                    cos = np.cos(theta)
                    sin = np.sin(theta)
                    s1 = (cos * e21 - e23) / (sin * r21)
                    s3 = (cos * e23 - e21) / (sin * r23)
                    BMatrix[i, i1*3:i1*3+3] = s1
                    BMatrix[i, i2*3:i2*3+3] = -s1 - s3
                    BMatrix[i, i3*3:i3*3+3] = s3
                # dihedrals
                else:
                    i1,i2,i3,i4 = zindex[0],zindex[1],zindex[2],zindex[3]
                    e12,r12 = Vector(qCoord[i2], qCoord[i1])
                    e23,r23 = Vector(qCoord[i3], qCoord[i2])
                    e32,r32 = Vector(qCoord[i2], qCoord[i3])
                    e43,r43 = Vector(qCoord[i3], qCoord[i4])
                    theta2 = Angle([qCoord[i1], qCoord[i2], qCoord[i3]])
                    sin2 = np.sin(theta2)
                    cos2 = np.cos(theta2)
                    theta3 = Angle([qCoord[i2], qCoord[i3], qCoord[i4]])
                    sin3 = np.sin(theta3)
                    cos3 = np.cos(theta3)
                    s1 = -np.cross(e12, e23) / (r12 * sin2 * sin2)
                    s4 = -np.cross(e43, e32) / (r43 * sin3 * sin3)
                    BMatrix[i, i1*3:i1*3+3] = s1
                    BMatrix[i, i2*3:i2*3+3] = -s1 * (r23 - r12 * cos2) / r23 - s4 * r43 * cos3 / r23
                    BMatrix[i, i3*3:i3*3+3] = -s4 * (r32 - r43 * cos3) / r32 - s1 * r12 * cos2 / r32
                    BMatrix[i, i4*3:i4*3+3] = s4

            return BMatrix

        B = WilsonBMatrix(qCoord, zIndices)
        G = np.matmul(B, np.matmul(np.diag(1.0 / Mass), B.T))

        e,a = np.linalg.eigh(G)
        r = 0
        for i in np.arange(len(e)):
            if (e[i] <= input.de):
                r += 1
        e = e[r:]
        a = a[:, r:]

        E = np.matmul(a, np.matmul(np.diag(1. / e), a.T))
        tmp = np.diag(1.0 / np.sqrt(Mass))
        return np.matmul(np.matmul(Mode.T, tmp), np.matmul(B.T, E))

    if (input.QCFlag.lower() == 'g16'):
        Name,InitCoord = g16.FchkRead('S1opt', 'coord')
        Name,FinalCoord = g16.FchkRead('S0freq', 'coord')
        Zindices = g16.FchkRead('S1opt', 'zindices')
        AtomMass,ModeFreq,ModeQ,ModeVect = g16.FchkRead('S0freq', 'vibration')
        ES0S0 = g16.FchkRead('S0opt', 'energy-gd')
        ES1S0 = g16.FchkRead('S1force', 'energy-ex')[0]
        ES0S1 = g16.FchkRead('S1opt', 'energy-gd')
        ES1S1 = g16.FchkRead('S1opt', 'energy-ex')[0]
    elif (input.QCFlag.lower() == 'orca'):
        Name,InitCoord = orca.LogRead('S1opt', 'coord')
        Name,FinalCoord = orca.LogRead('S0freq', 'coord')
        Zindices = orca.LogRead('S1opt', 'zindices')
        AtomMass,ModeFreq,ModeQ,ModeVect = orca.LogRead('S0freq', 'vibration')
        ES0S0 = orca.LogRead('S0opt', 'energy-gd')
        ES1S0 = orca.LogRead('S1force', 'energy-ex')[0] + ES0S0
        ES0S1 = orca.LogRead('S1opt', 'energy-gd')
        ES1S1 = orca.LogRead('S1opt', 'energy-ex')[0] + ES0S1
    else:
        print('"QCFlag" now has only two options, "g16" or "orca"')
        exit()

    Lambda10 = ES0S1 - ES0S0
    Lambda01 = ES1S0 - ES1S1

    # remove imaginary frequency modes
    ModeNumber = len(ModeFreq)
    for i in np.arange(ModeNumber - 1, -1, -1):
        if (ModeFreq[i] < 0.):
            ModeFreq = np.delete(ModeFreq, i)
            ModeQ = np.delete(ModeQ, i)
            ModeVect = np.delete(ModeVect, i, 0)
        
    R = WilsonMatrix(FinalCoord, Zindices, AtomMass, ModeVect.T)
    ModeQ0 = np.matmul(R, ZCoordDiff(FinalCoord, InitCoord, Zindices)) * ModeQ
    HR = 0.5 * ModeQ0 * ModeQ0
    LambdaDisp = ModeFreq * HR

    ModeFreq /= (2 * np.pi * input.c * input.au2fs)
    fout = open('../result_folder/HRfactor_disp.dat', 'w')
    fout.writelines(' freq/cm^-1        Q            HR        lambda/cm^-1\n')
    for i in np.arange(len(ModeFreq)):
        fout.writelines('{:10.4f}{:14.7f}{:14.7f}{:14.7f}\n'.format(ModeFreq[i],
            ModeQ0[i], HR[i], LambdaDisp[i] / (2 * np.pi * input.c * input.au2fs)))
    fout.close()

    return Lambda10, Lambda01, np.sum(LambdaDisp)

# calculate lambda/HR via force approach
def LambdaForce():
    if (input.QCFlag.lower() == 'g16'):
        Gradient = g16.FchkRead('S1force', 'gradient')
        AtomMass,ModeFreq,ModeQ,ModeVect = g16.FchkRead('S0freq', 'vibration')
    elif (input.QCFlag.lower() == 'orca'):
        Gradient = orca.LogRead('S1force', 'gradient')
        AtomMass,ModeFreq,ModeQ,ModeVect = orca.LogRead('S0freq', 'vibration')
    else:
        print('"QCFlag" now has only two options, "g16" or "orca"')
        exit()

    # remove imaginary frequency modes
    ModeNumber = len(ModeFreq)
    for i in np.arange(ModeNumber - 1, -1, -1):
        if (ModeFreq[i] < 0.):
            ModeFreq = np.delete(ModeFreq, i)
            ModeQ = np.delete(ModeQ, i)
            ModeVect = np.delete(ModeVect, i, 0)

    Gradient /= np.sqrt(AtomMass)
    ModeForce = np.dot(ModeVect, Gradient) / ModeQ
    ModeQ0 = ModeForce / (input.hbar * ModeFreq)
    HR = 0.5 * ModeQ0 * ModeQ0
    Lambda = ModeFreq * HR

    ModeFreq /= (2 * np.pi * input.c * input.au2fs)
    fout = open('../result_folder/HRfactor_force.dat', 'w')
    fout.writelines(' freq/cm^-1       Q             HR       lambda/cm^-1\n')
    for i in np.arange(len(ModeFreq)):
        fout.writelines('{:10.4f}{:14.7f}{:14.7f}{:14.7f}\n'.format(ModeFreq[i],
            ModeQ0[i], HR[i], Lambda[i] / (2 * np.pi * input.c * input.au2fs)))
    fout.close()

    return np.sum(Lambda)

# calculate bond order difference
def BOD():
    # Bond order
    # [Mayer, I., Chem. Phys. Lett., 1983, 97, 270.] for close shell
    # [Mayer, I. Int. J. Quantum Chem., 1986, XXIX, 73.] for open shell 
    def BondOrder(indices, S, Ca, Cb):
        Pa = np.matmul(Ca, Ca.T)  # density matrix of alpha electron
        Pb = np.matmul(Cb, Cb.T)  # density matrix of beta electron
        PSa = np.matmul(Pa, S)
        PSb = np.matmul(Pb, S)
        B = np.zeros((len(indices), len(indices)), dtype=float)
        for i,indexa in enumerate(indices):
            for j,indexb in enumerate(indices):
                PSai = PSa[indexa[0]:indexa[1], indexb[0]:indexb[1]]
                PSaj = PSa[indexb[0]:indexb[1], indexa[0]:indexa[1]]
                PSbi = PSb[indexa[0]:indexa[1], indexb[0]:indexb[1]]
                PSbj = PSb[indexb[0]:indexb[1], indexa[0]:indexa[1]]
                B[i, j] = 2. * np.trace(np.matmul(PSai, PSaj) + np.matmul(PSbi, PSbj))

        return B

    names,coords,charge,spin = prepare.FromSmiles(input.SystemSmiles)
    if (spin != 1):
        print('error in computation.BOD(): now only valid for close shell system')
        exit()

    for name in names:
        if (str(name) not in input.BasisFunc.keys()):
            prepare.QCinput('S0sp', name)
            QCRun('S0sp')
            if (input.QCFlag.lower() == 'g16'):
                input.BasisFunc[str(name)] = g16.FchkRead('S0sp', 'basis')
            elif (input.QCFlag.lower() == 'orca'):
                input.BasisFunc[str(name)] = orca.LogRead('S0sp', 'basis')
            else:
                print('"QCFlag" now has only two options, "g16" or "orca"')
                exit()
            os.system('rm S0sp*')

    if (input.QCFlag.lower() == 'g16'):
        AOIndices,S,COcc,CUnocc = g16.FchkRead('S0opt', 'orbital')
        Zindices = g16.FchkRead('S0opt', 'zindices')
    elif (input.QCFlag.lower() == 'orca'):
        AOIndices,S,COcc,CUnocc = orca.LogRead('S0opt', 'orbital')
        Zindices = orca.LogRead('S0opt', 'zindices')
    else:
        print('"QCFlag" now has only two options, "g16" or "orca"')
        exit()

    # calculate bond order for S0 and S1 (assume HOMO->LUMO excitation)
    CS0a = np.zeros_like(COcc)
    CS0a[:, :] = COcc
    CS0b = np.zeros_like(COcc)
    CS0b[:, :] = COcc
    BOS0 = BondOrder(AOIndices, S, CS0a, CS0b)
    CS1a = np.zeros((np.shape(COcc)[0], np.shape(COcc)[1] + 1), dtype=float)
    CS1a[:, :-1] = COcc
    CS1a[:, -2] *= np.sqrt(0.5)
    CS1a[:, -1] = np.sqrt(0.5) * CUnocc[:, 0]
    CS1b = np.zeros_like(CS1a)
    CS1b[:, :] = CS1a
    BOS1 = BondOrder(AOIndices, S, CS1a, CS1b)

    BODT = 0.0
    fout = open('../result_folder/BOD.dat', 'wt')
    fout.writelines('    a    b       BOD\n')
    for index in Zindices:
        if(index[2] == -1):
            a = index[0]
            b = index[1]
            BOD = BOS1[a, b] - BOS0[a, b]
            fout.writelines('{:5d}{:5d}{:14.5e}\n'.format(a + 1, b + 1, BOD))
            BODT += BOD * BOD
    fout.close()

    return BODT

def SOC():
    if (input.QCFlag.lower() == 'g16'):
        print('"SOC" is only for ORCA now')
        exit()
    elif (input.QCFlag.lower() == 'orca'):
        SOC0 = orca.LogRead('soc', 'soc')
    else:
        print('"QCFlag" now has only two options, "g16" or "orca"')
        exit()

    fout = open('../result_folder/SOC.dat', 'wt')
    fout.writelines('    S    T            m=0/cm^-1                   m=-1/cm^-1                  m=1/cm^-1          total/cm^-1\n')
    for soc in SOC0:
        Result = '{:5d}{:5d}'.format(soc[0], soc[1])

        soc0 = soc[2] / (2 * np.pi * input.c * input.au2fs)
        soct = (soc0 * soc0.conj()).real
        Result += '{:14.5e}{:14.5e}'.format(soc0.real, soc0.imag)

        socd = soc[3] / (2 * np.pi * input.c * input.au2fs)
        soct += (socd * socd.conj()).real
        Result += '{:14.5e}{:14.5e}'.format(socd.real, socd.imag)

        socu = soc[4] / (2 * np.pi * input.c * input.au2fs)
        soct += (socu * socu.conj()).real
        Result += '{:14.5e}{:14.5e}'.format(socu.real, socu.imag)

        soct = np.sqrt(soct)
        soc.append(soct * 2 * np.pi * input.c * input.au2fs)
        Result += '{:14.5e}{:14.5e}'.format(soc0.real, soc0.imag)

        fout.writelines('{}\n'.format(Result))
    fout.close()

    return SOC0

#### need to be checked
def NAC():
    if (input.QCFlag.lower() == 'g16'):
        NACs0 = g16.FchkRead('S1nac', 'nac')
        AtomMass,ModeFreq,ModeQ,ModeVect = g16.FchkRead('S0freq', 'vibration')
    elif (input.QCFlag.lower() == 'orca'):
        NACs0 = orca.LogRead('S1nac', 'nac')
        AtomMass,ModeFreq,ModeQ,ModeVect = orca.LogRead('S0freq', 'vibration')
    else:
        print('"QCFlag" now has only two options, "g16" or "orca"')
        exit()

    # remove imaginary frequency modes
    ModeNumber = len(ModeFreq)
    for i in np.arange(ModeNumber - 1, -1, -1):
        if (ModeFreq[i] < 0.):
            ModeFreq = np.delete(ModeFreq, i)
            ModeQ = np.delete(ModeQ, i)
            ModeVect = np.delete(ModeVect, i, 0)

    NACs0 /= np.sqrt(AtomMass)
    NACs = input.hbar * ModeFreq * np.dot(ModeVect, NACs0) / ModeQ

    ModeFreq /= (2 * np.pi * input.c * input.au2fs)
    fout = open('../result_folder/NAC.dat', 'w')
    fout.writelines(' freq/cm^-1    NAC/cm^-1\n')
    for i in np.arange(len(ModeFreq)):
        fout.writelines('{:10.4f}{:14.7f}\n'.format(ModeFreq[i],
            NACs[i] / (2 * np.pi * input.c * input.au2fs)))
    fout.close()

    return NACs                                                                 # unit: Hartree

# Get S0 geometry via QC package
def GeomCalculate():
    # S0opt
    QCCalculate('S0opt')

    # S1freq
    QCCalculate('S0freq')

    # check if S0 opt need to be rerun
    while (QCCalculate('CheckFreq') == -1):
        QCCalculate('S0freq')

    return

# Calculate properties via QC package and homemake script, including
# 'lambda_4p': reorganization energy via 4-point approach
#              and reorganization energy/HR factor via displacement approach
#              (J. Chem. Phys., 2001, 115, 9103.)
# 'lambda_force': reorganization energy/HR factor via force approach
#                 (https://doi.org/10.1021/acs.jpclett.3c00749)
# 'BOD': bond order difference under HOM->LUMO excitation assumption
# 'SOC': spin-orbit coupling (only in ORCA)
# 'MAC': nonadiabatic coupling
def PropCalculate():
    results = {}
    # 'lambda_4p': reorganization energy via 4-point/displacement approach
    if ('lambda_4p' in input.Properties):
        #QCCalculate('S1opt')
        #QCCalculate('S1force')

        Lambda10,Lambda01,LambdaDisp = Lambda4p()
        results['lambda_10'] = float(Lambda10)                                  # unit: Hartree
        results['lambda_01'] = float(Lambda01)                                  # unit: Hartree
        results['lambda_4p'] = float((Lambda01 + Lambda10) * 0.5)               # unit: Hartree
        results['lambda_disp'] = float(LambdaDisp)                              # unit: Hartree

        # 'lambda_force': reorganization energy/HR factor via force approach
        if ('lambda_force' in input.Properties):
            results['lambda_force'] = float(LambdaForce())                      # unit: Hartree

            os.system('mv S1force* ../result_folder/')

        os.system('mv S1opt* ../result_folder/')

    # 'lambda_force': reorganization energy/HR factor via force approach
    elif ('lambda_force' in input.Properties):
        #QCCalculate('S1force')

        results['lambdaForce'] = float(LambdaForce())                           # unit: Hartree

        os.system('mv S1force* ../result_folder/')

    if ('NAC' in input.Properties):
        #QCCalculate('S1nac')

        NACs = NAC()
        results['NAC'] = float(np.sqrt(np.sum(NACs * NACs)))                    # unit: Hartree

        os.system('mv S1nac* ../result_folder/')

    if ('SOC' in input.Properties):
        #QCCalculate('soc')

        SOC0 = SOC()
        results['SOC'] = float(SOC0[input.IRoot - 1][5])                        # unit: Hartree

        os.system('mv soc* ../result_folder/')

    # 'BOD': bond order difference under HOM->LUMO excitation assumption
    if ('BOD' in input.Properties):
        results['BOD'] = float(BOD())

    if (len(results) == 0):
        print('no property is calculated.')
        exit()
    else:
        fout = open('../result_folder/result.yml', 'wt')
        yaml.dump(results, fout)
        fout.close()

    os.system('mv S0opt* ../result_folder/')
    os.system('mv S0freq* ../result_folder/')