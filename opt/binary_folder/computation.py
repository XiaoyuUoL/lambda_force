import numpy as np
import os
import yaml

import input
import g16
import orca
import prepare

# prepare QC input and run QC calculation
def QCCalculate(task):
    if (task in ['S0opt', 'S0freq', 'S1force', 'S1opt']):
        prepare.QCinput(task)
        os.system('./{}_{}.sh'.format(input.QCFlag.lower(), task))     #### check the command of G16/ORCA
    elif (task == 'CheckFreq'):
        result = prepare.CheckFreq()
        if (result == -1):
            os.system('./{}_{}.sh'.format(input.QCFlag.lower(), 'S0opt'))     #### check the command of G16/ORCA
        return result
    else:
        print('error: task is not valid for computation.QCcalculate')
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
    zCoordDiff = np.zeros(zCoordNum, dtype = float)

    for i, zIndex in enumerate(zIndices):
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
            zCoordDiff[i] = (Angle(qcoord) - Angle(qcoord0))
        # dihedral
        else:
            zCoordDiff[i] = (Dihedral(qcoord) - Dihedral(qcoord0))
    
    return zCoordDiff

# transfer Cartesian coordinates to internal coordinates according to the 4-number ZIndices
# using notation in JCP, 2001, 115, 9103
# return R = cpp^T * Gp^-1/2 * a^T in equation (20)
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
        ## numberial derivative Bij = partail zi / partial qj (dq = 0.1)
        #if (mode == 'n'):
        #    qCoordP = np.zeros_like(qCoord)
        #    qCoordM = np.zeros_like(qCoord)
        #    for i in np.arange(AtomNum):
        #        for j in np.arange(3):
        #            qCoordP[:] = qCoord[:]
        #            qCoordP[i, j] += dq
        #            qCoordM[:] = qCoord[:]
        #            qCoordM[i, j] -= dq
        #            BMatrix[:, i * 3 + j] = ZCoordDiff(qCoordP, qCoordM, zIndices) / (2.0 * dq)
        #else:
        # analytical results (Molecular Vibrations, E. Bright Wilson, etc.)
        for i,zindex in enumerate(zIndices):
            # length
            if (zindex[2] == -1):
                i1,i2 = zindex[0],zindex[1]
                e12,r12 = Vector(qCoord[i2], qCoord[i1])
                BMatrix[i, i1 * 3: i1 * 3 + 3] = -e12
                BMatrix[i, i2 * 3: i2 * 3 + 3] = e12
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
                BMatrix[i, i1 * 3: i1 * 3 + 3] = s1
                BMatrix[i, i2 * 3: i2 * 3 + 3] = -s1 - s3
                BMatrix[i, i3 * 3: i3 * 3 + 3] = s3
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
                BMatrix[i, i1 * 3: i1 * 3 + 3] = s1
                BMatrix[i, i2 * 3: i2 * 3 + 3] = -s1 * (r23 - r12 * cos2) / r23 - s4 * r43 * cos3 / r23
                BMatrix[i, i3 * 3: i3 * 3 + 3] = -s4 * (r32 - r43 * cos3) / r32 - s1 * r12 * cos2 / r32
                BMatrix[i, i4 * 3: i4 * 3 + 3] = s4
        
        return BMatrix
    
    BMatrix = WilsonBMatrix(qCoord, zIndices)
    GMatrix = np.matmul(BMatrix, np.matmul(np.diag(1.0 / Mass), BMatrix.T))

    e, a = np.linalg.eigh(GMatrix)
    r = 0
    for i in np.arange(len(e)):
        if (e[i] <= input.de):
            r += 1
    e = e[r:]
    a = a[:, r:]

    EMatrix = np.matmul(a, np.matmul(np.diag(1. / e), a.T))
    return np.matmul(np.matmul(Mode.T, np.diag(1.0 / np.sqrt(Mass))), np.matmul(BMatrix.T, EMatrix))

# calculate HR via opt geometries (dushin)
def Lambda4pCal():
    if (input.QCFlag.lower() == 'g16'):
        Name, InitCoord = g16.FchkRead('S1opt', 'coord')
        Name, FinalCoord = g16.FchkRead('S0freq', 'coord')
        Zindices = g16.FchkRead('S1opt', 'zindices')
        AtomMass, ModeFreq, ModeQ, ModeVect = g16.FchkRead('S0freq', 'vibration')
        ES0S0 = g16.FchkRead('S0opt', 'energy-gd')
        ES1S0 = g16.FchkRead('S1force', 'energy-ex')[0]
        ES0S1 = g16.FchkRead('S1opt', 'energy-gd')
        ES1S1 = g16.FchkRead('S1opt', 'energy-ex')[0]
    elif (input.QCFlag.lower() == 'orca'):
        Name, InitCoord = orca.LogRead('S1opt', 'coord')
        Name, FinalCoord = orca.LogRead('S0freq', 'coord')
        Zindices = orca.LogRead('S1opt', 'zindices')
        AtomMass, ModeFreq, ModeQ, ModeVect = orca.LogRead('S0freq', 'vibration')
        ########### NEED TO BE UPDATED TO ADD READING OPTION FOR EXCITED STATE ENERGY IN ORCA
        ES0S0 = orca.LogRead('S0opt', 'energy-gd')
        ES1S0 = orca.LogRead('S1force', 'energy-ex')[0] + ES0S0
        ES0S1 = orca.LogRead('S1opt', 'energy-gd')
        ES1S1 = orca.LogRead('S1opt', 'energy-ex')[0] + ES0S1
    else:
        print('"QCFlag" now has only two options: g16 or orca')
        exit()
    
    # remove imaginary frequency modes
    ModeNumber = len(ModeFreq)
    for i in np.arange(ModeNumber - 1, -1, -1):
        if (ModeFreq[i] < 0.):
            ModeFreq = np.delete(ModeFreq, i)
            ModeQ = np.delete(ModeQ, i)
            ModeVect = np.delete(ModeVect, i, 0)
        
    R = WilsonMatrix(FinalCoord, Zindices, AtomMass, ModeVect.T)
    ModeQ0 = np.matmul(R, ZCoordDiff(FinalCoord, InitCoord, Zindices)) * ModeQ

    return ES0S1 - ES0S0, ES1S0 - ES1S1, ModeFreq, ModeQ0
    
# calculate HR via force/gradient
def LambdaForceCal():
    if (input.QCFlag.lower() == 'g16'):
        Gradient = g16.FchkRead('S1force', 'gradient')
        AtomMass, ModeFreq, ModeQ, ModeVect = g16.FchkRead('S0freq', 'vibration')
    elif (input.QCFlag.lower() == 'orca'):
        Gradient = orca.LogRead('S1force', 'gradient')
        AtomMass, ModeFreq, ModeQ, ModeVect = orca.LogRead('S0freq', 'vibration')
    else:
        print('"QCFlag" now has only two options: g16 or orca')
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

    return ModeFreq, ModeQ0

def HRCalculate():
    results = {}
    if ('force' in input.Properties and '4p' in input.Properties):
        Lambda10,Lambda01,Freq,QDisp = Lambda4pCal()
        Lambda10 *= input.au2eV
        Lambda01 *= input.au2eV
        HRDisp = 0.5 * QDisp * QDisp
        LambdaDisp = Freq * HRDisp
        Freq,QForce = LambdaForceCal()
        HRForce = 0.5 * QForce * QForce
        LambdaForce = Freq * HRForce

        Freq /= (2 * np.pi * input.c * input.au2fs)
        fout = open('../result_folder/HRfactor.dat', 'w')
        fout.writelines(' freq/cm^-1       Qd            HRd      lambdad/cm^-1      Qf            HRf      lambdaf/cm^-1\n')
        for i in np.arange(len(Freq)):
            fout.writelines('{:10.4f}{:14.7f}{:14.7f}{:14.7f}{:14.7f}{:14.7f}{:14.7f}\n'.format(Freq[i],
                            QDisp[i], HRDisp[i], LambdaDisp[i] / (2 * np.pi * input.c * input.au2fs),
                            QForce[i], HRForce[i], LambdaForce[i] / (2 * np.pi * input.c * input.au2fs)))
        fout.close()

        results['lambda10'] = float(Lambda10)
        results['lambda01'] = float(Lambda01)
        results['lambda4p'] = float((Lambda01 + Lambda10) * 0.5)
        results['lambdaDisp'] = float(np.sum(LambdaDisp) * input.au2eV)
        results['lambdaForce'] = float(np.sum(LambdaForce) * input.au2eV)
    
    elif ('force' in input.Properties):
        Freq,QForce = LambdaForceCal()
        HRForce = 0.5 * QForce * QForce
        LambdaForce = Freq * HRForce

        Freq /= (2 * np.pi * input.c * input.au2fs)
        fout = open('../result_folder/HRfactor.dat', 'w')
        fout.writelines(' freq/cm^-1      Qf            HRf      lambdaf/cm^-1\n')
        for i in np.arange(len(Freq)):
            fout.writelines('{:10.4f}{:14.7f}{:14.7f}{:14.7f}\n'.format(Freq[i],
                            QForce[i], HRForce[i], LambdaForce[i] / (2 * np.pi * input.c * input.au2fs)))
        fout.close()
        
        results['lambdaForce'] = float(np.sum(LambdaForce) * input.au2eV)
    
    elif ('4p' in input.Properties):
        Lambda10,Lambda01,Freq,QDisp = Lambda4pCal()
        Lambda10 *= input.au2eV
        Lambda01 *= input.au2eV
        HRDisp = 0.5 * QDisp * QDisp
        LambdaDisp = Freq * HRDisp

        Freq /= (2 * np.pi * input.c * input.au2fs)
        fout = open('../result_folder/HRfactor.dat', 'w')
        fout.writelines(' freq/cm^-1       Qd            HRd      lambdad/cm^-1\n')
        for i in np.arange(len(Freq)):
            fout.writelines('{:10.4f}{:14.7f}{:14.7f}{:14.7f}\n'.format(Freq[i],
                            QDisp[i], HRDisp[i], LambdaDisp[i] / (2 * np.pi * input.c * input.au2fs)))
        fout.close()
        
        results['lambda10'] = float(Lambda10)
        results['lambda01'] = float(Lambda01)
        results['lambda4p'] = float((Lambda01 + Lambda10) * 0.5)
        results['lambdaDisp'] = float(np.sum(LambdaDisp) * input.au2eV)
    
    else:
        print('no property is calculated')
        exit()
    
    os.system('mv S0opt* ../result_folder/')
    os.system('mv S0freq* ../result_folder/')
    if ('force' in input.Properties):
        os.system('mv S1force* ../result_folder/')
    if ('4p' in input.Properties):
        os.system('mv S1opt* ../result_folder/')

    fout = open('../result_folder/result.yml', 'wt')
    yaml.dump(results, fout)
    fout.close()