import os
import sys

import computation
import input
import prepare

WorkDir = sys.argv[1]

fin = open('{}/molecules.dat'.format(WorkDir), 'r')
lines = fin.readlines()
SystemNames = []
SystemSmiles = []
for line in lines:
    data = line.rstrip().split()
    SystemNames.append(data[-1])
    SystemSmiles.append(data[0])

fout = open('{}/run.sh'.format(WorkDir), 'w')
for name,smiles in zip(SystemNames,SystemSmiles):
    #print(name)
    input.SystemName = name
    input.SystemSmiles = smiles
    #print(name)
    input.FilePath = WorkDir + input.SystemName
    if (not os.path.isdir(input.FilePath)):
        print('no workdir')
        exit()
        os.makedirs(input.FilePath)
    
    if (len(sys.argv) != 3):
        print('usage: python main.py [WorkDir] 0/1/2')
        print('0: prepare S0opt files based on molecules.dat file')
        print('1: prepare S0freq files based on S0 geometry')
        print('2: check if S0 opt need to be rerun  based on S0freq results')
        print('3: prepare S1force/S1opt files based on S0 geometry')
        print('4: calculation HR factors')
        exit()

    if (sys.argv[2] == '0'):        # prepare S0opt files
        prepare.Prepare1()
        fout.writelines('cd {}\n'.format(name))
        fout.writelines('sbatch S0opt.sh\n')
        fout.writelines('cd ..\n')
    elif (sys.argv[2] == '1'):      # prepare S0freq files
        prepare.Prepare2()
        fout.writelines('cd {}\n'.format(name))
        fout.writelines('sbatch S0freq.sh\n')
        fout.writelines('cd ..\n')
    elif (sys.argv[2] == '2'):      # check if S0 opt need to be rerun (if not, then S1force)
        if (prepare.CheckFreq() == -1):
        #if (False):
            print(name)
            fout.writelines('cd {}\n'.format(name))
            fout.writelines('sbatch S0opt.sh\n')
            fout.writelines('cd ..\n')
        else:
            prepare.Prepare3()      # prepare S1force files
            fout.writelines('cd {}\n'.format(name))
            fout.writelines('sbatch S1force.sh\n')
            fout.writelines('cd ..\n')
    elif (sys.argv[2] == '3'):      # prepare S1opt files
        prepare.Prepare4()
        fout.writelines('cd {}\n'.format(name))
        fout.writelines('sbatch S1opt.sh\n')
        fout.writelines('cd ..\n')
    elif (sys.argv[2] == '4'):      # calculate HR factors
        computation.HRcalculate()
    else:
        print('usage: python main.py [WorkDir] 0/1/2/3/4')
        print('0: prepare S0opt files based on molecules.dat file')
        print('1: prepare S0freq files based on S0 geometry')
        print('2: check if S0 opt need to be rerun  based on S0freq results')
        print('3: prepare S1force/S1opt files based on S0 geometry')
        print('4: calculation HR factors')
        exit()
fout.close()
