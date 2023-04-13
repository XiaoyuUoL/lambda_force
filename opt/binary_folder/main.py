import computation
import input

# BOD calculation
if ('BOD' in input.Properties):
    computation.BODCalculation()
exit()

## S0opt and S0freq calculations
# S0opt
computation.QCCalculate('S0opt')

# S1freq
computation.QCCalculate('S0freq')
# check if S0 opt need to be rerun
while (computation.QCCalculate('CheckFreq') == -1):
    computation.QCCalculate('S0freq')

# S1force calculation
if ('force' in input.Properties):
    computation.QCCalculate('S1force')

# S1opt calculation
if ('4p' in input.Properties):
    computation.QCCalculate('S1opt')

# HR factor and lambda calculation
computation.HRCalculate()

# write results into result.yml file
computation.WriteOutput()