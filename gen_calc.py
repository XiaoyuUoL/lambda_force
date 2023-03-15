#!/usr/bin/env python3
import yaml

sample_calc = {
    "calculator_id": "XTB_example",
    "specifications": { 
        "numsteps": 200,
    },
    "provides": [
      "HOMO",
      "LUMO",
    ]
}
with open("xtb_example_homolumo.yml",'wt') as outfile:
    yaml.dump(sample_calc, outfile)

inputfile = open("xtb_example_homolumo.yml")
print(yaml.load(inputfile, Loader=yaml.BaseLoader))
inputfile.close()