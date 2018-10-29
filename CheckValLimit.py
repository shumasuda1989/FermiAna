#!/usr/bin/env python

import sys
import os
import xml.etree.ElementTree as ET
import numpy as np
import glob

RED = '\033[31m'
END = '\033[0m'

xml = sys.argv[1]

tree = ET.parse(xml)
root = tree.getroot()

for src in root:
    for param in src.findall("./spectrum/parameter[@free='1']"):
        print src.get('name'), param.get('name')
        val=float(param.get('value'))
        min=float(param.get('min'))
        max=float(param.get('max'))
        if val <= min:
            print RED+'  Fit value is below the lower limit.'+END
        elif val >=max:
            print RED+'  Fit value exceeds the upper limit.'+END

