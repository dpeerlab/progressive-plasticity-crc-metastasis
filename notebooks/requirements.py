'''
This file contains standardized imports and variables used throughout the 
analysis notebooks for "Progressive Plasticity in Colorectal Cancer Metastasis"
by Moorman et al (DOI: )
'''
# Packages used throughout
import scanpy as sc
import anndata
import numpy as np
import scipy as sp
import pandas as pd
from matplotlib import pyplot as plt
import os
import sys
import json
import tqdm

# Local modules used throughout
module_path = os.path.abspath('../src')
if module_path not in sys.path:
    sys.path.append(module_path)

# Standard variables referenced throughout
data_dir = f"{os.getcwd()}/../data"
media_dir = f"{os.getcwd()}/../media"

# Plotting styles used throughout
module_path = os.path.abspath('../src')
stylesheet = f'{module_path}/utils/pl/assets/default.mplstyle'
plt.style.use(stylesheet)
with open(f'{module_path}/utils/pl/assets/named_colors.json', 'r') as f:
    named_colors = json.load(f)