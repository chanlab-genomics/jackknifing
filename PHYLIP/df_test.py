#!/usr/bin/env python3

__author__ = 'Michael Ciccotosto-Camp'
__version__ = ''

import os
import sys
import pandas as pd

new_df = pd.DataFrame({'A': [1.0, .22352, 3.0], 'B': [4.0, 5.0, 6.0]})

# Write the remaining matrix, omit the column (header) names
new_df.to_csv(sys.stdout, sep='\t', header=False,
              index=True, index_label=False, float_format="%.4f")
