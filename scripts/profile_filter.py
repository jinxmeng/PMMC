#!/share/data1/software/miniconda3/envs/jinxin/bin/python3

import sys
import pandas as pd

profile = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
profile = profile.loc[(profile!=0).any(axis=1), (profile!=0).any(axis=0)]
profile.to_csv(sys.argv[2], sep="\t", index=True)
