import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
df = pd.read_csv("~/MondalLab/NCC_RNAseq/BGI.diff.expr.tables/case1.1_vs_case1.3.xls",sep="\t")
df.head()
df.dtypes
df.columns.to_list()
plt.scatter(x=df.iloc[:,3],y=df.iloc[:,4])
plt.plot
plt.show()

from typing_extensions import TypeAlias