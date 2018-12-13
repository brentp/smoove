import sys
import os
import numpy as np
import json
import pandas as pd
import tabulate


js = [json.load(open(f)) for f in sys.argv[1:]]

methods = [os.path.dirname(f).split("/")[1] for f in sys.argv[1:]]

df = pd.DataFrame.from_records(js)
df['method'] = methods
df['FDR'] = 1 - df['precision']

sdf = df['method FDR FN   FP  TP-call precision recall f1'.split()]
m = pd.melt(sdf, id_vars=['method'], value_vars=['recall', 'FDR'])

print sdf.to_csv(sys.stdout, index=False, sep="\t", float_format="%.3f")

sdf = sdf.drop('f1', axis=1)
sdf['recall-%'] = 100.0 * sdf['recall'] / float(sdf['recall'][0])
sdf['FP-%'] = 100.0 * sdf['FP'] / float(sdf['FP'][0])
print(tabulate.tabulate(sdf.values, sdf.columns, tablefmt="pipe", floatfmt=".3f"))
