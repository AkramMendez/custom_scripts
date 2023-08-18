import sys
import pandas as pd

if len(sys.argv) < 3:
    print("Usage: python script.py file1.narrowPeak file2.xls")
    sys.exit(1)

file1_path = sys.argv[1]
file2_path = sys.argv[2]

# Read the first file into a pandas DataFrame
df1 = pd.read_csv(file1_path, sep='\t')

# Read the second file (xls) into a pandas DataFrame
df2 = pd.read_csv(file2_path, sep='\t', header=None)

# Join the DataFrames based on the common column
result = pd.merge(df1, df2, left_on=df1.columns[3], right_on=df2.columns[9])

selected_columns = result.iloc[:,0:10]

# Print the result
print(selected_columns.to_csv(header=None,sep='\t',index=False))


