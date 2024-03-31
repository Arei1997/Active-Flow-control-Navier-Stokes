import pandas as pd
from matplotlib import pyplot as plt
%matplotlib inline
import seaborn as sns  
df = pd.read_csv('malefemale column 1', index_col=0)
df.head()
sns.countplot (x='Male', data=df)