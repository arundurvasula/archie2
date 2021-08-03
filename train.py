import statsmodels.api as sm
import pandas as pd

training_data = pd.read_csv('train.txt', sep='\t', header=0)

X = training_data[['Dist.Mean','Dist.Var','Dist.Skew','Dist.Kurtosis','Min.Dist','Num.Private']]
y = training_data['Archaic']

log_reg = sm.Logit(y, X).fit(maxiter=1000, method='nm') #note, if Y doesn't contain any archaic samples, you will get errors
