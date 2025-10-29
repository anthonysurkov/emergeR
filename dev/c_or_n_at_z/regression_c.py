import pandas as pd
import statsmodels.api as sm

df = pd.read_csv("r270x_cstatus.csv")
X = df['c_status'].to_numpy().reshape(-1,1).astype(float)
X = sm.add_constant(X)

df_h = df[(df['mle'] > 0.30)]
X_h = df_h['c_status'].to_numpy().reshape(-1,1).astype(float)
X_h = sm.add_constant(X_h)

y_mle = df['mle']
y_map = df['map']
y_high = df_h['mle']

model_mle = sm.OLS(y_mle, X).fit()
model_map = sm.OLS(y_map, X).fit()
model_high = sm.OLS(y_high, X_h).fit()

print("Regression over MLE:")
print(model_mle.summary())
print()
print("Regression over MAP:")
print(model_map.summary())
print()
print("Regression over MLE for winner N9s:")
print(model_high.summary())
print()

