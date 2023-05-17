# kedharnath1992@gmail.com
# Kindly cite if you use the script ""

#############################################
# Import necessary modules
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt
import xgboost as xgb
from sklearn.metrics import mean_squared_error
import shap as shap

# Read the input data and select data for analysis 
df = pd.read_excel("data_ML.xlsx")
df.head()
X = df.drop(df.columns[[2,3,4,5,6,8,10,12,13, 16,17,18,20]], axis ='columns')   # Drop in stress = 7; yield stress = 8; critical distance = 9; Drop in PE = 11
# Rotation axis = 0; slip plane = 1
y = df.Yieldstress
list(X.columns)

# Separate training and testing data sets
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=10)   # 25% data for testing
list(X_test)

#############################################
cormat = pearcol.corr()
round(cormat,2)
fig, ax = plt.subplots(figsize=(15,15))         # Sample figsize in inches
ax.figure.axes[-1].yaxis.label.set_size(20)

sns.set(font_scale=2)
plt.xticks(weight = 'bold', fontsize= 15)
plt.yticks(weight = 'bold', fontsize= 15)
plt.rcParams["axes.linewidth"] = 1.5
sns.heatmap(cormat, annot=True, linewidth=.5, ax = ax, fmt=".2f", cmap='coolwarm', cbar_kws={'shrink': 0.8}, annot_kws={'size': 15}, square=True);
for text in ax.texts:
    text.set_size(14)
    if text.get_text() >= '0.7':
        text.set_size(18)
        text.set_weight('bold')
        text.set_style('italic')
plt.savefig('pearson.png',dpi = 300, bbox_inches = "tight")

#############################################




