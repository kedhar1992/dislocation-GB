# Mail id: kedharnath1992@gmail.com 
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
!pip install shap
import shap as shap

# Read the input data and select data for analysis 
df = pd.read_excel("data_ML.xlsx")
df.head()
pearcol = df.drop(df.columns[[2,3,4,5,6,10,12,13]], axis ='columns') 
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
xg_reg = xgb.XGBRegressor(objective ='reg:squarederror', use_rmm = 'True', 
                          max_depth = 3,
                          gamma = 85,
                          colsample_bytree = 0.191421, 
                          subsample = 0.284894,
                          learning_rate = 0.843006, 
                          max_delta_step = 67, 
                          random_state = 1,
                          n_estimators = 99,
                          alpha = 6,
                          seed = 19
                          )

xg_reg.fit(X_train,y_train)
preds = xg_reg.predict(X_test)     # predict y data using X_test 
# Low RMSE is desired 
rmse = np.sqrt(mean_squared_error(y_test, preds))  # calc RMSE bet y_test & y_predict 
print("RMSE: %f" % (rmse))

# Plotting 
xgb.plot_importance(xg_reg)
plt.savefig('feature imp.png',dpi = 300)
plt.rcParams['figure.figsize'] = [5,5]
plt.show()

# sort = xg_reg.feature_importances_.argsort()
# plt.barh(X.columns[sort], xg_reg.feature_importances_[sort])
# plt.xlabel("Feature Importance")

pred_ytrain = xg_reg.predict(X_train)
pred_ytest = xg_reg.predict(X_test)

plt.scatter(y_train, pred_ytrain, c='b', label='Training data', s = 70)
plt.scatter(y_test, pred_ytest, c='r', label='Test data', s = 70)
plt.plot([0, 1000], [0, 1000], color = 'black', linewidth = 3, linestyle='dashed')
plt.xlim(0, 400)
plt.ylim(0, 400)

plt.legend(loc="upper left", fontsize = 17, prop={'weight': 'bold'})
plt.xticks(weight = 'bold', fontsize= 13)
plt.yticks(weight = 'bold', fontsize= 13)
plt.rcParams["axes.linewidth"] = 1.5
plt.tick_params(direction='out', length=6, width=2, grid_alpha=0.5)


plt.xlabel('Calculated values', fontsize= 20, fontweight='bold')
plt.ylabel('Predicted values', fontsize= 20, fontweight='bold')
plt.savefig('scatter.png',dpi = 300, bbox_inches = "tight")

#############################################
shap.initjs()
explainer = shap.Explainer(xg_reg, X) # (xg_reg, X_test), (xg_reg, X_train), & (xg_reg) shows DIFFERENCE   
shap_values = explainer(X)   # if X_test is changed to X_train then shows DIFFERENCE; but mostly X_test is used https://www.kaggle.com/code/dansbecker/shap-values
feature_names = [ a + ": " + str(b) for a,b in zip(X.columns, np.abs(shap_values.values).mean(0).round(1))]
shap.summary_plot(shap_values, plot_type='violin', feature_names = feature_names, color_bar = False, show = False) 

cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15)
cbar.ax.tick_params(direction='out', length=6, width=2,  grid_alpha=0.5)

plt.xticks(weight = 'bold', fontsize= 13)
plt.yticks(weight = 'bold', fontsize= 13)
plt.rcParams["axes.linewidth"] = 1.5
plt.tick_params(direction='out', length=6, width=2, grid_alpha=0.5)
plt.xlabel('SHAP value', fontsize= 15, fontweight='bold')

plt.savefig("summary_plot.png", dpi = 300, bbox_inches = "tight")
plt.close()

shap.summary_plot(shap_values, plot_type='violin', feature_names=feature_names)

shap.plots.scatter(shap_values[:,"Slipplane"], color=shap_values[:,"mprime"], dot_size=40, cmap='rainbow') 
plt.xticks(weight = 'bold', fontsize= 13)
plt.yticks(weight = 'bold', fontsize= 13)
plt.rcParams["axes.linewidth"] = 1.5
plt.tick_params(direction='out', length=6, width=2, grid_alpha=0.5)
plt.savefig("Slipplane.png", dpi = 300, bbox_inches = "tight")
plt.close()

shap.plots.scatter(shap_values[:,"Slipplane"], color=shap_values[:,"Drop in shear stress"], dot_size=40, cmap='rainbow') 
plt.xticks(weight = 'bold', fontsize= 13)
plt.yticks(weight = 'bold', fontsize= 13)
plt.rcParams["axes.linewidth"] = 1.5
plt.tick_params(direction='out', length=6, width=2, grid_alpha=0.5)
plt.savefig("Slipplane.png", dpi = 300, bbox_inches = "tight")
plt.close()

shap.plots.scatter(shap_values[:,"Rotationaxis"], color=shap_values[:,"mprime"], dot_size=40, cmap='rainbow') 
plt.xticks(weight = 'bold', fontsize= 13)
plt.yticks(weight = 'bold', fontsize= 13)
plt.rcParams["axes.linewidth"] = 1.5
plt.tick_params(direction='out', length=6, width=2, grid_alpha=0.5)
plt.savefig("Rotationaxis.png", dpi = 300, bbox_inches = "tight")
plt.close()

shap.plots.scatter(shap_values[:,"Rotationaxis"], color=shap_values[:,"Drop in shear stress"], dot_size=40, cmap='rainbow') 
plt.xticks(weight = 'bold', fontsize= 13)
plt.yticks(weight = 'bold', fontsize= 13)
plt.rcParams["axes.linewidth"] = 1.5
plt.tick_params(direction='out', length=6, width=2, grid_alpha=0.5)
plt.savefig("Rotationaxis.png", dpi = 300, bbox_inches = "tight")
plt.close()
