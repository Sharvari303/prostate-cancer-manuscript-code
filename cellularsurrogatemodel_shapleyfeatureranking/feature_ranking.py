#conda create Shapley_env 
#conda activate Shapley_env

#!pip install shap
#conda install -c conda-forge shap
import shap
#from tabulate import tabulate
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
#import torch.nn.functional as F
from scipy import stats
from sklearn.model_selection import train_test_split
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import balanced_accuracy_score
from sklearn.preprocessing import StandardScaler
#pip install imbalanced-learn
from imblearn.over_sampling import SMOTE
from sklearn.metrics import balanced_accuracy_score
import tensorflow as tf
from tensorflow import keras
from keras import Sequential
from keras.layers import Dense
from sklearn.metrics import balanced_accuracy_score

#Collating all data of boolean, erk-akt, pten,egf for each cohort (df_total has data for all cohorts)

#FIRST INPUT CSV -> boolean + ncg data
df_boolean_and_ncg_ip = pd.read_csv('~/Shapley_PCa/boolean_input_output_default_run_Testo0.csv')
df_boolean_and_ncg_ip.head() #this csv has time series data for boolean+ncg value at all iterations
result_df = df_boolean_and_ncg_ip[df_boolean_and_ncg_ip['step'] == 'SS'].reset_index(drop=True)
#result_df is a df size 1920 x 26 -> after exractting SS rows for each iteration
#print(result_df.head())
value_counts=result_df['cohort'].value_counts()
# Convert the result to a DataFrame for better presentation
count_df = pd.DataFrame({'Cohort': value_counts.index, 'Count': value_counts.values})
#count_df - we have 640 data points for each of the three cohorts - CNT, BR, TF
#print(count_df)

#we now extract 640 points in each cohort to a separate df
df_boolean_and_ncg_CNT= result_df[result_df['cohort']=='CNT'].reset_index(drop=True)
df_boolean_and_ncg_BR= result_df[result_df['cohort']=='BR'].reset_index(drop=True)
df_boolean_and_ncg_TF= result_df[result_df['cohort']=='TF'].reset_index(drop=True)
df_boolean_and_ncg_TF.head()

#SECOND INPUT CSV -> PTEN, EGF analog data
df_pten_egf_ip=pd.read_csv('~/Shapley_PCa/bool_input_threshold_case_spec_modified.csv')
#df_pten_egf_ip.head()

##we now extract 640 points in each cohort to a separate df
df_pten_egf_CNT= df_pten_egf_ip[df_pten_egf_ip['cohort']=='CNT'].reset_index(drop=True)
df_pten_egf_BR= df_pten_egf_ip[df_pten_egf_ip['cohort']=='BR'].reset_index(drop=True)
df_pten_egf_TF= df_pten_egf_ip[df_pten_egf_ip['cohort']=='TF'].reset_index(drop=True)
df_pten_egf_TF.head()

#THIRD INPUT CSV
df_erkakt_ip= pd.read_csv('~/Shapley_PCa/pp_SS_boolean_input_threshold_0.05_validation_v1.csv')
df_erkakt_ip.head()

#we now extract 640 points in each cohort to a separate df
df_erkakt_CNT= df_erkakt_ip[df_erkakt_ip['cohort']=='CNT'].reset_index(drop=True)
df_erkakt_BR= df_erkakt_ip[df_erkakt_ip['cohort']=='BR'].reset_index(drop=True)
df_erkakt_TF= df_erkakt_ip[df_erkakt_ip['cohort']=='TF'].reset_index(drop=True)
df_erkakt_BR.head()

#now combine entries for all inputs in a single df (one for each cohort)
df_CNT = pd.concat([df_boolean_and_ncg_CNT,df_erkakt_CNT[['ERK_PP_CONC','AKT_PP_CONC']],df_pten_egf_CNT[['PTEN','EGF_nM']] ], axis=1)
df_BR = pd.concat([df_boolean_and_ncg_BR,df_erkakt_BR[['ERK_PP_CONC','AKT_PP_CONC']],df_pten_egf_BR[['PTEN','EGF_nM']] ], axis=1)
df_TF = pd.concat([df_boolean_and_ncg_TF,df_erkakt_TF[['ERK_PP_CONC','AKT_PP_CONC']],df_pten_egf_TF[['PTEN','EGF_nM']] ], axis=1)
#print(df_TF)

# Concatenate row-wise
df_total = pd.concat([df_CNT, df_BR, df_TF], axis=0)

# Reset index
df_total.reset_index(drop=True, inplace=True)
#print(df_total)

# defining inputs as X, binary output as y -> output decided using set threshold
threshold = 0.15
X = df_total[["AKT_PP","TP53","ERK_PP","E2F1","PTEN","WIP1","CDKN1A","ATM","BAX","CASP9","MDM4","CDK1","CDK2","AR","RB1","ARF","Raf_P","MDM2","CCNG1","BCL2","ERK_PP_CONC","AKT_PP_CONC"]] #,"PTEN", "EGF_nM"]] #[["ERK_PP_CONC", "AKT_PP_CONC"]] #,"PTEN_analog", "EGF"
# Add an output column to binary based on the ncg threshold
df_total['binary_class'] = (df_total['NCG'] > threshold).astype(int)
y = df_total['binary_class']
#storing all input featires in a list 
features = X.columns.values
#print(features)

#NN Model

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
# Feature scaling
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Apply SMOTE to the training set only
smote = SMOTE(random_state=42)
X_train_resampled, y_train_resampled = smote.fit_resample(X_train_scaled, y_train)

nn_model = Sequential([
            Dense(10, input_shape=(X_train.shape[1],), activation='relu'),
            Dense(20, activation='relu'),
            Dense(40, activation='relu'),
            Dense(20, activation='relu'),
            Dense(1, activation='sigmoid')]) # No activation function here for logits (raw scores)
nn_model.compile(optimizer='adam', loss='binary_crossentropy',metrics=['accuracy'])

nn_model.fit(X_train_resampled, y_train_resampled, epochs=10)

#SHAPLEY ANALYSIS

explainer = shap.KernelExplainer(nn_model.predict,X_train_resampled,silent=True)
shap_values = explainer.shap_values(X_test_scaled,nsamples=500)
print("shap_values")
shap.summary_plot(shap_values,X_test_scaled,feature_names=features,plot_type='bar')

plt.savefig('PCa_NN_shap_summaryplot.png')

"""SVM"""

from sklearn.svm import SVC
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
# Feature scaling
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Apply SMOTE to the training set only
smote = SMOTE(random_state=42)
X_train_resampled, y_train_resampled = smote.fit_resample(X_train_scaled, y_train)

svm_model = SVC(kernel='rbf', C=10, gamma='scale')

svm_model.fit(X_train_resampled, y_train_resampled)

"""SVM- Shapley"""

explainer = shap.KernelExplainer(svm_model.predict,X_train_resampled,silent=True)

shap_values = explainer.shap_values(X_test_scaled,nsamples=500)
print("shap_values")
shap.summary_plot(shap_values,X_test_scaled,feature_names=features, plot_type='bar')
plt.savefig('PCa_SVM_shap_summaryplot.png')



"""Random Forest"""

from sklearn.ensemble import RandomForestClassifier
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
# Feature scaling
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Apply SMOTE to the training set only
smote = SMOTE(random_state=42)
X_train_resampled, y_train_resampled = smote.fit_resample(X_train_scaled, y_train)

rf_model = RandomForestClassifier(n_estimators=200, random_state=42)

rf_model.fit(X_train_resampled, y_train_resampled)

"""RF - Shapley Analysis"""

explainer = shap.KernelExplainer(rf_model.predict,X_train_resampled,silent=True)

shap_values = explainer.shap_values(X_test_scaled,nsamples=500)
print("shap_values")
shap.summary_plot(shap_values,X_test_scaled,feature_names=features, plot_type='bar')

plt.savefig('PCa_RF_shap_summaryplot.png')