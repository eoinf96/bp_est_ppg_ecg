import pandas as pd
from classes.MOLLIE_session_class import MOLLIE_session
from pickle_functions import *
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from classes.regression_class import Regression_model

MOLLIE = MOLLIE_session()

device = 'Cardioscreen'

df_ECG = pickle_load('pickles/ECG/df_ECG_set_delays_no_filter_' + device)
# df_ECG = pickle_load('pickles/ECG/df_ECG_' + device)

# df_PPG = pickle_load('pickles/PPG/df_PPG_set_delays_' + device)
df_PPG = pickle_load('pickles/PPG/df_PPG_set_delays_no_filter_' + device)
# df_PPG = pickle_load('pickles/PPG/df_PPG_' + device)

# df_PAT = pickle_load('pickles/PAT/df_PAT_set_delays_' + device)
df_PAT = pickle_load('pickles/PAT/df_PAT_set_delays_no_filter_' + device)
# df_PAT = pickle_load('pickles/PAT/df_PAT_' + device)

df = pd.concat([df_ECG, df_PPG, df_PAT], axis=1, join='inner')
df = df.drop('SBP', axis=1)
df['BP'] = df_ECG.SBP

mdl = Regression_model(MOLLIE = MOLLIE, df=df, method = 'Individual', eval_method='cv', model_name='SVM', add_demographics=True)
mdl.z_norm_df()
mdl.run_model(run_hyper_param_tuning = False, plot_flag=True, num_features = 30, do_feat_optimisation=False)
# mdl.save_model(file_name = 'All_feat_delay_SVM_30_z_no_filter_'+device)



plt.show()