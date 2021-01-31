from funcs.regression_funcs import linear_regression, RF_regression, SVM_regression, feature_selection
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import copy
from sklearn.model_selection import KFold


def LOSOCV_num_features(mdl, vars=None, max_num_features=30,  run_hyper_param_tuning=False):

    if mdl.model_name is None:
        #Then loop through all models as weve done the LASSO anyway
        mdl.model_name = ['Linear', 'RF', 'SVM']
    else:
        mdl.model_name = [mdl.model_name]


    #Initialise fold info
    mdl.fold_info = {model_name: {} for model_name in mdl.model_name}
    for model_name in mdl.model_name:
        for num_feat in range(1,max_num_features):
            mdl.fold_info[model_name][num_feat] = {Name: [] for Name in ['RMSE','Correlation','R_squared']}
            mdl.fold_info[model_name][num_feat]['model_params'] = {}
            mdl.fold_info[model_name][num_feat]['feature_list'] = {}

    feature_selection.initialise_feature_selection(mdl=mdl)


    for volunteer_num in mdl.MOLLIE.volunteer_list:
        mdl.train_Y, mdl.train_X, mdl.test_Y, mdl.test_X, mdl._msk, mdl.train_groups = train_test_split_LOSOCV(
            mdl.df, mdl.MOLLIE.get_volunteer_id_name(volunteer_num))

        df = feature_selection.LASSO_feature_selection(mdl = mdl, max_iter = 30, num_features=max_num_features)

        for num_feat in range(1,max_num_features):
            df_subset = df.drop(list(df.columns)[(num_feat+1):], axis =1)

            for model_name in mdl.model_name:
                mdl.train_Y, mdl.train_X, mdl.test_Y, mdl.test_X, mdl._msk, mdl.train_groups = train_test_split_LOSOCV(
                    df_subset, mdl.MOLLIE.get_volunteer_id_name(volunteer_num))

                mdl.fold_info[model_name][num_feat]['feature_list'][volunteer_num] = (list(df_subset.drop('BP', axis =1).columns))


                if model_name.upper() == 'LINEAR':
                    linear_regression(mdl, method=None)
                elif model_name.upper() == 'RF':
                    RF_regression(mdl, vars=vars, run_hyper_param_tuning=run_hyper_param_tuning)
                    if run_hyper_param_tuning: mdl.fold_info[model_name][num_feat]['model_params'][volunteer_num] = mdl.model.best_params_
                elif model_name.upper() == 'SVM':
                    SVM_regression(mdl, vars=vars, run_hyper_param_tuning=run_hyper_param_tuning)
                    if run_hyper_param_tuning: mdl.fold_info[model_name][num_feat]['model_params'][volunteer_num] = mdl.model.best_params_
                else:
                    raise ValueError('Model name not known')
                mdl.fold_info[model_name][num_feat]['RMSE'].append(mdl.test_RMSE)
                mdl.fold_info[model_name][num_feat]['Correlation'].append(mdl.get_coefficient_of_determination())
                mdl.fold_info[model_name][num_feat]['R_squared'].append(mdl.get_correlation_coefficient())



#For effective RMSE estimation when features are z normalised, this needs to update the mask
def train_test_split_LOSOCV(df, volunteer_num):

    #make volunteer_id
    volunteer_id = '00' + str(volunteer_num)
    volunteer_id = volunteer_id[-3:]

    msk = ~np.isin(df.index, volunteer_id)

    train_df = df[msk]
    test_df = df[~msk]

    train_Y = train_df['BP'].values
    train_X = train_df.drop(columns='BP', axis=1).values

    train_groups = train_df.index.values.astype('int')

    test_Y = test_df['BP'].values
    test_X = test_df.drop(columns='BP', axis=1).values


    return train_Y, train_X, test_Y, test_X, msk, train_groups