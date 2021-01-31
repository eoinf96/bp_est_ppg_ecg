from funcs.regression_funcs import linear_regression, RF_regression, SVM_regression, feature_selection
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import copy
from sklearn.model_selection import KFold


def LOSOCV(mdl, plot_flag=False, vars=None, num_features=np.inf, do_feat_optimisation= True,  run_hyper_param_tuning=False, return_fig=False):
    mdl.method = 'LOSOCV'
    mdl.fold_info ={}
    mdl.fold_info['model_params'] = {}
    train_RMSE_store = []
    test_RMSE_store = []
    R_squared_store = []
    correlation_store = []

    y_store = []
    yest_store = []

    if do_feat_optimisation:
        feature_selection.initialise_feature_selection(mdl=mdl)

    if plot_flag:
        # Initialise figure to store subplot of all LOSOCV
        n_rows = np.ceil(np.sqrt(len(mdl.MOLLIE.volunteer_list)))
        n_cols = np.floor(np.sqrt(len(mdl.MOLLIE.volunteer_list)))
        fig, ax = plt.subplots(nrows=int(n_rows), ncols=int(n_cols))
        plot_row_idx = 0
        plot_col_idx = 0

    for volunteer_num in mdl.MOLLIE.volunteer_list:
        mdl.train_Y, mdl.train_X, mdl.test_Y, mdl.test_X, mdl._msk, mdl.train_groups = train_test_split_LOSOCV(
            mdl.df, mdl.MOLLIE.get_volunteer_id_name(volunteer_num))

        if do_feat_optimisation:
            df = feature_selection.LASSO_feature_selection(mdl = mdl, max_iter = 30, num_features=num_features)
            mdl.train_Y, mdl.train_X, mdl.test_Y, mdl.test_X, mdl._msk, mdl.train_groups = train_test_split_LOSOCV(
                df, mdl.MOLLIE.get_volunteer_id_name(volunteer_num))
            # mdl.fold_info[volunteer_num]['feature_list_rank'] = (list(df.drop('BP', axis =1).columns))
            mdl.fold_info['feature_list_rank'][volunteer_num] = {df.columns.get_loc(feature_name) : feature_name  for feature_name in list(df.drop('BP', axis =1).columns)}

        if mdl.model_name is None:
            pass
            #This feature is added in order to allow for feature selection only analysis
        else:
            if mdl.model_name.upper() == 'LINEAR':
                linear_regression(mdl, method=None)
            elif mdl.model_name.upper() == 'RF':
                RF_regression(mdl, vars=vars, run_hyper_param_tuning=run_hyper_param_tuning)
                if run_hyper_param_tuning: mdl.fold_info['model_params'][volunteer_num] = mdl.model.best_params_
            elif mdl.model_name.upper() == 'SVM':
                SVM_regression(mdl, vars=vars, run_hyper_param_tuning=run_hyper_param_tuning)
                if run_hyper_param_tuning: mdl.fold_info['model_params'][volunteer_num] = mdl.model.best_params_
            else:
                raise ValueError('Model name not known')
        train_RMSE_store.append(mdl.train_RMSE)
        test_RMSE_store.append(mdl.test_RMSE)
        R_squared_store.append(mdl.get_coefficient_of_determination())
        correlation_store.append(mdl.get_correlation_coefficient())


        y_store.append(mdl.test_Y)
        yest_store.append(mdl.y_est)
        if plot_flag:
            axes = ax[plot_row_idx, plot_col_idx]
            axes.plot(mdl.test_Y, '-o', label='Ground Truth', markersize=2)
            axes.plot(mdl.y_est, '-o', label='Estimated', markersize=2)
            axes.set_xticklabels([])
            axes.set_title(mdl.MOLLIE.get_volunteer_id_name(volunteer_num))
            axes.grid()
            plot_col_idx = plot_col_idx + 1

            if plot_col_idx == n_cols:
                plot_col_idx = 0
                plot_row_idx = plot_row_idx + 1

            if plot_row_idx == 0 and plot_col_idx == 1:
                axes.legend()

    mdl.fold_info['RMSE'] = test_RMSE_store
    mdl.fold_info['Correlation'] = correlation_store
    mdl.fold_info['R_squared'] = R_squared_store
    mdl.train_RMSE = np.median(train_RMSE_store)
    mdl.test_RMSE = np.median(test_RMSE_store)
    mdl.R_squared = np.median(R_squared_store)
    mdl.correlation_coefficient = np.median(correlation_store)

    mdl.test_Y = np.hstack(y_store)
    mdl.y_est = np.hstack(yest_store)
    mdl._msk = np.ones((len(mdl.y_est),))

   #Store values
    if not return_fig:
        fig = 0

    return fig

#Function to get performance statistics for all Regression algorithms and save the results
# @todo just implement this function into LOSCV
def LOSOCV_all_models(mdl, vars=None,  run_hyper_param_tuning=False, do_feat_optimisation= True):

    if mdl.model_name is None:
        #Then loop through all models as weve done the LASSO anyway
        mdl.model_name = ['Linear', 'RF', 'SVM']
    else:
        mdl.model_name = [mdl.model_name]


    #Initialise fold info
    mdl.fold_info = {model_name: {} for model_name in mdl.model_name}
    for model_name in mdl.model_name:
        mdl.fold_info[model_name] = {Name: [] for Name in ['RMSE','Correlation','R_squared']}
        if do_feat_optimisation:
            mdl.fold_info[model_name]['model_params'] = {}
            mdl.fold_info[model_name]['feature_list_rank'] = {}



    for volunteer_num in mdl.MOLLIE.volunteer_list:
        mdl.train_Y, mdl.train_X, mdl.test_Y, mdl.test_X, mdl._msk, mdl.train_groups = train_test_split_LOSOCV(
            mdl.df, mdl.MOLLIE.get_volunteer_id_name(volunteer_num))

        if do_feat_optimisation:
            df = feature_selection.LASSO_feature_selection(mdl = mdl, max_iter = 30)
        else:
            df = copy.deepcopy(mdl.df)


        for model_name in mdl.model_name:
            mdl.train_Y, mdl.train_X, mdl.test_Y, mdl.test_X, mdl._msk, mdl.train_groups = train_test_split_LOSOCV(
                df, mdl.MOLLIE.get_volunteer_id_name(volunteer_num))
            if do_feat_optimisation:
                mdl.fold_info[model_name]['feature_list_rank'][volunteer_num] = {df.columns.get_loc(feature_name): feature_name for
                                                                     feature_name in list(df.drop('BP', axis=1).columns)}

            if model_name.upper() == 'LINEAR':
                linear_regression(mdl, method=None)
            elif model_name.upper() == 'RF':
                RF_regression(mdl, vars=vars, run_hyper_param_tuning=run_hyper_param_tuning)
                if run_hyper_param_tuning: mdl.fold_info[model_name]['model_params'][volunteer_num] = mdl.model.best_params_
            elif model_name.upper() == 'SVM':
                SVM_regression(mdl, vars=vars, run_hyper_param_tuning=run_hyper_param_tuning)
                if run_hyper_param_tuning: mdl.fold_info[model_name]['model_params'][volunteer_num] = mdl.model.best_params_
            else:
                raise ValueError('Model name not known')
            mdl.fold_info[model_name]['RMSE'].append(mdl.test_RMSE)
            mdl.fold_info[model_name]['Correlation'].append(mdl.get_coefficient_of_determination())
            mdl.fold_info[model_name]['R_squared'].append(mdl.get_correlation_coefficient())


def ten_fold_CV(mdl,  do_feat_optimisation= True, plot_flag= False, vars = vars, run_hyper_param_tuning = False, return_fig = False):
    mdl.method = '10FCV'

    kf = KFold(n_splits=5)
    X = mdl.df.drop(columns='BP').values
    Y= mdl.df['BP'].values

    train_RMSE_store = []
    train_correlation_store = []
    test_RMSE_store = []
    test_correlation_store = []
    R_squared_store = []
    correlation_store = []

    y_store = copy.copy(Y)
    yest_store = copy.copy(Y)

    for train_index, test_index in kf.split(X):
        mdl._msk = train_index
        # mdl.msk = test_index
        mdl.train_X, mdl.test_X = X[train_index], X[test_index]
        mdl.train_Y, mdl.test_Y = Y[train_index], Y[test_index]

        #Run the model on the fold and get the performance metrics
        if mdl.model_name.upper() == 'LINEAR':
            linear_regression(mdl)
        elif mdl.model_name.upper() == 'RF':
            RF_regression(mdl, run_hyper_param_tuning = run_hyper_param_tuning)
        elif mdl.model_name.upper() == 'SVM':
            SVM_regression(mdl, run_hyper_param_tuning=run_hyper_param_tuning)
        else:
            raise ValueError('Model name not known')

        train_RMSE_store.append(mdl.train_RMSE)
        train_correlation_store.append(mdl.train_corr)
        test_RMSE_store.append(mdl.test_RMSE)
        test_correlation_store.append(mdl.test_corr)
        R_squared_store.append(mdl.get_coefficient_of_determination())
        correlation_store.append(mdl.get_correlation_coefficient())

        y_store[test_index] = mdl.test_Y
        yest_store[test_index] = mdl.y_est

    mdl.train_RMSE = np.median(train_RMSE_store)
    mdl.train_corr = np.median(train_correlation_store)
    mdl.test_RMSE = np.median(test_RMSE_store)
    mdl.test_corr = np.median(test_correlation_store)
    mdl.R_squared = np.median(R_squared_store)
    mdl.correlation_coefficient = np.median(correlation_store)

    mdl.test_Y = y_store
    mdl.y_est = yest_store
    mdl._msk = np.ones((len(mdl.y_est),))

    if plot_flag:
        mdl._msk = np.zeros((len(Y),)) > 1
        fig = mdl.plot_individual_response(return_fig=return_fig)

        if return_fig:
            return fig





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