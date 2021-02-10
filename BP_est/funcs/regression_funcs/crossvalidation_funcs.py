'''

Group of functions to handle nested cross validation regression.
Can implement 10 fold CV, LOSOCV (Example shown)
--
 Released under the GNU General Public License
 Copyright (C) 2021  Eoin Finnegan
 eoin.finnegan@eng.ox.ac.uk

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''


from funcs.regression_funcs import (
    linear_regression,
    RF_regression,
    SVM_regression,
    feature_selection,
)
import numpy as np
import matplotlib.pyplot as plt
def LOSOCV(
    mdl,
    run_hyper_param_tuning=False,
    model_vars=None,
    do_feat_selection=True,
    num_features=np.inf,
    plot_flag=False,
    return_fig=False,
):
    '''
      Arguments
        mdl - regression model of type RegressionModel
        run_hyper_param_tuning -- boolean to indicate whether to run hyperparameter tuning
        model_vars -- dict containing the hyperparameters of the regression model  (None if run_hyper_param_tuning = True)
        do_feat_optimisation -- boolean to indicate whether to run feature selection algorithm
        num_features -- The maximum number of features for the feature selection algorithm
        plot_flag -- boolean to plot data or not
        return_flag -- boolean to return figure for saving
      '''

    #Initialise variables
    mdl.fold_info = {}
    if do_feat_selection:
        mdl.fold_info["model_params"] = {}
        mdl.fold_info["feature_list_rank"] ={}
    rmse_store = []
    r_squared_store = []
    correlation_store = []
    y_store = []
    yest_store = []

    if plot_flag:
        # Initialise figure to store subplot of all LOSOCV
        n_rows = np.ceil(np.sqrt(len(mdl.MOLLIE.volunteer_list)))
        n_cols = np.ceil(np.sqrt(len(mdl.MOLLIE.volunteer_list)))
        fig, ax = plt.subplots(nrows=int(n_rows), ncols=int(n_cols))
        plot_row_idx = 0
        plot_col_idx = 0
    else:
        return_fig = False

    # LOSOCV loop -- loop through all subjects in the dataset
    for volunteer_num in mdl.MOLLIE.volunteer_list:
        #Split dataset based on which volunteer is in the out of bag sample
        (
            mdl.train_y,
            mdl.train_x,
            mdl.test_y,
            mdl.test_x,
            mdl._msk,
            mdl.train_groups,
        ) = train_test_split_LOSOCV(
            mdl.df, mdl.MOLLIE.get_volunteer_id_name(volunteer_num)
        )

        #Run LASSO feature selection algorithm based
        if do_feat_selection:
            df = feature_selection.LASSO_feature_selection(
                mdl=mdl, max_iter=30, num_features=num_features
            )
            (
                mdl.train_y,
                mdl.train_x,
                mdl.test_y,
                mdl.test_x,
                mdl._msk,
                mdl.train_groups,
            ) = train_test_split_LOSOCV(
                df, mdl.MOLLIE.get_volunteer_id_name(volunteer_num)
            )
            mdl.fold_info["feature_list_rank"][volunteer_num] = {
                df.columns.get_loc(feature_name): feature_name
                for feature_name in list(df.drop("BP", axis=1).columns)
            }

        # Run regression model on selected data
        if mdl.model_name.upper() == "LINEAR":
            linear_regression(mdl, method=None)
        elif mdl.model_name.upper() == "RF":
            RF_regression(mdl, model_vars=model_vars, run_hyper_param_tuning=run_hyper_param_tuning)
            if run_hyper_param_tuning:
                mdl.fold_info["model_params"][volunteer_num] = mdl.model.best_params_
        elif mdl.model_name.upper() == "SVM":
            SVM_regression(
                mdl, model_vars=model_vars, run_hyper_param_tuning=run_hyper_param_tuning
            )
            if run_hyper_param_tuning:
                mdl.fold_info["model_params"][volunteer_num] = mdl.model.best_params_
        else:
            raise ValueError("Model name not known")


        ##### Save the results of the regression
        rmse_store.append(mdl.rmse)
        r_squared_store.append(mdl.r_squared)
        correlation_store.append(mdl.correlation_coefficient)

        y_store.append(mdl.test_y)
        yest_store.append(mdl.y_est)



        if plot_flag:
            axes = ax[plot_row_idx, plot_col_idx]
            axes.plot(mdl.test_y, "-o", label="Ground Truth", markersize=2)
            axes.plot(mdl.y_est, "-o", label="Estimated", markersize=2)
            axes.set_xticklabels([])
            axes.set_title(mdl.MOLLIE.get_volunteer_id_name(volunteer_num))
            axes.grid()
            plot_col_idx = plot_col_idx + 1

            if plot_col_idx == n_cols:
                plot_col_idx = 0
                plot_row_idx = plot_row_idx + 1

            if plot_row_idx == 0 and plot_col_idx == 1:
                axes.legend()

    mdl.fold_info["RMSE"] = rmse_store
    mdl.fold_info["Correlation"] = correlation_store
    mdl.fold_info["R_squared"] = r_squared_store
    mdl.rmse = np.median(rmse_store)
    mdl.r_squared = np.median(r_squared_store)
    mdl.correlation_coefficient = np.median(correlation_store)

    mdl.test_y = np.hstack(y_store)
    mdl.y_est = np.hstack(yest_store)

    return fig if return_fig else 0


def train_test_split_LOSOCV(df, volunteer_num):
    '''
    Function generate train and test data based on which volunteer makes up the out of bag sample

    Arguments
        df -- dataframe being used for regression
        volunteer_num --  index of volunteer in test set.
    '''

    # make volunteer_id
    volunteer_id = "00" + str(volunteer_num)
    volunteer_id = volunteer_id[-3:]

    msk = ~np.isin(df.index, volunteer_id)
    train_df = df[msk]
    test_df = df[~msk]

    train_y = train_df["BP"].values
    train_x = train_df.drop(columns="BP", axis=1).values

    train_groups = train_df.index.values.astype("int")

    test_y = test_df["BP"].values
    test_x = test_df.drop(columns="BP", axis=1).values
    return train_x, train_x, test_y, test_x, msk, train_groups
