'''
Class to handle regression
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

import pandas as pd
from joblib import Parallel, delayed
from statsmodels.stats.outliers_influence import variance_inflation_factor
from funcs.regression_funcs import *
import numpy as np
import matplotlib.pyplot as plt
from pickle_functions import pickle_load
from sklearn.model_selection import LeaveOneGroupOut

###################################################################
### class to perform BP estimation regression
###################################################################
class RegressionModel:
    '''
      Arguments
          dataset -- dataset class containing information on the train and test population
          df -- dataframe of feature values
          model_name -- name of the regression model: Linear, RF and SVM supported
          eval_method -- Set to 'CV' to perform Leave One Subject Out Cross Validation. Else Train Test split is used
          add_demographics -- Flag on whether to add demographics to the dataset. 1 for yes. 0 for no
      '''

    # pylint: disable=too-many-instance-attributes
    # pylint: disable=too-many-arguments
    def __init__(
        self,
        dataset=None,
        df=None,
        model_name=None,
        eval_method=None,
        add_demographics=True,
    ):
        self.model_name = model_name
        self.demo_added = add_demographics
        self.eval_method = eval_method  # Eval method can be train/test or CV
        self.model = model_name

        # Add df and dataset
        #This defines self.df, self.feature_names, self.train_ids, self.test_ids
        if dataset is not None:
            self.upload_dataset(dataset=dataset)
            if df is not None:
                self.upload_df(df=df)


        # Performance metrics that we want to track
        self.rmse = []
        self.correlation_coefficient = []
        self.r_squared = []

        #Arrays that store train and test info
        self.train_y = None
        self.train_x = None
        self.test_y = None
        self.test_x = None

        #Estimated Y to be compared to test_y
        self.y_est = []

        # Hidden variables
        self._msk = None
        self._mean_norm_flag = False
        self._z_norm_flag = False

    ##############
    # Model Initialisation
    ##############
    def upload_df(self, df):
        '''Upload dataframe containing features here.'''
        if self.demo_added:
            df_demo = self.dataset.get_demographics_df()
            # Inner join the df with df_demo
            self.df = pd.merge(df, df_demo, left_index=True, right_index=True)
        else:
            self.df = df

        self.df = self.df.dropna()
        self.feature_names = list(self.df.drop(["BP"], axis=1).columns)

    def upload_dataset(self, dataset):
        '''Upload dataset containing train and test ids and demographics.'''
        self.dataset = dataset
        self.train_ids = dataset.TrainID
        self.test_ids = dataset.TestID

    def mean_norm_df(self):
        '''Subtract norm per feature per subject'''
        _norm_df(self, std_flag=0)
        self._mean_norm_flag = True

    def z_norm_df(self, do_calibration=False, num_cali=5):
        '''Z- normalise per feature per subject'''
        _norm_df(self, std_flag=1, do_calibration=do_calibration, num_cali=num_cali)
        self._z_norm_flag = True

    def generate_loso(self):
        '''Initialise Leave One Subject Out Generator using the dataset -- used if CV is eval_method'''
        loso = LeaveOneGroupOut().split(
            X=StandardScaler().fit_transform(self.train_x),
            y=self.train_y,
            groups=self.train_groups,
        )
        self._loso = loso  # loso is required for the nested CV in hyperparameter tuning

    ##############
    # Run model
    ##############
    def run_model(
        self,
        vars=None,
        run_hyper_param_tuning=True,
        do_feat_optimisation=True,
        num_features=np.inf,
        plot_flag=False,
        return_fig=False,
    ):
        ''' Run regression models for BP estimation using features

        Inputs:
        Vars -- hyperparameters for models if neccesary -> if None then default will be used
        run_hyper_param_tuning -- tune hyperparameters
        do_feature_optimisation -- run LASSO feature selection
        num_features -- Maximum number of features to accept by LASSO feature selection
        plot_flag -- Show model responses
        return_fig -- for saving

        '''
        if self.eval_method.upper() == 'CV':
            # Run Leave one subject out cross validation
            if self.model_name is None:
                crossvalidation_funcs.LOSOCV_all_models(
                    self,
                    vars=vars,
                    do_feat_optimisation=do_feat_optimisation,
                    run_hyper_param_tuning=run_hyper_param_tuning,
                )
                fig = 0
            else:
                fig = crossvalidation_funcs.LOSOCV(
                    self,
                    model_vars=vars,
                    num_features=num_features,
                    do_feat_selection=do_feat_optimisation,
                    run_hyper_param_tuning=run_hyper_param_tuning,
                    plot_flag=plot_flag,
                    return_fig=False,
                )
        else: #Else we use train test split of dataset
            #Set train and test values
            create_train_test_data(mdl=self)

            if self.model_name.upper() == "LINEAR":
                linear_regression(self, method="LASSO")
            elif self.model_name.upper() == "RF":
                rf_regression(
                    self, model_vars=vars, run_hyper_param_tuning=run_hyper_param_tuning
                )
            elif self.model_name.upper() == "SVM":
                svm_regression(
                    self, model_vars=vars, run_hyper_param_tuning=run_hyper_param_tuning
                )
            else:
                raise ValueError("Model name not known")

            if plot_flag:
                fig = self.plot_individual_response(return_fig=return_fig)


        return fig

    def save_model(self, file_name, folder_name="/"):
        ''' Save output of model '''
        from pickle_functions import pickle_save

        a = self.__dict__
        if '_loso' in a:
            a.pop("_loso")
        a.pop("model")  # As also generator
        pickle_save(
            folder_loc="/pickles/models/" + folder_name + "/",
            file_name=file_name,
            var=a,
        )

    ##############
    # Model statistics
    ##############
    def get_coefficient_of_determination(self):
        ''' Compute and return coefficient of determination '''
        if self.y_est.size == 0:
            ValueError('Run model first')
        y_dash = np.mean(self.test_y)
        sum_sqared_errors = np.sum((self.test_y - self.y_est) ** 2)
        total_sum_of_squares = np.sum((self.test_y - y_dash) ** 2)
        self.r_squared = 1 - sum_sqared_errors / total_sum_of_squares

        return self.r_squared

    def get_correlation_coefficient(self):
        ''' Compute and return correlation coefficient '''
        if self.y_est.size == 0:
            ValueError('Run model first')
        # Correlation coefficient between yest and testy
        self.correlation_coefficient = np.corrcoef(self.test_y, self.y_est)[0, 1]

        return self.correlation_coefficient


    ##############
    # Plots
    ##############
    def plot_correlation(self, return_fig = False):
        ''' Plot scatter of y test and y est'''
        fig = plt.figure()
        plt.scatter(self.test_y, self.y_est)
        plt.xlabel("BP cuff (mmHg)")
        plt.ylabel("BP est (mmHg)")
        plt.grid()

        return fig if return_fig else 0

    def plot_individual_response(self, return_fig=False):
        ''' Plot Individual response'''
        fig = plotting_functions.plot_individual_response(self, return_fig=return_fig)
        return fig

    def show_feature_contributions(self, return_fig=False):

        if self.eval_method == "CV":
            return 0
        else:
            if self.model_name.upper() == "LINEAR":
                feature_importance = abs(self.model.coef_)
            elif self.model_name.upper() == "RF":
                feature_importance = self.model.feature_importances_
            elif self.model_name.upper() == "SVM":
                feature_importance = abs(self.model.coef_)

        # feature_weights = (feature_weights - np.mean(feature_weights))/np.std(feature_weights)
        parameter_names = self.df.drop("BP", axis=1).columns

        plt.rc("axes", axisbelow=True)
        fig = plt.figure()
        plt.bar(parameter_names, feature_importance)
        plt.xticks(rotation="vertical")
        plt.grid()
        plt.ylabel("Feature importance")
        plt.xlabel("Feature")

        return fig if return_fig else 0


def _norm_df(self, std_flag=0, do_calibration=False, num_cali=5):
    df = self.df
    if self.demo_added:
        demo_names = self.dataset.get_demographics_names()
        # Drop demo for now and then re apply them
        # not_demo = np.setdiff1d(list(df.columns), demo_names)
        not_demo = [feat for feat in list(df.columns) if feat not in demo_names]
        df = df.reindex(columns=not_demo)
    if do_calibration:
        id_values = np.array(df.index).astype("int")
        unique_ids = np.unique(id_values)
        loc_cali = np.zeros((len(id_values),))
        for id in unique_ids:
            loc_cali[np.where(id_values == id)[0][:num_cali]] = 1
        a = df.loc[loc_cali == 1].groupby("ID").mean()
        means = a.reindex(a.index.repeat(df["BP"].groupby("ID").count()))
        a = df.loc[loc_cali == 1].groupby("ID").std()
        std = (
            a.reindex(a.index.repeat(df["BP"].groupby("ID").count())) if std_flag else 1
        )

    else:
        means = df.groupby("ID").transform("mean")
        std = df.groupby("ID").transform("std") if std_flag else 1

    df = (df - means) / std
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    if self.demo_added:
        self.upload_df(df=df)
    self._means = means
    if std_flag:
        self._std = std


# The mask indicates which datapoints are in the training set
def create_train_test_data(mdl):

    if mdl.dataset == None:
        raise ValueError("Need to load dataset!")
    mdl._msk = np.isin(
        mdl.df.index, mdl.dataset.get_volunteer_id_name(mdl.dataset.TrainID)
    )

    train_df = mdl.df[mdl._msk]
    test_df = mdl.df[~mdl._msk]

    mdl.train_y = train_df["BP"].values
    mdl.train_x = train_df.drop(columns="BP", axis=1).values
    mdl.test_y = test_df["BP"].values
    mdl.test_x = test_df.drop(columns="BP", axis=1).values

    # mdl.df_test_Y = test_df["BP"]

def drop_collinear_features(self, thresh=10, plot_flag=False):

    X = self.df
    if self.demo_added:
        demo_names = self.dataset.get_demographics_names()
        # Drop demo for now and then re apply them
        not_demo = np.setdiff1d(list(X.columns), demo_names)
        X = X.reindex(columns=not_demo)

    if plot_flag:
        print("Warning the plots can take a little while to run")
        axs = pd.plotting.scatter_matrix(X)
        n = len(X.columns)
        for x in range(n):
            for y in range(n):
                ax = axs[x, y]
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                ax.xaxis.label.set_rotation(90)
                ax.yaxis.label.set_rotation(0)
                ax.yaxis.labelpad = 50

    variables = [X.columns[i] for i in range(X.shape[1])]
    variables_drop = []
    dropped = True
    while dropped:
        dropped = False
        vif = Parallel(n_jobs=-1)(
            delayed(variance_inflation_factor)(X[variables].values, ix)
            for ix in range(len(variables))
        )

        maxloc = vif.index(max(vif))
        if max(vif) > thresh:
            variables_drop.append(variables[maxloc])
            variables.pop(maxloc)
            dropped = True

    if plot_flag:
        X = X.drop(variables_drop, axis=1)
        axs = pd.plotting.scatter_matrix(X)
        n = len(X.columns)
        for x in range(n):
            for y in range(n):
                ax = axs[x, y]
                ax.xaxis.set_ticks([])
                ax.yaxis.set_ticks([])
                ax.xaxis.label.set_rotation(90)
                ax.yaxis.label.set_rotation(0)
                ax.yaxis.labelpad = 50

    # Drop the high VIF varaibles
    self.df = self.df.drop(variables_drop, axis=1)

# This will run much faster than the collinear features
# for a dataframe with a large number of features
def drop_correlated_features(self, thresh=0.8):

    X = self.df
    if self.demo_added:
        demo_names = self.dataset.get_demographics_names()
        # Drop demo for now and then re apply them
        not_demo = np.setdiff1d(list(X.columns), demo_names)
        X = X.reindex(columns=not_demo)

    corr_matrix = X.corr().abs()
    upper = corr_matrix.where(
        np.triu(np.ones(corr_matrix.shape), k=1).astype(np.bool)
    )
    to_drop = [column for column in upper.columns if any(upper[column] > thresh)]
    self.df = self.df.drop(to_drop, axis=1)




if __name__ == "__main__":

    #### MOLLIE session cannot be provided as it contains
    #### healthy volunteer demographic data that is confidential
    MOLLIE = MOLLIE_session()

    for device in ["Stowood", "Cardioscreen", "Philips"]:

        df_ECG = pickle_load("pickles/ECG/df_ECG_set_delays_no_filter_" + device)
        df_PPG = pickle_load("pickles/PPG/df_PPG_set_delays_no_filter_" + device)
        df_PAT = pickle_load("pickles/PAT/df_PAT_set_delays_no_filter_" + device)

        df = pd.concat([df_ECG, df_PPG, df_PAT], axis=1, join="inner")
        df = df.drop("SBP", axis=1)
        df["BP"] = df_ECG.SBP

        mdl = RegressionModel(
            dataset=MOLLIE,
            df=df,
            eval_method="cv",
            model_name="SVM",
            add_demographics=True,
        )
        mdl.z_norm_df()
        mdl.run_model(
            run_hyper_param_tuning=True,
            plot_flag=False,
            num_features=30,
            do_feat_optimisation=False,
        )
        mdl.save_model(file_name="All_feat_delay_SVM_30_z_no_filter_" + device)


    plt.show()
