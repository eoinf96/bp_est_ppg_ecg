import pandas as pd
from joblib import Parallel, delayed
from statsmodels.stats.outliers_influence import variance_inflation_factor
from funcs import remove_empty_subplots
from funcs.regression_funcs import *
import numpy as np
from sklearn.model_selection import LeaveOneGroupOut

class Regression_model():
    def __init__(self, MOLLIE=None, df = None, model_name = None, method = None, eval_method = None, msk = None, add_demographics = True):
        self.model_name = model_name
        self.method = method
        self.demo_added = add_demographics
        self.eval_method = eval_method # Eval method can be train/test or CV
        self._msk= msk # needed so that we can ensure the same mask between tests for comparison
        self.model = model_name

        #Add df and MOLLIE
        if MOLLIE is not None:
            self.upload_MOLLIE_dataset(MOLLIE = MOLLIE)
            if df is not None:
                self.upload_df(df = df)

        self.train_RMSE = []
        self.train_corr = []
        self.test_RMSE = []
        self.test_corr = []

        self.train_Y = None
        self.train_X = None
        self.test_Y = None
        self.test_X = None

        self.y_est = []


        self._mean_norm_flag = False
        self._z_norm_flag = False


    def upload_df(self, df):
        if self.demo_added:
            df_demo = self.MOLLIE.get_demographics_df()
            # Inner join the df with df_demo
            self.df = pd.merge(df, df_demo, left_index=True, right_index=True)
        else:
            self.df = df

        self.df = self.df.dropna()
        self.feature_names = list(self.df.drop(['BP'], axis=1).columns)

    def set_mask(self, msk):
        self._msk = msk

    def upload_MOLLIE_dataset(self, MOLLIE):
        self.MOLLIE = MOLLIE
        self.train_ids = MOLLIE.TrainID
        self.test_ids = MOLLIE.TestID


    def mean_norm_df(self):
        _norm_df(self, std_flag=0)
        self._mean_norm_flag = True

    def z_norm_df(self, do_calibration = False, num_cali=5):
        _norm_df(self, std_flag=1, do_calibration=do_calibration, num_cali=num_cali)
        self._z_norm_flag = True

    def drop_collinear_features(self, thresh = 10, plot_flag = False):
        X = self.df
        if self.demo_added:
            demo_names = self.MOLLIE.get_demographics_names()
            # Drop demo for now and then re apply them
            not_demo = np.setdiff1d(list(X.columns), demo_names)
            X = X.reindex(columns=not_demo)

        if plot_flag:
            print('Warning the plots can take a little while to run')
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
                delayed(variance_inflation_factor)(X[variables].values, ix) for ix in range(len(variables)))

            maxloc = vif.index(max(vif))
            if max(vif) > thresh:
                variables_drop.append(variables[maxloc])
                variables.pop(maxloc)
                dropped = True

        if plot_flag:
            X = X.drop(variables_drop, axis = 1)
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


        #Drop the high VIF varaibles
        self.df = self.df.drop(variables_drop, axis = 1)

    #This will run much faster than the collinear features for a dataframe with a large number of features
    def drop_correlated_features(self, thresh = 0.8):

        X = self.df
        if self.demo_added:
            demo_names = self.MOLLIE.get_demographics_names()
            # Drop demo for now and then re apply them
            not_demo = np.setdiff1d(list(X.columns), demo_names)
            X = X.reindex(columns=not_demo)

        corr_matrix = X.corr().abs()
        upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(np.bool))
        to_drop = [column for column in upper.columns if any(upper[column] > thresh)]
        self.df = self.df.drop(to_drop, axis=1)

    def generate_loso(self):
        loso = LeaveOneGroupOut().split(X=StandardScaler().fit_transform(self.train_X), y=self.train_Y, groups=self.train_groups)
        self._loso = loso  # loso is required for the nested CV in hyperparameter tuning



    def run_model(self,vars =None, run_hyper_param_tuning = True, num_features = np.inf,
                  plot_flag= False, return_fig= False, do_feat_optimisation=True):
        if not self.eval_method.upper() == 'CV':
            create_train_test_data(self=self)

            if self.model_name.upper() == 'LINEAR':
                linear_regression(self, method= 'LASSO')
            elif self.model_name.upper() == 'RF':
                RF_regression(self, vars = vars, run_hyper_param_tuning = run_hyper_param_tuning)
            elif self.model_name.upper() =='SVM':
                SVM_regression(self,vars=vars,  run_hyper_param_tuning = run_hyper_param_tuning)
            else:
                raise ValueError('Model name not known')

            self.get_coefficient_of_determination()
            if plot_flag:
                fig = self.plot_individual_response(return_fig=return_fig)
        else:
            #Add in code here for Cross validation of each method
            if self.method is 'Individual':
                if self.model_name is None:
                    crossvalidation_funcs.LOSOCV_all_models(self, vars = vars,do_feat_optimisation=do_feat_optimisation,  run_hyper_param_tuning = run_hyper_param_tuning)
                    fig = 0
                else:
                    fig = crossvalidation_funcs.LOSOCV(self, vars = vars,num_features=num_features,do_feat_optimisation=do_feat_optimisation,  run_hyper_param_tuning = run_hyper_param_tuning,  plot_flag= plot_flag, return_fig= False)
            elif self.method is 'Mix':
                fig = crossvalidation_funcs.ten_fold_CV(self, do_feat_optimisation=do_feat_optimisation, plot_flag = plot_flag, return_fig= False)


        return fig

    def save_model(self, file_name, folder_name = '/'):
        from pickle_functions import pickle_save
        a = self.__dict__
        if hasattr(a, 'loso'):
            a.pop('_loso')
        a.pop('model') # As also generator
        pickle_save(folder_loc= '/pickles/models/'+folder_name + '/',file_name=file_name, var=a)


    def get_coefficient_of_determination(self):
        y_dash = np.mean(self.test_Y)
        sum_sqared_errors = np.sum((self.test_Y - self.y_est) ** 2)
        total_sum_of_squares = np.sum((self.test_Y - y_dash) ** 2)
        self.R_squared = 1 - sum_sqared_errors/total_sum_of_squares

        return self.R_squared


    def get_correlation_coefficient(self):
        #Correlation coefficient between yest and testy
        self.correlation_coefficient = np.corrcoef(self.test_Y , self.y_est)[0,1]

        return self.correlation_coefficient

    def plot_correlation(self):
        plt.figure()
        plt.scatter(self.test_Y, self.y_est)
        plt.xlabel('BP cuff (mmHg)')
        plt.ylabel('BP est (mmHg)')
        plt.grid()

    def plot_individual_response(self, return_fig=False):
        fig =  plotting_functions.plot_individual_response(self, return_fig=return_fig)
        return fig

    def plot_test_response(self, return_fig=False):
        fig = plotting_functions.plot_test_response(self, return_fig=return_fig)
        return fig


    def show_feature_contributions(self, return_fig = False):
        #Note that this only works if the model parameters have been normalised prior to being trained

        if self.eval_method =='CV':
            return 0
        else:
            if self.model_name.upper()  == 'LINEAR':
               feature_importance = abs(self.model.coef_)
            elif self.model_name.upper() == 'RF':
                feature_importance = self.model.feature_importances_
            elif self.model_name.upper() == 'SVM':
                feature_importance = abs(self.model.coef_)

        # feature_weights = (feature_weights - np.mean(feature_weights))/np.std(feature_weights)
        parameter_names = self.df.drop('BP', axis=1).columns

        plt.rc('axes', axisbelow=True)
        fig = plt.figure()
        plt.bar(parameter_names, feature_importance)
        plt.xticks(rotation='vertical')
        plt.grid()
        plt.ylabel('Feature importance')
        plt.xlabel('Feature')

        if return_fig: return fig



def _norm_df(self, std_flag = 0, do_calibration = False, num_cali=5):
    df = self.df
    if self.demo_added:
        demo_names = self.MOLLIE.get_demographics_names()
        # Drop demo for now and then re apply them
        # not_demo = np.setdiff1d(list(df.columns), demo_names)
        not_demo = [feat for feat in list(df.columns) if feat not in demo_names]
        df = df.reindex(columns=not_demo)
    if do_calibration:
        id_values =np.array(df.index).astype('int')
        unique_ids = np.unique(id_values)
        loc_cali = np.zeros((len(id_values),))
        for id in unique_ids:
            loc_cali[np.where(id_values == id)[0][:num_cali]] = 1
        a = df.loc[loc_cali == 1].groupby('ID').mean()
        means = a.reindex(a.index.repeat(df['BP'].groupby('ID').count()))
        a = df.loc[loc_cali == 1].groupby('ID').std()
        std = a.reindex(a.index.repeat(df['BP'].groupby('ID').count()))if std_flag else 1


    else:
        means = df.groupby('ID').transform('mean')
        std = df.groupby('ID').transform('std') if std_flag else 1

    df = (df - means) / std
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    if self.demo_added:
        self.upload_df(df=df)
    self._means = means
    if std_flag : self._std = std

def intialise_random_mask(df, frac_train = 0.5):
    #This requires df to be uploaded
    if df is None:
        raise ValueError('Need to load df!')
    msk = np.random.rand(len(df)) < frac_train
    return msk


#The mask indicates which datapoints are in the training set
def create_train_test_data(self):
    #Method can either be mix or individual
    if self.method.upper() == 'MIX':
        if self._msk is None:
            print('Initialising mask with 50/50 split')
            self._msk = intialise_random_mask(self.df)

    elif self.method.upper() == 'INDIVIDUAL':
        #Make sure that MOLLIE has been uploaded as this contains the train and test IDs
        if self.MOLLIE == None:
            raise ValueError('Need to load MOLLIE!')
        self._msk = np.isin(self.df.index, self.MOLLIE.get_volunteer_id_name(self.MOLLIE.TrainID))


    train_df = self.df[self._msk]
    test_df = self.df[~self._msk]


    self.train_Y = train_df['BP'].values
    self.train_X = train_df.drop(columns='BP', axis =1).values
    self.test_Y = test_df['BP'].values
    self.test_X = test_df.drop(columns='BP', axis =1).values

    self.df_test_Y = test_df['BP']
