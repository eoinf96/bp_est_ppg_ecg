'''
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
import numpy as np
from sklearn.model_selection import LeaveOneGroupOut, StratifiedKFold
from sklearn.preprocessing import StandardScaler
import pickle
import joblib
import copy

class DatasetHandler:
    '''
    Class to handle dataset preprocessing dand manipulation for BP estimation regression.
    '''

    def __init__(self, df, feature_names, target_name):
        '''
        Initializes the dataset handler object.

        Args:
        - df (pandas.DataFrame): A dataframe of feature values.
        - feature_names (list of str): A list of names of the columns in the dataframe that contain the feature values.
        - target_name (str): The name of the column in the dataframe that contains the target values.
        '''
        self.df = df
        self.feature_names = feature_names
        self.target_name = target_name

    def __eq__(self, other):
        '''
        Compares this dataset handler object to another dataset handler object to check if they are equal.

        Args:
        - other (DatasetHandler): Another dataset handler object.

        Raises:
        - ValueError: If the two dataset handler objects are not equal.
        '''
        if len(self.df.columns.intersection(other.df.columns)) != self.df.shape[1]:
            raise ValueError('The two dataframes are not the same.')
        if self.feature_names != other.feature_names:
            raise ValueError('The two feature_names are not the same.')
        if self.target_name != other.target_name:
            raise ValueError('The two target_name are not the same.')

    def return_X_Y(self, to_numpy=True, relevant_IDs=None):
        '''
        Returns the feature values and target values as numpy arrays.

        Args:
        - to_numpy (bool): Whether to convert the feature values and target values to numpy arrays.
        - relevant_IDs (int or list of int): The ID or IDs of the rows to return. If None, returns all rows.

        Returns:
        - X (pandas.DataFrame or numpy.ndarray): The feature values.
        - y (pandas.Series or numpy.ndarray): The target values.
        '''
        if relevant_IDs is None:
            relevant_IDs = self.df['ID'].unique()
        if isinstance(relevant_IDs, np.integer) or isinstance(relevant_IDs, int):
            relevant_IDs = np.array([relevant_IDs])

        X = self.df.loc[self.df['ID'].isin(relevant_IDs)][self.feature_names]
        y = self.df.loc[self.df['ID'].isin(relevant_IDs)][self.target_name]

        if to_numpy:
            X = X.to_numpy()
            y = y.to_numpy()

        return X, y

    def drop_collinear_features(self, VIF_max=10):
        '''
        Removes collinear features from the dataset.

        Args:
        - VIF_max (float): The maximum acceptable variance inflation factor. Features with VIFs greater than this value will be removed.

        Returns:
        - None
        '''
        (X, _) = self.return_X_Y(to_numpy=True)
        dropping = True

        # Center X
        X = StandardScaler().fit_transform(X)

        while dropping:
            R0 = np.corrcoef(X.T)
            VIF = np.diag(np.linalg.inv(R0))
            if max(VIF) > VIF_max:
                ii = np.argmax(VIF)
                X = np.delete(X, ii, 1)
                self.feature_names.remove(self.feature_names[ii])
            else:
                dropping = False



class BP_est_dataset:
    """
    This class implements a dataset for blood pressure (BP) estimation.
    It includes functions to upload data, calibrate data, drop collinear features, and iterate through IDs.
    """

    def __init__(self, study=None, df=None, df_aug=None, BP_name='SBP', add_demographics=True, all_model_names=None):
        """
        Constructor for the BP_est_dataset class.

        :param study: A Study object
        :param df: A Pandas DataFrame containing the original dataset
        :param df_aug: A Pandas DataFrame containing the augmented dataset
        :param BP_name: A string representing the name of the BP column (default: 'SBP')
        :param add_demographics: A boolean indicating whether to add demographics (default: True)
        :param all_model_names: A list of strings containing all the model names to use (default: ['LASSO', 'RF'])
        """

        # Disable pylint warnings for too many instance attributes and too many arguments
        # pylint: disable=too-many-instance-attributes
        # pylint: disable=too-many-arguments
        self.BP_name = BP_name
        self.demo_added = add_demographics

        # Add df and dataset
        self.Original = None
        self.Augmented = None
        self.waveform_feature_names = None
        if study is not None:
            self.upload_study(study=study)
            if df is not None:
                if df_aug is None:
                    df_aug = df
                self.upload_df(df=df, df_aug=df_aug)

        self.ID_i = None

        if all_model_names is None:
            all_model_names = ['LASSO', 'RF']
        self.results = {}
        for model_name in all_model_names:
            self.results[model_name] = model_results(model_name=model_name)

    def upload_df(self, df, df_aug):
        """
        This function uploads the original and augmented datasets.

        :param df: A Pandas DataFrame containing the original dataset
        :param df_aug: A Pandas DataFrame containing the augmented dataset
        """

        df_attr_name = {0: 'Original', 1: 'Augmented'}
        for (ii, curr_df) in enumerate([df, df_aug]):
            if self.demo_added:
                # df_demo = self.dataset.get_demographics_df()
                df_demo = self.study.df.copy()
                df_demo['Sex'] = df_demo['Sex'].map({'F':1, 'M':0})
                # Inner join the df with df_demo
                curr_df = pd.merge(curr_df, df_demo, on='ID').dropna()
            else:
                curr_df = curr_df.dropna()

            feature_names = [x for x in list(curr_df.columns) if x not in ['BP', 'SBP', 'MAP', 'DBP', 'ID', 'calibration']]
            curr_df = DatasetHandler(df=curr_df, feature_names=feature_names, target_name=self.BP_name)
            setattr(self, df_attr_name[ii], curr_df)
        self.feature_names = self.Original.feature_names
        self.waveform_feature_names = [x for x in feature_names if x not in self.study.get_demographics_names()]

        # Check that the two datasets are for the same thing
        self.Original == self.Augmented

    def upload_study(self, study):
        """
        This function uploads the study.

        :param study: A Study object
        """
        self.study = study

    def calibrate_dataset(self, std_flag=False, num_cali=5):
        """
        Calibrates the dataset by applying normalization to the waveform and BP features.

        Args:
            - std_flag (bool): Whether to standardize the features by dividing by their standard deviation instead of their mean.
            - num_cali (int): Number of calibration readings to use per subject. If None, all readings for calibration will be used.

        Returns:
            - None
        """

        for (ii, curr_df) in enumerate([self.Original.df, self.Augmented.df]):
            if self.demo_added:
                demo_names = self.study.get_demographics_names()

                # Drop demo for now and then reapply them
                df_demo = curr_df.reindex(columns=['ID'] + demo_names)

                # Select non-demo features
                not_demo = [feat for feat in list(curr_df.columns) if feat not in demo_names]
                curr_df = curr_df.reindex(columns=not_demo)

            unique_ids = curr_df['ID'].unique()

            # Add calibration column if it doesn't exist
            if 'calibration' not in list(curr_df.columns):
                if num_cali is not None:
                    curr_df['calibration'] = np.zeros((len(curr_df),))
                    for id in unique_ids:
                        curr_df['calibration'][np.where(curr_df['ID'] == id)[0][:num_cali]] = 1
                else:
                    raise ValueError('Unable to get the calibration indices. Pass a calibration column or num_cali')

            # Calculate means and standard deviations for calibration readings
            means = curr_df.loc[curr_df['calibration'] == 1].groupby("ID", as_index=False).mean()
            means = means.reindex(means.index.repeat(curr_df.groupby('ID')['SBP'].count())).reset_index(drop=True)
            stds = curr_df.loc[curr_df['calibration'] == 1].groupby("ID").std()
            stds = stds.reindex(stds.index.repeat(curr_df.groupby('ID')['SBP'].count()))

            # Replace 0 with small positive number to avoid division by zero
            means = means.replace(0, np.finfo(float).eps)
            stds = stds.replace(0, np.finfo(float).eps)

            # Normalize features by means or standard deviations
            if std_flag:
                curr_df[self.waveform_feature_names] = (curr_df[self.waveform_feature_names] - means[
                    self.waveform_feature_names]) / stds[self.waveform_feature_names]
                curr_df[self.BP_name] = (curr_df[self.BP_name] - means[self.BP_name]) / stds[self.BP_name]
            else:
                curr_df[self.waveform_feature_names] = (curr_df[self.waveform_feature_names] - means[
                    self.waveform_feature_names]) / means[self.waveform_feature_names]
                curr_df[self.BP_name] = curr_df[self.BP_name] - means[self.BP_name]

            # Reapply demographic features
            if self.demo_added:
                curr_df = pd.merge(curr_df, df_demo, on='ID')

            # Update Original or Augmented dataset
            if ii == 0:
                self.Original.df = curr_df
            else:
                self.Augmented.df = curr_df

    def drop_collinear_features(self, VIF_max=10):
        '''
        Removes collinear features from the total dataset.
        Args:
        - VIF_max (float): The maximum acceptable variance inflation factor. Features with VIFs greater than this value will be removed.

        '''
        self.Original.drop_collinear_features(VIF_max=VIF_max)
        self.feature_names = self.Original.feature_names
        self.Augmented.feature_names = self.Original.feature_names

        # iterator over unique IDs in dataset
        def ID_iterator(self, val=None, restart=False):
            """
            Iterator over unique IDs in the dataset.

            Args:
            - val (int): Value to set the iterator to. If None, will increment the iterator to the next ID in the list.
            - restart (bool): Whether to restart the iterator from the beginning.

            Returns:
            - None
            """
            if restart:
                self.ID_i = None
            if val is None or restart:
                unique_ids = self.Original.df['ID'].unique()
                if self.ID_i is None:
                    self.ID_i = unique_ids[0]
                else:
                    ii = np.where(unique_ids == self.ID_i)[0][0]
                    self.ID_i = unique_ids[ii + 1]
            else:
                self.ID_i = val

    # iterator over unique IDs in dataset
    def ID_iterator(self, val=None, restart=False):
        """
        Iterator over unique IDs in the dataset.

        Args:
        - val (int): Value to set the iterator to. If None, will increment the iterator to the next ID in the list.
        - restart (bool): Whether to restart the iterator from the beginning.

        Returns:
        - None
        """
        if restart:
            self.ID_i = None
        if val is None or restart:
            unique_ids = self.Original.df['ID'].unique()
            if self.ID_i is None:
                self.ID_i = unique_ids[0]
            else:
                ii = np.where(unique_ids == self.ID_i)[0][0]
                self.ID_i = unique_ids[ii + 1]
        else:
            self.ID_i = val

    # return train-test data for leave-one-subject-out (LOSO) cross-validation
    def return_LOSO_train_test_data(self):
        """
        Return train-test data for leave-one-subject-out (LOSO) cross-validation.

        Args:
        - None

        Returns:
        - tuple: (X_train, X_test, y_train, y_test) - numpy arrays of training and testing data.
        """
        unique_ids = self.Original.df['ID'].unique()
        (X_train, y_train) = self.Augmented.return_X_Y(to_numpy=True, relevant_IDs=np.delete(unique_ids,
                                                                                             np.argwhere(
                                                                                                 unique_ids == self.ID_i)))
        (X_test, y_test) = self.Original.return_X_Y(to_numpy=True, relevant_IDs=self.ID_i)
        return (X_train, X_test, y_train, y_test)

    # return inner cross-validation folds
    def return_inner_cv(self, eval_method='LOSO', nFolds=5):
        """
        Return inner cross-validation folds.

        Args:
        - eval_method (str): Cross-validation method to use. Either 'LOSO' or 'CV' (k-fold cross-validation).
        - nFolds (int): Number of folds to use in k-fold cross-validation.

        Returns:
        - list: List of indices to split the data into train/test sets.
        """
        unique_ids = self.Original.df['ID'].unique()
        (X, y) = self.Augmented.return_X_Y(to_numpy=True,
                                           relevant_IDs=np.delete(unique_ids, np.argwhere(unique_ids == self.ID_i)))

        if eval_method == 'LOSO':
            cv = LeaveOneGroupOut().split(
                X=X,
                y=y,
                groups=self.Augmented.df['ID'][(self.Augmented.df['ID'] != self.ID_i)]
            )
        elif eval_method == 'CV':
            cv = StratifiedKFold(n_splits=nFolds).split(
                X=X,
                y=y,
                groups=self.Augmented.df['ID'][(self.Augmented.df['ID'] != self.ID_i)]
            )
        return cv

    # pickle the class object to a file
    def pickle_class(self, file_name):
        """
        Pickle the class object to a file.

        Args:
        - file_name (str): Name of the file to save the object to.

        Returns:
        - None
        """
        with open(file_name, 'wb') as f:
        # joblib.dump(self, file_name)
            pickle.dump(self, f, 2)
            # joblib.dump(self, f)
            f.close()

    @classmethod
    def load_pickled_class(cls, file_name):
        """
        Loads a pickled class from a file.

        Args:
            - file_name (str): The name of the file.

        Returns:
            - The loaded class object.
        """
        with open(file_name, 'rb') as f:
            return pickle.load(f)
            # return joblib.load(f)


class model_results():
    def __init__(self, model_name):
        self.model_name = model_name
        self.correlation = []
        self.RMSE = []
        self.MAE = []

        self.mdl_store = []
        self.data_store = {
            'X_test_store': [],
            'y_test_store': [],
            'X_train_store': [],
            'y_train_store': []
        }

    def update_results(self, mdl, X_test, y_test, X_train=None, y_train=None):
        y_est = mdl.predict(X_test)

        self.correlation.append(np.corrcoef(y_est, y_test)[0, 1])
        self.RMSE.append(np.sqrt(np.mean((y_est - y_test) ** 2)))
        self.MAE.append(np.mean(np.abs((y_est - y_test))))

        # Store model and data for later reference
        # Have to remove any CV generators so that pickle works
        if hasattr(mdl, 'cv'):
            mdl.cv = None
        self.mdl_store.append(mdl)
        self.data_store['X_test_store'].append(X_test)
        self.data_store['y_test_store'].append(y_test)
        self.data_store['X_train_store'].append(X_train)
        self.data_store['y_train_store'].append(y_train)


