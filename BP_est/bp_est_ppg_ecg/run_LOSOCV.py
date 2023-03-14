'''
This script runs the regression model for the features of the ECG and PPG
--
 Released under the GNU General Public License
 Copyright (C) 2022  Eoin Finnegan
 eoin.finnegan@eng.ox.ac.uk

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
w
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import pandas as pd
from bp_est_ppg_ecg import MOLLIE_session, BP_est_dataset
import regression_model_funcs
from sklearn.preprocessing import StandardScaler

def load_df(data_loc='../data', is_aug = False):
    suffix = '.csv' if is_aug else '_aug.csv'

    df_PPG = pd.read_csv(data_loc+'/PPG_features'+suffix)
    df_ECG = pd.read_csv(data_loc+'/ECG_features'+suffix)
    df_PAT = pd.read_csv(data_loc+'/PAT_features'+suffix)

    df = pd.concat([df_PPG, df_ECG], axis=1, join="inner")
    df = pd.concat([df, df_PAT], axis=1, join="inner")
    df = df.loc[:, ~df.columns.duplicated()].copy()
    return df


if __name__ == "__main__":


    model_name = 'LASSO'

    df = load_df()
    df_aug = load_df(is_aug=True)

    BP_name = 'SBP'

    MOLLIE = MOLLIE_session()
    BP_est = BP_est_dataset(df=df, df_aug=df_aug, study=MOLLIE, BP_name=BP_name)
    BP_est.calibrate_dataset()
    BP_est.drop_collinear_features(VIF_max=10)

    BP_est.ID_iterator(restart=True)

    ### PERFORM LOSOCV

    for ii in range(len(BP_est.Original.df['ID'].unique())):
        print(ii)
        BP_est.ID_iterator(val=ii + 1)
        (X_train, X_test, y_train, y_test) = BP_est.return_LOSO_train_test_data()
        ### Standardise
        sc = StandardScaler()
        X_train = sc.fit_transform(X_train)
        X_test = sc.transform(X_test)
        ### LASSO + OLS
        cv = BP_est.return_inner_cv()
        mdl = regression_model_funcs.linear_regression(X_train=X_train, y_train=y_train, method='LASSO', cv=cv)
        BP_est.results['LASSO'].update_results(mdl, X_test=X_test, y_test=y_test, X_train=X_train, y_train=y_train)

        # ### RF
        cv = BP_est.return_inner_cv()
        mdl = regression_model_funcs.RF_regression(X_train=X_train, y_train=y_train, run_hyper_param_tuning=True, cv=cv)
        BP_est.results['RF'].update_results(mdl, X_test=X_test, y_test=y_test, X_train=X_train, y_train=y_train)

    BP_est.pickle_class(file_name='../data/BP_est_results')

