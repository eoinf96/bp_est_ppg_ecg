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

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

# To do -- Add augmented dataset. Add SHAP. Add feature ranking coefficient

import pandas as pd
import matplotlib.pyplot as plt
from classes import MOLLIE_session, BP_est_dataset


def load_df(data_loc='data'):
    df_PPG = pd.read_csv(data_loc+'/PPG_features.csv')
    df_ECG = pd.read_csv(data_loc+'/ECG_features.csv')
    df_PAT = pd.read_csv(data_loc+'/PAT_features.csv')

    df = pd.concat([df_PPG, df_ECG], axis=1, join="inner")
    df = pd.concat([df, df_PAT], axis=1, join="inner")
    df = df.loc[:, ~df.columns.duplicated()].copy()
    return df


if __name__ == "__main__":

    df = load_df()

    MOLLIE = MOLLIE_session()
    BP_est = BP_est_dataset(df=df, study=MOLLIE)
    BP_est.calibrate_dataset()
    BP_est.drop_collinear_features(VIF_max=10)

    # for ii in range(len(BP_est.df['ID'].unique())):
