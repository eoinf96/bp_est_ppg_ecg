'''
This script runs the regression model for the features of the ECG and PPG
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
import matplotlib.pyplot as plt
# from classes import MOLLIE_session, RegressionModel

if __name__ == "__main__":

    #### MOLLIE session cannot be provided as it contains
    #### healthy volunteer demographic data that is confidential
    # MOLLIE = MOLLIE_session()


    df_PPG = pd.read_csv('data/PPG_features.csv')
    df_ECG = pd.read_csv('data/ECG_features.csv')
    df_PAT = pd.read_csv('data/PAT_features.csv')

    df = pd.concat([df_ECG, df_PPG], axis=1, join="inner")
    df = pd.concat([df, df_PAT], axis=1, join="inner")
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


plt.show()
