"""
Group of functions for handling feature selection algorithms

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
"""

from sklearn import linear_model
import numpy as np
from sklearn.preprocessing import StandardScaler
import copy


def LASSO_feature_selection(mdl, max_iter=30, num_features=np.inf):
    '''

    :param mdl: mdl of type RegressionModel
    :param max_iter: Maximum number of iterations for LASSO
    :param num_features: Maximum number of features with non zero coefficinets to allow (default is all)
    :return: dataframe
    '''

    #Initalise a LOSO generator for nested CV
    mdl.generate_loso()
    feat_model = linear_model.LassoCV(
        max_iter=max_iter, cv=mdl._loso
    )

    feat_model.fit(X=StandardScaler().fit_transform(mdl.train_x), y=mdl.train_y)

    #Return a dataframe containing only features with non zero coefficients.
    df = copy.deepcopy(mdl.df)
    df = df.drop(np.setdiff1d(list(df.columns), ["BP"]), axis=1)

    best_cont_features = np.argsort(abs(feat_model.coef_))[::-1]
    idx = -1
    # Store the features that are in the top 25 and top 5 of features
    # for idx in range(num_features):
    while True:
        idx += 1
        # Exit conditions
        if idx == len(best_cont_features):
            break
        if feat_model.coef_[best_cont_features[idx]] == 0.0:
            break
        if idx == num_features:
            break
        feature_name = mdl.feature_names[best_cont_features[idx]]
        df[feature_name] = mdl.df[feature_name]

    return df

