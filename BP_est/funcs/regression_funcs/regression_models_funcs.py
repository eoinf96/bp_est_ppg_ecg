"""

Group of functions to run and evaluate regression models for BP estimation.
Contains functions for Linear Regression, RF regression and SVM regression


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

import numpy as np
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import RandomizedSearchCV
from sklearn import linear_model


SVM_CONST_VARS = {
    "kernel": "linear",
    "C": 10,
    "epsilon": 0.1,
}

RF_CONST_VARS = {
    "n_estimators": 128,
    "max_features": "sqrt",
    "max_depth": 10,
    "bootstrap": True,
}


def svm_regression(mdl, model_vars=SVM_CONST_VARS, run_hyper_param_tuning=True):
    """

    :param mdl: mdl of type RegressionMOdel
    :param model_vars: Model hyperparamters
    :param run_hyper_param_tuning: Optimise run hyperparamter optimisation
    """

    if model_vars is None:
        # Use default values
        model_vars = SVM_CONST_VARS
    sc = StandardScaler()
    train_df_x = sc.fit_transform(mdl.train_x)
    test_df_x = sc.transform(mdl.test_x)

    # Set up model
    mdl.model = SVR(
        kernel=model_vars["kernel"], C=model_vars["C"], epsilon=model_vars["epsilon"]
    )

    if run_hyper_param_tuning:
        print("Running Hyperparameter tuning ... ")
        # kernels = ['linear', 'RBF', 'poly']
        kernels = ["linear"]
        # degree = [2, 3]
        # gammas =  [1e-4, 1e-3, 0.01, 0.1, 0.2, 0.5, 0.6, 0.9]
        C = [0.1, 1, 10]
        epsilon = np.logspace(-3, 3, num=7)

        grid = {"kernel": kernels, "C": C, "epsilon": epsilon}
        run_hyperparamter_optimising(grid, mdl)

    mdl.model.fit(X=train_df_x, y=mdl.train_y)

    if run_hyper_param_tuning:
        print(mdl.model.best_params_)

    # Evaluate model and save model statistics
    evaluate_model(mdl, train_df_x, test_df_x)


def rf_regression(mdl, model_vars=RF_CONST_VARS, run_hyper_param_tuning=True):
    """

    :param mdl: mdl of type RegressionMOdel
    :param model_vars: Model hyperparamters
    :param run_hyper_param_tuning: Optimise run hyperparamter optimisation
    """
    if model_vars is None:
        # Use default values
        model_vars = RF_CONST_VARS

    sc = StandardScaler()
    train_df_x = sc.fit_transform(mdl.train_x)
    test_df_x = sc.transform(mdl.test_x)

    mdl.model = RandomForestRegressor(
        n_estimators=model_vars["n_estimators"],
        max_features=model_vars["max_features"],
        max_depth=model_vars["max_depth"],
        bootstrap=model_vars["bootstrap"],
    )

    if run_hyper_param_tuning:
        print("Running Hyperparameter tuning ... ")
        # Number of trees in random forest -- Dont need to tune-- just set to a large enough number. There is no chance of overfitting based on the number of trees
        # n_estimators = [int(x) for x in np.linspace(start=500, stop=1000, num=2)]
        n_estimators = [128]
        # n_estimators = 128
        # Number of features to consider at every split
        max_features = ["auto"]
        # Maximum number of levels in tree -- This one is more important
        max_depth = [int(x) for x in np.linspace(5, 60, num=3)]
        max_depth.append(None)
        # Method of selecting samples for training each tree
        # bootstrap = [True, False]
        bootstrap = [True]
        # Create the random grid
        random_grid = {
            "n_estimators": n_estimators,
            "max_features": max_features,
            "max_depth": max_depth,
            "bootstrap": bootstrap,
        }

        run_hyperparamter_optimising(random_grid, mdl)

    mdl.model.fit(X=train_df_x, y=mdl.train_y)

    if run_hyper_param_tuning:
        print(mdl.model.best_params_)

    evaluate_model(mdl, train_df_x, test_df_x)


# Run hyperparameter tuning here
def run_hyperparamter_optimising(grid, mdl):
    """
    :param grid: dict defining the range of variables to optimise over
    :param mdl: mld of type RegressionModel
    :return:
    """
    mdl.generate_loso()
    mdl.model = RandomizedSearchCV(
        estimator=mdl.model,
        param_distributions=grid,
        n_iter=10,
        cv=mdl._loso,
        verbose=0,
        random_state=42,
    )


def linear_regression(mdl, method=None):
    """
    :param mdl: mdl of type RegressionModel
    :param method: regualrisation method to run
    """
    if method is not None:
        mdl.model_name = mdl.model_name + "_" + method

    sc = StandardScaler()
    train_df_x = sc.fit_transform(mdl.train_x)
    test_df_x = sc.transform(mdl.test_x)

    if method is None:
        mdl.model = linear_model.LinearRegression()
    elif method == "LASSO":
        mdl.model = linear_model.LassoCV(max_iter=100)
    elif method == "RIDGE":
        mdl.model = linear_model.RidgeCV(max_iter=100)
    elif method == "ELASTIC":
        mdl.model = linear_model.ElasticNetCV()
    else:
        ValueError("Unknown Linear method")

    mdl.model.fit(X=train_df_x, y=mdl.train_y)

    evaluate_model(mdl, train_df_x, test_df_x)


def evaluate_model(
    mdl,
    train_df_x,
    test_df_x,
):
    """
    Computes model statistics after regregression model is run
    :param mdl:  mdl of class RegressionModel
    :param train_df_x: Training features to evaluate model with
    :param test_df_x: Testing features to evaluate model with
    """

    sigma = np.array(mdl._std["BP"][mdl._msk]) if mdl._z_norm_flag else 1
    y_est_train = mdl.model.predict(train_df_x) * sigma
    mdl.train_y *= sigma
    mdl.train_rmse = np.sqrt(np.mean((y_est_train - mdl.train_y * sigma) ** 2))
    mdl.train_corr = np.corrcoef(y_est_train, mdl.train_y)[0, 1]

    sigma = np.array(mdl._std["BP"][~mdl._msk]) if mdl._z_norm_flag else 1
    mdl.y_est = mdl.model.predict(X=test_df_x) * sigma
    mdl.test_y *= sigma
    mdl.rmse = np.sqrt(np.mean((mdl.y_est - mdl.test_y) ** 2))

    mdl.correlation_coefficient = mdl.get_correlation_coefficient()
    mdl.r_squared = mdl.get_coefficient_of_determination()
