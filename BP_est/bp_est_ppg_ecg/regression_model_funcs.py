"""
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
from sklearn.model_selection import RandomizedSearchCV
from sklearn import linear_model

# Constants for SVM and Random Forest regressors
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


def SVM_regression(X_train, y_train, model_vars=None, run_hyper_param_tuning=True, cv=None):
    """
    Train a Support Vector Machine regression model.

    Args:
        X_train (array-like): Training data features.
        y_train (array-like): Training data target values.
        model_vars (dict, optional): Dictionary of model hyperparameters.
        run_hyper_param_tuning (bool, optional): Whether to run hyperparameter tuning.
        cv (int, cross-validation generator or an iterable, optional): Determines the cross-validation splitting strategy.

    Returns:
        A trained SVR model.
    """
    if model_vars is None:
        # Use default values
        model_vars = SVM_CONST_VARS

    # Set up model
    mdl = SVR(
        kernel=model_vars["kernel"], C=model_vars["C"], epsilon=model_vars["epsilon"]
    )

    if run_hyper_param_tuning:
        print("Running Hyperparameter tuning ... ")
        kernels = ['linear', 'RBF', 'poly']
        # kernels = ["linear"]
        degree = [2, 3]
        # gammas =  [1e-4, 1e-3, 0.01, 0.1, 0.2, 0.5, 0.6, 0.9]
        C = [0.1, 1, 10]
        epsilon = np.logspace(-3, 3, num=7)

        random_grid = {"kernel": kernels, "C": C, "epsilon": epsilon}
        run_hyperparamter_optimising(random_grid, mdl, cv)

    mdl.fit(X=X_train, y=y_train)
    return mdl


def RF_regression(X_train, y_train, model_vars=None, run_hyper_param_tuning=True, cv=None):
    """
    Train a Random Forest regression model.

    Args:
        X_train (array-like): Training data features.
        y_train (array-like): Training data target values.
        model_vars (dict, optional): Dictionary of model hyperparameters.
        run_hyper_param_tuning (bool, optional): Whether to run hyperparameter tuning.
        cv (int, cross-validation generator or an iterable, optional): Determines the cross-validation splitting strategy.

    Returns:
        A trained RandomForestRegressor model.
    """
    if model_vars is None:
        # Use default values
        model_vars = RF_CONST_VARS

    # Will need to have a check here that all values are accounted for
    mdl = RandomForestRegressor(
        n_estimators=model_vars["n_estimators"],
        max_features=model_vars["max_features"],
        max_depth=model_vars["max_depth"],
        bootstrap=model_vars["bootstrap"],
    )

    if run_hyper_param_tuning:
        print("Running Hyperparameter tuning ... ")
        # n_estimators = [int(x) for x in np.linspace(start=500, stop=1000, num=2)]
        n_estimators = [300]
        # Number of features to consider at every split
        num_input_features = X_train.shape[1]
        max_features = [int(x) for x in np.linspace(start=1, stop=num_input_features, num=2)]
        # Maximum number of levels in tree
        max_depth = [int(x) for x in np.linspace(5, 60, num=3)]
        max_depth.append(None)
        # Method of selecting samples for training each tree
        bootstrap = [True]
        # Create the random grid
        random_grid = {
            "n_estimators": n_estimators,
            "max_features": max_features,
            "max_depth": max_depth,
            "bootstrap": bootstrap,
        }

        mdl = run_hyperparamter_optimising(random_grid, mdl, cv)

    mdl.fit(X=X_train, y=y_train)
    return mdl


# Function for running hyperparameter optimization
def run_hyperparameter_optimizing(grid, mdl, cv):
    """
    Run hyperparameter optimization for a given model and hyperparameter grid.

    Args:
    - grid (dict): A dictionary of hyperparameters and their possible values to search over.
    - mdl: A machine learning model object.
    - cv (int): The number of cross-validation folds to use.

    Returns:
    - A RandomizedSearchCV object with the best hyperparameters found.
    """
    mdl = RandomizedSearchCV(
        estimator=mdl,
        param_distributions=grid,
        n_iter=10,  # Number of iterations to run the search for.
        cv=cv,  # Number of cross-validation folds to use.
        verbose=0,
        random_state=42,
    )
    return mdl


def linear_regression(X_train, y_train, method=None, cv=None):
    """
    Fit a linear regression model with or without regularization.

    Args:
    - X_train (array-like): The training data for the independent variables.
    - y_train (array-like): The training data for the dependent variable.
    - method (str): The regularization method to use. Must be one of: 'LASSO', 'RIDGE', 'ELASTIC', or None.
    - cv (int): The number of cross-validation folds to use.

    Returns:
    - A trained linear regression model object.
    """
    if method is None:
        mdl = linear_model.LinearRegression()
    elif method == "LASSO":
        mdl = linear_model.LassoCV(max_iter=100, cv=cv)
    elif method == "RIDGE":
        mdl = linear_model.RidgeCV(max_iter=100, cv=cv)
    elif method == "ELASTIC":
        mdl = linear_model.ElasticNetCV(max_iter=100, cv=cv)
    else:
        raise ValueError("Unknown Linear method")  # Raise an error for unknown method.

    mdl.fit(X=X_train, y=y_train)

    return mdl
