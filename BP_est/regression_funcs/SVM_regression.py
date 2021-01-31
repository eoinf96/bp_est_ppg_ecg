import numpy as np
from sklearn.svm import SVR
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from funcs.regression_funcs.evaluate_model import evaluate_model
from sklearn.model_selection import RandomizedSearchCV

CONST_VARS ={
    'kernel': 'linear',
    'C':10,
    'epsilon':0.1,
}


def SVM_regression(mdl,vars=CONST_VARS, run_hyper_param_tuning=True):
    # mdl.model_name = 'SVM'
    if vars is None:
        vars = CONST_VARS
    sc = StandardScaler()
    train_df_X = sc.fit_transform(mdl.train_X)
    test_df_X = sc.transform(mdl.test_X)

    mdl.model = SVR(kernel=vars['kernel'], C=vars['C'],  epsilon=vars['epsilon'])

    if run_hyper_param_tuning: run_hyperparamter_optimising(mdl)

    mdl.model.fit(X=train_df_X, y=mdl.train_Y)

    if run_hyper_param_tuning: print(mdl.model.best_params_)

    evaluate_model(mdl, train_df_X, test_df_X)


# Run hyperparameter tuning here
def run_hyperparamter_optimising(mdl):
    print('Running Hyperparameter tuning ... ')
    # kernels = ['linear', 'RBF', 'poly']
    kernels = ['linear']
    # degree = [2, 3]
    # gammas =  [1e-4, 1e-3, 0.01, 0.1, 0.2, 0.5, 0.6, 0.9]
    C = [0.1, 1, 10]
    epsilon = np.logspace(-3,3, num=7)

    random_grid = {'kernel': kernels,
                      'C': C,
                      'epsilon': epsilon}

    # random_grid = {'kernel': kernels,
    #                'C': C,
    #                # 'gamma': gammas,
    #                # 'degree': degree,
    #                'epsilon': epsilon}

    mdl.generate_loso()
    mdl.model = RandomizedSearchCV(estimator=mdl.model, param_distributions=random_grid,
                                    n_iter =10, cv=mdl._loso, verbose=1, random_state=42)