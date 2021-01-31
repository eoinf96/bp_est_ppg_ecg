from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
import numpy as np
from sklearn.model_selection import RandomizedSearchCV
from funcs.regression_funcs.evaluate_model import evaluate_model

CONST_VARS= {
'n_estimators': 128,
'max_features': 'sqrt',
'max_depth': 10,
'bootstrap': True
}

def RF_regression(mdl, vars = CONST_VARS, run_hyper_param_tuning=True):
    # mdl.model_name = 'RF'
    if vars is None:
        vars = CONST_VARS

    sc = StandardScaler()
    train_df_X = sc.fit_transform(mdl.train_X)
    test_df_X = sc.transform(mdl.test_X)

    mdl.model = RandomForestRegressor(n_estimators=vars['n_estimators'], max_features = vars['max_features'], max_depth=vars['max_depth'], bootstrap = vars['bootstrap'])

    if run_hyper_param_tuning: run_hyperparamter_optimising(mdl)

    mdl.model.fit(X=train_df_X, y=mdl.train_Y)

    if run_hyper_param_tuning: print(mdl.model.best_params_)

    evaluate_model(mdl, train_df_X, test_df_X)


# Run hyperparameter tuning here
def run_hyperparamter_optimising(mdl):
    print('Running Hyperparameter tuning ... ')
    # Number of trees in random forest -- Dont need to tune-- just set to a large enough number. There is no chance of overfitting based on the number of trees
    # n_estimators = [int(x) for x in np.linspace(start=500, stop=1000, num=2)]
    n_estimators = [128]
    # n_estimators = 128
    # Number of features to consider at every split
    max_features = ['auto']
    # Maximum number of levels in tree -- This one is more important
    max_depth = [int(x) for x in np.linspace(5, 60, num=3)]
    max_depth.append(None)
    # Method of selecting samples for training each tree
    # bootstrap = [True, False]
    bootstrap = [True]
    # Create the random grid
    random_grid = {'n_estimators': n_estimators,
                   'max_features': max_features,
                   'max_depth': max_depth,
                   'bootstrap': bootstrap}

    mdl.generate_loso()
    mdl.model = RandomizedSearchCV(estimator=mdl.model, param_distributions=random_grid,
                                    n_iter=10, cv=mdl._loso, verbose=0, random_state=42)