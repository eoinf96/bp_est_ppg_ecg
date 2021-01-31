from sklearn.preprocessing import StandardScaler
from sklearn import linear_model
import numpy as np
from funcs.regression_funcs.evaluate_model import evaluate_model


def linear_regression(mdl, method = None):
    # mdl.model_name = 'Linear'
    if method is not None:
        mdl.model_name = mdl.model_name + '_' + method

    sc = StandardScaler()
    train_df_X = sc.fit_transform(mdl.train_X)
    test_df_X = sc.transform(mdl.test_X)

    if method is None:
        mdl.model = linear_model.LinearRegression()
    elif method is 'LASSO':
        mdl.model = linear_model.LassoCV(max_iter=100)
    elif method is 'RIDGE':
        mdl.model = linear_model.RidgeCV(max_iter=100)
    elif method is 'ELASTIC':
        mdl.model = linear_model.ElasticNetCV()
    else:
        ValueError('Unknown Linear method')

    mdl.model.fit(X=train_df_X, y=mdl.train_Y)


    evaluate_model(mdl, train_df_X, test_df_X)