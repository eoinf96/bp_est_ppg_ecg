import numpy as np


def evaluate_model(mdl, train_df_X, test_df_X, ):

    sigma = np.array(mdl._std['BP'][mdl._msk]) if mdl._z_norm_flag else 1
    y_est_train = mdl.model.predict(train_df_X) * sigma
    mdl.train_Y *= sigma
    mdl.train_RMSE = np.sqrt(np.mean((y_est_train - mdl.train_Y* sigma) ** 2))
    mdl.train_corr = np.corrcoef(y_est_train, mdl.train_Y)[0, 1]

    sigma = np.array(mdl._std['BP'][~mdl._msk]) if mdl._z_norm_flag else 1
    mdl.y_est = mdl.model.predict(X=test_df_X) * sigma
    mdl.test_Y *= sigma
    mdl.test_RMSE = np.sqrt(np.mean((mdl.y_est - mdl.test_Y) ** 2))

    mdl.correlation_coefficient = mdl.get_correlation_coefficient()
    mdl.R_squared = mdl.get_coefficient_of_determination()