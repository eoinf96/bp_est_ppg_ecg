from sklearn import linear_model
import numpy as np
from sklearn.preprocessing import StandardScaler
import copy

def initialise_feature_selection(mdl):
    # mdl.feat_in_lasso = {feat: 0 for feat in mdl.feature_names}
    # mdl.top_5 = {feat: 0 for feat in mdl.feature_names}
    # mdl.top_15 = {feat: 0 for feat in mdl.feature_names}
    # mdl.fold_info['feature_list'] = {}
    mdl.fold_info['feature_list_rank'] = {}

def LASSO_feature_selection(mdl, max_iter= 30, num_features = np.inf):
    if not hasattr(mdl, 'fold_info'):
        initialise_feature_selection(mdl)

    ### Note attempt to use Least Angle Regression and see the results
    # feat_model = linear_model.LassoLarsCV(cv = loso)
    mdl.generate_loso()
    feat_model = linear_model.LassoCV(max_iter=max_iter, cv=mdl._loso) # Cant use the original loso or else itl get exhausted (bless it!)

    feat_model.fit(X=StandardScaler().fit_transform(mdl.train_X), y=mdl.train_Y)

    df = copy.deepcopy(mdl.df)
    df = df.drop(np.setdiff1d(list(df.columns), ['BP']), axis = 1)

    best_cont_features = np.argsort(abs(feat_model.coef_))[::-1]
    idx = -1
    #Store the features that are in the top 25 and top 5 of features
    # for idx in range(num_features):
    while True:
        idx +=1
        #Exit conditions
        if idx == len(best_cont_features):
            break
        if feat_model.coef_[best_cont_features[idx]] == 0.0:
            break
        if idx == num_features:
            break
        feature_name = mdl.feature_names[best_cont_features[idx]]
        # mdl.feat_in_lasso[feature_name] += 1
        df[feature_name] = mdl.df[feature_name]
        # if idx < 5:
        #     mdl.top_5[feature_name] += 1
        # if idx < 15:
        #     mdl.top_15[feature_name] += 1

    return df

def convert_to_prop_in_lasso(mdl):
    # mdl.prop_in_lasso = {feat_name: mdl.feat_in_lasso[feat_name] / len(mdl.MOLLIE.volunteer_list) for feat_name in
    #                  mdl.feat_in_lasso}
    # mdl.prop_in_5 = {feat_name: mdl.top_5[feat_name] / len(mdl.MOLLIE.volunteer_list) for feat_name in
    #              mdl.top_5}
    # mdl.prop_in_15 = {feat_name: mdl.top_15[feat_name] / len(mdl.MOLLIE.volunteer_list) for feat_name in
    #               mdl.top_15}
    pass

