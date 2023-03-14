import pandas as pd
import numpy as np
import unittest
from bp_est_ppg_ecg import DatasetHandler
from sklearn.datasets import make_regression

class TestDatasetHandler(unittest.TestCase):

    def setUp(self):
        # Set up a test dataset
        X, y = make_regression(n_samples=100, n_features=5, random_state=42)
        feature_names = ['feat_1', 'feat_2', 'feat_3', 'feat_4', 'feat_5']
        df = pd.DataFrame(X, columns=feature_names)
        df['target'] = y
        self.df = df
        self.feature_names = feature_names
        self.target_name = 'target'
        self.dh = DatasetHandler(self.df, self.feature_names, self.target_name)

    def test_return_X_Y(self):
        # Test the return_X_Y method
        X, y = self.dh.return_X_Y(to_numpy=True, relevant_IDs=[0, 2])
        self.assertTrue(np.testing.assert_array_almost_equal(X, np.array([[-0.93782504,  0.51504769,  0.51503527,  3.85273149,  0.51378595],
                                                    [-0.60170661, -1.05771093,  1.85227818,  0.82254491, -0.01349722]]), decimal=5) is None)
        self.assertTrue(np.testing.assert_array_almost_equal(y, np.array([271.31612081,  11.86102446]), decimal=5) is None)


if __name__ == '__main__':
    unittest.main()