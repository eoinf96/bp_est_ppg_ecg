import unittest
import pandas as pd
import numpy as np
from bp_est_ppg_ecg import BP_est_dataset
from sklearn.datasets import make_regression

class TestBPEstDataset(unittest.TestCase):

    def setUp(self):
        X, y = make_regression(n_samples=100, n_features=5, random_state=42)
        feature_names = ['feat_1', 'feat_2', 'feat_3', 'feat_4', 'feat_5']
        df = pd.DataFrame(X, columns=feature_names)
        df['SBP'] = y
        self.BP_est_dataset = BP_est_dataset(df=df, df_aug=df, study = None)

    def test_upload_study(self):
        self.assertIsNone(self.BP_est_dataset.study)
        self.assertFalse(self.BP_est_dataset.demo_added)

    def test_upload_df(self):
        self.assertEqual(self.BP_est_dataset.Original.df.shape, (100, 7))
        self.assertEqual(self.BP_est_dataset.Augmented.df.shape, (100, 7))

    def test_calibrate_dataset(self):
        self.BP_est_dataset.calibrate_dataset(std_flag=False, num_cali=5)

        # Check that the calibration column was added
        self.assertTrue('calibration' in self.BP_est_dataset.Original.df.columns)
        self.assertTrue('calibration' in self.BP_est_dataset.Augmented.df.columns)

    def test_drop_collinear_features(self):
        # Create a dataset with collinear features
        X = np.random.rand(100, 8)
        y = np.random.rand(100, 1)
        X[:, 1] = X[:, 0] * 2
        X[:, 2] = X[:, 0] * 3
        df = pd.DataFrame(X, columns=['feat1', 'feat2', 'feat3', 'feat4', 'feat5', 'feat6', 'feat7', 'feat8'])
        df['SBP'] = y
        dataset = BP_est_dataset(df=df, study=None)

        # Check that collinear features are dropped correctly
        dataset.drop_collinear_features(VIF_max=2)
        expected_features = ['feat3', 'feat4', 'feat5', 'feat6', 'feat7', 'feat8']

        assert dataset.feature_names == expected_features


if __name__ == '__main__':
    unittest.main()