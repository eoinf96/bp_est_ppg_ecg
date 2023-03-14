# Features from the photoplethysmogram and the electrocardiogram for estimating changes in blood pressure

This Python codebase performed leave-one-subject-out cross validation (LOSOCV) for BP estimation using the features extracted in [feature extraction](/Feature_extraction) across a number of participants. 

## Structure
``` 
├── /bp_est_ppg_ecg/                     # BP estimation using PPG and ECG features
│   ├── /__init__.py      
│   ├── /Dataset_utils.py      		 # Classes to handle datasets for blood pressure (BP) estimation.
│   ├── /MOLLIE_session.py               # Class specific to the MOLLIE clinical study -- NO participant identifying information is present in this file
│   ├── /regression_model_funcs.py       # Regression models
│   ├── /run_LOSOCV.py      		 # Run LOSOCV
├── /example_data 	                 # example data  
├── /tests                               # Unit tests
├── requirements.txt                     # configs file 
└── README.md
```
## bp_est_ppg_ecg

The below figure (from the initial publication) details the process implemented for BP estimation.

<p float="center">
  <img src="../figs/LOSOCV_framework.png" width="99%" alt>
  <b>Figure -</b> <em> Schematic of the Delta BP estimation pipeline for each of the proposed models. We extracted features from the PPG and ECG and averaged their values within a window of size 40s centred on times of cuff inflations. We then implemented a hybrid calibration approach such that the proposed models estimate Delta BP from a baseline calibration value determined during the rest period. Data augmentation was implemented to increase the training and validation set size by interpolating between cuff inflations. Models were trained and evaluated in a nested leave-one-subject-out cross-validation (LOSOCV) framework shown here by the iterator j which indicates the test participant for that iteration. Participant j was then removed from the training/ validation set (XAug) for that iteration. </em>
</p>

### Running the code

See **`run_LOSOCV.py`** to run BP estimation using the example data, this script performs the following taks:

- Loads dataset and augmented dataset.
- Initialises custom **`MOLLIE_session()`** class detailing study attributes.
- Initialises custom **`BP_est_dataset(df=df, df_aug=df_aug, study=MOLLIE, BP_name=BP_name)`** class.
  -  **`BP_est_dataset`** handles the backend of managing the main dataset and the augmented dataset, as well as defining the LOSOCV framework, and storing model results.
- Calibrate dataset using the calibration period.
- Remove collinear features.
- Loop through each ID and train LASSO+OLS and RF models for BP estimation.
- Store the models and model performance metrics.
- Pickle the results for later analysis of performance metrics.

## Example data

Some example data of ECG and PPG features have have been supplied to assist in understanding this codebase. This data is in .csv format. Additionally, example participant demographics.csv file has been provided. Due to data privacy, this data is NOT from the clinical trial we ran as part of our work - the demographics file does not represent any individuals.

## License

This project is licensed under the [MIT License](LICENSE).

## <a name="cite"/> :clipboard: Citation

If you use this code in your research, please consider citing our paper:
```
@article{Finnegan2023,
               title = {{Features from the photoplethysmogram and the electrocardiogram for estimating changes in blood pressure}},
               author = {Finnegan, Eoin and Davidson, Shaun and Harford, Mirae and Watkinson, Peter and Tarassenko, Lionel and Villarroel, Mauricio},
               journal = {Scientific Reports},
               month = {jan},
               number = {1}, 
               pages = {1--20},
               volume = {13},
               year = {2023}}

```

## Contact

If you have any questions or comments, please contact the authors at eoin.finnegan@stx.ox.ac.uk


