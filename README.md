# Features from the photoplethysmogram and the electrocardiogram for estimating changes in blood pressure

This repository contains the supplementary code for feature extraction and subsequent blood pressure (BP) estimation from the photoplethysmogram (PPG) and electrocardiogram (ECG) waveforms. The code can be used to replicate the experiments and further any analysis.

Note that due to data protection laws, data from the clinical trial used to develop our algorithms cannot be provided - any example data provided in this repository is either from independent signal recordings or an [open source ICU database.](https://mimic.mit.edu/)

## Structure

This repository is split into two parts:

1. [Feature extraction](/Feature_extraction) - MATLAB code to extract PPG features, ECG features, and the pulse arrival time (PAT).
2. [BP estimation](/BP_est) - Python code for BP estimation using the relevant features in a LOSOCV framework.


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


