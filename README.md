# Features from the photoplethysmogram and the electrocardiogram for estimating changes in blood pressure

This repository contains the supplementary code for feature extraction and subsequent blood pressure (BP) estimation from the photoplethysmogram (PPG) and electrocardiogram (ECG) waveforms. The code can be used to replicate the experiments and suggest further analysis.

Note that due to data protection laws, data from the clinical trial used to develop our algorithms cannot be provided - all example data provided in this repository is from an [open source ICU database.](https://mimic.mit.edu/)

This repository is split into two parts:

1. [Feature extraction](/Feature_extraction) - MATLAB code to extract PPG features, ECG features, and the pulse arrival time (PAT).
2. [BP estimation](/BP_est) - Python code for BP estimation using the relevant features in a LOSOCV framework.

## Contributing

If you want to contribute to this project, please create a pull request with your changes. Before submitting a pull request, ensure that all tests pass and the code follows the style guidelines.

## License

This project is licensed under the [MIT License](LICENSE).

## Contact

If you have any questions or comments, please contact the authors at eoin.finnegan@stx.ox.ac.uk
