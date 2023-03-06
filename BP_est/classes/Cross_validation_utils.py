





    ##############
    # Run model
    ##############
    def run_model(
        self,
        vars=None,
        run_hyper_param_tuning=True,
        do_feat_optimisation=True,
        num_features=np.inf,
        plot_flag=False,
        return_fig=False,
    ):
        ''' Run regression models for BP estimation using features

        Inputs:
        Vars -- hyperparameters for models if neccesary -> if None then default will be used
        run_hyper_param_tuning -- tune hyperparameters
        do_feature_optimisation -- run LASSO feature selection
        num_features -- Maximum number of features to accept by LASSO feature selection
        plot_flag -- Show model responses
        return_fig -- for saving

        '''
        # Run Leave one subject out cross validation
        fig = crossvalidation_funcs.LOSOCV(
            self,
            model_vars=vars,
            num_features=num_features,
            do_feat_optimisation=do_feat_optimisation,
            run_hyper_param_tuning=run_hyper_param_tuning,
            plot_flag=plot_flag,
            return_fig=False,
        )
        # if plot_flag:
        #     fig = self.plot_individual_response(return_fig=return_fig)


        return fig

    def save_model(self, file_name, folder_name="/"):
        ''' Save output of model '''
        from pickle_functions import pickle_save

        a = self.__dict__
        if '_loso' in a:
            a.pop("_loso")
        a.pop("model")  # As also generator
        pickle_save(
            folder_loc='/pickles/'+self.dataset.study_id+'/models/' + folder_name + '/',
            file_name=file_name,
            var=a,
        )

    ##############
    # Model statistics
    ##############
    def get_coefficient_of_determination(self):
        ''' Compute and return coefficient of determination '''
        if self.y_est.size == 0:
            ValueError('Run model first')
        y_dash = np.mean(self.test_y)
        sum_sqared_errors = np.sum((self.test_y - self.y_est) ** 2)
        total_sum_of_squares = np.sum((self.test_y - y_dash) ** 2)
        self.r_squared = 1 - sum_sqared_errors / total_sum_of_squares

        return self.r_squared

    def get_correlation_coefficient(self):
        ''' Compute and return correlation coefficient '''
        if self.y_est.size == 0:
            ValueError('Run model first')
        # Correlation coefficient between yest and testy
        self.correlation_coefficient = np.corrcoef(self.test_y, self.y_est)[0, 1]

        return self.correlation_coefficient


    ##############
    # Plots
    ##############
    def plot_correlation(self, return_fig = False):
        ''' Plot scatter of y test and y est'''
        fig = plt.figure()
        plt.scatter(self.test_y, self.y_est)
        plt.xlabel("BP cuff (mmHg)")
        plt.ylabel("BP est (mmHg)")
        plt.grid()

        return fig if return_fig else 0

    def plot_individual_response(self, return_fig=False):
        ''' Plot Individual response'''
        fig = plotting_functions.plot_individual_response(self, return_fig=return_fig)
        return fig

    def show_feature_contributions(self, return_fig=False):

        if self.eval_method == "CV":
            return 0
        else:
            if self.model_name.upper() == "LINEAR":
                feature_importance = abs(self.model.coef_)
            elif self.model_name.upper() == "RF":
                feature_importance = self.model.feature_importances_
            elif self.model_name.upper() == "SVM":
                feature_importance = abs(self.model.coef_)

        # feature_weights = (feature_weights - np.mean(feature_weights))/np.std(feature_weights)
        parameter_names = self.df.drop("BP", axis=1).columns

        plt.rc("axes", axisbelow=True)
        fig = plt.figure()
        plt.bar(parameter_names, feature_importance)
        plt.xticks(rotation="vertical")
        plt.grid()
        plt.ylabel("Feature importance")
        plt.xlabel("Feature")


        return fig if return_fig else 0