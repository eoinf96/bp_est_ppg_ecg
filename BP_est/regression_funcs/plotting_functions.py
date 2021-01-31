import numpy as np
import matplotlib.pyplot as plt
from funcs import remove_empty_subplots

def plot_individual_response(self, return_fig=False):
    n_rows = int(np.ceil(np.sqrt(len(self.MOLLIE.volunteer_list))))
    n_cols = int(np.floor(np.sqrt(len(self.MOLLIE.volunteer_list))))
    fig, ax = plt.subplots(nrows=int(n_rows), ncols=int(n_cols), figsize=(15, 8), sharex=True)
    gap = 0.071
    plt.subplots_adjust(left=gap, bottom=gap, right=1 - gap, top=1 - gap, wspace=0.345, hspace=0.405)
    plot_row_idx = 0
    plot_col_idx = 0

    ids = self.df.index
    y_est_marker = 0

    BP_list = self.df['BP'] * self._std['BP']

    for volunteer_num in self.MOLLIE.volunteer_list:
        axes = ax[plot_row_idx, plot_col_idx]
        curr_bps = BP_list.values[ids == self.MOLLIE.get_volunteer_id_name(volunteer_num)]
        curr_mask = self._msk[ids == self.MOLLIE.get_volunteer_id_name(volunteer_num)]
        curr_nans = np.any(np.isnan(self.df[ids == self.MOLLIE.get_volunteer_id_name(volunteer_num)].values), axis=1)
        index_test_and_not_nan = np.logical_and(~curr_nans, ~curr_mask)

        curr_yest = self.y_est[y_est_marker + np.arange(0, sum(index_test_and_not_nan))]

        y_est_marker += sum(index_test_and_not_nan)

        curr_idx = np.arange(0, len(curr_bps))
        axes.plot(curr_idx, curr_bps, '-o', label='Ground Truth', markersize=2)
        axes.plot(curr_idx[curr_mask], curr_bps[curr_mask], '-o', label='Train', markersize=2)
        axes.plot(curr_idx[index_test_and_not_nan], curr_yest, '-o', label='Test', markersize=2)
        axes.set_xticklabels([])
        axes.set_title(self.MOLLIE.get_volunteer_id_name(volunteer_num))
        if self._z_norm_flag:
            axes.set_ylabel('Norm BP')
        elif self._mean_norm_flag:
            axes.set_ylabel(r'$\Delta$ BP (mmHg)')
        else:
            axes.set_ylabel('BP (mmHg)')

        axes.grid()
        plot_col_idx = plot_col_idx + 1

        if plot_col_idx == n_cols:
            plot_col_idx = 0
            plot_row_idx = plot_row_idx + 1

        # if plot_row_idx == 0 and plot_col_idx == 1:
        if volunteer_num == 31:
            axes.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=3)
            # axes.legend(bbox_to_anchor=(1,0), loc="best",
            # bbox_transform=fig.transFigure, ncol=3)

    # Assuming we are on the last row here
    remove_empty_subplots.remove_empty_subplots(n_cols, plot_col_idx, ax)

    if not return_fig:
        fig = 0
    return fig

def plot_test_response(self, return_fig=False):
    n_rows = int(np.ceil(np.sqrt(len(self.MOLLIE.TestID))))
    n_cols = int(np.ceil(np.sqrt(len(self.MOLLIE.TestID))))
    fig, ax = plt.subplots(nrows=int(n_rows), ncols=int(n_cols))
    plot_row_idx = 0
    plot_col_idx = 0

    test_ids = list(self.df_test_Y.index.unique())

    for volunteer_id in test_ids:
        axes = ax[plot_row_idx, plot_col_idx]
        curr_bps = self.df_test_Y.values[self.df_test_Y.index == volunteer_id]
        curr_yest = self.y_est[self.df_test_Y.index == volunteer_id]

        curr_idx = np.arange(0, len(curr_bps))
        axes.plot(curr_idx, curr_bps, '-o', label='Ground Truth', markersize=2)
        axes.plot(curr_idx, curr_yest, '-o', label='Test', markersize=2)
        axes.set_xticklabels([])
        axes.set_title(volunteer_id)
        if self._z_norm_flag:
            axes.set_ylabel('Norm BP')
        elif self._mean_norm_flag:
            axes.set_ylabel(r'$\Delta$ BP (mmHg)')
        else:
            axes.set_ylabel('BP (mmHg)')
        axes.grid()
        plot_col_idx = plot_col_idx + 1

        if plot_col_idx == n_cols:
            plot_col_idx = 0
            plot_row_idx = plot_row_idx + 1

        if plot_row_idx == 0 and plot_col_idx == 1:
            axes.legend(bbox_to_anchor=(1, 0), loc="lower right",
                        bbox_transform=fig.transFigure, ncol=3)

    # Assuming we are on the last row here
    subplots_left = n_cols - plot_col_idx
    for s_idx in range(subplots_left):
        ax[-1, -(s_idx + 1)].axis('off')

    if not return_fig:
        fig = 0
    return fig
