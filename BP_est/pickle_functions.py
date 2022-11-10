import pickle
import _pickle as cPickle
import os
# from set_root import set_root


def pickle_save(folder_loc, file_name, var, overwrite_results_root = False):
    #Saving things in the results section unless overwrite_results_root = True
    if not overwrite_results_root:
        root = set_root()
        folder_loc = root + '/Studies/MOLLIE/results/' + folder_loc


    #Check if file_loc exists and if not make it
    if not os.path.isdir(folder_loc):
        print('Folder'+ folder_loc + ' does not exist. Making ....')
        os.mkdir(folder_loc)

    file_loc = folder_loc + '/' + file_name
    with open(file_loc, 'wb') as outf:
        cPickle.dump(var, outf, protocol=pickle.HIGHEST_PROTOCOL)
        # cPickle.dump(var, outf, protocol=pickle.HIGHEST_PROTOCOL)


def pickle_load(file_loc, overwrite_results_root = False):
    #loading things from results
    if not overwrite_results_root:
        root = set_root()
        file_loc = root + '/Studies/MOLLIE/results/' + file_loc

    if os.path.getsize(file_loc) > 0:
        with open(file_loc, 'rb') as inpf:
            var_out = cPickle.load(inpf)
    else:
        raise ValueError('No variable found: ' + file_loc)

    return var_out
