import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer
import os
from sklearn.preprocessing import StandardScaler

data_save_path = "/Users/raghavtandon/Documents/PhD/multi-modal/data"
# Identifies which biomarkers correspond to zeros in the probmat
def zero_in_probmat(prob_mat):
	# peps which show 0 probability for the no AD prob_mat
	zero_peps_no = np.sum(prob_mat[:,:,0] == 0, axis=0)
	# peps which show 0 probability for the yes AD prob_mat
	zero_peps_yes = np.sum(prob_mat[:,:,1] == 0, axis=0)
	return zero_peps_no, zero_peps_yes
	# Index of the no peptides
	# idx_peps_no = np.nonzero(zero_peps_no)[0]
	# # Index of the yes peptides
	# idx_peps_yes = np.nonzero(zero_peps_yes)[0]
	# # Concatenated list of all peptides which ever show a 0 probability
	# return 
	# bmname_no = [bmname[i] for i in idx_peps_no]
	# bmname_yes = [bmname[i] for i in idx_peps_yes]
	# bmname_to_drop = bmname_no + bmname_yes
	# print(zero_peps_no, "\n", zero_peps_yes)
	# print(bmname_to_drop)

def zero_bm(prob_mat, bmname, cr):
	zero_peps_no, zero_peps_yes = zero_in_probmat(prob_mat)
	idx_no = np.where(zero_peps_no > cr)[0]
	idx_yes = np.where(zero_peps_yes > cr)[0]
	bmname_no = [bmname[i] for i in idx_no]
	bmname_yes = [bmname[i] for i in idx_yes]
	bmname_to_drop = bmname_no + bmname_yes
	bmname_remove_indices = [bmname.index(i) for i in bmname_to_drop]; 
	bmname_remove_indices.sort()
	bmname_keep_indices = [i for i in range(len(bmname)) if i not in bmname_remove_indices]; 
	bmname_keep_indices.sort()
	bmname_keep_names = [bmname[i] for i in bmname_keep_indices]
	return bmname_keep_indices, bmname_remove_indices

def imptKNN(x, n_neighbors, weights="uniform"):
	imputer = KNNImputer(n_neighbors=n_neighbors, weights=weights)
	x_imputed = imputer.fit_transform(x)
	return x_imputed

def imputeData(x):
	x_imp = x.T.values.copy()
	n_neighbors=10
	x_imp = imptKNN(x_imp, n_neighbors)
	Ximp = pd.DataFrame(x_imp, index=x.columns, columns=x.index)
	return Ximp

def readLabels(peptide_df):
	label_df = peptide_df[["Replicate", "Condition"]]
	label_df = label_df[label_df["Condition"].isin(["AD", "Control", "AsymAD"])]
	label_dict = dict(zip(label_df["Replicate"], label_df["Condition"]))
	label_df = pd.DataFrame.from_dict(label_dict, orient="index").reset_index()
	label_df.columns = ["sbj", "DX"]
	label_df.set_index("sbj", inplace=True)
	return label_df

def readData():
	peptide_fname = os.path.join(data_save_path,"EHBS-2-210611", "SRM data from Caroline", "BSR2020-102", "Peptide Area Report_BSR2020-102_80pep.csv")
	skyline_fname = os.path.join(data_save_path,"EHBS-2-210611", "SRM data from Caroline", "BSR2020-102", "SkylineRatios-FullPercesion_2021_0608.csv")
	skyline_df = pd.read_csv(skyline_fname, index_col=0)
	peptide_df = pd.read_csv(peptide_fname)
	# Get the labels for all subjects in the data
	label_df = readLabels(peptide_df)
	# Impute missing values
	Ximp = imputeData(skyline_df)
	ss = StandardScaler()
	scaledX = ss.fit_transform(Ximp)
	scaledX = pd.DataFrame(scaledX, index=Ximp.index, columns=Ximp.columns)
	scaledX = scaledX.reindex(label_df.index)
	Ximp = Ximp.reindex(label_df.index)
	return scaledX, Ximp, label_df

def extract_classes(X, label_df, dx):
	merged = pd.merge(X, label_df, left_index=True, right_index=True)
	merged_subset = merged[merged["DX"].isin(dx)]
	# y = merged_subset["DX"].map({"Control":0, "AD":1}).to_numpy()
	y = merged_subset["DX"]
	x = merged_subset[merged_subset.columns.difference(["DX"])]
	x = x.to_numpy().astype(float)
	return x, y


def stage_histogram_asym(stages, y, max_stage=None, class_names=None):
	colors = ['C{}'.format(x) for x in range(10)]
	fig, ax = plt.subplots(figsize=(8, 4))
	hist_dat = [stages[y == 0],
                stages[y == 1],
                stages[y == 2]]
	if class_names is None:
		class_names = ['CN', 'AD', "AsymAD"]
	if max_stage is None:
		max_stage = stages.max()
	hist_c = colors[:3]
	n, bins, patch = ax.hist(hist_dat,
                             label=class_names,
                             density=True,
                             color=hist_c,
                             stacked=False,
                             bins=max_stage)
	ax.legend(loc=0, fontsize=20)

	idxs = np.arange(max_stage+2)
	bin_w = bins[1] - bins[0]
	ax.set_xticks(bins+bin_w/2)
	ax.set_xticklabels([str(x) for x in idxs], rotation=90)

	ax.set_ylabel('Fraction', fontsize=20)
	ax.set_xlabel('EBM Stage', fontsize=20)
	ax.tick_params(axis='both', which='major', labelsize=13)
	fig.tight_layout()
	return fig, ax



def sequence_error(order):
	clusters = len(order)
	ebm_order = []
	bm_arrays = []
	for _ in range(clusters):
		v = len(order[_])
		ebm_order.append(v)
		bm_arrays.append(np.array(range(v)))

	ebm_cumsum = np.cumsum(ebm_order)
	for i, j in enumerate(ebm_cumsum[:-1]):
		bm_arrays[i+1] += j

	overlap = []
	for _ in range(clusters):
		overlap.append(len(set(order[_]).intersection(set(bm_arrays[_]))))

	return overlap/np.array(ebm_order)





