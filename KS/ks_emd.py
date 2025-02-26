print("imports ...")

import argparse
import numpy as np
from sklearn.preprocessing import RobustScaler
import pandas as pd
import matplotlib.pyplot as plt
import sys





print("arg parser ...")


parser = argparse.ArgumentParser(
        description="Parse four file paths and store them."
    )

# Adding arguments for 4 file paths
parser.add_argument(
    "--patientid",
    type=str,
    required=True,
    help="patient ID file"
)
parser.add_argument(
    "--data",
    type=str,
    required=True,
    help="data file of shape marker x cells"
)
parser.add_argument(
    "--metaclustering",
    type=str,
    required=True,
    help="file storing clustering and manual gating for each cell"
)

parser.add_argument(
    "--output",
    type=str,
    required=True,
    help="output vector file"
)

# Parse the arguments
args = parser.parse_args()


print(f"   --patientid    {args.patientid}")
print(f"   --data    {args.data}")
print(f"   --metaclustering    {args.metaclustering}")
print(f"   --output    {args.output}")
print(" ")



print("load data ...")

patient_id = pd.read_csv(args.patientid , index_col = 0)
print("   patient_id:")
print(patient_id)


markerXcells = pd.read_csv(args.data , index_col= 0)
print(f"\n   markerXcells ({sys.getsizeof(markerXcells)} byte):")
print(markerXcells)

marker = list(markerXcells.index)
print(f"\n   marker ({sys.getsizeof(marker)} byte):")
print(marker)

markerXcells = markerXcells.transpose().reset_index()[marker]
print(f"\n   markerXcells ({sys.getsizeof(markerXcells)} byte):")
print(markerXcells)


meta8_merging1 = pd.read_csv(args.metaclustering , index_col = 0).reset_index()
meta8_merging1 = meta8_merging1[["meta50"]] 
#meta8_merging1["meta50"] = ["cluster_" + str(x) for x in list(meta8_merging1.meta50)]
print(f"\n   meta8_merging1 ({sys.getsizeof(meta8_merging1)} byte):")
print(meta8_merging1)

anno_matrix = pd.read_csv("data/with_negative_markers.csv")










### Robust scaler  -->  skiped



print("Robust Scaler ... skipped.")


data = markerXcells.copy()
data["patient_id"] = list(patient_id.patient_id)
data["cluster"] = list(meta8_merging1.meta50)
data
print(f"\n   data ({sys.getsizeof(data)} byte):")


patients = sorted(list(set(list(patient_id.patient_id))))

norm = data.copy()
norm.to_csv(f"data/{args.output}/norm.csv")
#"string".copy()

### Annotation

q = ["Cell type"] + marker
q
anno_matrix_small = anno_matrix[q]
anno_matrix_small



### EMD & KS

print("compute EMD / KS ...")

from scipy.stats import wasserstein_distance
from scipy.stats import ks_2samp , kstest
bins = 100


celltypes = sorted(list(anno_matrix_small["Cell type"]))
clusters = sorted(list(set(list(norm["cluster"]))))
#'T-Helferzellen' in celltypes





class distribution:

    def __init__(self , cluster_data , background_data , method = "ks" , celltype = None , cluster = None , marker = None , bins = 100):
        self.celltype = celltype
        self.cluster = cluster
        self.marker = marker
        self.bins = bins
        self.cluster_data = cluster_data
        self.background_data = background_data

        #self.cluster_hist = self.hist_data(self.cluster_data) #/ np.sum(self.hist_data(self.cluster_data))
        #self.background_hist = self.hist_data(self.background_data) #/ np.sum(self.hist_data(self.background_data))
        self.cluster_hist = self.cluster_data
        self.background_hist = self.background_data

        if method != "ks":
            self.emd = wasserstein_distance(self.cluster_hist, self.background_hist)
            self.emd = self.emd * self.mean_shift()
        else:
            self.emd = None

        if method != "emd":
            self.ks_statistic = self.ks()
        else:
            self.ks_statistic = None
    
    def clear(self):
        self.cluster_data = None
        self.background_data = None
        self.cluster_hist = None
        self.background_hist = None

    def ks(self):
        obj = kstest(np.array(self.cluster_hist) , np.array(self.background_hist))
        stat = obj.statistic
        sign = obj.statistic_sign
        #print(stat , sign)
        return stat * sign

    
    def mean_shift(self):
        if np.mean(self.cluster_data) < np.mean(self.background_data):
            return -1
        else:
            return 1

        


    def hist_data(self , data):
        temp , _ , _ = plt.hist(data, bins = self.bins)
        plt.close()
        return temp

    def plot_hist(self):
        plt.figure(figsize=(10, 6))

        plt.hist(self.cluster_data, bins=self.bins, alpha=0.6, label='cluster', edgecolor='black')
        plt.hist(self.background_data, bins=self.bins, alpha=0.6, label='background', edgecolor='black')

        plt.xlabel('Value')
        plt.ylabel('Frequency')
        plt.title(f'{self.cluster} - {self.marker} - {self.celltype}\nEMD = {round(self.emd , 2)}   KS stat. = {round(self.ks_statistic , 2)}')
        plt.legend(loc='upper right')
        plt.grid(axis='y', linestyle='--', alpha=0.7)

        plt.show()

    def query(self , cluster , marker , celltype):
        return (self.celltype == celltype) and (self.cluster == cluster) and (self.marker == marker)
    

    






import gc
del data
del markerXcells
gc.collect()





from tqdm import tqdm
#results = []

cluster_log = []
marker_log = []
score_log = []
score_emd_log = []

matrix2 = np.zeros((len(marker) , len(clusters)))
matrix_emd = np.zeros((len(marker) , len(clusters)))
#print(matrix2.shape)

for j , cluster in tqdm(enumerate(clusters)):
    cell_type_score = 0
    for m , mark in enumerate(marker):

        cluster_data = list(norm[ norm.cluster == cluster ][mark])
        background_data = list(norm[ norm.cluster != cluster ][mark])

        d = distribution(cluster_data=cluster_data , background_data=background_data , marker=mark , cluster=cluster , celltype= None , method="both")
        matrix2[m,j] = d.ks_statistic
        matrix_emd[m,j] = d.emd


        cluster_log.append(cluster)
        marker_log.append(mark)
        score_log.append(d.ks_statistic)
        score_emd_log.append( d.emd )

        d.clear()
        #results.append(d)



print("export matrices ...")

print("\n   matrix_ks:")
print(matrix2)
print(matrix2.shape)
df = pd.DataFrame(matrix2)
df.columns = clusters
df.index = marker
df.to_csv(f"data/{args.output}/matrix_ks.csv")


print("\n   matrix_emd:")
print(matrix_emd)
print(matrix_emd.shape)
df = pd.DataFrame(matrix_emd)
df.columns = clusters
df.index = marker
df.to_csv(f"data/{args.output}/matrix_emd.csv")



pd.DataFrame({"cluster" : cluster_log , "marker" : marker_log , "ks" : score_log , "emd" : score_emd_log}).to_csv(f"data/{args.output}/computation_log.csv")

