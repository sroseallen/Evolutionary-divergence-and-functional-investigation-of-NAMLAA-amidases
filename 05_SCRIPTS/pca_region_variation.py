from Bio import AlignIO
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from tqdm import tqdm

def pca_prep(coordinates:list, alignment_path:str = "./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa") -> pd.DataFrame:
    alignment = AlignIO.read(alignment_path, "fasta")

    # initialise dataframe with cols for each insertion region
    df = pd.DataFrame({
        "name": [],
        "I-1": [],
        "I-2": [],
        "I-3": [],
        "I-4": [],
        "I-5": [],
        "I-6": [],
        "I-7": [],
        "I-8": []
    })

    for seq in tqdm(alignment):
        seq_occupancy = 0
        # pulls out phylum only (labelling for the PCA later on)
        seq_add = [seq.id.split(":")[1]]
        # quantify the number of occuped (non-gap) columns in the region
        for coords in coordinates:
            insertion_region = list(range(coords[0],coords[1]))
            for pos in insertion_region:
                if seq[pos-1] != '-':
                    seq_occupancy += 1
            seq_add.append(seq_occupancy)
        # reads species identifier and the occupancy into the dataframe
        to_add = pd.DataFrame([seq_add], columns=df.columns)
        df = pd.concat([df, to_add], ignore_index=True)
    
    # make the dataframe coercible to a matrix (ie column 1 becomes row names)
    df = df.set_index("name")
    df.index.names = [None]
    return df

coordinates=[[1,6],[15,57],[81,85],[96,128],[154,175],[184,241],[252,303],[355,399]]
df = pca_prep(coordinates)

# Scale data within the dataframe
scaler = StandardScaler()
scaler.fit(df)
array_scaled = scaler.transform(df)

df_scaled = pd.DataFrame(data=array_scaled, 
                         columns=df.columns)

# Run initial PCA (to identify PCs that explain the majority of the variance)
pca = PCA(n_components=8) #one component for each insertion region
pca.fit_transform(df_scaled)
prop_var = pca.explained_variance_ratio_
#eigenvalues = pca.explained_variance_

# Scree plot
PC_numbers = np.arange(pca.n_components_) + 1
 
plt.plot(PC_numbers, 
         prop_var, 
         'ro-')
plt.title("Scree Plot", fontsize=8)
plt.ylabel("Proportion of Variance", fontsize=8)
plt.show()
plt.close()

# Run secondary PCA on PCs explaining most of the variance (PCs1-2 explan >90% of variance per Scree plot)
pca_2 = PCA(n_components=3)
pca_2_nums = pca_2.fit_transform(df_scaled)

pca_output = pd.DataFrame(data = pca_2_nums, 
                          columns = ["PC1","PC2","PC3"])

# add phylum as labels using the names assigned as rownames in the original df
df.index.name = "label"
df.reset_index(inplace=True)
df = df[["label"]]
pca_output["label"] = df

# Plot the output, colour by phylum
sns_ax = sns.scatterplot(x="PC1",y="PC2",hue="label",data=pca_output,palette="Accent")
sns.move_legend(sns_ax, "upper left", bbox_to_anchor=(1, 1))
plt.show()

sns_ax2 = sns.scatterplot(x="PC1",y="PC3",hue="label",data=pca_output,palette="Accent")
sns.move_legend(sns_ax2, "upper left", bbox_to_anchor=(1, 1))
plt.show()

sns_ax3 = sns.scatterplot(x="PC2",y="PC3",hue="label",data=pca_output,palette="Accent")
sns.move_legend(sns_ax3, "upper left", bbox_to_anchor=(1, 1))
plt.show()

# Plot the output, colour by cluster predicted with k-means
kmeans =KMeans(n_clusters=10).fit(df_scaled)
pca_output["cluster"] = pd.Categorical(kmeans.labels_)
sns.scatterplot(x="PC1",y="PC2",hue="cluster",data=pca_output)
plt.show()

# def biplot(score,coef,labels=None):
#     xs = score[:,0]
#     ys = score[:,1]
#     n = coef.shape[0]
#     scalex = 1.0/(xs.max() - xs.min())
#     scaley = 1.0/(ys.max() - ys.min())
#     plt.scatter(xs * scalex,ys * scaley,
#                 s=5, 
#                 color='orange')
#     for i in range(n):
#         plt.arrow(0, 0, coef[i,0], 
#                   coef[i,1],color = 'purple',
#                   alpha = 0.5)
#         plt.text(coef[i,0]* 1.15, 
#                  coef[i,1] * 1.15, 
#                  labels[i], 
#                  color = 'darkblue', 
#                  ha = 'center', 
#                  va = 'center')
 
#     plt.xlabel("PC{}".format(1))
#     plt.ylabel("PC{}".format(2))    

#     plt.figure()
#     plt.show()

# plt.title('Biplot of PCA')
# biplot(pca_2_nums, 
#        np.transpose(pca_2.components_), 
#        list(df.columns))