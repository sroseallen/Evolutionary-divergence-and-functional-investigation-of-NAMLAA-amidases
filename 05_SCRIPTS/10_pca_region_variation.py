from Bio import AlignIO
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from tqdm import tqdm


def pca_prep(
    coordinates: list,
    alignment_path: str = "./01_DATA/Amidase_3/04_Multiple_Alignment/alignment_3_thresh1.0.fa",
) -> pd.DataFrame:
    alignment = AlignIO.read(alignment_path, "fasta")

    # initialise dataframe with cols for each insertion region
    df = pd.DataFrame(
        {
            "name": [],
            "I-1": [],
            "I-2": [],
            "I-3": [],
            "I-4": [],
            "I-5": [],
            "I-6": [],
            "I-7": [],
            "I-8": [],
        }
    )

    for seq in tqdm(alignment):
        seq_occupancy = 0
        # pulls out phylum only (labelling for the PCA later on)
        seq_add = [seq.id.split(":")[1]]
        # quantify the number of occuped (non-gap) columns in the region
        for coords in coordinates:
            insertion_region = list(range(coords[0], coords[1]))
            for pos in insertion_region:
                if seq[pos - 1] != "-":
                    seq_occupancy += 1
            seq_add.append(seq_occupancy)
        # reads species identifier and the occupancy into the dataframe
        to_add = pd.DataFrame([seq_add], columns=df.columns)
        df = pd.concat([df, to_add])

    # make the dataframe coercible to a matrix (ie column 1 becomes row names)
    df = df.set_index("name")
    df.index.names = [None]
    return df


coordinates = [
    [1, 6],
    [15, 57],
    [81, 85],
    [96, 128],
    [154, 175],
    [184, 241],
    [252, 303],
    [355, 399],
]
df = pca_prep(coordinates)
print(df.head(10))

# Scale data within the dataframe
scaler = StandardScaler()
scaler.fit(df)
array_scaled = scaler.transform(df)

df_scaled = pd.DataFrame(data=array_scaled, columns=df.columns)

# Run initial PCA (to identify PCs that explain the majority of the variance)
pca = PCA(n_components=8)  # one component for each insertion region
pca.fit_transform(df_scaled)
prop_var = pca.explained_variance_ratio_
# eigenvalues = pca.explained_variance_

# scree plot
PC_numbers = np.arange(pca.n_components_) + 1
plt.bar(PC_numbers, prop_var, color="navy")
for i in range(len(PC_numbers)):
    plt.text(
        i + 1, prop_var[i] + 0.01, f"{prop_var[i]*100:.1f}%", ha="center", va="center"
    )
plt.title("Scree Plot", fontsize=8)
plt.ylabel("Proportion of Variance", fontsize=8)
plt.xlabel("Principal Component", fontsize=8)
plt.show()

# Run secondary PCA on PCs explaining most of the variance (PCs1-2 explan >90% of variance per Scree plot)
pca_2 = PCA(n_components=3)
pca_2_nums = pca_2.fit_transform(df_scaled)

pca_output = pd.DataFrame(data=pca_2_nums, columns=["PC1", "PC2", "PC3"])
loadings = pd.DataFrame(
    pca_2.components_.T, columns=["PC1", "PC2", "PC3"], index=df.columns
)
print(loadings)

# add phylum as labels using the names assigned as rownames in the original df
df.index.name = "label"
df.reset_index(inplace=True)
df = df[["label"]]
pca_output["label"] = df

# Plot the output, colour by phylum
sns_ax = sns.scatterplot(
    x="PC1", y="PC2", hue="label", data=pca_output, palette="Accent"
)
sns.move_legend(sns_ax, "upper left", bbox_to_anchor=(1, 1))
plt.rcParams["axes.labelsize"] = 15
plt.title("PCA plot: PC1, PC2")
plt.xlabel("PC1: 51.2% Variance")
plt.ylabel("PC2: 32.7% Variance")
plt.show()

sns_ax = sns.scatterplot(
    x="PC1", y="PC3", hue="label", data=pca_output, palette="Accent"
)
sns.move_legend(sns_ax, "upper left", bbox_to_anchor=(1, 1))
plt.rcParams["axes.labelsize"] = 15
plt.title("PCA plot: PC1, PC3")
plt.xlabel("PC1: 51.2% Variance")
plt.ylabel("PC3: 12.3% Variance")
plt.show()

sns_ax = sns.scatterplot(
    x="PC2", y="PC3", hue="label", data=pca_output, palette="Accent"
)
sns.move_legend(sns_ax, "upper left", bbox_to_anchor=(1, 1))
plt.rcParams["axes.labelsize"] = 15
plt.title("PCA plot: PC2, PC3")
plt.xlabel("PC2: 32.7% Variance")
plt.ylabel("PC3: 12.3% Variance")
plt.show()


# create biplot - overlay eigenvectors onto the plot
def biplot(score, components, labels=None):
    x_score = score[:, 0]
    y_score = score[:, 1]
    features = components.shape[0]  # number of features in the model (I-1 to I-8)
    scalex = 1.0 / (x_score.max() - x_score.min())  # range of values on the x axis
    scaley = 1.0 / (y_score.max() - y_score.min())  # range of values on the y axis
    plt.scatter(x_score * scalex, y_score * scaley, s=5, color="gray")
    for i in range(features):
        plt.arrow(
            0,  # start/end positions for the arrow
            0,
            components[i, 0],
            components[i, 1],
            head_width=0.05,
            head_length=0.05,
            overhang=1,
            color="red",
        )
        plt.text(
            components[i, 0]
            * 0.7,  # place label text slightly offset from the end of the arrow
            components[i, 1] * 0.7,
            labels[i],
            color="maroon",
        )

    plt.title("PCA Biplot: PC1, PC2")
    plt.xlabel("PC1: 51.2% Variance")
    plt.ylabel("PC2: 32.7% Variance")
    plt.show()


biplot(
    pca_2_nums,
    np.transpose(pca_2.components_),
    ["I-1", "I-2", "I-3", "I-4", "I-5", "I-6", "I-7", "I-8"],
)
