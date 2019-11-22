import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Number of clusters used in K-means
k = 3
numIterations = 100




def countOccurrences(k, realLabels, predictedLabels):

    occurrences = {}
    for idx, label in enumerate(realLabels):
        try:
            occurrences[label+ " - cluster: " + str(predictedLabels[idx])] += 1
        except:
            occurrences[label+ " - cluster: " + str(predictedLabels[idx])] = 1

    return occurrences


def printOccurrences(occurrences):
    print("\n\nOccurrences found:\n\n")
    for x, y in occurrences.items():
        print(x, "- occurrences:", y) 
    
    return



def convertLabelsTo01(realLabels, ALL, AML):

    realLabels_01 = []  
    for label in list(realLabels):
        if label == 'ALL':
            realLabels_01.append(ALL)
        elif label == 'AML':
            realLabels_01.append(AML)
    return realLabels_01


def calculateAccuracy(realLabels, predictedLabels):

    realLabels_01 = convertLabelsTo01(realLabels, 0, 1)
    realLabels_10 = convertLabelsTo01(realLabels, 1, 0)

    rigthGuesses01 = (np.array(realLabels_01) == predictedLabels)
    rigthGuesses10 = (np.array(realLabels_10) == predictedLabels)

    # Because there's no way to know which cluster will be assigned to each class
    rigthGuesses = max(np.sum(rigthGuesses01), np.sum(rigthGuesses10))

    numSamples = len(realLabels) 
    numRigthGuesses = np.sum(rigthGuesses)

    accuracy = numRigthGuesses / numSamples * 100

    print("Number of samples:", numSamples)
    print("Number of rigth classifications:", numRigthGuesses)
    print("Accuracy:", accuracy)
    return




# Loads the dataframe
df = pd.read_csv("/home/colombelli/Documents/datasets/leukemia_big.csv", header=None)

# Gets the transposed dataframe without the column that labels the samples
x = df.T.iloc[0:, 1:]

# Standardizes the features 
#x = StandardScaler().fit_transform(x)

# Clusters the data
kmeans = KMeans(n_clusters=k, n_init=numIterations)
kmeans.fit(x)

# Get information about the clusterization
predictedLabels = kmeans.predict(x)
centroids = kmeans.cluster_centers_

realLabels = df.iloc[0:1].values[0]

occurrences = countOccurrences(k, realLabels, predictedLabels)
printOccurrences(occurrences)
print("\n\n\n")

if k == 2:
    calculateAccuracy(realLabels, predictedLabels)
    print("\n\n\n")



# Plotting the results with PCA

# PCAing
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
pcaDf = pd.DataFrame(data = principalComponents, 
                           columns = ['principal component 1', 'principal component 2'])
pcaDf['cluster'] = predictedLabels

# Plotting
fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)

if k == 2:
    targets = [0, 1]
    colors = ['r', 'b']
    imageName = 'k2.png'

elif k == 3:
    targets = [0, 1, 2]
    colors = ['r', 'g', 'b']
    imageName = 'k3.png'

else:
    raise("k must be 2 or 3")

for target, color in zip(targets,colors):
    indicesToKeep = pcaDf['cluster'] == target
    ax.scatter(pcaDf.loc[indicesToKeep, 'principal component 1']
               , pcaDf.loc[indicesToKeep, 'principal component 2']
               , c = color
               , s = 50)
ax.legend(targets)
ax.grid()
plt.savefig(imageName)




