import sys
import pandas as pd
import numpy as np
import mykmeanssp as kmc

np.random.seed(1234)


def min_distance_from_centroids(data_point, centroids):
    min_distance = float("inf")
    for centroid in centroids:
        distance = np.linalg.norm(data_point - centroid)
        if distance < min_distance:
            min_distance = distance
    return min_distance


def kmeans_pp(K, file_name_1, file_name_2, iter, eps):
    # Choose one center uniformly at random from among the data points
    file_1 = pd.read_csv(file_name_1, header=None)
    file_2 = pd.read_csv(file_name_2, header=None)

    file_1.columns = ["key"] + [f"col_{i}" for i in range(1, len(file_1.columns))]
    file_2.columns = ["key"] + [f"col_{i}" for i in range(1, len(file_2.columns))]

    # Inner join the two files, based on key column
    merged_file = pd.merge(file_1, file_2, on="key", how="inner")

    # Sort the merged file by the key column in ascending order
    merged_file = merged_file.sort_values(by="key")

    # Convert the merged file to a numpy array
    numpy_data_points = merged_file.to_numpy()

    # remove first element from each data points, which is the key, as it is not needed for clustering
    numpy_data_points = numpy_data_points[:, 1:]

    n_samples = len(numpy_data_points)
    n_features = len(numpy_data_points[0])

    # choose the first center uniformly at random, and add it to the list of centers
    centroids = [numpy_data_points[np.random.choice(n_samples)]]

    for i in range(1, K):
        distances = []
        for data_point in numpy_data_points:
            if not np.any(np.all(data_point == centroids, axis=1)):
                min_distance = min_distance_from_centroids(data_point, centroids)
                distances.append(min_distance)
            else:
                distances.append(0)
        sum_distances = np.sum(distances)
        probabilities = [distance / sum_distances for distance in distances]
        new_center = numpy_data_points[np.random.choice(n_samples, p=probabilities)]
        centroids.append(new_center)

    centroids = np.array(centroids)
    centroid_indexes = []

    # Find the index of the centroids in the numpy_data_points array
    for centroid in centroids:
        index = np.where(numpy_data_points == centroid)
        centroid_indexes.append(index[0][0])

    # Run k-means clustering with the chosen initial centers
    centroids = kmc.fit(
        centroids.tolist(),
        numpy_data_points.tolist(),
        eps,
        K,
        iter,
        n_samples,
        n_features,
    )

    # Print first the index of the centroids
    print(",".join(f"{index}" for index in centroid_indexes))

    for centroid in centroids:
        # Print the elements of each centroids
        print(
            ",".join(
                "{:.4f}".format(element_of_centroid) for element_of_centroid in centroid
            )
        )

    # Now that we have chosen the initial centers, proceed using standard k-means clustering


def main():

    args = sys.argv[1:]
    num_of_args = len(args)
    k = args[0]
    if num_of_args == 4:
        iter = 300
        eps = args[1]
        input_text_1 = args[2]
        input_text_2 = args[3]
    else:
        iter = args[1]
        eps = args[2]
        input_text_1 = args[3]
        input_text_2 = args[4]

    # Check k integer validity
    try:
        k = float(k)
        if not k.is_integer():
            raise ValueError
        k = int(k)
    except ValueError:
        print("Invalid number of clusters!")
        sys.exit(1)

    # Check iteration validity
    try:
        iter = float(iter)
        if not iter.is_integer():
            raise ValueError
        iter = int(iter)
        if not 1 < iter < 1000:
            raise ValueError
    except ValueError:
        print("Invalid maximum iteration!")
        sys.exit(1)

    # Check epsilon validity
    try:
        if float(eps) < 0:
            raise ValueError
        eps = float(eps)
    except ValueError:
        print("Invalid epsilon!")
        sys.exit(1)

    kmeans_pp(k, input_text_1, input_text_2, iter, eps)


if __name__ == "__main__":
    main()
