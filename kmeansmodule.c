#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef INFINITY
#define INFINITY (1.0 / 0.0)
#endif

struct cord
{
    double value;
    struct cord *next;
};

struct vector
{
    struct vector *next;
    struct cord *cords;
};

void sum_clusters_reset(double **sum_clusters, int num_clusters, int vector_length)
{
    int i;
    int j;
    for (i = 0; i < num_clusters; i++)
    {
        for (j = 0; j < vector_length; j++)
        {
            sum_clusters[i][j] = 0;
        }
    }
}

void counters_reset(int *counters, int num_clusters)
{
    int i;
    for (i = 0; i < num_clusters; i++)
    {
        counters[i] = 0;
    }
}

int have_centroids_changed(double *delta_centroids, double eps, int num_clusters)
{
    int i;
    for (i = 0; i < num_clusters; i++)
    {
        if (delta_centroids[i] >= eps)
        {
            return 1;
        }
    }

    return 0;
}

double euclidean_distance(double *cord1, double *cord2, int vector_length)
{
    double sum = 0;
    int i;
    double sqrt_sum;
    double distance;
    double squared_cord_distance;
    for (i = 0; i < vector_length; i++)
    {
        distance = cord1[i] - cord2[i];
        squared_cord_distance = pow(distance, 2);
        sum += squared_cord_distance;
    }
    sqrt_sum = sqrt(sum);
    return sqrt_sum;
}

int find_closest_centroid(struct cord *curr_cord, double **centroids, int num_clusters, int vector_length)
{
    double *curr_cord_array = malloc(vector_length * sizeof(double));
    int cord_index = 0;
    double min_index;
    double min_distance;
    int i;

    while (curr_cord != NULL)
    {
        curr_cord_array[cord_index] = curr_cord->value;
        curr_cord = curr_cord->next;
        cord_index++;
    }
    min_distance = euclidean_distance(curr_cord_array, centroids[0], vector_length);
    min_index = 0;
    for (i = 1; i < num_clusters; i++)
    {
        double current_distance = euclidean_distance(curr_cord_array, centroids[i], vector_length);
        if (current_distance < min_distance)
        {
            min_distance = current_distance;
            min_index = i;
        }
    }
    free(curr_cord_array);
    return min_index;
}

static PyObject *fit(PyObject *self, PyObject *args)
{

    PyObject *datapoints;
    PyObject *datapoint;
    PyObject *numpy_centroids;
    PyObject *centroid;
    PyObject *returned_py_val;
    double **centroids;
    double eps;
    int num_clusters, iter, n_samples, n_features;
    int i, j;

    if (!PyArg_ParseTuple(args, "OOdiiii", &numpy_centroids, &datapoints, &eps, &num_clusters, &iter, &n_samples, &n_features))
    {
        return NULL;
    }

    centroids = (double **)malloc(k * sizeof(double *));
    for (i = 0; i < k; i++)
    {
        numpy_centroids = PyList_GetItem(numpy_centroids, i);
        centroids[i] = (double *)malloc(d * sizeof(double));
        if (centroids[i] == NULL)
        {
            for (j = 0; j < i; j++)
            {
                free(centroids[j]);
            }
            free(centroids);
            printf("An Error Has Occurred\n");
            return NULL;
        }
        for (j = 0; j < d; j++)
        {
            item = PyList_GetItem(center, j);
            entry = PyFloat_AsDouble(item);
            centroids[i][j] = entry;
        }
    }

    // implemntation of K-means clustering here
    int i;
    int j;
    double **sum_clusters;
    int *counters;
    double delta;
    int num_clusters;
    int vector_length;
    int curr_iter;
    struct vector *head_vec, *curr_vec, *next_vec;
    struct cord *head_cord, *curr_cord, *next_cord;
    double n;
    char c;
    double *delta_centroids;

    head_vec = malloc(sizeof(struct vector));
    curr_vec = head_vec;
    curr_vec->next = NULL;

    for (i = 0; i < n_samples; i++)
    {
        curr_vec->cords = head_cord;
        curr_vec->next = malloc(sizeof(struct vector));
        head_cord = malloc(sizeof(struct cord));
        curr_cord = head_cord;
        curr_cord->next = NULL;

        datapoint = PyList_GetItem(datapoints, i);

        for (j = 0; j < n_features; j++) /*iterating over a single vector's coords*/
        {
            item = PyList_GetItem(datapoint, j);
            cord = PyFloat_AsDouble(item);
            curr_cord->value = cord;
            curr_cord->next = malloc(sizeof(struct cord));
            curr_cord = curr_cord->next;
            curr_cord->next = NULL;
        }
        curr_vec = curr_vec->next;
    }
    curr_vec->next = NULL;

    sum_clusters = malloc(num_clusters * sizeof(double *));
    counters = malloc(num_clusters * sizeof(int));

    for (i = 0; i < num_clusters; i++)
    {
        sum_clusters[i] = malloc(vector_length * sizeof(double));
        for (j = 0; j < vector_length; j++)
        {
            sum_clusters[i][j] = 0;
        }
    }

    curr_iter = 0;

    delta_centroids = malloc(num_clusters * sizeof(double));

    for (i = 0; i < num_clusters; i++)
    {
        delta_centroids[i] = INFINITY;
    }

    while (curr_iter < num_iterations && have_centroids_changed(delta_centroids, eps, num_clusters))
    {
        sum_clusters_reset(sum_clusters, num_clusters, vector_length);
        counters_reset(counters, num_clusters);

        curr_vec = head_vec;
        while (curr_vec != NULL && curr_vec->cords != NULL)
        {
            struct cord *curr_cord = curr_vec->cords;

            int closest_centroid = find_closest_centroid(curr_cord, centroids, num_clusters, vector_length);

            for (j = 0; j < vector_length; j++)
            {
                sum_clusters[closest_centroid][j] += curr_cord->value;
                curr_cord = curr_cord->next;
            }
            counters[closest_centroid]++;
            curr_vec = curr_vec->next;
        }

        for (i = 0; i < num_clusters; i++)
        {
            double *temp_centroid = malloc(vector_length * sizeof(double));
            for (j = 0; j < vector_length; j++)
            {
                temp_centroid[j] = centroids[i][j];
            }

            for (j = 0; j < vector_length; j++)
            {
                centroids[i][j] = sum_clusters[i][j] / counters[i];
            }
            delta = euclidean_distance(temp_centroid, centroids[i], vector_length);
            delta_centroids[i] = delta;
            free(temp_centroid);
        }

        curr_iter++;
    }

    PyObject *py_centroids = PyList_New(k);
    for (i = 0; i < k; i++)
    {
        PyObject *py_centroid = PyList_New(d);
        for (j = 0; j < d; j++)
        {
            PyList_SetItem(py_centroid, j, PyFloat_FromDouble(centroids[i][j]));
        }
        PyList_SetItem(py_centroids, i, py_centroid);
    }
    returned_py_val = Py_BuildValue("O", py_centroids);

    for (i = 0; i < num_clusters; i++)
    {
        free(centroids[i]);
    }
    free(centroids);

    curr_vec = head_vec;
    while (curr_vec != NULL)
    {
        curr_cord = curr_vec->cords;
        while (curr_cord != NULL)
        {
            next_cord = curr_cord->next;
            free(curr_cord);
            curr_cord = next_cord;
        }
        next_vec = curr_vec->next;
        free(curr_vec);
        curr_vec = next_vec;
    }

    for (i = 0; i < num_clusters; i++)
    {
        free(sum_clusters[i]);
    }
    free(sum_clusters);
    free(counters);
    free(delta_centroids);

    return Py_BuildValue("O", returned_py_val);
}

static PyMethodDef Kmeans_FuncTable[] = {
    {"fit",
     fit,
     METH_VARARGS,
     "Calculates centroids for k means algorithm"},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef Kmeans_Mod = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    "Python wrapper for C kmeans algorithm",
    -1,
    Kmeans_FuncTable};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    return PyModule_Create(&Kmeans_Mod);
}
