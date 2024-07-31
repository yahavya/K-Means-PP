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

void sum_clusters_reset(double **sum_clusters, int num_clusters, int n_features)
{
    int i, j;
    for (i = 0; i < num_clusters; i++)
    {
        for (j = 0; j < n_features; j++)
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

double euclidean_distance(double *cord1, double *cord2, int n_features)
{
    double sum = 0;
    int i;
    for (i = 0; i < n_features; i++)
    {
        double distance = cord1[i] - cord2[i];
        sum += distance * distance;
    }
    return sqrt(sum);
}

int find_closest_centroid(struct cord *curr_cord, double **centroids, int num_clusters, int n_features)
{
    double *curr_cord_array = malloc(n_features * sizeof(double));
    if (curr_cord_array == NULL)
    {
        printf("Memory allocation failed for curr_cord_array\n");
        exit(1);
    }

    int cord_index = 0;
    while (curr_cord != NULL)
    {
        curr_cord_array[cord_index] = curr_cord->value;
        curr_cord = curr_cord->next;
        cord_index++;
    }

    double min_distance = euclidean_distance(curr_cord_array, centroids[0], n_features);
    int min_index = 0;
    for (int i = 1; i < num_clusters; i++)
    {
        double current_distance = euclidean_distance(curr_cord_array, centroids[i], n_features);
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
    PyObject *datapoints, *numpy_centroids;
    double eps;
    int num_clusters, iter, n_samples, n_features;

    if (!PyArg_ParseTuple(args, "OOdiiii", &numpy_centroids, &datapoints, &eps, &num_clusters, &iter, &n_samples, &n_features))
    {
        return NULL;
    }


    double **centroids = (double **)malloc(num_clusters * sizeof(double *));
    if (centroids == NULL)
    {
        printf("Memory allocation failed for centroids\n");
        return NULL;
    }

    for (int i = 0; i < num_clusters; i++)
    {
        PyObject *numpy_centroid = PyList_GetItem(numpy_centroids, i);
        centroids[i] = (double *)malloc(n_features * sizeof(double));
        if (centroids[i] == NULL)
        {
            for (int j = 0; j < i; j++)
            {
                free(centroids[j]);
            }
            free(centroids);
            printf("Memory allocation failed for centroid %d\n", i);
            return NULL;
        }
        for (int j = 0; j < n_features; j++)
        {
            PyObject *item = PyList_GetItem(numpy_centroid, j);
            centroids[i][j] = PyFloat_AsDouble(item);
        }
    
    }

    double **sum_clusters = (double **)malloc(num_clusters * sizeof(double *));
    if (sum_clusters == NULL)
    {
        for (int i = 0; i < num_clusters; i++)
        {
            free(centroids[i]);
        }
        free(centroids);
        printf("Memory allocation failed for sum_clusters\n");
        return NULL;
    }

    for (int i = 0; i < num_clusters; i++)
    {
        sum_clusters[i] = (double *)malloc(n_features * sizeof(double));
        if (sum_clusters[i] == NULL)
        {
            for (int j = 0; j < i; j++)
            {
                free(sum_clusters[j]);
            }
            for (int j = 0; j < num_clusters; j++)
            {
                free(centroids[j]);
            }
            free(centroids);
            free(sum_clusters);
            printf("Memory allocation failed for sum_clusters[%d]\n", i);
            return NULL;
        }
    }

    int *counters = (int *)malloc(num_clusters * sizeof(int));
    if (counters == NULL)
    {
        for (int i = 0; i < num_clusters; i++)
        {
            free(centroids[i]);
            free(sum_clusters[i]);
        }
        free(centroids);
        free(sum_clusters);
        printf("Memory allocation failed for counters\n");
        return NULL;
    }

    double *delta_centroids = (double *)malloc(num_clusters * sizeof(double));
    if (delta_centroids == NULL)
    {
        for (int i = 0; i < num_clusters; i++)
        {
            free(centroids[i]);
            free(sum_clusters[i]);
        }
        free(centroids);
        free(sum_clusters);
        free(counters);
        printf("Memory allocation failed for delta_centroids\n");
        return NULL;
    }

    for (int i = 0; i < num_clusters; i++)
    {
        delta_centroids[i] = INFINITY;
    }

    struct vector *head_vec = (struct vector *)malloc(sizeof(struct vector));
    if (head_vec == NULL)
    {
        for (int i = 0; i < num_clusters; i++)
        {
            free(centroids[i]);
            free(sum_clusters[i]);
        }
        free(centroids);
        free(sum_clusters);
        free(counters);
        free(delta_centroids);
        printf("Memory allocation failed for head_vec\n");
        return NULL;
    }
    head_vec->next = NULL;

    struct vector *curr_vec = head_vec;
    for (int i = 0; i < n_samples; i++)
    {
        curr_vec->cords = (struct cord *)malloc(sizeof(struct cord));
        if (curr_vec->cords == NULL)
        {
            printf("Memory allocation failed for curr_vec->cords\n");
            return NULL;
        }
        struct cord *curr_cord = curr_vec->cords;
        curr_cord->next = NULL;

        PyObject *datapoint = PyList_GetItem(datapoints, i);
       

        for (int j = 0; j < n_features; j++)
        {
            PyObject *item = PyList_GetItem(datapoint, j);
            curr_cord->value = PyFloat_AsDouble(item);
            if (j < n_features - 1)
            {
                curr_cord->next = (struct cord *)malloc(sizeof(struct cord));
                if (curr_cord->next == NULL)
                {
                    printf("Memory allocation failed for curr_cord->next\n");
                    return NULL;
                }
                curr_cord = curr_cord->next;
                curr_cord->next = NULL;
            }
        }
        if (i < n_samples - 1)
        {
            curr_vec->next = (struct vector *)malloc(sizeof(struct vector));
            if (curr_vec->next == NULL)
            {
                printf("Memory allocation failed for curr_vec->next\n");
                return NULL;
            }
            curr_vec = curr_vec->next;
            curr_vec->next = NULL;
        }
    }

    int curr_iter = 0;
    while (curr_iter < iter && have_centroids_changed(delta_centroids, eps, num_clusters))
    {
        sum_clusters_reset(sum_clusters, num_clusters, n_features);
        counters_reset(counters, num_clusters);

        curr_vec = head_vec;
        while (curr_vec != NULL && curr_vec->cords != NULL)
        {
            struct cord *curr_cord = curr_vec->cords;
            int closest_centroid = find_closest_centroid(curr_cord, centroids, num_clusters, n_features);

            for (int j = 0; j < n_features; j++)
            {
                sum_clusters[closest_centroid][j] += curr_cord->value;
                curr_cord = curr_cord->next;
            }
            counters[closest_centroid]++;
            curr_vec = curr_vec->next;
        }

        for (int i = 0; i < num_clusters; i++)
        {
            double *temp_centroid = (double *)malloc(n_features * sizeof(double));
            if (temp_centroid == NULL)
            {
                printf("Memory allocation failed for temp_centroid\n");
                return NULL;
            }
            for (int j = 0; j < n_features; j++)
            {
                temp_centroid[j] = centroids[i][j];
                centroids[i][j] = sum_clusters[i][j] / counters[i];
            }
            delta_centroids[i] = euclidean_distance(temp_centroid, centroids[i], n_features);
            free(temp_centroid);
        }
        curr_iter++;
    }

    PyObject *py_centroids = PyList_New(num_clusters);
    for (int i = 0; i < num_clusters; i++)
    {
        PyObject *py_centroid = PyList_New(n_features);
        for (int j = 0; j < n_features; j++)
        {
            PyList_SetItem(py_centroid, j, PyFloat_FromDouble(centroids[i][j]));
        }
        PyList_SetItem(py_centroids, i, py_centroid);
    }

    for (int i = 0; i < num_clusters; i++)
    {
        free(centroids[i]);
        free(sum_clusters[i]);
    }
    free(centroids);
    free(sum_clusters);
    free(counters);
    free(delta_centroids);

    curr_vec = head_vec;
    while (curr_vec != NULL)
    {
        struct cord *curr_cord = curr_vec->cords;
        while (curr_cord != NULL)
        {
            struct cord *next_cord = curr_cord->next;
            free(curr_cord);
            curr_cord = next_cord;
        }
        struct vector *next_vec = curr_vec->next;
        free(curr_vec);
        curr_vec = next_vec;
    }

    return Py_BuildValue("O", py_centroids);
}

static PyMethodDef Kmeans_FuncTable[] = {
    {"fit", fit, METH_VARARGS, "Calculates centroids for k means algorithm"},
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