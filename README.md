# Statistically Validated Hierarchical Clustering

## Installation

### On Ubuntu Linux

```
$ sudo apt-get install pip python-tk
$ sudo pip install svhc
```

### On Mac OsX

The easist way to install the library is with easy_install (or with port).

```
$ sudo easy_install install pip
$ sudo pip install svhc
```

if the last command raise conflicts, then you should try:

```
$ sudo pip install --ignore-installed svhc
```

Then you can try to call the script:
```
$ svhc
```

if it does not find the script, you can add the path temporaneously with the command:
```
export PATH="/Users/$USER/Library/Python/2.7/bin:$PATH"
```

or permanently with (da provare):
```
$ echo "export PATH=$PATH:/Users/$USER/Library/Python/2.7/bin" >>  ~/.bash_profile
```
### on Windows
Install python2.7 and pip by following the instruction at the [link](https://pip.pypa.io/en/stable/installing/#do-i-need-to-install-pip). Then from the prompt

```
pip install svhc
```



## Usage

### Generate a Benchmark
If you want generate a dataset starting from a factor model, then you need a factor loading matrix. It is possible to download an example of such matrix from  the link ([pattern_example.dat](https://github.com/cbongiorno/svhc/blob/master/svhc/example/pattern_example.dat)).
Then, to generate the data series run:

```
$ svhc_benchmark pattern.dat 500 test
```
where 500 is the lenght of the data series. Test is the output name.  The program will produce **test_dataSeries_benchmark.dat**, that is the data matrix (a matrix 100x500), and  **test_cluster_reference.dat**, that is the list of the nodes that belong to each cluster; each line is a different cluster, nodes are comma separated.

If you want to add noise to the data you can use the optional parameter:

```
$ svhc_benchmark pattern.dat 500 test --noise 0.3
```


### Evaluate Statistically Validated Hierachical Clusters

To estimate the validated clusters on the dataseries generated in the previous example, you can run:

```
$ svhc test_dataSeries_benchmark.dat 1000 test
```

where **1000** is the number of bootstrap copies, and **test** is the output name. If you want to use your own dataset, please remember to store your object by row, and your attribute by colomn tab separated and without header of index numbers. Then the algoithm will find the clusters of objects.

Few optional parameters are allowed:
```
$ svhc test_dataSeries_benchmark.dat 1000 test --alpha 0.5 --nan 0 --ncpu 1
```
where **alpha** is the confidence of the FDR multiple comparison correction (default 0.05); **nan** is a boolean entry that must be fixed to 1 if there are NaN entries in the dataset; **ncpu** is the number of core used for the evaluation of the bootstrap copies (default 1).


The program will produce three files: **test_Validated_Cluster.dat**, that is the list (by row) of validated clusters, each row is the list of nodes within the cluster (comma separated); **prova_pvalue.dat** that is the list of pvalues associated to each cluster; **test_dendrogram.dat** that contain the full information of the dendrogram.

### Plot the Dendrogram

To plot the dendrogram with the output of the previous example please run:

```
$ svhc_plot test_dataSeries_benchmark.dat test_Validated_Cluster.dat test_dendrogram.dat test_pic
```

the program will procude a **test_pic.pdf**






