# Statistically Validated Hierarchical Clustering

## Installation

### On Ubuntu Linux

```
$ sudo apt-get install pip python-tk
$ sudo pip install svhc
```

### On Mac OsX

The easist way to install the library is by using easy_install.

```
$ sudo easy_install install pip
$ echo "export PATH=$PATH:/Users/$USER/Library/Python/2.7/bin" >>  ~/.bash_profile
$ sudo pip install svhc
```

if the last command raise some conflict, you can try to run:

```
$ sudo pip install --ignore-installed svhc
```


```
export PATH="/Users/$USER/Library/Python/2.7/bin:$PATH"
```

### on Windows
Install python2.7 and pip by following the instruction at the [link](https://pip.pypa.io/en/stable/installing/#do-i-need-to-install-pip). Then from the prompt

```
pip install svhc
```



## Usage

### Generate a Benchmark
If you want generate a dataset starting from a factor model you need a factor loading matrix. It is possible to download an example of the matrix [pattern_example.dat](https://github.com/cbongiorno/svhc/blob/master/svhc/example/pattern_example.dat).
Then, to generate the data series run:

```
$ svhc_benchmark pattern.dat 500 test
```
where 500 is the lenght of the data series. Test is the output name.  The program will produce "test_dataSeries_benchmark.dat", that is the data matrix ( a matrix 100x500), and  "test_cluster_reference.dat", that is the list of the nodes that belong to each cluster; each line is a different cluster, nodes are comma separated.

If you want to add noise to the data you can use the optional parameter:

```
$ svhc_benchmark pattern.dat 500 test --noise 0.3
```


### Evaluate Statistically Validated Hierachical Clusters

To estimate the validated clusters on the dataseries generated in the previous example you can run:

```
$ svhc test_dataSeries_benchmark.dat 1000 test
```

where 1000 is the number of bootstrap copies, and test is the output name. If you want to use your own dataset, please remember to store your object by row, and your attribute by colomn tab separated and without header of index numbers. Then the algoithm will find the clusters of objects.


