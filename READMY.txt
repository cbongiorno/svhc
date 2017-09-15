
###############################################################
		INSTALLAZIONE
###########################################################

Entrare nella cartella:

	$ cd HClustValid-master/

Installare il software manager "pip" con:

	$ sudo easy_install pip ( oppure $ sudo brew install python ... o altri )


Quindi, installare tramite pip le seguenti librerie:

	$  sudo pip install scipy matplotlib numpy python-igraph pandas seaborn


(Spero che igraph non faccia errori .. in caso fatemi sapere)


Infine, installare la libreria hclustval:

	$ sudo python setup.py install

Ora il software è pronto per l'uso


########################################################
		DESCRIZIONE -- BENCHMARK FATTORIALE
###################################################

++++ Per generare i dati del benchmark:

	$ python hclustval/Benchmark.py examples/pattern_example.dat 1000 0

- il file "examples/pattern_example.dat" contiene la factor loading matrix del modello fattoriale. Ogni riga è un nodo, ogni colonna un fattore. I valori sono separati da tab

- 1000 è la lunghezza della serie temporale.

- 0 è il noise aggiunto, varia tra [0,1]. Come nell'equazione di sezione IV.A del report. (Questo parametro non è indispensabile, è stato introdotto per comodità, infatti si può ottenere lo stesso risultato modificando unicamente la factor loading matrix) 


Dopo avere eseguito il comando verranno creati 3 files:"dataSeries_benchmark.dat", "HyperNodeList_reference.dat", "CompressDendrogram_reference.dat".
1) Il file "dataSeries_benchmark.dat" contiene la serie di dati. Ogni riga è un nodo, ed i valori sono separati da tab.
2) Il file "HyperNodeList_reference.dat" contiene la partizione di riferimento. Ogni riga è un ipernodo, i valori separati da virgola sono i nodi-leave che l'ipernodo contiene.
3) Il file "CompressDendrogram_reference.dat" contiene un altro modo di rappresentare la partizione di riferimento. Ogni riga è un layer, ed ogni colonna un nodo. due nodi appartengono alla stessa comunità in un layer se hanno lo stesso valore.

##########################################################
		DESCRIZIONE -- BENCHMARK DI GARLASCHELLI
#########################################################

++++ Per generare i dati del benchmark:

	$ python hclustval/Benchmark_Garlaschelli.py examples/PartitionSize_example.dat 0.4 0.7 5000

- il file "examples/PartitionSize_example.dat" contiene le size delle partizioni. Sono separate da virgola.

- 0.4 è mu, il parametro di mercato

- 0.7 è ni, il parametro di noise

- 5000 è la lunghezza della serie temporale ( Il numero di nodi viene ricavato automaticamente dal partitionSize sommando i valori )


Dopo avere eseguito il comando verranno creati 3 files:"dataSeries_benchmark.dat", "HyperNodeList_reference.dat", "CompressDendrogram_reference.dat". La descrizione è la stessa di sopra.


#############################################################
		DESCRIZIONE -- CLUSTERING VALIDATO
############################################################

+++ per eseguire il metodo:

	$ python hclustval/MainLib.py dataSeries_benchmark.dat 0.05 1000 average


- dataSeries_benchmark.dat è la serie di dati del benchmark... o un altra.
- 0.05 è la soglia del single test.
- 1000 è il numero di copie di bootstrap. (per 100 nodi questi ultimi due valori funzionano molto bene)
- average è il metodo. Le alternative sono [ single, average, complete ]


Il metodo ritorna due file: "HyperNodeList_out.dat", "CompressDendrogram_out.dat". Sono spiegati in DESCRIZIONE-BENCHMARK.



#############################################################
		DESCRIZIONE -- METRICHE
############################################################

.. per calcolare la metriche:

	$ python hclustval/metric.py CompressDendrogram_out.dat CompressDendrogram_reference.dat 

Il programma stamperà su schermo i valori di HARI ed HAWI. In oltre, scriverà su file "metric_out.dat" con i valori. ATTENZIONE, l'ordine è importante perché l'HARI non è simmetrica.






