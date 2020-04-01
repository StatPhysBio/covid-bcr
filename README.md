# covid-bcr
Analysis of BCR repertoires from covid-19 infections

# Raw sequence processing
1. Download [pRESTO](https://presto.readthedocs.io/en/stable/overview.html) at this [link](https://bitbucket.org/kleinstein/presto/downloads/).
2. Use the '''presto_pipeline.sh''' script to process the data. Pairs of raw sequences are assembled and filtered with a QScore of 30, V primers are masked and the C primer is cut, and finally sequences are deduplicated.
