# covid-bcr
Analysis of BCR repertoires from covid-19 infections

Dependencies
---

Dependences are listed in [`env.yml`](). You can set up a [conda](https://docs.conda.io/en/latest/) environment with
```bash
$ conda env create -f env.yml
```
and activate the new environment (including all dependencies) with
```bash
$ conda activate covid-bcr
```

Pipeline
---

### Raw sequence processing

Use the [`presto_pipeline.sh`]() script to process the data with [pRESTO](https://presto.readthedocs.io/en/stable/overview.html). Pairs of raw sequences are assembled and filtered with a QScore of 30, V primers are masked and the C primer is cut, and finally sequences are deduplicated.

### Error correction

[Prototype notebook](error_correct.ipynb) for clustering singletons into large clones if similar in sequence
