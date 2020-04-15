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

[Prototype notebook](error_correct.ipynb) for clustering singletons (and other low-frequency sequences) into large clones if similar in sequence.
Command line usage can be displayed as follows:
```bash
$ python error_correct.py -h
usage: error_correct.py [-h] [--delta_r DELTA_R] [--delta_a DELTA_A]
                        [--passes PASSES] [--keep_singletons]
                        fasta outbase

FASTA error correction, streams to stdout

positional arguments:
  fasta              path to FASTA
  outbase            basename for output files

optional arguments:
  -h, --help         show this help message and exit
  --delta_r DELTA_R  marginal Hamming distance tolerance per decade in log
                     ratio abundances (default 1)
  --delta_a DELTA_A  marginal abundance tolerance of clusterable sequences per
                     decade in log ratio abundances (default 1)
  --passes PASSES    number of times to repeat greedy clustering (default 1)
  --keep_singletons  don't discard uncorrected singletons
```

- A marginal hamming threshold x means you need a clone 10^(d/x) times larger in abundance to absorb a sequence with at hamming distance d.
- Marginal abundance threshold y means you need a clone 10^(a/y) times larger in abundance to absorb a sequence with abundance a
- Both of the above conditions must be met for absorbtion to occur.
- Error correction is performed greedily (and thus approximately) by working down a list ranked by abundance.
