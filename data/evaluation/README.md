### Evaluation

Using WES estimates from Edgecopy's `agcn` module and WGS estimates from Parascopy, you can compare and evaluate the concordance, precision, and recall values. More specifically, we have the following commands: 

```
# Example command to create evaluation files
python evaluate.py \
--wes-dir wes-estimates \
--wgs-dir wgs-estimates \
--out-dir eval \
--refcn wes-estimates/summary.input.info.tsv \
--allexons out-depth/exons/allexons.parascopy.examine.output
```
where
* `--wes-dir` is the output directory from running `edgecopy agcn`
* `--wgs-dir` is the directory containing WGS estimates from Parascopy (e.g. benchmark-data/parascopy-wgs-calls)
* `--out-dir` is the directory you want to store the evaluation files.
* `--refcn` is a TSV file containing reference copy numbers of genes examined (located inside `--wes-dir`)
* `--allexons` is a comprehensive list of loci and reference copy numbers (located inside output directory from running `edgecopy depth`)


With the evaluation files, you can compute the performance by the following commands:
```
# Example command to compute concordance
python compute_concordance.py eval/common-cnvs/ wes-estimates/summary.input.info.tsv

# Example command to compute precision and recall 
python compute_precision_recall.py eval/rare-cnvs/ wes-estimates/summary.input.info.tsv
```

where you want to pass in (1) the directory containing genes with common/rare CNVs and (2) reference copy number TSV file (as previously described). 

