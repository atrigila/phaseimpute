# Development

## Features and tasks

- [] Add automatic detection of chromosome name to create a renaming file for the vcf
- [] Make the different tests workflows work
  - []  Simulation
  - []  Validation
  - []  Preprocessing
  - [X] Imputation
  - []  Validation
  - []  Postprocessing
- [] Add support of `anyOf()` or `oneOf()` in the nf-core schema for the map, panel and region files
- [] Add nf-test for all modules and subworkflows
- [] Remove all TODOs
- [] Check if panel is necessary depending on the tool selected

## Run tests

```bash
nextflow run main.nf -profile singularity,test --outdir results -resume
```

## Problematic

### Channel management and combination

If only one specie at a time, then only one fasta file and only one map file (normally ?)
Do we want to be able to compute multiple panel at the same time ?
If so we need to correctly combine the different channel depending on their meta map.

All channel need to be identified by a meta map as follow:

- I : individual id
- P : panel id
- R : region used
- M : map used
- T : tool used
- G : reference genome used (is it needed ?)
- D : depth

## Open questions

How to use different schema ?

- Use nf-validation
  For the moment use different input / step.
  In the futur, if/else logic will be added in the yml nf-core schema.

What's the use of dumpcustomsoftware ?
Will be deleted

How to add to multiQC ?
Take exemple on Sarek.
All report file are in a dedicated channel.

How to add nf-test ?
Add in `tests` folder and run with tag.
Add tags.yml

How to run stub tests ?
Use nf-test

How to run the tests ?
nf-test option tag

What's the use of the template branch ?
TEMPLATE branch have the skeleton for all common part of the pipeline.
Will be asked to be merged to dev from time to time.

When is it necessary to merge to master / main ?
First release, create a false PR to first commit that will be checked by whole community + 2 reviewers approval.

What should be the Github action ?
All GA come from the TEMPLATE branch.
