Meeting 30/09/2025

State of project: 
- Trying to make sense of snipar output (imputed genotypes)
- Half-read paper Anna sent about some ways to measure imputation accuracy
- Comparing real vs predicted genotypes, but anecdotally (i.e. looking at specific individuals w given family structure, & specific snp-s, but no overall metric of accuracy etc)
- Problem: I did find occasions where snipar predicted e.g. 1.0 (so no uncertainty) whereas the original genotype was reference allele homozygous which in snipar shoudl be 2.0. When I looked (at 2 snp-s only though), this happened where all relatives were homozygous for the ref allele so snipar had no way of telling based on snp-s only whether mother was 1.0 or 2.0. Maybe to do with the estimated ibd segments between sibs which imputation is based on? 
- Also, NA values in imputation : also in cases (from a sample size of 2) where all relatives homozygous for the ref allele.

Discussed during meeting // To do:
- A few suggestions for how to get a more systematic measure of imputation accuracy & errors: Per individual: confusion matrix for genotypes that are predicted 0.0/1.0/2.0; correlation matrix for all predictions ; Overall: correlation of imputed & real per snp for all individuals; filter based on number of children to see if there is difference
- Run snipar on 1 family with only supplying pedigree of this family so that it's sure that only these genotypes are used for imputation. Check if better/worse/no change for this family.
- Potentially add genetic map in a later run. Snipar can accept this file as input but I didn't know we had it so ran without, using a default one. Might improve imputation slightly
- Also get basic values like % of NA-s per individual (how many snps it didnt impute)
% of incorrectly imputed snp-s (where snipar says e.g. 1.0 while its actually 2.0 etc)

