---
'@platforma-open/milaboratories.star-read-mapping.workflow': patch
---

Fix non-deterministic R1/R2 order in STAR `--readFilesIn` for single-key (lane-less) input. Read-file args were added during a `for sKey, inputFile in inputData.inputs()` loop, which has non-deterministic iteration order. When STAR received `--readFilesIn input_R2.fq.gz input_R1.fq.gz`, it either failed with `wrong read ID line format` or produced a BAM that featureCounts then rejected with `No paired-end reads were detected in paired-end read library`. Args are now added in a separate loop in fixed `R1` then `R2` order, mirroring the multi-lane (keyLength==2) path.
