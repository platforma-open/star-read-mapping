---
'@platforma-open/milaboratories.star-read-mapping.workflow': patch
---

Bump STAR memory request from 4 GB to 32 GiB. Loading the human GRCh38 STAR index requires ~30 GB resident; the previous 4 GB allocation caused the STAR process (and on localfs deployments often the backend itself) to be OOM-killed mid-run, surfacing as `Emergency bootstrap mode enabled - command marked as failed to escape failure loop` after backend restart. 32 GiB matches what STAR typically needs for mammalian genomes; on smaller hosts the queue scheduler caps the request at the queue limit, and small genomes (yeast, test species) use only a fraction of the request.
