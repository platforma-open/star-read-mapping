---
'@platforma-open/milaboratories.star-read-mapping.workflow': patch
---

Set STAR memory request per genome instead of the previous hardcoded `mem("4Gb")`. STAR's memory footprint is dominated by mmap of the suffix-array index, which scales with genome size; the 4 GB request was below what mammalian genomes need (~30 GB resident for human GRCh38) and caused STAR (and on localfs deployments often the backend itself) to be OOM-killed mid-run, surfacing as `Emergency bootstrap mode enabled - command marked as failed to escape failure loop` after the next backend start.

Per-species sizing (chosen for ~1.5× genome size + overhead):

| genome | mem |
| --- | --- |
| `saccharomyces-cerevisiae`, `test-species` | 2 GiB |
| `drosophila-melanogaster`, `arabidopsis-thaliana`, `caenorhabditis-elegans` | 4 GiB |
| `danio-rerio`, `gallus-gallus` | 16 GiB |
| `homo-sapiens`, `mus-musculus`, `rattus-norvegicus`, `bos-taurus`, `sus-scrofa` (and unknown) | 32 GiB |

Pinning per genome rather than relying on heavy-queue defaults insulates the block from future changes to the queue's default memory allocation.
