//  import { PlTableState } from "@milaboratory/pl-table";
import {
  BlockModel,
  ImportFileHandle,
  InferOutputsType,
  Option,
  Ref,
  isPColumnSpec
} from "@platforma-sdk/model";
import { parseResourceMap } from "./helpers";

/**
 * Block arguments coming from the user interface
 */
export type BlockArgs = {
  /**
   * Reference to the fastq data
   */
  ref?: Ref;

  /**
   * Genome index upload
   */
  indexFile?: ImportFileHandle;

  /**
   * Genome index upload
   */
  genomeAnnFile?: ImportFileHandle;

  /**
   * Species settings
   */
  species?: "hsa" | "mmu";

  /**
   * Library settings
   */
  libraryType?: "SingleEnd" | "PairedEnd";

  /**
   * Strandness settings
   */
  strandness?: "0" | "1" | "2";
};

// /**
//  * UI state
//  */
// export type UiState = {
//   // selected chain
//   libraryTypeOptions?: string;
// };


export const model = BlockModel.create<BlockArgs>("Heavy")

  .initialArgs({
    species: "hsa",
    libraryType: "SingleEnd",
    strandness: "0",
  })

  /**
   * Find possible options for the fasta input
   */
  .output("dataOptions", (ctx) => {
    return ctx.resultPool
      .getSpecsFromResultPool()
      .entries.filter((v) => {
        if (!isPColumnSpec(v.obj)) return false;
        const domain = v.obj.domain;
        return (
          v.obj.name === 'pl7.app/sequencing/data' &&
          (v.obj.valueType as string) === 'File' &&
          domain !== undefined &&
          (domain['pl7.app/fileExtension'] === 'fastq' ||
            domain['pl7.app/fileExtension'] === 'fastq.gz')
        );
      })
      .map(
        (v) =>
          ({
            ref: v.ref,
            label: `${ctx.getBlockLabel(v.ref.blockId)} / ${
              v.obj.annotations?.['pl7.app/label'] ?? `unlabelled`
            }`
          } satisfies Option)
      );
  })

  .output("indexUploadProgress", (wf) => 
    wf.outputs?.resolve("indexImportHandle")?.getImportProgress()
  )

  .output("genomeAnnUploadProgress", (wf) => 
    wf.outputs?.resolve("genomeAnnImportHandle")?.getImportProgress()
  )

  /**
   * Preprocessing progress
   */
  //.output("starProgress", (wf) => wf.outputs?.resolve("starProgress")?.getLastLogs(100))
  .output("starProgress", (wf) => {
    return parseResourceMap(wf.outputs?.resolve({ field: 'starProgress', assertFieldType: 'Input' }), (acc) =>
      acc.getLastLogs(100),
      false
    )
  })
  .output("starQc", (wf) => wf.outputs?.resolve("starQc")?.getLastLogs(100)) // Does this work with this type of file?
  .output("featureCountsProgress", (wf) => wf.outputs?.resolve("featureCountsProgress")?.getLastLogs(100))
  .output("featureCountsQc", (wf) => wf.outputs?.resolve("featureCountsQc")?.getLastLogs(100)) // Does this work with this type of file?
  

  /**
   * P-frame with rawCounts
   */
  .output("pf", (wf) => {
    const pCols = wf.outputs?.resolve("pf")?.getPColumns();
    if (pCols === undefined) return undefined;

    return wf.createPFrame(pCols);
  })

  /**
   * P-table with counts output for all samples
   */
  .output("pt", (wf) => {
    const pCols = wf.outputs?.resolve("pf")?.getPColumns();
    if (pCols === undefined || pCols.length === 0) return undefined;

    return wf.createPTable({
      columns: pCols,
    });
  })

  .sections([
    { type: "link", href: "/", label: "Settings" },
    { type: "link", href: "/QC", label: "Sequence Data QC" },
    { type: "link", href: "/PCA", label: "Principal Component Analysis" },
    { type: "link", href: "/HC", label: "Hierarchical Clustering" },
  ])

  .done();

export type BlockOutputs = InferOutputsType<typeof model>;
