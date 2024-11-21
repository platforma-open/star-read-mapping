//  import { PlTableState } from "@milaboratory/pl-table";
import {
  BlockModel,
  InferOutputsType,
  isPColumn,
  isPColumnSpec,
  Ref,
  ValueType,
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
   * Species settings
   */
  species?: "homo-sapiens" | "mus-musculus" | "saccharomyces-cerevisiae";

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

export const model = BlockModel.create()

  .withArgs<BlockArgs>({
    species: "homo-sapiens",
    strandness: "0",
  })

  /**
   * Find possible options for the fastq input
   */
  .output("dataOptions", (ctx) => {
    return ctx.resultPool.getOptions((v) => {
      if (!isPColumnSpec(v)) return false;
      const domain = v.domain;
      return (
        v.name === "pl7.app/sequencing/data" &&
        (v.valueType as string) === "File" &&
        domain !== undefined &&
        (domain["pl7.app/fileExtension"] === "fastq" ||
          domain["pl7.app/fileExtension"] === "fastq.gz")
      );
    });
  })

  .output("labels", (ctx): Record<string, string> | undefined => {
    const inputRef = ctx.args.ref;
    if (inputRef === undefined) return undefined;

    const inputSpec = ctx.resultPool.getSpecByRef(inputRef);

    if (inputSpec === undefined || !isPColumnSpec(inputSpec)) return undefined;

    const sampleLabelsObj = ctx.resultPool.findDataWithCompatibleSpec({
      kind: "PColumn",
      name: "pl7.app/label",
      valueType: "String",
      axesSpec: [inputSpec.axesSpec[0]], // samplesId axis
      domain: inputSpec.domain,
    });

    if (sampleLabelsObj.length === 0) return undefined;

    // @TODO implement standard method for getting labels
    const labels = Object.fromEntries(
      Object.entries(
        sampleLabelsObj[0].data.getDataAsJson<{
          data: Record<string, string>;
        }>().data
      ).map((e) => [JSON.parse(e[0])[0], e[1]])
    ) satisfies Record<string, string>;

    return labels;
  })

  /**
   * Preprocessing progress
   */
  .output("starProgress", (wf) => {
    return parseResourceMap(
      wf.outputs?.resolve("starProgress"),
      (acc) => acc.getLogHandle(),
      false
    );
  })

  /**
   * Last line from StAR output
   */
  .output("starProgressLine", (wf) => {
    return parseResourceMap(
      wf.outputs?.resolve("starProgress"),
      (acc) => acc.getLastLogs(1),
      false
    );
  })

  .output("starQc", (wf) =>
    parseResourceMap(
      wf.outputs?.resolve("starQc"),
      (acc) => acc.getFileContentAsString(),
      false
    )
  )

  .output("alignedBAM", (wf) =>
    wf.outputs?.resolve("alignedBAM")?.getLastLogs(1)
  )

  .output("featureCountsProgress", (wf) => {
    return parseResourceMap(
      wf.outputs?.resolve("featureCountsProgress"),
      (acc) => acc.getLogHandle(),
      false
    );
  })
  .output("featureCountsQc", (wf) =>
    parseResourceMap(
      wf.outputs?.resolve("featureCountsQc"),
      (acc) => acc.getFileContentAsString(),
      false
    )
  ) // Does this work with this type of file?

  /**
   * P-frame with rawCounts
   */
  .output("pf", (wf) => {
    //return wf.outputs?.resolve("pf")?.resolve("rawCounts.data")?.listInputFields()
    const pCols = wf.outputs?.resolve("pf")?.getPColumns();
    if (pCols === undefined) return undefined;

    return wf.createPFrame(pCols);
  })

  .output("pcaPf", (wf) => {
    //return wf.outputs?.resolve("pf")?.resolve("rawCounts.data")?.listInputFields()
    const pCols = wf.outputs?.resolve("pcaComponents")?.getPColumns();
    if (pCols === undefined) return undefined;

    //return wf.createPFrame(pCols);
    // enriching with upstream data
    const valueTypes = [
      "Int",
      "Long",
      "Float",
      "Double",
      "String",
      "Bytes",
    ] as ValueType[];
    const upstream = wf.resultPool
      .getData()
      .entries.map((v) => v.obj)
      .filter(isPColumn)
      .filter((column) =>
        valueTypes.find((valueType) => valueType === column.spec.valueType)
      );

    return wf.createPFrame([...pCols, ...upstream]);
  })

  .sections([
    { type: "link", href: "/", label: "Settings" },
    //    { type: "link", href: "/QC", label: "Sequence Data QC" },
    { type: "link", href: "/PCA", label: "Principal Component Analysis" },
    //    { type: "link", href: "/HC", label: "Hierarchical Clustering" },
  ])

  .done();

export type BlockOutputs = InferOutputsType<typeof model>;
