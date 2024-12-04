//  import { PlTableState } from "@milaboratory/pl-table";
import {
  BlockModel,
  InferOutputsType,
  isPColumn,
  isPColumnSpec,
  parseResourceMap,
  PlRef,
  ValueType,
} from "@platforma-sdk/model";

import { GraphMakerState } from "@milaboratories/graph-maker";
/**
 * Block arguments coming from the user interface
 */
export type BlockArgs = {
  /**
   * Reference to the fastq data
   */
  ref?: PlRef;

  /**
   * Species settings
   */
  species?: "homo-sapiens" | "mus-musculus" | "saccharomyces-cerevisiae";

  /**
   * Strandness settings
   */
  strandness?: "0" | "1" | "2";

  /**
   * Block title
   */
  title?: string;
};

/**
 * UI state
 */
export type UiState = {
  pcaGraphState: GraphMakerState;
};

export const model = BlockModel.create()

  .withArgs<BlockArgs>({
    species: "homo-sapiens",
    strandness: "0",
  })

  .withUiState<UiState>({
    pcaGraphState: {
      template: "dots",
      title: "Principal Components Analysis",
    },
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

  .output("labels", (ctx) => {
    const inputRef = ctx.args.ref;
    if (inputRef === undefined) return undefined;

    const inputSpec = ctx.resultPool.getSpecByRef(inputRef); // @TODO use resultPool.getPColumnSpecByRef after updating SDK
    if (inputSpec === undefined || !isPColumnSpec(inputSpec)) return undefined;

    const labels = ctx.findLabels(inputSpec.axesSpec[0]);
    if (!labels) return undefined;

    // @TODO change to return labels after next SDK update.
    return Object.fromEntries(labels);
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
  )

  /**
   * P-frame with rawCounts
   */
  .output("pf", (wf) => {
    //return wf.outputs?.resolve("pf")?.resolve("rawCounts.data")?.listInputFields()
    const pCols = wf.outputs?.resolve("pf")?.getPColumns();
    if (pCols === undefined) return undefined;

    return wf.createPFrame(pCols);
  })

  /**
   * Returns true if the block is currently in "running" state
   */
  .output("isRunning", (ctx) => ctx.outputs?.getIsReadyOrError() === false)

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
    { type: "link", href: "/PCA", label: "Principal Component Analysis" },
  ])

  .title((ctx) =>
    ctx.args.title
      ? `STAR Read Mapping - ${ctx.args.title}`
      : "STAR Read Mapping"
  )

  .done();

export type BlockOutputs = InferOutputsType<typeof model>;
