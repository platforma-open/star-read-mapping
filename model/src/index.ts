import type {
  PColumnIdAndSpec,
} from '@platforma-sdk/model';
import {
  BlockModel,
  createPFrameForGraphs,
  type InferOutputsType,
  isPColumnSpec,
  parseResourceMap,
  type PlRef,
} from '@platforma-sdk/model';

import { type GraphMakerState } from '@milaboratories/graph-maker';
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
  species?: string;

  /**
   * Strandness settings
   */
  strandness?: string;

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
  sDistGraphState: GraphMakerState;
};

export const model = BlockModel.create()

  .withArgs<BlockArgs>({
    species: 'homo-sapiens',
    strandness: '0',
  })

  .withUiState<UiState>({
    pcaGraphState: {
      template: 'dots',
      title: 'Principal Components Analysis',
    },
    sDistGraphState: {
      template: 'heatmapClustered',
      title: 'Sample Distances',
    },
  })

  /**
   * Find possible options for the fastq input
   */
  .output('dataOptions', (ctx) => {
    return ctx.resultPool.getOptions((v) => {
      if (!isPColumnSpec(v)) return false;
      const domain = v.domain;
      return (
        v.name === 'pl7.app/sequencing/data'
        && (v.valueType as string) === 'File'
        && domain !== undefined
        && (domain['pl7.app/fileExtension'] === 'fastq'
          || domain['pl7.app/fileExtension'] === 'fastq.gz')
      );
    });
  })

  .output('labels', (ctx) => {
    const inputRef = ctx.args.ref;
    if (inputRef === undefined) return undefined;

    const inputSpec = ctx.resultPool.getSpecByRef(inputRef); // @TODO use resultPool.getPColumnSpecByRef after updating SDK
    if (inputSpec === undefined || !isPColumnSpec(inputSpec)) return undefined;

    const labels = ctx.findLabels(inputSpec.axesSpec[0]);
    if (!labels) return undefined;

    return labels;
  })

  /**
   * Preprocessing progress
   */
  .output('starProgress', (wf) => {
    return parseResourceMap(
      wf.outputs?.resolve('starProgress'),
      (acc) => acc.getLogHandle(),
      false,
    );
  })

  /**
   * Last line from StAR output
   */
  .output('starProgressLine', (wf) => {
    return parseResourceMap(
      wf.outputs?.resolve('starProgress'),
      (acc) => acc.getLastLogs(1),
      false,
    );
  })

  .output('starQc', (wf) =>
    parseResourceMap(
      wf.outputs?.resolve('starQc'),
      (acc) => acc.getFileContentAsString(),
      false,
    ),
  )

  .output('featureCountsProgress', (wf) => {
    return parseResourceMap(
      wf.outputs?.resolve('featureCountsProgress'),
      (acc) => acc.getLogHandle(),
      false,
    );
  })

  .output('featureCountsQc', (wf) =>
    parseResourceMap(
      wf.outputs?.resolve('featureCountsQc'),
      (acc) => acc.getFileContentAsString(),
      false,
    ),
  )

  /**
   * P-frame with rawCounts
   */
  .output('rawCountsPf', (wf) => {
    // return wf.outputs?.resolve("pf")?.resolve("rawCounts.data")?.listInputFields()
    const pCols = wf.outputs?.resolve('rawCountsPf')?.getPColumns();
    if (pCols === undefined) return undefined;

    return wf.createPFrame(pCols);
  })

/**
   * P-frame with normCounts
   */
  .output('normCountsPf', (wf) => {
    const pCols = wf.outputs?.resolve('normCountsPf')?.getPColumns();
    if (pCols === undefined) return undefined;

    return wf.createPFrame(pCols);
  })

  /**
   * Returns true if the block is currently in "running" state
   */
  .output('isRunning', (ctx) => ctx.outputs?.getIsReadyOrError() === false)

  .output('pcaPf', (wf) => {
    // return wf.outputs?.resolve("pf")?.resolve("rawCounts.data")?.listInputFields()
    const pCols = wf.outputs?.resolve('pcaComponents')?.getPColumns();
    if (pCols === undefined) return undefined;

    // Check for metadata in upstream data
    // const upstream = wf.resultPool
    //   .getData()
    //   .entries.map((v) => v.obj)
    //   .filter(isPColumn)
    //   .filter((column) =>
    //     column.spec.name === 'pl7.app/metadata' || column.spec.name === 'pl7.app/label',
    //   );

    return createPFrameForGraphs(wf, pCols);
  })

  // List of pcolumns used for PCA default options
  .output('pcaPcols', (wf) => {
    const pCols = wf.outputs?.resolve('pcaComponents')?.getPColumns();
    if (pCols === undefined) return undefined;

    return pCols.map(
      (c) =>
        ({
          columnId: c.id,
          spec: c.spec,
        } satisfies PColumnIdAndSpec),
    );
  })

  .output('sampleDistancesSpec', (wf) => {
    const pCols = wf.outputs?.resolve('sampleDistances')?.getPColumns();
    if (pCols === undefined) return undefined;
    return pCols[0].spec;
  })

  .output('sampleDistancesPf', (wf) => {
    const pCols = wf.outputs?.resolve('sampleDistances')?.getPColumns();
    if (pCols === undefined) return undefined;

    // Check for metadata in upstream data
    // const upstream = wf.resultPool
    //   .getData()
    //   .entries.map((v) => v.obj)
    //   .filter(isPColumn)
    //   .filter((column) =>
    //     column.spec.name === 'pl7.app/metadata' || column.spec.name === 'pl7.app/label',
    //   );

    // @TODO find why createPFrameForGraphs is not detecting the metadata without the previous check
    // return createPFrameForGraphs(wf, [...pCols, ...upstream]);
    return createPFrameForGraphs(wf, pCols);
  })

  .sections([
    { type: 'link', href: '/', label: 'Settings' },
    { type: 'link', href: '/PCA', label: 'Principal Component Analysis' },
    { type: 'link', href: '/SDist', label: 'Sample Distances' },
  ])

  .title((ctx) =>
    ctx.args.title
      ? `STAR Read Mapping - ${ctx.args.title}`
      : 'STAR Read Mapping',
  )

  .done(2);

export type BlockOutputs = InferOutputsType<typeof model>;
