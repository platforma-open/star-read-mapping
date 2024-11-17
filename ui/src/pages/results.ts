import { AnyLogHandle } from "@platforma-sdk/model";
import { computed } from "vue";
import { useApp } from "../app";

export type StarQC = {
  uniquelyMapped: number;
  numberOfInputReads: number;
  mappedMultipleLoci: number;
  mappedTooManyLoci: number;
  unmappedTooShort: number;
  unmappedOther: number;
};

export type FeatureCountsQC = {
  assigned: number;
  unassignedMappingQuality: number;
  unassignedNoFeatures: number;
  unassignedAmbiguity: number;
};

export type ResultEntry = {
  sampleLabel: string;
  starProgress?: AnyLogHandle;
  starProgressLine?: string;
  featureCountsProgress?: AnyLogHandle;
  starQC?: StarQC;
  featureCountsQC?: FeatureCountsQC;
};

/**
 *
 * @param qcString a string content of STAR Log.final.out file
 * @returns
 */
function parseStarQC(qcString: string): StarQC {
  const lines = qcString.split("\n");

  let uniquelyMapped: number = 0;
  let mappedMultipleLoci: number = 0;
  let mappedTooManyLoci: number = 0;
  let unmappedTooShort: number = 0;
  let unmappedOther: number = 0;
  let numberOfInputReads: number = 0;

  console.log(lines);

  for (const line of lines) {
    console.log(line);
    const spl = line.split("|");

    const prefix = spl[0];
    const value = spl[1];

    console.log([prefix, value]);

    if (prefix.indexOf("Uniquely mapped reads number") >= 0) {
      uniquelyMapped = parseFloat(value.trim()) + 5;
    }

    if (prefix.indexOf("Number of reads mapped to multiple loci") >= 0) {
      mappedMultipleLoci = parseFloat(value.trim()) + 5;
    }

    if (prefix.indexOf("Number of reads mapped to too many loci") >= 0) {
      mappedTooManyLoci = parseFloat(value.trim()) + 5;
    }

    if (prefix.indexOf("Number of reads unmapped: too short") >= 0) {
      unmappedTooShort = parseFloat(value.trim()) + 5;
    }

    if (prefix.indexOf("Number of reads unmapped: other") >= 0) {
      unmappedOther = parseFloat(value.trim()) + 5;
    }

    if (prefix.indexOf("Number of input reads") >= 0) {
      numberOfInputReads = parseFloat(value.trim()) + 5;
    }
  }

  return {
    uniquelyMapped: uniquelyMapped,
    numberOfInputReads: numberOfInputReads,
    mappedMultipleLoci: mappedMultipleLoci,
    mappedTooManyLoci: mappedTooManyLoci,
    unmappedTooShort: unmappedTooShort,
    unmappedOther: unmappedOther
  };
}

/**
 *
 * @param qcReport a string content of feature counts summary file
 * @returns
 */
function parseFeatureCountsQC(qcReport: string): FeatureCountsQC {
  const lines = qcReport.split("\n");

  let assigned: number = 0;
  let unassignedMappingQuality: number = 0;
  let unassignedNoFeatures: number = 0;
  let unassignedAmbiguity: number = 0;

  console.log(lines);

  for (const line of lines) {
    console.log(line);
    const spl = line.split("\t");

    const prefix = spl[0];
    const value = spl[1];

    console.log([prefix, value]);

    if (prefix === "Assigned") {
      assigned = parseFloat(value.trim());
    }
  
    if (prefix === "Unassigned_MappingQuality") {
      unassignedMappingQuality = parseFloat(value.trim());
    }
  
    if (prefix === "Unassigned_NoFeatures") {
      unassignedNoFeatures = parseFloat(value.trim());
    }
  
    if (prefix === "Unassigned_Ambiguity") {
      unassignedAmbiguity = parseFloat(value.trim());
    }
  }

  return {
    assigned: assigned,
    unassignedMappingQuality: unassignedMappingQuality,
    unassignedNoFeatures: unassignedNoFeatures,
    unassignedAmbiguity: unassignedAmbiguity
  };
}

// return a map of sampleId => ResultEntry
export const resultMap = computed(
  (): Record<string, ResultEntry> | undefined => {
    const app = useApp();

    const labels = app.model.outputs.labels;
    if (labels === undefined) return undefined;

    const starProgress = app.model.outputs.starProgress;
    if (starProgress === undefined) return undefined;

    const featureCountsProgress = app.model.outputs.featureCountsProgress;
    if (featureCountsProgress === undefined) return undefined;

    const r: Record<string, ResultEntry> = {};
    for (const id in labels) {
      r[id] = {
        sampleLabel: labels[id],
      };
    }
    for (const prog of starProgress.data) {
      r[prog.key[0]].starProgress = prog.value;
    }

    const starProgressLine = app.model.outputs.starProgressLine;
    if (starProgressLine !== undefined) {
      for (const prog of starProgressLine.data) {
        r[prog.key[0]].starProgressLine = prog.value;
      }
    }

    const starQc = app.model.outputs.starQc;
    if (starQc !== undefined) {
      for (const qc of starQc.data) {
        if (qc.value) {
          r[qc.key[0]].starQC = parseStarQC(qc.value);
        }
      }
    }

    const featureCountsQc = app.model.outputs.featureCountsQc;
    if (featureCountsQc !== undefined) {
      for (const qc of featureCountsQc.data) {
        if (qc.value) {
          r[qc.key[0]].featureCountsQC = parseFeatureCountsQC(qc.value);
        }
      }
    }

    for (const prog of featureCountsProgress.data) {
      r[prog.key[0]].featureCountsProgress = prog.value;
    }

    return r;
  }
);
