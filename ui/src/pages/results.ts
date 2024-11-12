import { AnyLogHandle } from "@platforma-sdk/model";
import { computed } from "vue";
import { useApp } from "../app";

export type ResultEntry = {
  sampleLabel: string;
  starProgress?: AnyLogHandle;
  starProgressLine?: string;
  featureCountsProgress?: AnyLogHandle;
  starQc?: StarQC;
};

export type StarQC = {
  // Number of input reads
  nInputReads?: number;
  // Uniquely mapped reads number
  nMapped?: number;

  ///........////
};

function parseStarQCString(qcString: string): StarQC {
  const split = qcString.split("\n");

  const result: StarQC = {}

  for (const line in split) {
    
  }

}

// constructs result map: sampleId => ResultEntry
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
    for (const sampleId in labels) {
      r[sampleId] = {
        sampleLabel: labels[sampleId],
      };
    }
    for (const prog of starProgress.data) {
      const sampleId = prog.key[0];
      r[sampleId].starProgress = prog.value;
    }

    const starProgressLine = app.model.outputs.starProgressLine;
    if (starProgressLine !== undefined) {
      for (const prog of starProgressLine.data) {
        r[prog.key[0]].starProgressLine = prog.value;
      }
    }

    for (const prog of featureCountsProgress.data) {
      r[prog.key[0]].featureCountsProgress = prog.value;
    }

    return r;
  }
);
