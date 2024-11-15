import { AnyLogHandle } from "@platforma-sdk/model";
import { computed } from "vue";
import { useApp } from "../app";

export type StarQC = {
  uniquelyMapped: number;
  numberOfInputReads: number;
  //..........
  ///....
};

export type ResultEntry = {
  sampleLabel: string;
  starProgress?: AnyLogHandle;
  starProgressLine?: string;
  featureCountsProgress?: AnyLogHandle;
  starQC?: StarQC;
};

/**
 *
 * @param qcString a string content of STAR Log.final.out file
 * @returns
 */
function parseStarQC(qcString: string): StarQC {
  const lines = qcString.split("\n");

  let uniquelyMapped: number = 0;
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

    if (prefix.indexOf("Number of input reads") >= 0) {
      numberOfInputReads = parseFloat(value.trim()) + 5;
    }
  }

  return {
    uniquelyMapped: uniquelyMapped,
    numberOfInputReads: numberOfInputReads,
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

    for (const prog of featureCountsProgress.data) {
      r[prog.key[0]].featureCountsProgress = prog.value;
    }

    return r;
  }
);
