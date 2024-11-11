import { AnyLogHandle } from "@platforma-sdk/model";
import { computed } from "vue";
import { useApp } from "../app";

export type ResultEntry = {
  sampleLabel: string;
  starProgress?: AnyLogHandle;
  starProgressLine?: string;
  featureCountsProgress?: AnyLogHandle;
};

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

    for (const prog of featureCountsProgress.data) {
      r[prog.key[0]].featureCountsProgress = prog.value;
    }

    return r;
  }
);
