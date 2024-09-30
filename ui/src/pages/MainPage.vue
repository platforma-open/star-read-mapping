<script setup lang="ts">
import { 
  PlDropdown, 
  PlTextArea,
  FileInput,
  PlProgressBar 
} from "@milaboratory/sdk-vue";
// import { PlBtnGroup } from '@milaboratory/platforma-uikit';
import { computed 
 // , reactive 
} from "vue";
import { useApp } from "../app";

const app = useApp();

const args = app.createArgsModel();

const dataOptions = computed(() =>
  app.outputValues.dataOptions?.map((v) => ({
    text: v.label,
    value: v.ref,
  }))
);

const inputOptions = [
  { text: "Single-end", value: "SingleEnd" },
  { text: "Paired-end", value: "PairedEnd" },
];

const inputOptionsStr = [
  { text: "Unstranded", value: "Unstranded" },
  { text: "Same strand", value: "SameStrand" },
  { text: "Opposite strand", value: "OppositeStrand" },
];

const speciesOptions = [
  { text: "Homo sapiens", value: "hsa" },
  { text: "Mus musculus", value: "mmu" },
];

// index file upload progress
const indexUploadProgress = computed(() => {
  const p = app.getOutputFieldOkOptional("indexUploadProgress");
  if (p?.done) {
    return 100;
  } else {
    return 100 * (p?.status?.progress ?? 0.0);
  }
});

// genome annotation file upload progress
const genomeAnnUploadProgress = computed(() => {
  const p = app.getOutputFieldOkOptional("genomeAnnUploadProgress");
  if (p?.done) {
    return 100;
  } else {
    return 100 * (p?.status?.progress ?? 0.0);
  }
});

</script>

<template>
  <div class="container">
    <PlDropdown
      v-if="dataOptions"
      :options="dataOptions"
      v-model="args.model.ref"
      label="Select dataset"
      clearable
    />

    <PlDropdown
      :options="speciesOptions"
      v-model="args.model.species"
      label="Select species"
    />
    <PlDropdown
      :options="inputOptions"
      v-model="args.model.libraryType"
      label="Select library type"
    />

    <PlDropdown
      :options="inputOptionsStr"
      v-model="args.model.strandness"
      label="Select strandness"
    />

    <FileInput
      v-model="args.model.indexFile"
      placeholder="Drag .tar index file"
      file-dialog-title="Select index tar file"
      clearable
    />
    <br />
    <PlProgressBar
      :progress="indexUploadProgress"
      loading
    />
    <FileInput
      v-model="args.model.genomeAnnFile"
      placeholder="Drag .gtf genome annoration file"
      file-dialog-title="Select genome annotation gtf file"
      clearable
    />
    <br />
    <PlProgressBar
      :progress="genomeAnnUploadProgress"
      loading
    />

    <PlTextArea
      v-if="app.outputs.starProgress && app.outputs.starProgress.ok"
      :model-value="app.outputs.starProgress.value"
      :readonly="true"
      :rows="3"
      label="STAR alignment"
    />

    <PlTextArea
      v-if="app.outputs.featureCountsProgress && app.outputs.featureCountsProgress.ok"
      :model-value="app.outputs.featureCountsProgress.value"
      :readonly="true"
      :rows="10"
      label="Feature counts"
    />
  </div>
</template>

<style lang="css">
.container {
  margin-top: 12px;
  display: flex;
  flex-direction: column;
  gap: 12px;
  width: 70%;
}
</style>
