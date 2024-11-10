<script setup lang="ts">
import {
  PlBlockPage,
  PlDropdown,
  PlDropdownRef,
  PlFileInput,
  PlLogView
} from "@platforma-sdk/ui-vue";
import { useApp } from "../app";

const app = useApp();

const inputOptions = [
  { text: "Single-end", value: "SingleEnd" },
  { text: "Paired-end", value: "PairedEnd" },
];

const inputOptionsStr = [
  { text: "Unstranded", value: "0" },
  { text: "Stranded", value: "1" },
  { text: "Reverse stranded", value: "2" },
];

const speciesOptions = [
  { text: "Homo sapiens", value: "hsa" },
  { text: "Mus musculus", value: "mmu" },
];
</script>

<template>
  <PlBlockPage>
    <template #title>STAR Read Mapping</template>

    <PlDropdownRef :options="app.model.outputs.dataOptions ?? []" v-model="app.model.args.ref" label="Select dataset"
      clearable />

    <PlDropdown :options="speciesOptions" v-model="app.model.args.species" label="Select species" />
    <PlDropdown :options="inputOptions" v-model="app.model.args.libraryType" label="Select library type" />

    <PlDropdown :options="inputOptionsStr" v-model="app.model.args.strandness" label="Select strandness" />

    <PlFileInput v-model="app.model.args.indexFile" :progress="app.model.outputs.indexUploadProgress"
      placeholder="Select .tar index file" fileDialogTitle="Select index tar file" clearable />

    <PlFileInput v-model="app.model.args.genomeAnnFile" :progress="app.model.outputs.genomeAnnUploadProgress"
      placeholder="Drag .gtf genome annoration file" file-dialog-title="Select genome annotation gtf file" clearable />

    <PlLogView :model-value="app.model.outputs.starProgress" label="STAR alignment log" />

    <PlLogView :model-value="app.model.outputs.featureCountsProgress" label="Feature counts log" />

  </PlBlockPage>
</template>
