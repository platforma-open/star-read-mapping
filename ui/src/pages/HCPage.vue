<script setup lang="ts">
import { PlTextArea } from '@milaboratories/uikit';
import { ref } from 'vue';
import { GraphMaker } from "@milaboratories/graph-maker";
import { GraphMakerSettings } from "@milaboratories/graph-maker/dist/GraphMaker/types";
import { model } from "@platforma-open/milaboratories.star-read-mapping.model";
import "@milaboratories/graph-maker/dist/style.css";
import { useApp } from "../app";

const placeholder = ref('Here should go the Hierarchical Clustering plot');

const app = useApp();

const settings = ref({
  chartType: "discrete",
  template: null,
  optionsState: null,
  statisticsSettings: null,
  axesSettings: null,
  layersSettings: null,
  dataBindAes: null,
} satisfies GraphMakerSettings);
</script>

<template>
  <PlTextArea
    v-model="placeholder"
    :rows="3"
    label="MultiQC"
  />

  <GraphMaker
    v-if="app.outputs.pf?.ok && app.outputs.pf.value"
    :p-frame-handle="app.outputs.pf.value"
    :settings="settings"
    :p-frame-driver="model.pFrameDriver"
    graph-title="Title"
  />
  <div v-else>"Nothing"</div>
</template>


<style lang="css">
.container {
  margin-top: 12px;
  display: flex;
  flex-direction: column;
  gap: 12px;
}
</style>
