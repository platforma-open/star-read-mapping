<script setup lang="ts">
import { PlTextArea } from '@milaboratory/platforma-uikit';
import { ref } from 'vue';
import { GraphMaker } from "@milaboratory/graph-maker";
import { GraphMakerSettings } from "@milaboratory/graph-maker/dist/GraphMaker/types";
import "@milaboratory/platforma-uikit/lib/dist/style.css";
import { model } from "@milaboratory/milaboratories.star-read-mapping.model";
import "@milaboratory/graph-maker/dist/style.css";
import { useApp } from "../app";

const placeholder = ref('Here should go the PCA plot');

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
