<script setup lang="ts">
import { GraphMaker, GraphMakerSettings } from "@milaboratories/graph-maker";
import '@milaboratories/graph-maker/styles';
import { ref } from 'vue';
import { useApp } from "../app";

const app = useApp();

const settings = ref({
  chartType: "discrete",
  template: "bar",
  title: "Tittle",
  defaultOptions: undefined,
  optionsState: undefined,
  statisticsSettings: undefined,
  axesSettings: undefined,
  layersSettings: undefined,
  dataBindAes: undefined,
} satisfies GraphMakerSettings);
</script>

<template>
  <GraphMaker :p-frame="app.model.outputs.pf" v-model="settings" graph-title="Title" />
</template>


<style lang="css">
.container {
  margin-top: 12px;
  display: flex;
  flex-direction: column;
  gap: 12px;
}
</style>
