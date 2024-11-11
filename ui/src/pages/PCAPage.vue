<script setup lang="ts">
// import { ref } from 'vue';
import { GraphMakerSettings, GraphMaker } from '@milaboratories/graph-maker';
//import { model } from "@platforma-open/milaboratories.star-read-mapping.model";
import "@milaboratories/graph-maker/dist/style.css";
import { useApp } from "../app";

const app = useApp();

const settings = {
  chartType: 'scatterplot',
  template: 'dots',
  title: '',
  defaultOptions: [{
    inputName: 'y',
    selectedSource: {
      kind: 'PColumn',
      name: 'countMatrix',
      valueType: 'Int'
    }
  }]
   } as GraphMakerSettings;
</script>

<template>
  <div class="container_graph_page" :key="app.queryParams.id">
    <graph-maker
      v-model="settings"
      :pFrame="app.model.outputs.pf"
    />
  </div>
</template>


<style lang="css">
.container {
  margin-top: 12px;
  display: flex;
  flex-direction: column;
  gap: 12px;
}
</style>
