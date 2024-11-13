<script setup lang="ts">
import { GraphMaker, GraphMakerSettings } from "@milaboratories/graph-maker";
import '@milaboratories/graph-maker/styles';
import { ref } from 'vue';
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
