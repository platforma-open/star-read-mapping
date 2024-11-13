<script setup lang="ts">
import { GraphMaker, GraphMakerSettings } from "@milaboratories/graph-maker";
import '@milaboratories/graph-maker/styles';
//import { ref } from 'vue';
import { useApp } from "../app";

const app = useApp();

const settings = {
  chartType: 'scatterplot',
  template: 'dots',
  title: '',
  defaultOptions: [
    {
    inputName: 'x',
    selectedSource: {
      kind: 'PColumn',
      name: 'Principal Component 1',
      valueType: 'Int'
    }
  },
  {
    inputName: 'y',
    selectedSource: {
      kind: 'PColumn',
      name: 'Principal Component 2',
      valueType: 'Int'
    }
  }
]
   } as GraphMakerSettings;
</script>

<template>
  <div>
    <graph-maker
      v-model="settings"
      :pFrame="app.model.outputs.pcaPf"
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
