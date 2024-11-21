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
        name: "pl7.app/rna-seq/pc1",
        valueType: "Double",
        axesSpec: [
          {
            name: "pl7.app/sampleId",
            type: "String"
          }
        ]
      }
    },
    {
      inputName: 'y',
      selectedSource: {
        kind: 'PColumn',
        name: "pl7.app/rna-seq/pc2",
        valueType: "Double",
        axesSpec: [
          {
            name: "pl7.app/sampleId",
            type: "String"
          }
        ]
      }
    }
  ]
} as GraphMakerSettings;
</script>

<template>
  <div>
    <graph-maker v-model="settings" :pFrame="app.model.outputs.pcaPf" />
  </div>
</template>
