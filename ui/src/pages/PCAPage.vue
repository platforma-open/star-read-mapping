<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
//import { ref } from 'vue';
import { GraphMaker, GraphMakerProps } from '@milaboratories/graph-maker';
import { useApp } from "../app";

const app = useApp();

const defaultOptions: GraphMakerProps['defaultOptions'] = [
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
</script>

<template>

  <GraphMaker chartType="scatterplot" v-model="app.model.ui.pcaGraphState" :p-frame="app.model.outputs.pcaPf"
    :defaultOptions="defaultOptions" />

</template>
