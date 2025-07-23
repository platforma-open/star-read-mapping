<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
// import { ref } from 'vue';
import type { PredefinedGraphOption } from '@milaboratories/graph-maker';
import { GraphMaker } from '@milaboratories/graph-maker';
import { useApp } from '../app';

const app = useApp();

const defaultOptions: PredefinedGraphOption<'scatterplot'>[] | undefined = [
  {
    inputName: 'x',
    selectedSource: {
      kind: 'PColumn',
      name: 'pl7.app/rna-seq/pc1',
      valueType: 'Double',
      axesSpec: [
        {
          name: 'pl7.app/sampleId',
          type: 'String',
        },
      ],
    },
  },
  {
    inputName: 'y',
    selectedSource: {
      kind: 'PColumn',
      name: 'pl7.app/rna-seq/pc2',
      valueType: 'Double',
      axesSpec: [
        {
          name: 'pl7.app/sampleId',
          type: 'String',
        },
      ],
    },
  },
];
</script>

<template>
  <GraphMaker
    v-model="app.model.ui.pcaGraphState" chartType="scatterplot" :p-frame="app.model.outputs.pcaPf"
    :defaultOptions="defaultOptions"
  />
</template>
