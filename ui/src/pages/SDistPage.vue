<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
// import { ref } from 'vue';
import type { GraphMakerProps } from '@milaboratories/graph-maker';
import { GraphMaker } from '@milaboratories/graph-maker';
import { useApp } from '../app';
import { computed } from 'vue';

const app = useApp();

// This ensures that we don't need to re-open the block to see the heatmap if
// the ui instance was created after opening it
// if (app.model.ui.sDistGraphState === undefined) {
//   app.model.ui.sDistGraphState = {template: "heatmap", title: "Sample Distances"}
// }

const defaultOptions = computed((): GraphMakerProps['defaultOptions'] => {
  const distanceSpec = app.model.outputs.sampleDistancesSpec;

  if (!distanceSpec) {
    return undefined;
  }

  const defaults: GraphMakerProps['defaultOptions'] = [
    {
      inputName: 'value',
      selectedSource: distanceSpec,
    },
    {
      inputName: 'x',
      selectedSource: distanceSpec.axesSpec[0],
    },
    {
      inputName: 'y',
      selectedSource: distanceSpec.axesSpec[1],
    },
  ];

  return defaults;
});

</script>

<template>
  <GraphMaker
    v-model="app.model.ui.sDistGraphState" chartType="heatmap"
    :p-frame="app.model.outputs.sampleDistancesPf"
    :defaultOptions="defaultOptions"
  />
</template>
