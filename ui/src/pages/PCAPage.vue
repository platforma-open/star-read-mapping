<script setup lang="ts">
import '@milaboratories/graph-maker/styles';
// import { ref } from 'vue';
import type { PredefinedGraphOption } from '@milaboratories/graph-maker';
import { GraphMaker } from '@milaboratories/graph-maker';
import type { PColumnIdAndSpec } from '@platforma-sdk/model';
import { computed } from 'vue';
import { useApp } from '../app';

const app = useApp();

const defaultOptions = computed((): PredefinedGraphOption<'scatterplot'>[] | undefined => {
  if (!app.model.outputs.pcaPcols)
    return undefined;

  const pcaPcols = app.model.outputs.pcaPcols;
  function getIndex(name: string, pcols: PColumnIdAndSpec[]): number {
    return pcols.findIndex((p) => (p.spec.name === name));
  }

  const defaults: PredefinedGraphOption<'scatterplot'>[] = [
    {
      inputName: 'x',
      selectedSource: pcaPcols[getIndex('pl7.app/rna-seq/pc1', pcaPcols)].spec,
    },
    {
      inputName: 'y',
      selectedSource: pcaPcols[getIndex('pl7.app/rna-seq/pc2', pcaPcols)].spec,
    },
  ];
  return defaults;
});

</script>

<template>
  <GraphMaker
    v-model="app.model.ui.pcaGraphState" chartType="scatterplot" :p-frame="app.model.outputs.pcaPf"
    :defaultOptions="defaultOptions"
  />
</template>
