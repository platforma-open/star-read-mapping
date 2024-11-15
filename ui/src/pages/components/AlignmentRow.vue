<script lang="ts" setup>
import { computed } from 'vue';
import { StarQC } from '../results';
import AlignmentLegend from './AlignmentLegend.vue';
import StackedRow from './StackedRow.vue';

const props = defineProps<{
  alignReport?: StarQC;
  size?: 'large';
  showFractionInLabel?: boolean;
}>();


const parts = computed(() => {
  const report = props.alignReport;
  console.log(report);
  if (!report) return [];

  return [
    {
      color: 'red',
      value: report?.numberOfInputReads - report?.uniquelyMapped,
      label: "Not mapped",
    },
    {
      color: 'green',
      value: report?.uniquelyMapped,
      label: "Mapped",
    }

  ];
});
const legends = computed(() => parts.value.map(p => ({
  color: p.color,
  text: p.label
})));
</script>

<template>
  <StackedRow :size="size" :value="parts" />
  <AlignmentLegend v-if="size === 'large' && legends.length" :legends="legends" />
</template>
