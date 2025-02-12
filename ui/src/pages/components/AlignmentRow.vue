<script lang="ts" setup>
import { computed } from 'vue';
import type { StarQC } from '../results';
import AlignmentLegend from './AlignmentLegend.vue';
import StackedRow from './StackedRow.vue';

const props = defineProps<{
  alignReport?: StarQC;
  size?: 'large';
  showFractionInLabel?: boolean;
}>();

const parts = computed(() => {
  const report = props.alignReport;
  // console.log(report);
  if (!report) return [];

  return [
    {
      color: '#6BD67D',
      value: report?.uniquelyMapped,
      label: 'Uniquely mapped',
    },
    {
      color: '#B8397A',
      value: report?.mappedMultipleLoci,
      label: 'Mapped to multiple loci',
    },
    {
      color: '#FEBF51',
      value: report?.mappedTooManyLoci,
      label: 'Mapped to too many loci',
    },
    {
      color: '#E75B64',
      value: report?.unmappedTooShort,
      label: 'Unmapped: too short',
    },
    {
      color: '#FB9361',
      value: report?.unmappedOther,
      label: 'Unmapped: other',
    },
    {
      color: 'grey',
      value: report?.numberOfInputReads
        - report?.uniquelyMapped
        - report?.mappedMultipleLoci
        - report?.mappedTooManyLoci
        - report?.unmappedTooShort
        - report?.unmappedOther,
      label: 'Other',
    },
  ];
});
const legends = computed(() => parts.value.map((p) => ({
  color: p.color,
  text: p.label,
})));
</script>

<template>
  <StackedRow :size="size" :value="parts" />
  <AlignmentLegend v-if="size === 'large' && legends.length" :legends="legends" />
</template>
