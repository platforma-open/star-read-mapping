<script lang="ts" setup>
import { computed } from 'vue';
import type { FeatureCountsQC } from '../results';
import AlignmentLegend from './AlignmentLegend.vue';
import StackedRow from './StackedRow.vue';

const props = defineProps<{
  alignReport?: FeatureCountsQC;
  size?: 'large';
  showFractionInLabel?: boolean;
}>();

const parts = computed(() => {
  const report = props.alignReport;
  console.log(report);
  if (!report) return [];

  return [
    {
      color: '#6BD67D',
      value: report?.assigned,
      label: 'Assigned',
    },
    {
      color: '#B8397A',
      value: report?.unassignedMappingQuality,
      label: 'Unassigned: Mapping Quality',
    },
    {
      color: '#E75B64',
      value: report?.unassignedNoFeatures,
      label: 'Unassigned: No Features',
    },
    {
      color: '#FB9361',
      value: report?.unassignedAmbiguity,
      label: 'Unassigned: Ambiguity',
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
