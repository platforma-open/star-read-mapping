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
      color: 'grey',
      value: report?.numberOfInputReads - 
              report?.uniquelyMapped -
              report?.mappedMultipleLoci -
              report?.mappedTooManyLoci -
              report?.unmappedTooShort - 
              report?.unmappedOther
              ,
      label: "Other",
    },
    {
      color: 'green',
      value: report?.uniquelyMapped,
      label: "Uniquely mapped",
    },
    {
      color: 'blue',
      value: report?.mappedMultipleLoci,
      label: "Mapped to multiple loci",
    },
    {
      color: 'orange',
      value: report?.mappedTooManyLoci,
      label: "Mapped to too many loci",
    },
    {
      color: 'red',
      value: report?.unmappedTooShort,
      label: "Unmapped: too short",
    },
    {
      color: 'purple',
      value: report?.unmappedOther,
      label: "Unmapped: other",
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
