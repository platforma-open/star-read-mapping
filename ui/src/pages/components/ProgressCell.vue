<script setup lang="ts">
import { ICellRendererParams } from 'ag-grid-enterprise';
import { computed, unref } from 'vue';
import { PlAgCellProgress, PlProgressCellProps } from '@platforma-sdk/ui-vue';

const props = defineProps<{
  params: ICellRendererParams;
}>();

const progressString = computed(() => {
  return props.params.value ?? 'Unknown';
});

type Parsed = {
  raw?: string;
  stage?: string;
  time?: string;
  percentage?: string;
  percentageLabel?: string;
};

const parsed = computed<Parsed>(() => {
  const raw = unref(progressString);

  const res: Parsed = {
    raw
  };

  if (!raw) {
    return res;
  }

  console.log(raw);
  if (raw.indexOf('.....') < 0) {
    return res;
  }

  const parts = raw.split('.....');

  res.time = parts[0].trim();
  res.stage = parts[1].trim();

  switch (res.stage) {
    case 'loading genome':
      res.percentage = '0';
      break;
    case 'started mapping':
      res.percentage = '33';
      break;
    case 'started sorting BAM':
      res.percentage = '66';
      break;
    case 'finished successfully':
      res.percentage = '100';
      break;
  }

  if (res.percentage) {
    res.percentageLabel = res.percentage + '%';
  }

  return res;
});

const ProgressProps = computed<PlProgressCellProps>(() => {
  return {
    stage: parsed.value.stage === 'Queued' ? 'not_started' : 'running',
    step: parsed.value.stage || '',
    progress: parsed.value.percentage ? +parsed.value.percentage : 0,
    progressString: parsed.value.percentageLabel || ''
  };
});
</script>

<template>
  <PlAgCellProgress v-bind="{ params: { ...props.params, ...ProgressProps } }" />
</template>
