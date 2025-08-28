<script setup lang="ts">
import type { ICellRendererParams } from 'ag-grid-enterprise';
import FeatureCountsRow from './components/FeatureCountsRow.vue';
import type { FeatureCountsQC } from './results';
import { PlAgCellProgress, type PlProgressCellProps } from '@platforma-sdk/ui-vue';
import { computed } from 'vue';

const props = defineProps<{
  params: ICellRendererParams<{ featureCountsQc: FeatureCountsQC }>;
}>();

const progressProps = computed<PlProgressCellProps>(() => ({
  stage: 'running',
  step: '',
  progressString: '',
}));
</script>

<template>
  <PlAgCellProgress
    v-if="!props.params.data?.featureCountsQc"
    v-bind="{ params: { ...props.params as any, ...progressProps } }"
  />
  <FeatureCountsRow v-else style="height: 100%" :align-report="props.params.data?.featureCountsQc" />
</template>
