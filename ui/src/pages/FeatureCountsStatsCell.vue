<script setup lang="ts">
import { ICellRendererParams } from 'ag-grid-enterprise';
import FeatureCountsRow from './components/FeatureCountsRow.vue';
import { FeatureCountsQC } from './results';
import { PlAgCellProgress } from '@platforma-sdk/ui-vue';

defineProps<{
  params: ICellRendererParams<{ featureCountsQc: FeatureCountsQC }>;
}>();
</script>

<template>
  <PlAgCellProgress
    v-if="!params.data?.featureCountsQc"
    v-bind="{ params: { ...params, stage: 'running', step: '', progressString: '' } }"
  />
  <FeatureCountsRow v-else style="height: 100%" :align-report="params.data?.featureCountsQc" />
</template>
