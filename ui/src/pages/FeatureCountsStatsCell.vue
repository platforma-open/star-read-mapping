<script setup lang="ts">
import { ICellRendererParams } from '@ag-grid-community/core';
import FeatureCountsRow from './components/FeatureCountsRow.vue';
import { FeatureCountsQC } from './results';

defineProps<{
    params: ICellRendererParams<{featureCountsQc: FeatureCountsQC}>;
}>();
</script>

<template>
    <FeatureCountsRow style="height: 100%;" :align-report="params.data?.featureCountsQc" />
</template>
