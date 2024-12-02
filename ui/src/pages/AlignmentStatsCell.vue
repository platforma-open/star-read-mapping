<script setup lang="ts">
import { ICellRendererParams } from '@ag-grid-community/core';
import AlignmentRow from './components/AlignmentRow.vue';
import { StarQC } from './results';

defineProps<{
    params: ICellRendererParams<{ starQc: StarQC }>;
}>();
</script>

<template>
    <AlignmentRow style="height: 100%;" :align-report="params.data?.starQc" />
</template>
