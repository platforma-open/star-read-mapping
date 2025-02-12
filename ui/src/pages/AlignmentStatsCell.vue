<script setup lang="ts">
import type { ICellRendererParams } from 'ag-grid-enterprise';
import AlignmentRow from './components/AlignmentRow.vue';
import type { StarQC } from './results';

defineProps<{
  params: ICellRendererParams<{ starQc: StarQC }>;
}>();
</script>

<template>
  <AlignmentRow style="height: 100%;" :align-report="params.data?.starQc" />
</template>
