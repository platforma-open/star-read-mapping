<script setup lang="ts">
import { PlBtnGroup, PlLogView } from '@platforma-sdk/ui-vue';
import { computed, ref } from 'vue';
import AlignmentRow from './components/AlignmentRow.vue';
import FeatureCountsRow from './components/FeatureCountsRow.vue';
import { resultMap } from './results';

const sampleId = defineModel<string | undefined>();

const starProgress = computed(() => {
  if (sampleId.value === undefined || resultMap.value === undefined)
    return undefined;
  else {
    return resultMap.value[sampleId.value].starProgress;
  }
});

const starQc = computed(() => {
  if (sampleId.value === undefined || resultMap.value === undefined)
    return undefined;
  else {
    return resultMap.value[sampleId.value].starQC;
  }
});

const featureCountsQc = computed(() => {
  if (sampleId.value === undefined || resultMap.value === undefined)
    return undefined;
  else {
    return resultMap.value[sampleId.value].featureCountsQC;
  }
});

const featureCountsProgress = computed(() => {
  if (sampleId.value === undefined || resultMap.value === undefined)
    return undefined;
  else {
    return resultMap.value[sampleId.value].featureCountsProgress;
  }
});

const options = [{
  label: 'Quality Controls',
  value: 'qc',
},
{
  label: 'Analysis Logs',
  value: 'logs',
}];

const currentView = ref('qc');

const logOptions = [{
  label: 'Mapping log',
  value: 'star',
},
{
  label: 'Feature counts log',
  value: 'subread',
}];

const currentLogView = ref('star');

</script>

<template>
  <PlBtnGroup v-model="currentView" :options="options" />
  <template v-if="currentView === 'qc'">
    <h3>Mapping</h3>
    <AlignmentRow :alignReport="starQc" size="large" showFractionInLabel />
    <h3>Feature assignment</h3>
    <FeatureCountsRow :alignReport="featureCountsQc" size="large" showFractionInLabel />
  </template>

  <template v-if="currentView === 'logs'">
    <PlBtnGroup v-model="currentLogView" :options="logOptions" />
    <PlLogView v-if="currentLogView === 'star'" :logHandle="starProgress" />
    <PlLogView v-if="currentLogView === 'subread'" :logHandle="featureCountsProgress" />
  </template>
</template>
