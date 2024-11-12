<script setup lang="ts">
import { PlBtnGroup, PlLogView } from '@platforma-sdk/ui-vue';
import { computed, ref } from 'vue';
import { resultMap } from './results';

const sampleId = defineModel<string | undefined>()

const starProgress = computed(() => {
    
    if (sampleId.value === undefined || resultMap.value === undefined)
        return undefined;
    else {
        console.log(resultMap.value[sampleId.value].starProgress);
        return resultMap.value[sampleId.value].starProgress;
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
    label: 'Mapping log',
    value: 'star'
},
{
    label: 'Feature counts log',
    value: 'subread'
}]

const currentView = ref('star')

</script>

<template>
    <PlBtnGroup v-model="currentView" :options="options" />
    <PlLogView v-if="currentView === 'star'" :logHandle="starProgress" />
    <PlLogView v-if="currentView === 'subread'" :logHandle="featureCountsProgress" />
</template>