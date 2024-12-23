<script setup lang="ts">
import { ICellRendererParams } from 'ag-grid-enterprise';
import FeatureCountsRow from './components/FeatureCountsRow.vue';
import { FeatureCountsQC } from './results';

defineProps<{
    params: ICellRendererParams<{ featureCountsQc: FeatureCountsQC }>;
}>();
</script>

<template>
    <div v-if="!params.data?.featureCountsQc" :class="{ 'spinner-cell__white-bg': true }" class="spinner-cell">
        <div class="spinner-cell__indicator" :style="{ width: '100%' }"></div>
    </div>
    <FeatureCountsRow v-else style="height: 100%;" :align-report="params.data?.featureCountsQc" />
</template>

<style lang="css">
.spinner-cell {
    background-color: transparent;
    height: 100%;
    position: relative;
    width: 100%;
    overflow: hidden;
    border-radius: 2px;
    /* border-left: 1px solid var(--border-color-div-grey);
    border-right: 1px solid var(--border-color-div-grey); */
}

.spinner-cell__white-bg {
    background-color: #fff;
}

.spinner-cell__indicator {
    position: absolute;
    height: 100%;
    transition: width 0.4s ease-in-out;
    /* background: linear-gradient(90deg, #FFF 0%, #D8FAC8 100%); */
    --c: no-repeat linear-gradient(90deg, #FFF 0%, #D8FAC8 100%);
    background: var(--c), var(--c), #FFF;
    background-size: 60% 100%;

    animation: l16 5s infinite;
}

@keyframes l16 {
    0% {
        background-position: -150% 0, -150% 0
    }

    66% {
        background-position: 250% 0, -150% 0
    }

    100% {
        background-position: 250% 0, 250% 0
    }
}
</style>