<script setup lang="ts">
import { ICellRendererParams } from 'ag-grid-enterprise';
import { computed, unref } from 'vue';

const props = defineProps<{
  params: ICellRendererParams
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
  if (raw.indexOf(".....") < 0) {
    return res;
  }

  const parts = raw.split(".....");

  res.time = parts[0].trim();
  res.stage = parts[1].trim();

  switch (res.stage) {
    case 'loading genome':
      res.percentage = '0';
      break;
    case 'started mapping':
      res.percentage = '33'
      break;
    case 'started sorting BAM':
      res.percentage = '66'
      break;
    case 'finished successfully':
      res.percentage = '100'
      break;
  }

  if (res.percentage) {
    res.percentageLabel = res.percentage + '%';
  }

  return res;
});

const canShowBackground = computed(
  () => parsed.value.stage !== 'Queued'
);
</script>

<template>
  <div :class="{ 'progress-cell__white-bg': canShowBackground }" class="progress-cell">
    <div class="progress-cell__indicator" :style="{ width: parsed.percentage + '%' }"></div>
    <div class="progress-cell__body">
      <div :class="{ 'progress-cell__stage--queued': parsed.stage === 'Queued' }" class="progress-cell__stage">
        {{ parsed.stage }}
      </div>
      <div class="progress-cell__percentage">{{ parsed.percentageLabel }}</div>
    </div>
  </div>
</template>

<style lang="css">
.progress-cell {
  background-color: transparent;
  height: 100%;
  position: relative;
  width: 100%;
  overflow: hidden;
  border-radius: 2px;
  /* border-left: 1px solid var(--border-color-div-grey);
    border-right: 1px solid var(--border-color-div-grey); */
}

.progress-cell__white-bg {
  background-color: #fff;
}

.progress-cell__indicator {
  position: absolute;
  height: 100%;
  transition: width 0.4s ease-in-out;
  /* background: linear-gradient(90deg, #FFF 0%, #D8FAC8 100%); */
  background: linear-gradient(90deg, #fff 0%, #d8fac8 100%);
}

.progress-cell__body {
  padding: 0 15px;
  display: flex;
  gap: 12px;
  align-items: center;
  height: 100%;
  width: 100%;
  position: absolute;
  z-index: 1;
}

.progress-cell__stage {
  overflow: hidden;
  text-overflow: ellipsis;
  flex-shrink: 1;
}

.progress-cell__percentage {
  flex-grow: 1;
  flex-shrink: 0;
  text-align: right;
}

.progress-cell__stage--queued {
  color: var(--txt-03);
}
</style>
