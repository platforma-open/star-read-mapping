<script setup lang="ts">
import { PlTextArea } from '@milaboratories/uikit';
import { ref } from 'vue';

const placeholder = ref('Here should go the html report');
</script>

<template>
  <PlTextArea
    v-model="placeholder"
    :rows="30"
    label="MultiQC"
  />
</template>


<style lang="css">
.container {
  margin-top: 12px;
  display: flex;
  flex-direction: column;
  gap: 12px;
}
</style>