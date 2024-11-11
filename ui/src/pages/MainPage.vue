<script setup lang="ts">
import { AgGridVue } from '@ag-grid-community/vue3';
import {
  AgGridTheme,
  PlAgOverlayLoading,
  PlAgOverlayNoRows,
  PlBlockPage,
  PlBtnPrimary,
  PlDropdown,
  PlDropdownRef,
  PlSlideModal
} from "@platforma-sdk/ui-vue";

import { ClientSideRowModelModule } from '@ag-grid-community/client-side-row-model';
import { GridApi, GridOptions, GridReadyEvent, ModuleRegistry } from '@ag-grid-community/core';
import { computed, reactive, shallowRef } from "vue";
import { useApp } from "../app";
import ReportPanel from './ReportPanel.vue';
import { resultMap } from './results';

const app = useApp();

const data = reactive<{
  settingsOpen: boolean,
  sampleReportOpen: boolean,
  selectedSample: string | undefined
}>({
  settingsOpen: app.args.ref === undefined,
  sampleReportOpen: false,
  selectedSample: undefined,
})

const inputOptions = [
  { text: "Single-end", value: "SingleEnd" },
  { text: "Paired-end", value: "PairedEnd" },
];

const inputOptionsStr = [
  { text: "Unstranded", value: "0" },
  { text: "Stranded", value: "1" },
  { text: "Reverse stranded", value: "2" },
];

const speciesOptions = [
  { text: "Homo sapiens", value: "hsa" },
  { text: "Mus musculus", value: "mmu" },
];


/** Rows for ag-table */
const results = computed<any[] | undefined>(() => {

  if (resultMap.value === undefined) return undefined;
  const rows = []
  for (const id in resultMap.value) {
    rows.push({
      "sampleId": id,
      "sampleLabel": resultMap.value[id].sampleLabel,
      "star": resultMap.value[id].starProgressLine, // @TODO status?
      "subread": "Running", // @TODO status?
    });
  }

  return rows;
});

ModuleRegistry.registerModules([ClientSideRowModelModule]);

const gridApi = shallowRef<GridApi<any>>();
const onGridReady = (params: GridReadyEvent) => {
  gridApi.value = params.api;
};

const columnDefs = [
  {
    colId: 'label',
    field: 'sampleLabel',
    headerName: "Sample"
  },
  {
    colId: 'star',
    field: 'star',
    // cellRenderer: 'ProgressCell',
    headerName: "Star progress"
  },
  {
    colId: 'subread',
    field: 'subread',
    // cellRenderer: 'AlignmentStatsCell',
    headerName: "Feature counts",
  }
];


const gridOptions: GridOptions = {
  getRowId: (row) => row.data.sampleId,
  onRowDoubleClicked: (e) => {
    data.selectedSample = e.data?.sampleId
    data.sampleReportOpen = data.selectedSample !== undefined;
  },
  // components: {
  //     AlignmentStatsCell,
  //     ProgressCell,
  //     ChainsStatsCell
  // }
};

</script>

<template>
  <PlBlockPage>
    <template #title>STAR Read Mapping</template>
    <template #append>
      <PlBtnPrimary :icon="'settings-2'" @click.stop="() => data.settingsOpen = true">Settings</PlBtnPrimary>
    </template>

    <AgGridVue :theme="AgGridTheme" :style="{ height: '100%' }" @grid-ready="onGridReady" :rowData="results"
      :columnDefs="columnDefs" :grid-options="gridOptions" :loadingOverlayComponentParams="{ notReady: true }"
      :loadingOverlayComponent=PlAgOverlayLoading :noRowsOverlayComponent=PlAgOverlayNoRows />


  </PlBlockPage>

  <PlSlideModal v-model="data.settingsOpen">
    <template #title>Settings</template>
    <PlDropdownRef :options="app.model.outputs.dataOptions ?? []" v-model="app.model.args.ref" label="Select dataset"
      clearable />

    <PlDropdown :options="speciesOptions" v-model="app.model.args.species" label="Select species" />
    <PlDropdown :options="inputOptions" v-model="app.model.args.libraryType" label="Select library type" />
    <PlDropdown :options="inputOptionsStr" v-model="app.model.args.strandness" label="Select strandness" />

  </PlSlideModal>

  <PlSlideModal v-model="data.sampleReportOpen" width="80%">
    <template #title>Results for {{ (data.selectedSample ? app.model.outputs.labels?.[data.selectedSample] :
      undefined) ?? "..." }}</template>
    <ReportPanel v-model="data.selectedSample" />
  </PlSlideModal>
</template>
