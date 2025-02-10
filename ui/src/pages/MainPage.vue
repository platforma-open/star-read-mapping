<script setup lang="ts">
import {
  AgGridTheme,
  makeRowNumberColDef,
  PlAgOverlayLoading,
  PlAgOverlayNoRows,
  PlAgTextAndButtonCell,
  PlBlockPage,
  PlBtnGhost,
  PlDropdown,
  PlDropdownRef,
  PlMaskIcon24,
  PlSlideModal
} from "@platforma-sdk/ui-vue";
import { AgGridVue } from 'ag-grid-vue3';

import { PlRef, plRefsEqual } from '@platforma-sdk/model';
import { ClientSideRowModelModule, ColDef, GridApi, GridOptions, GridReadyEvent, ModuleRegistry } from 'ag-grid-enterprise';
import { computed, reactive, shallowRef } from "vue";
import { useApp } from "../app";
import AlignmentStatsCell from './AlignmentStatsCell.vue';
import ProgressCell from './components/ProgressCell.vue';
import FeatureCountsStatsCell from './FeatureCountsStatsCell.vue';
import ReportPanel from './ReportPanel.vue';
import { resultMap } from './results';

const app = useApp();

const data = reactive<{
  settingsOpen: boolean,
  sampleReportOpen: boolean,
  selectedSample: string | undefined
}>({
  settingsOpen: app.model.args.ref === undefined,
  sampleReportOpen: false,
  selectedSample: undefined,
})

const inputOptionsStr = [
  { text: "Unstranded", value: "0" },
  { text: "Stranded", value: "1" },
  { text: "Reverse stranded", value: "2" },
];

const speciesOptions = [
  { text: "Homo sapiens (GRCh38)", value: "homo-sapiens" },
  { text: "Mus musculus (GRCm39)", value: "mus-musculus" },
  { text: "Saccharomyces cerevisiae (R64-1-1)", value: "saccharomyces-cerevisiae" },
  { text: "Rattus norvegicus (mRatBN7.2)", value: "rattus-norvegicus" },
  { text: "Danio rerio (GRCz11)", value: "danio-rerio" },
  { text: "Drosophila Melanogaster (BDGP6.46)", value: "drosophila-melanogaster" },
  { text: "Arabidopsis Thaliana (TAIR10)", value: "arabidopsis-thaliana" },
  { text: "Caenorhabditis Elegans (WBcel235)", value: "caenorhabditis-elegans" },
  { text: "Gallus Gallus (GRCg7b)", value: "gallus-gallus" },
  { text: "Bos Taurus (ARS-UCD1.3)", value: "bos-taurus" },
  { text: "Sus Scrofa (Sscrofa11.1)", value: "sus-scrofa" },
  { text: "Test genome (v1)", value: "test-species" },
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
      "starQc": resultMap.value[id].starQC,
      "featureCountsQc": resultMap.value[id].featureCountsQC
    });
  }

  return rows;
});

ModuleRegistry.registerModules([ClientSideRowModelModule]);

const gridApi = shallowRef<GridApi<any>>();
const onGridReady = (params: GridReadyEvent) => {
  gridApi.value = params.api;
};

const defaultColDef: ColDef = {
  suppressHeaderMenuButton: true,
  lockPinned: true,
  sortable: false
};

const columnDefs: ColDef[] = [
  makeRowNumberColDef(),
  {
    colId: 'label',
    field: 'sampleLabel',
    headerName: 'Sample',
    pinned: 'left',
    lockPinned: true,
    sortable: true,
    cellRenderer: PlAgTextAndButtonCell,
    cellRendererParams: {
      invokeRowsOnDoubleClick: true
    }
  },
  {
    colId: 'star',
    field: 'star',
    cellRenderer: ProgressCell,
    headerName: 'STAR Progress',
    cellStyle: {
      '--ag-cell-horizontal-padding': '0px',
      '--ag-cell-vertical-padding': '0px'
    },
  },
  {
    colId: 'starQc',
    field: 'starQc',
    headerName: 'Read alignment',
    cellRenderer: 'AlignmentStatsCell',
    cellStyle: {
      '--ag-cell-horizontal-padding': '0px',
      '--ag-cell-vertical-padding': '0px'
    },
    flex: 1,
  },
  {
    colId: 'featureCountsQc',
    field: 'featureCountsQc',
    headerName: 'Features assigned',
    cellRenderer: FeatureCountsStatsCell,
    cellStyle: {
      '--ag-cell-horizontal-padding': '0px',
      '--ag-cell-vertical-padding': '0px'
    },
    flex: 1,
  }
];

const gridOptions: GridOptions = {
  getRowId: (row) => row.data.sampleId,
  onRowDoubleClicked: (e) => {
    data.selectedSample = e.data?.sampleId
    data.sampleReportOpen = data.selectedSample !== undefined;
  },
  components: {
    AlignmentStatsCell,
    FeatureCountsStatsCell,
    PlAgTextAndButtonCell,
    ProgressCell
    //     ProgressCell,
    //     ChainsStatsCell
  }
};

function setInput(inputRef?: PlRef) {
  app.model.args.ref = inputRef;
  if (inputRef)
    app.model.args.title = app.model.outputs.dataOptions?.find(o => plRefsEqual(o.ref, inputRef))?.label
  else
    app.model.args.title = undefined;
}

</script>

<template>
  <PlBlockPage>
    <template #title>STAR Read Mapping</template>
    <template #append>
      <PlBtnGhost @click.stop="() => data.settingsOpen = true">
        Settings
        <template #append>
          <PlMaskIcon24 name="settings" />
        </template>
      </PlBtnGhost>
    </template>

    <AgGridVue :theme="AgGridTheme" :style="{ height: '100%' }" @grid-ready="onGridReady" :rowData="results"
      :columnDefs="columnDefs" :grid-options="gridOptions" :loadingOverlayComponentParams="{ notReady: true }"
      :defaultColDef="defaultColDef" :loadingOverlayComponent=PlAgOverlayLoading
      :noRowsOverlayComponent=PlAgOverlayNoRows />

  </PlBlockPage>

  <PlSlideModal v-model="data.settingsOpen">
    <template #title>Settings</template>
    <PlDropdownRef :options="app.model.outputs.dataOptions" v-model="app.model.args.ref" @update:model-value="setInput"
      label="Select dataset" clearable />

    <PlDropdown :options="speciesOptions" v-model="app.model.args.species" label="Select species" />
    <PlDropdown :options="inputOptionsStr" v-model="app.model.args.strandness" label="Select strandness" />

  </PlSlideModal>

  <PlSlideModal v-model="data.sampleReportOpen" width="80%">
    <template #title>Results for {{ (data.selectedSample ? app.model.outputs.labels?.[data.selectedSample] :
      undefined) ?? "..." }}</template>
    <ReportPanel v-model="data.selectedSample" />
  </PlSlideModal>
</template>
