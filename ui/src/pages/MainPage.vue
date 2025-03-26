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
  PlSlideModal,
} from '@platforma-sdk/ui-vue';
import { AgGridVue } from 'ag-grid-vue3';

import type { PlRef } from '@platforma-sdk/model';
import { plRefsEqual } from '@platforma-sdk/model';
import type { ColDef, GridApi, GridOptions, GridReadyEvent } from 'ag-grid-enterprise';
import { ClientSideRowModelModule, ModuleRegistry } from 'ag-grid-enterprise';
import { computed, reactive, shallowRef } from 'vue';
import { useApp } from '../app';
import AlignmentStatsCell from './AlignmentStatsCell.vue';
import ProgressCell from './components/ProgressCell.vue';
import FeatureCountsStatsCell from './FeatureCountsStatsCell.vue';
import ReportPanel from './ReportPanel.vue';
import { resultMap } from './results';

const app = useApp();

const data = reactive<{
  settingsOpen: boolean;
  sampleReportOpen: boolean;
  selectedSample: string | undefined;
}>({
  settingsOpen: app.model.args.ref === undefined,
  sampleReportOpen: false,
  selectedSample: undefined,
});

const inputOptionsStr = [
  { text: 'Unstranded', value: '0' },
  { text: 'Stranded', value: '1' },
  { text: 'Reverse stranded', value: '2' },
];

const speciesOptions = [
  { text: 'Homo sapiens (GRCh38)', value: 'homo-sapiens' },
  { text: 'Mus musculus (GRCm39)', value: 'mus-musculus' },
  { text: 'Saccharomyces cerevisiae (R64-1-1)', value: 'saccharomyces-cerevisiae' },
  { text: 'Rattus norvegicus (mRatBN7.2)', value: 'rattus-norvegicus' },
  { text: 'Danio rerio (GRCz11)', value: 'danio-rerio' },
  { text: 'Drosophila Melanogaster (BDGP6.46)', value: 'drosophila-melanogaster' },
  { text: 'Arabidopsis Thaliana (TAIR10)', value: 'arabidopsis-thaliana' },
  { text: 'Caenorhabditis Elegans (WBcel235)', value: 'caenorhabditis-elegans' },
  { text: 'Gallus Gallus (GRCg7b)', value: 'gallus-gallus' },
  { text: 'Bos Taurus (ARS-UCD1.3)', value: 'bos-taurus' },
  { text: 'Sus Scrofa (Sscrofa11.1)', value: 'sus-scrofa' },
  { text: 'Test genome (v1)', value: 'test-species' },
];

/** Rows for ag-table */
const results = computed(() => {
  if (resultMap.value === undefined) return undefined;
  const rows = [];
  for (const id in resultMap.value) {
    rows.push({
      sampleId: id,
      sampleLabel: resultMap.value[id].sampleLabel,
      star: resultMap.value[id].starProgressLine, // @TODO status?
      subread: 'Running', // @TODO status?
      starQc: resultMap.value[id].starQC,
      featureCountsQc: resultMap.value[id].featureCountsQC,
    });
  }

  return rows;
});

ModuleRegistry.registerModules([ClientSideRowModelModule]);

const gridApi = shallowRef<GridApi>();
const onGridReady = (params: GridReadyEvent) => {
  gridApi.value = params.api;
};

const defaultColDef: ColDef = {
  suppressHeaderMenuButton: true,
  lockPinned: true,
  sortable: false,
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
      invokeRowsOnDoubleClick: true,
    },
  },
  {
    colId: 'star',
    field: 'star',
    cellRenderer: ProgressCell,
    headerName: 'STAR Progress',
    cellStyle: {
      '--ag-cell-horizontal-padding': '0px',
      '--ag-cell-vertical-padding': '0px',
    },
  },
  {
    colId: 'starQc',
    field: 'starQc',
    headerName: 'Read alignment',
    cellRenderer: 'AlignmentStatsCell',
    cellStyle: {
      '--ag-cell-horizontal-padding': '0px',
      '--ag-cell-vertical-padding': '0px',
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
      '--ag-cell-vertical-padding': '0px',
    },
    flex: 1,
  },
];

const gridOptions: GridOptions = {
  getRowId: (row) => row.data.sampleId,
  onRowDoubleClicked: (e) => {
    data.selectedSample = e.data?.sampleId;
    data.sampleReportOpen = data.selectedSample !== undefined;
  },
  components: {
    AlignmentStatsCell,
    FeatureCountsStatsCell,
    PlAgTextAndButtonCell,
    ProgressCell,
    //     ProgressCell,
    //     ChainsStatsCell
  },
};

function setInput(inputRef?: PlRef) {
  app.model.args.ref = inputRef;
  if (inputRef)
    app.model.args.title = app.model.outputs.dataOptions?.find((o) => plRefsEqual(o.ref, inputRef))?.label;
  else
    app.model.args.title = undefined;
}

// // @TODO: Uncoment when user defined nCPU is ok with deduplication
// // variables for CPU usage
// type LocalState = {
//   tab: 'percent' | 'number' | undefined;
// };

// const state = reactive<LocalState>({
//   tab: undefined,
// });

// const cpuSelectionOptions = [
//   { text: '% from all available', value: 'percent' },
//   { text: 'User defined number', value: 'number' },
// ] as const satisfies ListOption [];

// const computedTab = computed({
//   // return tab or first tab assignment given cpuNumber value
//   get() {
//     return state.tab ?? (app.model.args.cpuNumber == -1 ? 'percent' : 'number');
//   },
//   set(tab) {
//     state.tab = tab;
//   },
// });

// // Make user see that selection in one tab overwrites value in the other one
// watch(() => app.model.args.cpuNumber, (_) => {
//   if (app.model.args.cpuNumber !== -1) {
//     app.model.args.cpuPercent = 0;
//   }
// });

// watch(() => app.model.args.cpuPercent, (_) => {
//   if (app.model.args.cpuPercent !== 0) {
//     app.model.args.cpuNumber = 0;
//   }
// });

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

    <AgGridVue
      :theme="AgGridTheme" :style="{ height: '100%' }" :rowData="results" :columnDefs="columnDefs"
      :grid-options="gridOptions" :loadingOverlayComponentParams="{ notReady: true }" :defaultColDef="defaultColDef"
      :loadingOverlayComponent="PlAgOverlayLoading" :noRowsOverlayComponent="PlAgOverlayNoRows"
      @grid-ready="onGridReady"
    />
  </PlBlockPage>

  <PlSlideModal v-model="data.settingsOpen">
    <template #title>Settings</template>
    <PlDropdownRef
      v-model="app.model.args.ref" :options="app.model.outputs.dataOptions" label="Select dataset"
      clearable @update:model-value="setInput"
    />

    <PlDropdown v-model="app.model.args.species" :options="speciesOptions" label="Select species" />
    <PlDropdown v-model="app.model.args.strandness" :options="inputOptionsStr" label="Select strandness" />

    <!-- Content hidden until you click ADDITIONAL SETTINGS -->
    <!-- // @TODO: Uncoment when user defined nCPU is ok with deduplication
    <PlAccordionSection label="ADDITIONAL SETTINGS">
      <PlBtnGroup v-model="computedTab" :options="cpuSelectionOptions" label="CPU allocation criteria">
        <template #tooltip>
          Select the CPU resources to allocate for the task. You can either specify a percentage of all available CPUs or assign a specific CPU number
        </template>
      </PlBtnGroup>
      <Slider
        v-if="computedTab === 'percent'"
        v-model="app.model.args.cpuPercent" :min="10" :max="100" :step="10"
        :breakpoints="true" label="CPU %"
      />
      <PlNumberField
        v-if="computedTab === 'number'"
        v-model="app.model.args.cpuNumber" label="Nº CPU" :minValue="1"
      />
    </PlAccordionSection> -->
  </PlSlideModal>

  <PlSlideModal v-model="data.sampleReportOpen" width="80%">
    <template #title>
      Results for {{ (data.selectedSample ? app.model.outputs.labels?.[data.selectedSample] :
        undefined) ?? "..." }}
    </template>
    <ReportPanel v-model="data.selectedSample" />
  </PlSlideModal>
</template>
