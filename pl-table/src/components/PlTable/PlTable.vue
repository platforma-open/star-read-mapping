<script lang="ts" setup>
import "./ag-theme.css";

import { ClientSideRowModelModule } from "@ag-grid-community/client-side-row-model";

import canonicalize from "canonicalize";

import {
  GetRowIdParams,
  GridApi,
  GridOptions,
  GridReadyEvent,
  GridState,
  ModuleRegistry,
  StateUpdatedEvent,
} from "@ag-grid-community/core";
import { InfiniteRowModelModule } from "@ag-grid-community/infinite-row-model";
import { AgGridVue } from "@ag-grid-community/vue3";
import { ClipboardModule } from "@ag-grid-enterprise/clipboard";
import { RangeSelectionModule } from "@ag-grid-enterprise/range-selection";
import { deepClone } from "@milaboratories/helpers";
import { PlDropdownLine } from "@milaboratories/uikit";
import { AxisId, PTableRecordFilter } from "@platforma-sdk/model";
import { computed, ref, shallowRef, watch } from "vue";
import OverlayLoading from "./OverlayLoading.vue";
import OverlayNoData from "./OverlayNoData.vue";
import { xsvGridOptions } from "./sources/file-source";
import { parseColId, pFrameGridOptions } from "./sources/pframe-source";
import { PlTableSettings, PlTableState } from "./types";

ModuleRegistry.registerModules([
  ClientSideRowModelModule,
  ClipboardModule,
  RangeSelectionModule,
  InfiniteRowModelModule,
]);

const props = defineProps<{
  settings?: Readonly<PlTableSettings>;
}>();

const model = defineModel<PlTableState>();

//////////////////////////////////////////

// inner states
const innerGridState = ref<GridState>();
// sheetId -> sheet selection
const innerSheetDeltasState = ref<Record<string, string>>();

// get ag row model type based settings
const rowModelType = computed(() => {
  switch (props.settings?.sourceType) {
    case "pframe":
      return "infinite";
    default:
      return "infinite"; // @TODO
  }
});

/** A ref storing calculated grid options with datasource etc */
const gridOptions = ref<GridOptions>();

/** Ref to the grid api instance */
const gridApi = shallowRef<GridApi<any>>();

/**
 * Utility to get axis id string from object
 */
const axisIdFn = (axis: AxisId) => canonicalize(axis)!;

/**
 * Save grid API instance callback
 */
const onGridReady = (params: GridReadyEvent) => {
  gridApi.value = params.api;
};

/**
 * Callback to persist updated Grid state
 */
const onStateUpdated = (params: StateUpdatedEvent<any>) => {
  innerGridState.value = params.state;
};

/**
 * Save sheet information and update filters
 */
const onSheetChanged = (sheetId: string, newValue: string) => {
  if (innerSheetDeltasState.value?.[sheetId] !== newValue) {
    // show loading spinner as we will wait for data update
    gridApi.value?.setGridOption("loading", true);

    const newState = deepClone(innerSheetDeltasState.value) ?? {};
    newState[sheetId] = newValue;
    innerSheetDeltasState.value = newState;
  }
};

const updateSheets = (sheets: Record<string, string> | undefined) => {
  if (!sheets) {
    return model.value?.gridState.sheets;
  }
  const r = model.value?.gridState.sheets ?? {};
  for (const [k, v] of Object.entries(sheets)) {
    r[k] = v;
  }
  return r;
};

// calculate & update state
watch(
  [innerGridState, innerSheetDeltasState],
  ([gridState, sheetStateDelta]) => {
    const sheetState = updateSheets(sheetStateDelta);
    // persist only the desired part of the state
    const data: PlTableState = {
      gridState: {
        columnOrder: gridState?.columnOrder,
        sort: gridState?.sort,
        sheets: sheetState,
      },
    };

    const source = props.settings?.sourceType;
    if (source === "pframe") {
      const sorting = gridState?.sort?.sortModel?.map((v) => ({
        column: parseColId(v.colId).id,
        ascending: v.sort == "asc",
        naAndAbsentAreLeastValues: true,
      }));

      const filters: PTableRecordFilter[] | undefined =
        props.settings?.sheets?.map((sheet) => {
          return {
            type: "bySingleColumn",
            column: {
              type: "axis",
              id: sheet.axis,
            },
            predicate: {
              operator: "Equal",
              reference:
                sheetState?.[axisIdFn(sheet.axis)] ??
                sheet.defaultValue ??
                sheet.options[0].value,
            },
          };
        });

      // calculate sort
      data.pTableParams = {
        sorting: sorting,

        filters: filters,
      };
    }

    model.value = data;
  },
  { immediate: true }
);

watch(
  [() => canonicalize(props.settings)],
  async ([_]) => {
    const { settings } = props;
    const source = settings?.sourceType;

    let gOps: GridOptions | undefined = undefined;
    switch (source) {
      case "pframe":
        const pfDriver = window.platforma?.pFrameDriver;
        if (!pfDriver) {
          throw Error("no p-frames driver found");
        }
        if (settings?.pTable && settings.pTable.ok && settings.pTable.value)
          gOps = await pFrameGridOptions(
            pfDriver,
            settings.pTable.value,
            gridApi.value
          );

        // set ag-grid data
        gridApi.value?.setGridOption("columnDefs", gOps?.columnDefs);
        gridApi.value?.setGridOption("datasource", gOps?.datasource);
        break;

      case "xsv":
        const blobDriver = window.platforma?.blobDriver;
        if (!blobDriver) {
          throw Error("no blob driver found");
        }
        if (
          settings?.xsvFile &&
          settings.xsvFile.ok &&
          settings.xsvFile.value
        ) {
          gOps = await xsvGridOptions(
            blobDriver,
            settings.xsvFile.value,
            gridApi.value
          );
          // @TODO
        }
    }

    gridOptions.value = gOps;
  },
  { immediate: true }
);

/**
 * A row id function. All sources must set "id" property of the rows.
 */
const getRowId = (params: GetRowIdParams) => {
  if (!params.data.id) {
    throw Error("id must be defined");
  }
  return params.data.id;
};
</script>

<template>
  <div v-if="settings?.sourceType === 'pframe' && settings.sheets">
    <div v-for="sheet in settings.sheets">
      <PlDropdownLine
        :model-value="
          model?.gridState?.sheets?.[axisIdFn(sheet.axis)] ??
          sheet.defaultValue ??
          sheet.options[0].value
        "
        @update:model-value="
          (newValue) => onSheetChanged(axisIdFn(sheet.axis), newValue)
        "
        :options="sheet.options"
        style="z-index: 3"
      />
    </div>
  </div>

  <div :style="{ height: '700px' }">
    <AgGridVue
      @grid-ready="onGridReady"
      @state-updated="onStateUpdated"
      :enableRangeSelection="true"
      :style="{ height: '100%', width: '100%' }"
      :noRowsOverlayComponent="OverlayNoData"
      :loadingOverlayComponent="OverlayLoading"
      :loadingOverlayComponentParams="{ noData: true }"
      :initialState="model?.gridState"
      :getRowId="getRowId"
      :rowModelType="rowModelType"
      :maxBlocksInCache="gridOptions?.maxBlocksInCache"
      :cacheBlockSize="gridOptions?.cacheBlockSize"
    />
  </div>
</template>
