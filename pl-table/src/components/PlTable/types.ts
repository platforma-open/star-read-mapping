import {
  AxisId,
  LocalBlobHandleAndSize,
  PTableHandle,
  PTableRecordFilter,
  PTableSorting,
  RemoteBlobHandleAndSize,
  ValueOrErrors,
} from "@milaboratory/sdk-ui";

/**
 * Data table settings
 */
export type PlTableSettings =
  | {
      /**
       * The type of the source to feed the data into the table.
       */
      sourceType: "pframe";
      /**
       * PTable handle output
       */
      pTable?: ValueOrErrors<PTableHandle | undefined> | undefined;
      /**
       * Sheets that we want to show in our table
       */
      sheets?: {
        /** id of the axis to use */
        axis: AxisId;
        /** options to show in the filter tan */
        options: {
          /** value of the option (should be one of the values in the axis) */
          value: any;
          /** corresponding label */
          text: string;
        }[];
        /** default (selected) value */
        defaultValue?: any;
      }[];
    }
  | {
      sourceType: "xsv";
      xsvFile?:
        | ValueOrErrors<RemoteBlobHandleAndSize | undefined>
        | ValueOrErrors<LocalBlobHandleAndSize | undefined>
        | undefined;
    };

/** A part of internal ag-grid GridState type */
export type PlTableGridState = {
  /** Includes column ordering */
  columnOrder?: {
    /** All colIds in order */
    orderedColIds: string[];
  };
  /** Includes current sort columns and direction */
  sort?: {
    /** Sorted columns and directions in order */
    sortModel: {
      /** Column Id to apply the sort to. */
      colId: string;
      /** Sort direction */
      sort: "asc" | "desc";
    }[];
  };

  /** current sheet selections */
  sheets?: Record<string, string>;
};

/**
 * Params used to get p-table handle from the driver
 */
export type PTableParams = {
  sorting?: PTableSorting[];
  filters?: PTableRecordFilter[];
};

export type PlTableState = {
  // internal ag-grid state
  gridState: PlTableGridState;

  // mapping of gridState onto the p-table data structures
  pTableParams?: PTableParams;
};
