import { model } from "@platforma-open/milaboratories.star-read-mapping.model";
import { defineApp } from "@platforma-sdk/ui-vue";
import HierarchicalClustering from "./pages/HCPage.vue";
import Settings from "./pages/MainPage.vue";
import PrincipalComponentAnalysis from "./pages/PCAPage.vue";
import SequenceDataQC from "./pages/QCPage.vue";

export const sdkPlugin = defineApp(model, () => {
  return {
    showErrorsNotification: true,
    routes: {
      "/": () => Settings,
      "/QC": () => SequenceDataQC,
      "/PCA": () => PrincipalComponentAnalysis,
      "/HC": () => HierarchicalClustering,
    },
  };
});

export const useApp = sdkPlugin.useApp;
