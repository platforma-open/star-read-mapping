import { model } from '@milaboratory/milaboratories.star-read-mapping.model';
import { defineApp } from '@milaboratory/sdk-vue';
import PrincipalComponentAnalysis from './pages/PCAPage.vue';
import Settings from './pages/MainPage.vue';
import SequenceDataQC from './pages/QCPage.vue';
import HierarchicalClustering from './pages/HCPage.vue';


export const sdkPlugin = defineApp(model, () => {
  return {
    routes: {
      '/': Settings,
      '/QC': SequenceDataQC,
      '/PCA': PrincipalComponentAnalysis,
      '/HC': HierarchicalClustering
    }
  };
});

export const useApp = sdkPlugin.useApp;