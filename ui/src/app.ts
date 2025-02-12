import { model } from '@platforma-open/milaboratories.star-read-mapping.model';
import { defineApp } from '@platforma-sdk/ui-vue';
import Settings from './pages/MainPage.vue';
import PrincipalComponentAnalysis from './pages/PCAPage.vue';
import SampleDistancesHeatmap from './pages/SDistPage.vue';

export const sdkPlugin = defineApp(model, (app) => {
  return {
    progress: () => {
      return app.model.outputs.isRunning;
    },
    showErrorsNotification: true,
    routes: {
      '/': () => Settings,
      '/PCA': () => PrincipalComponentAnalysis,
      '/SDist': () => SampleDistancesHeatmap,
    },
  };
});

export const useApp = sdkPlugin.useApp;
