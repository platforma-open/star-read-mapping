import { model } from '@platforma-open/milaboratories.star-read-mapping.model';
import { defineApp } from '@platforma-sdk/ui-vue';
import PrincipalComponentAnalysis from './pages/PCAPage.vue';
import Settings from './pages/MainPage.vue';


export const sdkPlugin = defineApp(model, () => {
  return {
    showErrorsNotification: true,
    routes: {
      '/': Settings,
      '/PCA': PrincipalComponentAnalysis
    }
  };
});

export const useApp = sdkPlugin.useApp;
