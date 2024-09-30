import '@milaboratories/uikit/styles';
import { BlockLayout } from '@platforma-sdk/ui-vue';
import '@platforma-sdk/ui-vue/dist/style.css'; // @todo (will also be sdk-vue/styles)
import { createApp } from 'vue';
import { sdkPlugin } from './app';

createApp(BlockLayout).use(sdkPlugin).mount('#app');