import '@milaboratory/platforma-uikit/styles';
import { BlockLayout } from '@milaboratory/sdk-vue';
import '@milaboratory/sdk-vue/lib/dist/style.css'; // @todo (will also be sdk-vue/styles)
import { createApp } from 'vue';
import { sdkPlugin } from './app';

createApp(BlockLayout).use(sdkPlugin).mount('#app');