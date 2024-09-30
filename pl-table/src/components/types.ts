import { Platforma } from "@milaboratory/sdk-ui";

declare global {
    const platforma: Platforma | undefined;
    interface Window {
      platforma: Platforma | undefined;
    }
  }