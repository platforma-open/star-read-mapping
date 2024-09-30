import { Platforma } from "@platforma-sdk/model";

declare global {
    const platforma: Platforma | undefined;
    interface Window {
      platforma: Platforma | undefined;
    }
  }