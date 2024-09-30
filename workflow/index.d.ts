declare type TemplateFromFile = { readonly type: "from-file"; readonly path: string; };
declare type TplName = "star-alignment" | "main";
declare const Templates: Record<TplName, TemplateFromFile>;
export { Templates };
