const blockTools = require("@milaboratory/pl-block-tools");

async function loadBlockDescription() {
  return await blockTools.tryLoadPackDescription(__dirname);
}

const blockSpec = {
  type: "dev-v2",
  folder: __dirname,
};

module.exports = {
  blockSpec,
  loadBlockDescription,
};
