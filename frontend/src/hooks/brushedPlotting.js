import { useBrushing } from './brushing';
import { usePlotting } from './plotting';

export function useBrushedPlotting(conf, layersData, remoteData) {
  const { covs, discreteCovs, continuousCovs } = layersData[conf.layer];

  let xys;
  switch (conf.projection) {
    case "umap":
      xys = "UMAP_1" in covs ? [covs.UMAP_1, covs.UMAP_2] : [covs["umap_1"], covs["umap_2"]];
      break;
    case "tsne":
      xys = "tSNE_1" in covs ? [covs.tSNE_1, covs.tSNE_2] : [covs["tsne_1"], covs["tsne_2"]];
      break;
    default:
      xys = [null, null];
  }

  const brushing = useBrushing(layersData[conf.layer]);

  const plot = usePlotting({
    xys, 
    colorBy: brushing.colorBy, 
    selected: brushing.selected,
    covs, 
    discreteCovs, 
    continuousCovs, 
    contrast: brushing.contrast, 
    zTransform: brushing.zTransform,
    remoteData
  });

  return [plot, brushing];
}
