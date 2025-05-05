
import React, { useState, useEffect, useMemo } from "react";
import { FormControl, Grid, Paper, IconButton, Icon } from "@material-ui/core";
import { ToggleButtonGroup, ToggleButton, Box } from "@mui/material";
import Plot from "react-plotly.js";
import HeatmapComponent from "./ScHeatmap";
import MarkersTableComponent from "./ScMarkersTable";
// import axios from "axios";
// import * as Plotly from "plotly.js";
// import Demo from "./ScVega";

import CircularProgress from "@mui/material/CircularProgress";

import { SearchAndHighlight, useSNHControls } from "./SearchAndHighlight";
import { useRemoteData } from "../hooks/remoteData";
import { FeatureTable } from "./ScFeatureTable";
import { useBrushedPlotting } from "../hooks/brushedPlotting";
import { PredictiveInput, CreatablePredictiveInput } from "./PredictiveInput";
// import { QueryClient, QueryClientProvider } from "react-query";
// Each viewer contains one or more panels
// each panel can display a plot of any kind, using the basedata associated with the choice
// brushing action can be applied to plots on each panel, using either basedata or remote data.
// At the viewer level, we can choose number of panels, hide or show controls
// at the panel level, we can choose the layer, projection, and feature to plot.
// Each panel can independently be associated with a dataset.

// const queryClient = new QueryClient();
export function ScPanel(props) {
  const {
    classes,
    settings,
    panelDataState,
    panelIdx,
    openChooseDataset,
    loadDataPanel,
    maximized,
    toggleMaximize,
  } = props;
  
  const datasetTerms = useMemo(
    () => settings.datasets.map((ds, idx) => ({ value: idx, label: ds.name })),
    [settings]
  );
  
  const { status, error, data: panelData } = panelDataState;
  const { data, dataAPI, genesets, name } = panelData;
  
  const [showSelectDataset, setShowSelectDataset] = useState(false);
  
  const availableProjectionsAllLayers = {};
  
  data.layers.forEach((layer) => {
    const tc = data.layersData[layer.id].covs;
    let ret = [];
    if ("UMAP_1" in tc || "umap_1" in tc) ret.push("umap");
    if ("TSNE_1" in tc || "tsne_1" in tc || "tSNE_1" in tc) ret.push("tsne");
    availableProjectionsAllLayers[layer.id] = ret;
  });
  
  function getDefaultClustering(tcs) {
    let nv=null;
    if(tcs!=null &&tcs.length > 0) {
      const tcsDefault = tcs.filter(tc=>tc.toLowerCase().includes("default"));
      const tcsCurated = tcs.filter(tc=>tc.toLowerCase().includes("curated"));
      const tcsPrelim = tcs.filter(tc=>tc.toLowerCase().includes("prelim"));
      const tcsMain = tcs.filter(tc=>tc.toLowerCase().includes("main"));
  
      if(tcsDefault.length>0) nv=tcsDefault[0];
      else if(tcsCurated.length>0) nv=tcsCurated[0];
      else if(tcsPrelim.length>0) nv=tcsPrelim[0];
      else if(tcsMain.length>0) nv=tcsMain[0];
      else if(tcs.length>0) nv=tcs[0];
    }
    return {value:nv,label:nv};
  }
  
  const panelDefault = {
    projection: availableProjectionsAllLayers[data.layers[0].id][0],
    feature: null,
    layer: data.layers[0].id,
    mode: ["plot"],
    modeAssay: data.meta.assays[0],
    heatmapMode: 'heatmap',
    clustering: getDefaultClustering(data.layersData[data.layers[0].id].clusterings),
    selectedHeatmapAssay: null
  };
  
  const [config, updateConfig] = useState(panelDefault);
  const clusterings = data.layersData[config.layer].clusterings;
  
  const legendsOptions=[
    { name: "plot", icon: "scatter_plot" }
  ];
  
  if(clusterings != null && Object.keys(clusterings).length > 0) {
    legendsOptions.push({ name: "heatmap", icon: "grid_view"})
  }
  
  legendsOptions.push({ name: "legends", icon: "table" });
  
  const [remoteData, loadRemoteData] = useRemoteData({
    dataAPI,
    dataConfig: config,
    covs: data.layersData[config.layer].covs,
  });
  
  const [plot, brushing] = useBrushedPlotting(
    config,
    data.layersData,
    remoteData
  );
  
  const [sNHControls, toggleSNHControls] = useSNHControls();
  
  const availableProjections = availableProjectionsAllLayers[config.layer];
  
  // const [plotDiv, setPlotDiv] = useState(null);
  
  plot.layout.showlegend = maximized;
  plot.layout.margin.r = maximized ? 0 : 0;
  
  function toggleShowWhatEmbedding(nv) {
    if (nv == null) return;
    let newConfig = { ...config };
    newConfig.projection = nv;
    updateConfig(newConfig);
  }
  
  function toggleLayer(nv) {
    if (nv == null) return;
    let newConfig = { ...config };
    newConfig.layer = nv;
    newConfig.clustering = getDefaultClustering(data.layersData[nv].clusterings);
    setClusteringData(null);
    updateConfig(newConfig);
  }
  
  function toggleZTransform(nv) {
    if(nv===null) return;
    console.log('toggling ztransform');
    console.log(nv);
    console.log(`setting brushing.zTransform to ${nv}`);
    brushing.setZTransform(nv);

  }

  function toggleMode(nv) {
    if (nv == null) return;
    let newConfig = { ...config };
    newConfig.mode = nv;
    updateConfig(newConfig);
  }
  
  function toggleModeAssay(nv) {
    if (nv == null) return;
    let newConfig = { ...config };
    newConfig.modeAssay = nv;
    updateConfig(newConfig);
  }
  
  function toggleHeatmapMode(nv) {
    if (nv == null) return;
    let newConfig = { ...config };
    newConfig.heatmapMode = nv;
    updateConfig(newConfig);
  }
  
  function renderSNHControlsToggleIcon() {
    switch (sNHControls) {
      case "none":
        return "fullscreen_exit";
      case "with-select-paint":
        return "fullscreen_exit";
      case "with-select-paint-geneset":
        return "fullscreen";
      default:
        return null;
    }
  }
  
  function downloadReports() {
    window.open(`report/${data.zip}`);
  }
  
  function handleSelectedDatasetChange(e, nv) {
    setShowSelectDataset(false);
    if (e.value === panelDataState.data.datasetIdx) return;
    loadDataPanel(e.value, panelIdx);
  }
  
  function handleClusteringChange(nv) {
    if (nv == null) return;
    if (nv.value == config.clustering.value) return;
    let newConfig = { ...config };
    newConfig.clustering = nv;
    setClusteringData(null);
    updateConfig(newConfig);
  }
  

  const [clusteringData, setClusteringData] = useState(null);
  const [isLoading, setIsLoading] = useState(true);
  const [usingFullMatrix, setFullMatrix] = useState(true);
  
  useEffect( () => {
    if (config.clustering == null || config.clustering.value === null) return;
    const url = `${data.url}/${config.layer}/clustering/${config.clustering.value}${usingFullMatrix ? '/true' : ''}`;
    setIsLoading(true);  
    // setClusteringData(null);
    const fetchData = async () => {
      const clusteringData = await (await fetch(url)).json();
      setClusteringData(clusteringData);
      setIsLoading(false);
    };
    fetchData();
  }, [config.clustering]);
  
  let heatmapAssayInEffect = config.selectedHeatmapAssay;
  if(clusteringData) {
    const keys = Object.keys(clusteringData);
    if(keys.includes('markers')) 
      heatmapAssayInEffect = null;
    else if(heatmapAssayInEffect==null || !keys.includes(heatmapAssayInEffect)) {
      heatmapAssayInEffect = keys.pop();
    } 
  }
  // console.log(`heatmapAssayInEffect ${heatmapAssayInEffect}`);
  
  function updateSelectedHeatmapAssay(nv) {
    let newConfig = { ...config, selectedHeatmapAssay: nv };
    updateConfig(newConfig);
  }
  
  function toggleSelectedHeatmapAssay(nv) {
    if (nv == null) return;
    updateSelectedHeatmapAssay(nv);
  }
  
  function activeModes() {
    return legendsOptions.filter(r => config.mode.includes(r.name)).map(r => r.name);
  }
  
  function DatasetControls() {
    return showSelectDataset ? (
      <FormControl className={classes.chooseDatasetControl} variant="outlined">
        <PredictiveInput
          classes={classes}
          options={datasetTerms}
          value={{
            value: panelDataState.data.datasetIdx,
            label: panelDataState.data.name,
          }}
          onChange={handleSelectedDatasetChange}
          defaultMenuIsOpen={true}
          onMenuClose={() => {
            setShowSelectDataset(false);
          }}
          placeholder="Select "
        />
      </FormControl>
    ) : (
      <>
        <IconButton color="primary" onClick={() => setShowSelectDataset(true)}>
          <Icon>folder</Icon>
        </IconButton>
        <span>
          <i>{name}</i>
        </span>
      </>
    );
  }
  
  
  function HeatmapContainer() {
    // const { heatmapMode, clustering, selectedHeatmapAssay, layer } = config;
    const { heatmapMode, clustering, layer } = config;
    const isActiveHeatmap = activeModes().includes("heatmap");
  
    if (!isActiveHeatmap) return null;
  
    const handleToggleHeatmapMode = (_, newValue) => toggleHeatmapMode(newValue);
    const handleToggleSelectedHeatmapAssay = (_, newValue) => toggleSelectedHeatmapAssay(newValue);
  
    const style1 = { margin: "0px", top: 0, marginLeft: '100px' };
    const style2 = { margin: "0px", top: 0, marginLeft: '20px' };
    const style3 = { margin: "0px", top: 0, marginLeft: activeModes().includes('plot')|!maximized?'0px':'300px' };

    const isToggleGroupVisible = !isLoading && isActiveHeatmap;
    const isClusteringDataValid = !isLoading && sNHControls !== 'none' && clusteringData && !Object.keys(clusteringData).includes('markers');
    
    return (
      <Grid item xs={gridItemSize} className={classes.plot} style={{overflowY: 'scroll' }}>
        {sNHControls !== "none"  && (
          <div className={maximized?classes.floatingControls:classes.fixedControls} style={{ top: 0 }}>
            <FormControl className={classes.parintControl} style={style3} variant="outlined">
              <PredictiveInput
                classes={classes}
                options={clusterings.map(c => ({ value: c, label: c }))}
                value={clustering}
                onChange={handleClusteringChange}
                placeholder="Select"
              />
            </FormControl>
    
            {isToggleGroupVisible && (
              <ToggleButtonGroup style={style1} color="primary" value={heatmapMode} exclusive onChange={handleToggleHeatmapMode}>
                <ToggleButton key="toggleHeatmapMode-plot" value="heatmap">HEATMAP</ToggleButton>
                <ToggleButton key="toggleHeatmapMode-markers-table" value="markers-table">MARKERS</ToggleButton>
              </ToggleButtonGroup>
            )}
    
            {isClusteringDataValid && (
              // <ToggleButtonGroup style={style2} color="primary" value={selectedHeatmapAssay} exclusive onChange={handleToggleSelectedHeatmapAssay}>
              <ToggleButtonGroup style={style2} color="primary" value={heatmapAssayInEffect} exclusive onChange={handleToggleSelectedHeatmapAssay}>
              {Object.keys(clusteringData).map((heatmapAssay, index) => (
                  <ToggleButton key={`toggleSelectedHeatmapAssay-${index}`} value={heatmapAssay}>
                    {heatmapAssay}
                  </ToggleButton>
                ))}
              </ToggleButtonGroup>
            )}
          </div>
        )}
    
        <div style={{ marginTop: 50, width: "100%"}}>
          {isLoading && (
            <Box sx={{ display: "flex", justifyContent: "center" }}>
              <CircularProgress />
            </Box>
          )}
    
          {!isLoading && clusteringData && heatmapMode === "markers-table" && (
            <MarkersTableComponent
              clusteringData={clusteringData}
              // selectedHeatmapAssay={selectedHeatmapAssay}
              selectedHeatmapAssay={heatmapAssayInEffect}
            />
          )}
    
    {!isLoading && clusteringData && heatmapMode === "heatmap" && (
          // {heatmapMode === "heatmap" && (
            <HeatmapComponent
              data={clusteringData}
              usingFullMatrix={usingFullMatrix}
              clustering={clustering.value}
              layer={layer}
              // selectedHeatmapAssay={selectedHeatmapAssay}
              heatmapAssayInEffect={heatmapAssayInEffect}
              sNHControls={sNHControls}
            />
          )}
        </div>
      </Grid>
    );
  }
  
  



  const isInPlotMode = config.mode.includes("plot");
  const isInLegendsMode = config.mode.includes("legends");
  const gridItemSize = 12 / activeModes().length;
  const isPlotOrHeatmap = isInPlotMode || config.mode.includes("heatmap");
  const style = { margin: "5px" };


  function PlotAndLayerControls() {
    if (sNHControls === "none" ) return null;
  
    const handleToggleMode = (event, newValue) => toggleMode(newValue);
    const handleShowWhatEmbedding = (event, newValue) => toggleShowWhatEmbedding(newValue);
    const handleToggleLayer = (event, newValue) => toggleLayer(newValue);
    const handleToggleZTransform = (event, newValue) => toggleZTransform(newValue);
    return (
      <Grid container className={classes.toggleProjection}>
        <ToggleButtonGroup style={style} color="primary" value={config.mode} exclusive={!maximized} onChange={handleToggleMode}>
          {legendsOptions.map(m => (
            <ToggleButton data-test-id="mode-type" key={`mode_${m.name}`} value={m.name}>
              <Icon>{m.icon}</Icon>
            </ToggleButton>
          ))}
        </ToggleButtonGroup>
        {isInPlotMode && maximized && false ? 
          <ToggleButtonGroup style={style} color="primary" value={brushing.zTransform} exclusive onChange={handleToggleZTransform}>
            {[{value:false,name:"default"}, {value:true,name:"contrast"}].map(
              l => <ToggleButton key={`toggleZTransform_${l.name}`} value={l.value}>{l.name}</ToggleButton>)}
          </ToggleButtonGroup> : null}

  
        {isInPlotMode && (
          <ToggleButtonGroup style={style} color="primary" value={config.projection} exclusive onChange={handleShowWhatEmbedding}>
            {availableProjections.map(p => (
              <ToggleButton data-test-id="projection-type" key={`toggleProjection${p}`} value={p}>
                {p}
              </ToggleButton>
            ))}
          </ToggleButtonGroup>
        )}
  
        {isPlotOrHeatmap && (
          <ToggleButtonGroup style={style} color="primary" value={config.layer} exclusive onChange={handleToggleLayer}>
            {data.layers.map(l => (
              <ToggleButton data-test-id="layer-type" key={`layer_${l.id}`} value={l.id}>
                {l.name}
              </ToggleButton>
            ))}
          </ToggleButtonGroup>
        )}


      </Grid>
    );
  }
  
  const handleToggleModeAssay = (event, newValue) => toggleModeAssay(newValue);
  const style1 = { margin: "0px", top: 0, marginLeft: '100px' };


  return (
    <Paper className={classes.paper}>
      <Grid container className={classes.headerPanel} alignItems="center">
        <DatasetControls />
      </Grid>
  
      {isInPlotMode && (
        <div className={classes.floatingControls}>
          <SearchAndHighlight
            classes={classes}
            covs={data.layersData[config.layer].covs}
            features={data.meta.features}
            continuousCovs={data.layersData[config.layer].continuousCovs}
            discreteCovs={data.layersData[config.layer].discreteCovs}
            brushing={brushing}
            genesets={genesets}
            remoteData={remoteData}
            loadRemoteData={loadRemoteData}
            sNHControls={sNHControls}
          />
        </div>
      )}
  
      <Grid container className={classes.plotMain} spacing={0}>
        {isInPlotMode && (
          <Grid item xs={gridItemSize} className={classes.plot}>
            <Plot
              className={classes.plot}
              data={plot.data}
              layout={plot.layout}
              config={plot.config}
              useResizeHandler={plot.useResizeHandler}
              // onInitialized={(fig, div) => setPlotDiv(div)}
              // onUpdate={() => Plotly.newPlot(plotDiv, plot.data, plot.layout, plot.config)}
              // onWebGlContextLost={() => Plotly.newPlot(plotDiv, plot.data, plot.layout, plot.config)}
            />
          </Grid>
        )}
  
        <PlotAndLayerControls />
        <HeatmapContainer />
  
        {isInLegendsMode && (
          <Grid item xs={gridItemSize} className={classes.plot} style={{ overflowY: 'scroll' }}>
            {sNHControls !== "none"  && (
              <div className={maximized?classes.floatingControls:classes.fixedControls}  style={{ top: 0 }}>

                <ToggleButtonGroup style={style1} color="primary" value={config.modeAssay} exclusive onChange={handleToggleModeAssay}>
                  {data.meta.assays.concat("genesets").map(p => (
                    <ToggleButton data-test-id="assay-type" key={`toggleModeAssay${p}`} value={p}>
                      {p}
                    </ToggleButton>
                  ))}
                </ToggleButtonGroup>
              </div>)}
            <FeatureTable
              // className={classes.plot}
              style={{ overflow: 'scroll' }}
              features={config.modeAssay === "genesets" ? genesets : data.meta.features[config.modeAssay]}
              assay={config.modeAssay}
            />
          </Grid>
        )}
      </Grid>
  
      <Grid container className={classes.footerPanel} alignItems="center">
        <div style={{ margin: "auto" }}>
          <IconButton onClick={toggleMaximize} color="primary">
            <Icon>view_column</Icon>
          </IconButton>
          <IconButton onClick={toggleSNHControls} color="primary">
            <Icon>{renderSNHControlsToggleIcon()}</Icon>
          </IconButton>
  
          {data.zip && (
            <IconButton style={{ marginLeft: "20px" }} onClick={downloadReports} color="primary">
              <Icon>download</Icon>
            </IconButton>
          )}
        </div>
      </Grid>
    </Paper>
  );
  
}
