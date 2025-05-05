diff --git a/src/components/ScPanel.js b/src/components/ScPanel.js
index a4b0b5d..e4bd6fa 100644
--- a/src/components/ScPanel.js
+++ b/src/components/ScPanel.js
@@ -6,7 +6,7 @@ import Plot from "react-plotly.js";
import HeatmapComponent from "./ScHeatmap";
import MarkersTableComponent from "./ScMarkersTable";
// import axios from "axios";
-// import * as Plotly from "plotly.js";
+import * as Plotly from "plotly.js";
// import Demo from "./ScVega";
 
import CircularProgress from "@mui/material/CircularProgress";
@@ -16,6 +16,7 @@ import { useRemoteData } from "../hooks/remoteData";
import { FeatureTable } from "./ScFeatureTable";
import { useBrushedPlotting } from "../hooks/brushedPlotting";
import { PredictiveInput, CreatablePredictiveInput } from "./PredictiveInput";
+// import { isArray } from "util";
// import { QueryClient, QueryClientProvider } from "react-query";
// Each viewer contains one or more panels
// each panel can display a plot of any kind, using the basedata associated with the choice
@@ -79,10 +80,11 @@ export function ScPanel(props) {
     feature: null,
     layer: data.layers[0].id,
     mode: ["plot"],
-    modeAssay: data.meta.assays[0],
+    // modeAssay: data.meta.assays[0],
+    modeAssay: Object.keys(genesets)[0],
     heatmapMode: 'heatmap',
     clustering: getDefaultClustering(data.layersData[data.layers[0].id].clusterings),
-    selectedHeatmapAssay: null
+    selectedHeatmapAssay: [null]
   };
   
   const [config, updateConfig] = useState(panelDefault);
@@ -114,7 +116,7 @@ export function ScPanel(props) {
   
   const availableProjections = availableProjectionsAllLayers[config.layer];
   
-  // const [plotDiv, setPlotDiv] = useState(null);
+  const [plotDiv, setPlotDiv] = useState(null);
   
   plot.layout.showlegend = maximized;
   plot.layout.margin.r = maximized ? 0 : 0;
@@ -208,22 +210,47 @@ export function ScPanel(props) {
     setIsLoading(true);  
     // setClusteringData(null);
     const fetchData = async () => {
-      const clusteringData = await (await fetch(url)).json();
-      setClusteringData(clusteringData);
-      setIsLoading(false);
+      try {
+        const clusteringData = await (await fetch(url)).json()
+        if(!clusteringData || Object.keys(clusteringData).length===0)
+          throw(new Error("Empty response from server."));
+        setClusteringData(clusteringData);  
+      } catch (e) {
+        console.log(`Failed fetching data from ${url}. ${e}`);
+        alert(`Failed fetching data from ${url}. ${e}`);
+      } finally {
+        console.log(`finally`);
+        setIsLoading(false);
+      }        
     };
     fetchData();
   }, [config.clustering]);
   
   let heatmapAssayInEffect = config.selectedHeatmapAssay;
-  if(clusteringData) {
+
+  heatmapAssayInEffect = Array.isArray(heatmapAssayInEffect) &&  heatmapAssayInEffect.length >0
+  ? heatmapAssayInEffect
+  : [heatmapAssayInEffect];
+
+
+  if (clusteringData) {
     const keys = Object.keys(clusteringData);
-    if(keys.includes('markers'))
-      heatmapAssayInEffect = null;
-    else if(heatmapAssayInEffect==null || !keys.includes(heatmapAssayInEffect)) {
-      heatmapAssayInEffect = keys.pop();
-    }
+    const defaultKeys = maximized ? keys : [keys.pop()];
+  
+    const isAssayInEffectInvalid = (aie) => {
+      return maximized ? aie.some(i => !keys.includes(i)) : aie.length > 1 || !keys.includes(aie[0]);
+    };
+  
+    if (keys.includes('markers')) {
+      heatmapAssayInEffect = [null];
+    } else {
+      heatmapAssayInEffect = isAssayInEffectInvalid(heatmapAssayInEffect)
+        ? defaultKeys
+        : heatmapAssayInEffect;
+    }
   }
+
+  
   // console.log(`heatmapAssayInEffect ${heatmapAssayInEffect}`);
   
   function updateSelectedHeatmapAssay(nv) {
@@ -289,10 +316,11 @@ export function ScPanel(props) {
     const isClusteringDataValid = !isLoading && sNHControls !== 'none' && clusteringData && !Object.keys(clusteringData).includes('markers');
     
     return (
-      <Grid item xs={gridItemSize} className={classes.plot} style={{overflowY: 'scroll' }}>
+      // <Grid item xs={gridItemSize} className={classes.plot} style={{overflowY: 'scroll' }}>
+      <Grid item xs={gridItemSize*heatmapAssayInEffect.length} className={classes.plot} style={{overflowY: 'scroll' }}>
         {sNHControls !== "none"  && (
-          <div className={maximized?classes.floatingControls:classes.fixedControls} style={{ top: 0 }}>
-            <FormControl className={classes.parintControl} style={style3} variant="outlined">
+          <div className={maximized?classes.floatingControls:classes.fixedControls} style={style3}>
+            <FormControl className={classes.parintControl} variant="outlined">
               <PredictiveInput
                 classes={classes}
                 options={clusterings.map(c => ({ value: c, label: c }))}
@@ -310,8 +338,8 @@ export function ScPanel(props) {
             )}
     
             {isClusteringDataValid && (
-              // <ToggleButtonGroup style={style2} color="primary" value={selectedHeatmapAssay} exclusive onChange={handleToggleSelectedHeatmapAssay}>
-              <ToggleButtonGroup style={style2} color="primary" value={heatmapAssayInEffect} exclusive onChange={handleToggleSelectedHeatmapAssay}>
+              // <ToggleButtonGroup style={style2} color="primary" value={heatmapAssayInEffect} exclusive onChange={handleToggleSelectedHeatmapAssay}>
+              <ToggleButtonGroup style={style2} color="primary" value={heatmapAssayInEffect} exclusive={!maximized} onChange={handleToggleSelectedHeatmapAssay}>
               {Object.keys(clusteringData).map((heatmapAssay, index) => (
                   <ToggleButton key={`toggleSelectedHeatmapAssay-${index}`} value={heatmapAssay}>
                     {heatmapAssay}
@@ -329,26 +357,46 @@ export function ScPanel(props) {
             </Box>
           )}
     
-          {!isLoading && clusteringData && heatmapMode === "markers-table" && (
+          {/* {!isLoading && clusteringData && heatmapMode === "markers-table" && (
             <MarkersTableComponent
               clusteringData={clusteringData}
               // selectedHeatmapAssay={selectedHeatmapAssay}
               selectedHeatmapAssay={heatmapAssayInEffect}
+              maxWidth={1500/activeModes().length/(maximized?1:3)/heatmapAssayInEffect.length}
             />
-          )}
+          )} */}
     
-    {!isLoading && clusteringData && heatmapMode === "heatmap" && (
-          // {heatmapMode === "heatmap" && (
-            <HeatmapComponent
-              data={clusteringData}
-              usingFullMatrix={usingFullMatrix}
-              clustering={clustering.value}
-              layer={layer}
-              // selectedHeatmapAssay={selectedHeatmapAssay}
-              heatmapAssayInEffect={heatmapAssayInEffect}
-              sNHControls={sNHControls}
-            />
-          )}
+          {!isLoading && clusteringData && (
+                // {heatmapMode === "heatmap" && (
+                <Grid container >
+                  {(heatmapAssayInEffect?heatmapAssayInEffect:[null]).map( (hmai, idx) =>
+                    <Grid item key={`heatmapContainer_${idx}`} xs={(12/heatmapAssayInEffect.length)}>
+                    {heatmapMode === "heatmap" &&
+                          <HeatmapComponent
+                            data={clusteringData}
+                            usingFullMatrix={usingFullMatrix}
+                            clustering={clustering.value}
+                            layer={layer}
+                            // selectedHeatmapAssay={selectedHeatmapAssay}
+                            heatmapAssayInEffect={hmai}
+                            sNHControls={sNHControls}
+                            // maxWidth={Math.floor(1600/activeModes().length/(maximized?1:3)/heatmapAssayInEffect.length)}
+                            maxWidth={Math.floor(1600/12*gridItemSize/(maximized?1:3))}
+                            // maxWidth={`100%`}
+                          />}
+                        {heatmapMode === "markers-table" &&  
+                          <MarkersTableComponent
+                            clusteringData={clusteringData}
+                            // selectedHeatmapAssay={selectedHeatmapAssay}
+                            selectedHeatmapAssay={hmai}
+                            // maxWidth={Math.floor(1600/activeModes().length/(maximized?1:3)/heatmapAssayInEffect.length)}
+                            maxWidth={Math.floor(1600/12*gridItemSize/(maximized?1:3))}
+
+                          />}                
+                    </Grid>
+                  )}
+                </Grid>
+                )}
         </div>
       </Grid>
     );
@@ -360,7 +408,14 @@ export function ScPanel(props) {
 
   const isInPlotMode = config.mode.includes("plot");
   const isInLegendsMode = config.mode.includes("legends");
-  const gridItemSize = 12 / activeModes().length;
+
+  
+  const ams = activeModes();
+  const gridItemSize = ams.includes("heatmap")?
+    12/(ams.length+heatmapAssayInEffect.length-1) :
+    12/ams.length;
+
+
   const isPlotOrHeatmap = isInPlotMode || config.mode.includes("heatmap");
   const style = { margin: "5px" };
 
@@ -449,9 +504,9 @@ export function ScPanel(props) {
               layout={plot.layout}
               config={plot.config}
               useResizeHandler={plot.useResizeHandler}
-              // onInitialized={(fig, div) => setPlotDiv(div)}
-              // onUpdate={() => Plotly.newPlot(plotDiv, plot.data, plot.layout, plot.config)}
-              // onWebGlContextLost={() => Plotly.newPlot(plotDiv, plot.data, plot.layout, plot.config)}
+              onInitialized={(fig, div) => setPlotDiv(div)}
+              onUpdate={() => Plotly.newPlot(plotDiv, plot.data, plot.layout, plot.config)}
+              onWebGlContextLost={() => Plotly.newPlot(plotDiv, plot.data, plot.layout, plot.config)}
             />
           </Grid>
         )}
@@ -461,21 +516,11 @@ export function ScPanel(props) {
   
         {isInLegendsMode && (
           <Grid item xs={gridItemSize} className={classes.plot} style={{ overflowY: 'scroll' }}>
-            {sNHControls !== "none"  && (
-              <div className={maximized?classes.floatingControls:classes.fixedControls}  style={{ top: 0 }}>
-
-                <ToggleButtonGroup style={style1} color="primary" value={config.modeAssay} exclusive onChange={handleToggleModeAssay}>
-                  {data.meta.assays.concat("genesets").map(p => (
-                    <ToggleButton data-test-id="assay-type" key={`toggleModeAssay${p}`} value={p}>
-                      {p}
-                    </ToggleButton>
-                  ))}
-                </ToggleButtonGroup>
-              </div>)}
             <FeatureTable
-              // className={classes.plot}
+              className={classes.plot}
               style={{ overflow: 'scroll' }}
-              features={config.modeAssay === "genesets" ? genesets : data.meta.features[config.modeAssay]}
+              features={genesets}
+              maxWidth={Math.floor(1600/12*gridItemSize/(maximized?1:3))}
               assay={config.modeAssay}
             />
           </Grid>