import React, { useState, useEffect, useMemo } from "react";
import Plot from "react-plotly.js";
import {
  ToggleButtonGroup,
  ToggleButton,
  Box,
  TextField,
} from "@mui/material";
import { transpose } from "../utils";

const HeatmapComponent = (({
  data,
  clustering,
  heatmapAssayInEffect,
  // selectedHeatmapAssay,
  usingFullMatrix,
  sNHControls,
  layer
}) => {
  // const [heatmapData, setHeatmapData] = useState({});
  // const {clustering, selectedHeatmapAssay, clusteringData: data} = config;
  const selectedHeatmapAssay = heatmapAssayInEffect;
  // console.log('rendering heatmap');
  // console.log(data);
  // console.log(clustering);
  // console.log(selectedHeatmapAssay);
  
  const [heatmapMode, setHeatmapMode] = useState(usingFullMatrix ? "allAvgsContrast" : "avgsContrast");
  const [heatmapHeight, setHeatmapHeight] = useState(1000);
  const [heatmapWidth, setHeatmapWidth] = useState(1000);
  const [heatmapMarkers, setHeatmapMarkers] = useState(null);
  const [selectedMarkersTestInput, setSelectedMarkersTestInput] = useState("");
  
  const heatmapTitle = `${layer} - ${clustering} - ${selectedHeatmapAssay} -${heatmapMode.includes("Contrast") ? "Contrast" : "LogNorm"}`;

  const toggleMode = (newValue) => {
    if(newValue && newValue !== heatmapMode){
      setHeatmapMode(newValue);
    }
  };

  const handleTextFieldChange = (event) => {
    setSelectedMarkersTestInput(event.target.value);
  };

  const handleSelectGenesButtonClick = () => {
    setHeatmapMarkers(selectedMarkersTestInput.split("\n"));
  };

  const handleResetGenesButtonClick = () => {
    setHeatmapMarkers(null);
    setSelectedMarkersTestInput("");
  };

  const normalizeData = () => {
    // console.log('normalizing data');
    // console.log(data);
    // console.log(selectedHeatmapAssay);
    // console.log(heatmapMarkers)
    const _data =
      selectedHeatmapAssay == null 
        ? data[heatmapMode]
        : data[selectedHeatmapAssay][heatmapMode]

    let _selectedMarkers =
      heatmapMarkers !== null
        ? heatmapMarkers
        : selectedHeatmapAssay !== null
        ? data[selectedHeatmapAssay]["heatmapMarkers"]
        : data["heatmapMarkers"];

    _selectedMarkers = _selectedMarkers.map((marker) => marker.toLowerCase());

    const _clusterNames =
      selectedHeatmapAssay !== null
        ? data[selectedHeatmapAssay]["clusterNames"]
        : data["clusterNames"];

    const selectedMarkersSet = new Set(_selectedMarkers);
    const selectedMarkersIndex = Object.fromEntries(_selectedMarkers.map((marker, index) => [marker, index]));

    const labels = _data.find((element) => Array.isArray(element));

    const seenLabels = new Set();
    const sorteredColumns = labels
      .reduce((result, item, index) => {
        if (!seenLabels.has(item) && selectedMarkersSet.has(item.toLowerCase())) {
          result.push({ item, index });
          seenLabels.add(item);
        }
        return result;
      }, [])
      .sort((a, b) => {
        return (
          selectedMarkersIndex[a.item.toLowerCase()] -
          selectedMarkersIndex[b.item.toLowerCase()]
        );
      });

    const heatmapLabels = sorteredColumns.map((column) => column.item);

    const filteredColumnsMatrix = _data.slice(1).map( (c,i) => 
      sorteredColumns.map((element) => c[element.index]));

    const _x = _clusterNames;
    const _y = heatmapLabels;
    const _z = transpose(filteredColumnsMatrix);

    const height =
      _z.length < 30
        ? _z.length < 3
          ? _z.length * 100 + 300
          : _z.length * 50 + 300
        : _z.length * 17;
    const width =
        (_x.length < 20
          ? _x.length < 3
            ? _x.length * 100
            : _x.length * 36
          : _x.length * 20)+300;

    setHeatmapHeight(height);
    setHeatmapWidth(width);
    // console.log('done normalizing data');
    // console.log({_x,_y,_z});

    return { x: _x, y: _y, z: _z };
  };
  
  const heatmapData = useMemo(()=> {
    //     console.log('in use memo setting heatmap data');
    // console.log(heatmapMode);
    // console.log(data);
    // console.log(selectedHeatmapAssay);
    // console.log(heatmapMarkers);

    if (data && heatmapMode && ( (selectedHeatmapAssay === null  && Object.keys(data).includes("markers"))  || 
                               (selectedHeatmapAssay !== null && Object.keys(data).includes(selectedHeatmapAssay)))) {
      return normalizeData();
    }
  }, [heatmapMode, selectedHeatmapAssay, heatmapMarkers]);



  // console.log(data);
  // console.log(heatmapData);
  return (
    <>
      {data && heatmapData && (
        <Box style={{ width: "100%" }}>
          {sNHControls !== 'none' && (
            <>
              <ToggleButtonGroup
                style={{ display: "block" }}
                color="primary"
                onChange={(e, nv) => toggleMode(nv)}
                value={heatmapMode}
                exclusive
              >
                <ToggleButton data-test-id="avgsContrast" key="avgsContrast-heatmap" value={usingFullMatrix ? "allAvgsContrast" : "avgsContrast"}>
                  CONTRAST
                </ToggleButton>
                <ToggleButton data-test-id="avgsLogNorm" key="toggleavgsLogNorm-heatmap" value={!usingFullMatrix ? "avgsLogNorm" : "allAvgsLogNorm"}>
                  LOGNORM
                </ToggleButton>
              </ToggleButtonGroup>

              <Box>
                <TextField
                  sx={{ margin: "20px", width: "400px" }}
                  placeholder="Enter gene labels..."
                  multiline
                  minRows={4}
                  maxRows={10}
                  onChange={handleTextFieldChange}
                  value={selectedMarkersTestInput}
                />
                <br />
                <ToggleButtonGroup style={{ display: "block" }} color="primary">
                  <ToggleButton data-test-id="selectGenesButton" onClick={handleSelectGenesButtonClick} value="A">
                    SELECT GENES
                  </ToggleButton>
                  <ToggleButton data-test-id="resetGenesButton" onClick={handleResetGenesButtonClick} value="B">
                    RESET
                  </ToggleButton>
                </ToggleButtonGroup>
              </Box>
            </>
          )}

          <Box sx={{ margin: "10px" }}>{heatmapTitle}</Box>
          <Box sx={{ minHeight: "1000px", width: "100%" }}>
            <Plot
              data={[
                {
                  x: heatmapData.x,
                  y: heatmapData.y,
                  z: heatmapData.z,
                  type: "heatmap",
                  colorscale: "Viridis",
                },
              ]}
              config={{
                displaylogo: false,
                modeBarButtonsToRemove: [
                  'sendDataToCloud',
                  'hoverCompareCartesian',
                  'zoom2d',
                  'pan2d',
                  'select2d',
                  'lasso2d',
                  'zoomIn2d',
                  'zoomOut2d',
                  'autoScale2d',
                  'resetScale2d', 
                  'hoverClosestCartesian',
                  'toggleSpikelines'
                ]
              }}
              layout={{
                height: heatmapHeight,
                width: heatmapWidth,
                margin: "auto",
                xaxis: {
                  side: "top",
                  tickangle: 35,
                  showticklabels: true,
                  automargin: true,
                },
                yaxis: {
                  autotick: true,
                  autorange: "reversed",
                  automargin: true,
                  autoshift: true,
                },
              }}
            />
          </Box>
        </Box>
      )}
    </>
  );
})


export default HeatmapComponent;