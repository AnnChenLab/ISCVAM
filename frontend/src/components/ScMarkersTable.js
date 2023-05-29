import React, { useState, useEffect } from "react";
import {
  DataGrid,
  GridToolbarContainer,
  GridToolbarExport,
} from "@mui/x-data-grid";
import { Button } from "@mui/material";
import { transpose2 } from "../utils";

export default function MarkersTableComponent({ clusteringData, selectedHeatmapAssay }) {
  const { heatmapMarkers, markers } = selectedHeatmapAssay === null 
    ? clusteringData 
    : clusteringData[selectedHeatmapAssay];

  const [data, setData] = useState(null);
  const [initialRowsState, setInitialRowsState] = useState(null);

  const getUniqueMarkers = (allMarkers) => {
    const hash = Object.create(null);

    return allMarkers.reduce((accumulator, row) => {
      if (!(row.gene in hash)) {
        hash[row.gene] = accumulator.push(row) - 1;
      } else if (row.diff > accumulator[hash[row.gene]].diff) {
        accumulator[hash[row.gene]] = row;
      }
      
      return accumulator;
    }, []);
  }

  const CustomToolbar = () => {
    const getHeatmapRows = () => {
      if (!heatmapMarkers) return;
      
      const filteredMarkers = data.rows.filter((row) => heatmapMarkers.includes(row.gene));
      const uniqueMarkers = getUniqueMarkers(filteredMarkers);
      
      uniqueMarkers.sort(
        (a, b) => heatmapMarkers.indexOf(a.gene.toLowerCase()) - heatmapMarkers.indexOf(b.gene.toLowerCase())
      );
      
      setData({
        columns: data.columns,
        rows: uniqueMarkers,
      });
    };

    const getAllMarkers = () => {
      setData({
        columns: data.columns,
        rows: initialRowsState,
      });
    };

    return (
      <GridToolbarContainer>
        <GridToolbarExport printOptions={{ disableToolbarButton: true }} />
        <Button onClick={getAllMarkers}>All Markers</Button>
        <Button onClick={getHeatmapRows}>Heatmap Markers</Button>
      </GridToolbarContainer>
    );
  }

  useEffect(() => {
    const headers = Object.keys(markers).filter(m=>m!=="gene");

    const columns = headers.map((header) => ({
        field: header,
        headerName: header.toUpperCase(),
        type: header.toUpperCase() !== "CLUSTER" ? "number" : "string",
        width: header.toUpperCase() !== "CLUSTER" ? 70 : 120
      }));

    columns.unshift({
      field: "gene",
      headerName: "GENE",
      width: 250,
    });

    const columnsValues = columns.map(c => markers[c.field]);
    columnsValues.forEach(c=> delete c.name);
    const rowsValues = transpose2(columnsValues);
    
    const rows = rowsValues.map((row, index) => {
      const result = {
        id: index,
        gene: row[0],
      };
      headers.forEach((header, index) => {
        result[header] = row[index + 1];
      });

      return result;
    });

    setInitialRowsState(rows);
    setData({
      rows: rows,
      columns: columns,
    });
  }, [selectedHeatmapAssay]);

  return (
    <>
      {data !== null && heatmapMarkers !== null && (
        <DataGrid
          sx={{ width: "100%" }}
          slots={{ toolbar: CustomToolbar }}
          {...data}
          pagination={{
            paginationModel: {
              pageSize: 25,
            },
          }}
        />
      )}
    </>
  );
}
