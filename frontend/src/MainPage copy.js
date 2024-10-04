import React, { useState, useEffect } from 'react';
import { withStyles } from '@material-ui/core';

import ScDataViewer from './components/ScDataViewer';
import ChooseData from './ChooseData';
import { versions } from './versions';

const title = 'ISCVAM: Interactive Single Cell Visual Analytics for Multiomics';
// // Add the hyperlink directly after the title constant
// const tutorialLink = (
//   <a href="/tutorial" style={{ display: 'block', marginTop: '10px', textAlign: 'center', textDecoration: 'none', color: 'blue', fontWeight: 'bold' }}>
//     View Tutorial
//   </a>
// );
const settingsURL = 'app-settings.json';

const styles = theme => ({
  fab: {
    margin: theme.spacing.unit,
    position: 'absolute',
    right: '10px',
    bottom: '10px',
    zIndex: 100
  },
  fab2: {
    margin: theme.spacing.unit,
    position: 'absolute',
    right: '100px',
    bottom: '10px',
    zIndex: 200
  }
});

const groupDatasets = (datasets) => {
  const grouped = {};

  datasets.forEach((dataset, idx) => {
    const name = dataset.name || ''; // Access the 'name' property
    const prefix = name.split('_')[0]; // Group by the prefix before the first underscore

    if (!grouped[prefix]) {
      grouped[prefix] = [];
    }
    grouped[prefix].push({ ...dataset, originalIndex: idx }); // Add the entire dataset object along with the original index
  });

  return grouped;
};

// Organ mapping
const organMapping = {
  "Bladder": ["BLCA"],
  "Blood": ["AEL", "AML", "ALL", "CLL", "PBMC"],
  "Bone": ["MM", "GCTB", "OS"],
  "Brain": ["Glioma", "MB"],
  "Breast": ["BRCA"],
  "Colorectum": ["CRC"],
  "Esophagus": ["ESCA"],
  "Eye": ["RB", "UVM"],
  "Head & Neck": ["HNSC", "THCA", "LSCC", "NPC", "OSCC"],
  "Kidney": ["KIRC", "KICH", "KIPAN"],
  "Liver": ["CHOL", "LIHC", "HB"],
  "Lung": ["NSCLC", "SCLC"],
  "Lymph node": ["NHL", "DLBC"],
  "Nervous system": ["MPNST", "NET", "NB", "Neurofibroma"],
  "Pancreas": ["PAAD"],
  "Pelvic cavity": ["CESC", "OV", "UCEC", "multiome discovery", "multiome int.validation"],
  "Prostate": ["PRAD"],
  "Skin": ["BCC", "MCC", "MF", "PCFCL", "SCC", "SKCM", "multiome PBMC", "GSE189341-acral-sc-seuratV5", "Acral melanoma CCR 2022", ], 
  "Soft tissue": ["GIST", "SARC", "PPB", "SS"],
  "Stomach": ["STAD"]
};

// Reverse mapping to get the organ by gene
const geneToOrgan = Object.entries(organMapping).reduce((acc, [organ, genes]) => {
  genes.forEach(gene => {
    acc[gene] = organ;
  });
  return acc;
}, {});

const MainPage = (props) => {
  const { classes } = props;
  const [state, setState] = useState({
    settings: null,
    dragged: false,
    openChooseData: true,
    pristineChooseData: true,
    initialDatasetIdx: null,
    groupedDatasets: {},
    expandedGroups: {}, // To manage which groups are expanded
  });

  useEffect(() => {
    fetch(settingsURL)
      .then(response => response.json())
      .then(settings => {
        const groupedDatasets = groupDatasets(settings.datasets);
        setState(prevState => ({ ...prevState, settings, groupedDatasets }));
      });
  }, []);

  const closeChooseData = (initialDatasetIdx) => {
    setState(prevState => ({
      ...prevState,
      openChooseData: false,
      initialDatasetIdx,
      pristineChooseData: initialDatasetIdx === null
    }));
  };

  const toggleGroupExpansion = (prefix) => {
    setState(prevState => ({
      ...prevState,
      expandedGroups: {
        ...prevState.expandedGroups,
        [prefix]: !prevState.expandedGroups[prefix],
      }
    }));
  };

  return (
    <div className="App">
      {state.settings !== null && (
        <ChooseData
          onClose={closeChooseData}
          datasets={Object.keys(state.groupedDatasets).map(prefix => ({
            displayName: `${prefix} (${state.groupedDatasets[prefix].length})`, // Display name with gene and count
            organ: geneToOrgan[prefix], // Organ information
            name: prefix, // Original prefix for toggling
            layers: [],
            modalities: [],
            reference: null,
            isGroup: true, // To differentiate group headers from dataset entries
            datasets: state.groupedDatasets[prefix], // Store the datasets under this group
            expanded: state.expandedGroups[prefix] || false, // Whether this group is expanded
          }))}
          open={state.openChooseData}
          pristine={state.pristineChooseData}
          title={title}
          versions={versions}
          onGroupClick={toggleGroupExpansion} // Handle group selection
        />
      )}
      {state.initialDatasetIdx !== null && (
        <ScDataViewer
          settings={state.settings}
          initialDatasetIdx={state.initialDatasetIdx}
        />
      )}
    </div>
  );
};

export default withStyles(styles)(MainPage);
