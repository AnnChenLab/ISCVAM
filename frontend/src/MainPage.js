import React, { useState, useEffect } from 'react';
import { withStyles } from '@material-ui/core';

import ScDataViewer from './components/ScDataViewer';
import sidemenu from './SideMenu';
import ChooseData from './ChooseData';
import { versions } from './versions';

//const title = 'ISCVAM: Interactive Single Cell Visual Analytics for Multiomics';
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
  "Blood": ["AEL", "AML", "ALL", "CLL", "PBMC", "multiome PBMC"],
  "Bone": ["MM", "GCTB", "OS"],
  "Brain": ["Glioma", "MB", "multiome human brain 3k"],
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
  "Skin": ["BCC", "MCC", "MF", "PCFCL", "SCC", "SKCM", "GSE189341-acral-sc-seuratV5", "Acral melanoma CCR 2022"],
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

const MainPage = ({ onDatasetClick, onShowHeader, classes }) => {
  //const { classes } = props;
  
  const [state, setState] = useState({
    settings: null,
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

    // Listen for browser back button events
    window.addEventListener('popstate', handleBackButton);

    return () => {
      // Cleanup event listener on component unmount
      window.removeEventListener('popstate', handleBackButton);
    };
  }, []);

  const handleBackButton = () => {
    onShowHeader(); // Show the header when user navigates back
    setState(prevState => ({
      ...prevState,
      initialDatasetIdx: null,
    }));
  };

  // Handle dataset selection and switch to the analysis view
  const handleDatasetSelect = (initialDatasetIdx) => {
    onDatasetClick();  // Call this function to hide the header
    setState(prevState => ({
      ...prevState,
      initialDatasetIdx,
    }));
    window.history.pushState(null, '', ''); // Push a new state to the history stack
  };

  // Handle returning to the dataset list
  const handleBackToList = () => {
    onShowHeader(); // Show the header when returning to the dataset list
    setState(prevState => ({
      ...prevState,
      initialDatasetIdx: null,
    }));
    window.history.back(); // Use the browser's back functionality
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
      {state.initialDatasetIdx === null ? (
        <div>
          <p style={{ fontSize: '16px', margin: '20px 240px', textAlign: 'center' }}>
            <strong>ISCVAM</strong> is a fast, interactive tool for analyzing scRNA-seq and multimodal datasets, <br></br>  
            including <strong>6,300,000+</strong> single cells across <strong>198</strong> studies (190 from TISCH).  
            It’s built in JavaScript for smooth visualization and analysis.
          </p>
          <ChooseData
            onDatasetSelect={handleDatasetSelect}
            datasets={Object.keys(state.groupedDatasets).map(prefix => ({
              displayName: `${prefix} (${state.groupedDatasets[prefix].length})`,
              organ: geneToOrgan[prefix],
              name: prefix,
              layers: [],
              modalities: [],
              reference: null,
              isGroup: true,
              datasets: state.groupedDatasets[prefix],
              expanded: state.expandedGroups[prefix] || false,
            }))}
            versions={versions}
            onGroupClick={toggleGroupExpansion}
          />
        </div>
      ) : (
        <div>
          <ScDataViewer
            settings={state.settings}
            initialDatasetIdx={state.initialDatasetIdx}
          />
        </div>
      )}
    </div>
  );
};

export default withStyles(styles)(MainPage);
