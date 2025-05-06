import React, { useState } from 'react';
import { Checkbox, FormGroup, FormControlLabel, FormLabel, FormControl, Slider } from '@mui/material';

const Sidebar = ({ organs, onFilterChange }) => {
  // State to hold the filter values
  const [selectedOrgans, setSelectedOrgans] = useState([]);
  const [layerRange, setLayerRange] = useState([1, 3]);
  const [selectedModalities, setSelectedModalities] = useState([1]);

  // Handle organ checkbox change
  const handleOrganChange = (organ) => {
    const newSelectedOrgans = selectedOrgans.includes(organ)
      ? selectedOrgans.filter((o) => o !== organ)
      : [...selectedOrgans, organ];
    setSelectedOrgans(newSelectedOrgans);
    onFilterChange({ organs: newSelectedOrgans, layers: layerRange, modalities: selectedModalities });
  };

  // Handle layer range change (Slider)
  const handleLayerChange = (event, newValue) => {
    setLayerRange(newValue);
    onFilterChange({ organs: selectedOrgans, layers: newValue, modalities: selectedModalities });
  };

  // Handle modalities checkbox change
  const handleModalitiesChange = (modality) => {
    const newSelectedModalities = selectedModalities.includes(modality)
      ? selectedModalities.filter((m) => m !== modality)
      : [...selectedModalities, modality];
    setSelectedModalities(newSelectedModalities);
    onFilterChange({ organs: selectedOrgans, layers: layerRange, modalities: newSelectedModalities });
  };

  return (
    <div style={{ padding: '20px', width: '200px', borderRight: '1px solid lightgray' }}>
      <FormControl component="fieldset">
        <FormLabel component="legend">Filter by Organ</FormLabel>
        <FormGroup>
          {organs.map((organ) => (
            <FormControlLabel
              key={organ}
              control={
                <Checkbox
                  checked={selectedOrgans.includes(organ)}
                  onChange={() => handleOrganChange(organ)}
                  name={organ}
                />
              }
              label={organ}
            />
          ))}
        </FormGroup>
      </FormControl>

      <div style={{ marginTop: '20px' }}>
        <FormLabel component="legend">Filter by Layers</FormLabel>
        <Slider
          value={layerRange}
          onChange={handleLayerChange}
          valueLabelDisplay="auto"
          min={1}
          max={3}
          step={1}
        />
      </div>

      <FormControl component="fieldset" style={{ marginTop: '20px' }}>
        <FormLabel component="legend">Filter by Modalities</FormLabel>
        <FormGroup>
          {[1, 2, 3, 4].map((modality) => (
            <FormControlLabel
              key={modality}
              control={
                <Checkbox
                  checked={selectedModalities.includes(modality)}
                  onChange={() => handleModalitiesChange(modality)}
                  name={`modality-${modality}`}
                />
              }
              label={`Modality ${modality}`}
            />
          ))}
        </FormGroup>
      </FormControl>
    </div>
  );
};

export default Sidebar;
