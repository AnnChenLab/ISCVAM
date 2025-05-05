import React, { Fragment, useState } from 'react';
import Divider from '@material-ui/core/Divider';
import DialogTitle from '@material-ui/core/DialogTitle';
import { withStyles } from '@material-ui/core/styles';
import {
  Table,
  TableHead,
  TableRow,
  TableCell,
  TableBody,
  TextField // Import TextField for the search input
} from '@material-ui/core';

const styles = (theme) => ({
  root: {
    width: '70%', // Reduce the width of the entire content
    backgroundColor: theme.palette.background.paper,
    paddingLeft: '50px',
    margin: '0 auto', // Center the root content
  },
  searchInput: {
    marginBottom: '16px', // Add margin below the search input
    width: '100%',
  },
  tableRow: {
    '&:hover': {
      backgroundColor: theme.palette.background.paper,
    },
    cursor: 'pointer',
  },
  stickyHeader: {
    position: 'sticky',
    top: 0,
    backgroundColor: theme.palette.background.paper,
    zIndex: 1,
  },
  table: {
    width: '70%', // Shrink the table width further to 70%
    margin: '0 auto', // Center the table
  },
  tableCell: {
    padding: '4px 8px', // Reduce the padding to make the content more compact
    wordWrap: 'break-word', // Ensure text wraps inside the cells
  },
});

const ChooseData = ({
  classes,
  title,
  datasets,
  versions,
  onGroupClick,
  onDatasetSelect = () => {} // Provide a default no-op function to avoid undefined errors
}) => {
  const [cursor, setCursor] = useState('grab');
  const [searchTerm, setSearchTerm] = useState(''); // Search term state
  const [expandedGroups, setExpandedGroups] = useState({}); // State to track expanded groups

  const handleDatasetSelect = (datasetIndex) => {
    setCursor('wait');
    onDatasetSelect(datasetIndex); // Call the parent function or trigger the analysis
  };

  // Handle the search term input change
  const handleSearchChange = (event) => {
    setSearchTerm(event.target.value.toLowerCase()); // Store the lowercase version of the search term
  };

  // Toggle the expand/collapse state of a group
  const handleGroupClick = (groupName) => {
    setExpandedGroups((prevState) => ({
      ...prevState,
      [groupName]: !prevState[groupName], // Toggle expand/collapse for the clicked group
    }));
  };

  const datasetsWithOrgan = datasets.map((dataset) => ({
    ...dataset,
    organ: dataset.organ || '',
  }));

  const sortedDatasets = datasetsWithOrgan.sort((a, b) => {
    if (a.name === 'multiome') return -1;
    if (b.name === 'multiome') return 1;
    if (!a.organ) return 1;
    if (!b.organ) return -1;
    if (a.organ < b.organ) return -1;
    if (a.organ > b.organ) return 1;
    return 0;
  });

  const filteredDatasets = sortedDatasets.reduce((result, { organ, displayName, isGroup, datasets: groupDatasets }) => {
    const searchTermRegex = new RegExp(searchTerm, 'i'); // Case-insensitive search term
    const organMatch = searchTermRegex.test(organ); // Match for organ name
    const groupNameMatch = searchTermRegex.test(displayName); // Match for group name
  
    if (isGroup) {
      if (searchTerm) {
        // If there's a search term, filter datasets within the group based on individual dataset names
        const matchingDatasets = groupDatasets.filter((dataset) =>
          searchTermRegex.test(dataset.name) // Match individual dataset names
        );
  
        // If there's a match in the organ name, group name, or datasets, include this group
        if (organMatch || groupNameMatch || matchingDatasets.length > 0) {
          result.push({
            organ,
            displayName, // Group name
            isGroup,
            datasets: matchingDatasets.length > 0 ? matchingDatasets : groupDatasets, // Show matching datasets or all if organ/group matches
            expanded: true, // Auto-expand group when filtering
          });
        }
      } else {
        // No search term, allow the group to be expanded/collapsed based on the state
        result.push({
          organ,
          displayName, // Group name
          isGroup,
          datasets: groupDatasets, // Show all datasets when no search term
          expanded: expandedGroups[displayName] || false, // Keep current expand/collapse state from the state
        });
      }
    } else if (organMatch || groupNameMatch) {
      // If there's a match on the organ or group name, show the entire group
      result.push({
        organ,
        displayName,
        isGroup,
        datasets: groupDatasets,
        expanded: expandedGroups[displayName] || false,
      });
    }
  
    return result;
  }, []);
  
  

  return (
    <div className={classes.root}>
      <DialogTitle id="draggable-dialog-title">{title}</DialogTitle>

      {/* Add Search Input */}
      <TextField
        className={classes.searchInput}
        label="Search by organ or dataset name"
        variant="outlined"
        value={searchTerm}
        onChange={handleSearchChange}
      />

      <Table className={classes.table}>
        <TableHead>
          <TableRow className={classes.stickyHeader}>
            <TableCell>organ/site</TableCell>
            <TableCell>dataset</TableCell>
            <TableCell>layers</TableCell>
            <TableCell>modalities</TableCell>
          </TableRow>
        </TableHead>
        <TableBody>
          {filteredDatasets.map(
            (
              {
                displayName,
                organ,
                name,
                layers,
                modalities,
                isGroup,
                datasets: groupDatasets,
                expanded
              },
              index
            ) => (
              <Fragment key={index}>
                <TableRow
                  className={classes.tableRow}
                  hover
                  onClick={() => {
                    if (isGroup) {
                      handleGroupClick(displayName); // Toggle expand/collapse state on group click
                    } else {
                      handleDatasetSelect(index);
                    }
                  }}
                >
                  <TableCell>{isGroup ? organ : ''}</TableCell>
                  <TableCell>
                    {isGroup
                      ? `${expanded ? '-' : '+'} ${displayName}`
                      : name}
                  </TableCell>
                  <TableCell>{isGroup ? '' : layers.length}</TableCell>
                  <TableCell>{isGroup ? '' : modalities.length}</TableCell>
                </TableRow>

                {isGroup &&
                  expanded &&
                  groupDatasets.map((dataset, subIndex) => (
                    <TableRow
                      key={`${index}-${subIndex}`}
                      className={classes.tableRow}
                      hover
                      onClick={() => handleDatasetSelect(dataset.originalIndex)}
                    >
                      <TableCell>{''}</TableCell> {/* Hide organ for individual datasets */}
                      <TableCell style={{ paddingLeft: '30px' }}>
                        {dataset.name}
                      </TableCell>
                      <TableCell>{dataset.layers.length}</TableCell>
                      <TableCell>{dataset.modalities.length}</TableCell>
                    </TableRow>
                  ))}
              </Fragment>
            )
          )}
        </TableBody>
      </Table>

      <br />
      <br />
      <Divider />
      <br />
      <span className={classes.buildInfo}>
        ISCVAM v{versions.version}, &nbsp;&nbsp; revision {versions.revision}
        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
        <a href="https://medicine.utah.edu/internal-medicine/epidemiology/chen-lab">
          <i>Created by Chen Lab</i>
        </a>
      </span>
    </div>
  );
};

export default withStyles(styles)(ChooseData);
