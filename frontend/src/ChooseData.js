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
  onDatasetSelect = () => {}  // Provide a default no-op function to avoid undefined errors
}) => {
  const [cursor, setCursor] = useState('grab');
  const [searchTerm, setSearchTerm] = useState(''); // Search term state

  const handleDatasetSelect = (datasetIndex) => {
    setCursor('wait');
    onDatasetSelect(datasetIndex); // Call the parent function or trigger the analysis
  };

  // Handle the search term input change
  const handleSearchChange = (event) => {
    setSearchTerm(event.target.value.toLowerCase()); // Store the lowercase version of the search term
  };

  const datasetsWithOrgan = datasets.map(dataset => ({
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

  // Filter datasets based on search term
  const filteredDatasets = sortedDatasets.filter(({ organ, name }) => {
    const organMatch = organ.toLowerCase().includes(searchTerm);
    const nameMatch = name.toLowerCase().includes(searchTerm);
    return organMatch || nameMatch;
  });

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
            <TableCell>organ</TableCell> 
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
                      onGroupClick(name);
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
                  <TableCell>
                    {isGroup ? '' : layers.length}
                  </TableCell>
                  <TableCell>
                    {isGroup ? '' : modalities.length}
                  </TableCell>
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
                      <TableCell>
                        {dataset.modalities.length}
                      </TableCell>
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
        ISCVAM v{versions.version}, &nbsp;&nbsp; revision{' '}
        {versions.revision}
        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
        <a href="https://chenlab.utah.edu/iscvam/about-iscva/">
          <i>About</i>
        </a>
      </span>
    </div>
  );
};

export default withStyles(styles)(ChooseData);
