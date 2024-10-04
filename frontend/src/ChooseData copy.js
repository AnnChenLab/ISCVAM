import React, { Fragment, useState, useEffect } from 'react';
import Dialog from '@material-ui/core/Dialog';
import Divider from '@material-ui/core/Divider';
import DialogTitle from '@material-ui/core/DialogTitle';
import Paper from '@material-ui/core/Paper';
import DialogContent from '@material-ui/core/DialogContent';
import { withStyles } from '@material-ui/core/styles';
import Draggable from 'react-draggable';
import {
  Table,
  TableHead,
  TableRow,
  TableCell,
  TableBody
} from '@material-ui/core';

function PaperComponent(props) {
  return (
    <Draggable>
      <Paper {...props} />
    </Draggable>
  );
}

const styles = (theme) => ({
  root: {
    width: '100%',
    backgroundColor: theme.palette.background.paper
  },
  subtitle: {
    marginLeft: '10px',
    backgroundColor: theme.palette.background.accent
  },
  tableRow: {
    '&:hover': {
      backgroundColor: theme.palette.background.paper
    },
    cursor: 'pointer'
  },
  buildInfo: {
    float: 'right',
    fontStyle: 'italic',
    fontSize: 'x-small'
  },
  stickyHeader: {
    position: 'sticky',
    top: 0,
    backgroundColor: theme.palette.background.paper,
    zIndex: 1,
  },
  // New style for the tutorial link
  tutorialLink: {  
    display: 'block',
    marginTop: '0.001px',  // Adjust the top margin to bring it closer
    marginBottom: '8px',  // Adjust the bottom margin
    textAlign: 'left',
    textDecoration: 'none',
    color: theme.palette.primary.main,
    fontSize: 'small',  // Smaller font size
    fontWeight: 'normal'  // Normal weight
  }
});

const DraggableDialog = ({
  classes,
  open,
  onClose,
  pristine,
  title,
  datasets,
  versions,
  onGroupClick
}) => {
  const [cursor, setCursor] = useState('grab');

  const handleClose = () => {
    setCursor('grab');
    onClose(null);
  };

  // Ensure all datasets have an organ field, even if empty
  const datasetsWithOrgan = datasets.map(dataset => ({
    ...dataset,
    organ: dataset.organ || '',
  }));

  // Sort datasets with "multiome" group at the top
  const sortedDatasets = datasetsWithOrgan.sort((a, b) => {
    if (a.name === 'multiome') return -1; // Always prioritize the "multiome" group
    if (b.name === 'multiome') return 1; // Always prioritize the "multiome" group
    if (!a.organ) return 1; // Push null/empty organ to the bottom
    if (!b.organ) return -1; // Push null/empty organ to the bottom
    if (a.organ < b.organ) return -1;
    if (a.organ > b.organ) return 1;
    return 0;
  });

  return (
    <div>
      <Dialog
        open={open}
        onClose={handleClose}
        PaperComponent={PaperComponent}
        className={classes.root}
        disableBackdropClick={pristine}
        disableEscapeKeyDown={pristine}
        aria-labelledby="draggable-dialog-title"
        maxWidth="md"
      >
        <DialogTitle id="draggable-dialog-title">{title}</DialogTitle>
        <DialogContent>

          {/* Added the tutorial link above "Select data" */}
          <a href="/iscvam/tutorial/" className={classes.tutorialLink}>
            View Tutorial
          </a>

          <h5 className={classes.subtitle}>
            {' '}
            <small>Select data</small>{' '}
          </h5>

          <Table>
            <TableHead>
              <TableRow className={classes.stickyHeader}>
                {/* Switched "organ" and "dataset" columns */}
                <TableCell>organ</TableCell> 
                <TableCell>dataset</TableCell>
                <TableCell>layers</TableCell>
                <TableCell>modalities</TableCell>
                <TableCell>reference</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {sortedDatasets.map(
                (
                  {
                    displayName,
                    organ,
                    name,
                    layers,
                    modalities,
                    reference,
                    isGroup,
                    datasets: groupDatasets,
                    expanded
                  },
                  index
                ) => (
                  <Fragment key={index}>
                    <TableRow
                      className={classes.tableRow}
                      style={{ cursor }}
                      hover
                      onClick={() => {
                        if (isGroup) {
                          onGroupClick(name); // Toggle group expansion
                        } else {
                          setCursor('wait');
                          onClose(index);
                        }
                      }}
                    >
                      {/* Switched "organ" and "dataset" columns */}
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
                      <TableCell>
                        {isGroup && !reference
                          ? ''
                          : reference !== null
                          ? (
                            <a
                              href={reference}
                              target="_blank"
                              rel="noreferrer noopener"
                            >
                              {reference}
                            </a>
                            )
                          : ''}
                      </TableCell>
                    </TableRow>

                    {isGroup &&
                      expanded &&
                      groupDatasets.map((dataset, subIndex) => (
                        <TableRow
                          key={`${index}-${subIndex}`}
                          className={classes.tableRow}
                          style={{ cursor }}
                          hover
                          onClick={() => {
                            setCursor('wait');
                            onClose(dataset.originalIndex); // Pass the original dataset index to onClose
                          }}
                        >
                          <TableCell>{''}</TableCell> {/* Hide organ for individual datasets */}
                          <TableCell style={{ paddingLeft: '30px' }}>
                            {dataset.name}
                          </TableCell>
                          <TableCell>{dataset.layers.length}</TableCell>
                          <TableCell>
                            {dataset.modalities.length}
                          </TableCell>
                          <TableCell>
                            {dataset.reference !== null ? (
                              <a
                                href={dataset.reference}
                                target="_blank"
                                rel="noreferrer noopener"
                              >
                                {dataset.reference}
                              </a>
                            ) : (
                              ''
                            )}
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
            ISCVA v{versions.version}, &nbsp;&nbsp; revision{' '}
            {versions.revision}
            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
            <a href="https://chenlab.utah.edu/iscvam/about-iscva/">
              <i>About</i>
            </a>
          </span>
        </DialogContent>
      </Dialog>
    </div>
  );
};

export default withStyles(styles)(DraggableDialog);
