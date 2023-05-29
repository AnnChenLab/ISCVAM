import React, { Fragment, useState } from 'react';
import Dialog from '@material-ui/core/Dialog';
import Divider from '@material-ui/core/Divider';
import DialogTitle from '@material-ui/core/DialogTitle';
import Paper from '@material-ui/core/Paper';
import DialogContent from '@material-ui/core/DialogContent';
import {withStyles} from '@material-ui/core/styles';
import Draggable from 'react-draggable';
import {Table, TableHead, TableRow, TableCell, TableBody} from '@material-ui/core';

function PaperComponent(props) {
  return (
    <Draggable>
      <Paper {...props} />
    </Draggable>
  );
}

const styles = theme => ({
  root: {
    width: '100%',
    backgroundColor: theme.palette.background.paper,
  },
  subtitle: {
    marginLeft: '10px',
    backgroundColor: theme.palette.background.accent
  },
  tableRow: {
    "&:hover": {
      backgroundColor: theme.palette.background.paper
    },
    cursor: 'pointer'
  },
  buildInfo: {
    float: "right",
    fontStyle: "italic",
    fontSize: "x-small",
  },
});

const DraggableDialog = ({classes, open, onClose, pristine, title, datasets, versions}) => {
  const [cursor, setCursor] = useState('grab');

  const handleClose = () => {
    setCursor('grab');
    onClose(null);
  };

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
          <h5 className={classes.subtitle}> <small>Select data</small> </h5>

          <Table>
            <TableHead>
              <TableRow>
              <TableCell>dataset</TableCell>
              <TableCell>layers</TableCell>
              <TableCell>modalities</TableCell>
              <TableCell>reference</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
            {datasets.map( ({name, layers, modalities, reference}, datasetIdx) => (
                <Fragment key={name}>
                <TableRow className={classes.tableRow} style={{cursor}} hover onClick={()=>{
                          setCursor('wait');
                          onClose(datasetIdx);
                          }}> 
                  <TableCell>{name}</TableCell>
                  <TableCell>{layers.length}</TableCell>
                  <TableCell>{modalities.length}</TableCell>
                  <TableCell>{reference !== null ? <a href={reference} target="_blank" rel="noreferrer noopener">{reference}</a> : ''}</TableCell>
                </TableRow>
                </Fragment>
                ))}
            </TableBody>
          </Table>

          <br></br>
          <br></br>
          <Divider></Divider>
          <br></br>
          <span  className={classes.buildInfo}>
            ISCVA  
            v{versions.version}, &nbsp;&nbsp; 
            revision {versions.revision}
            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<a href="https://lab.moffitt.org/chen/software/about-iscva/"><i>About</i></a>
          </span>
        </DialogContent>
      </Dialog>
    </div>
  );
}

export default withStyles(styles)(DraggableDialog);
