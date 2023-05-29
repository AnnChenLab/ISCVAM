import React, { useState } from 'react';
import { Grid, withStyles } from '@material-ui/core';
// import { useAutoFade } from '../hooks/autoFade';
import { useDataPanels } from '../hooks/dataPanels';
import { ScPanel } from './ScPanel';

const styles = theme => ({
  root: {
      flexGrow: 1,
      height: '100%',
      width: '100%'
  },
  input: {
      display: 'flex',
      padding: 0,
  },
  title: {
      marginLeft: theme.spacing.unit,
      marginRight: theme.spacing.unit,
      marginTop: 0,
      marginBottom: 0,
      wordWrap: 'break-word'
  },
  card: {
      marginLeft: theme.spacing.unit,
      marginRight: theme.spacing.unit,
  },
  paper: {
      height: '100%',
      width: '100%',
  },
  textField: {
      marginLeft: theme.spacing.unit,
      marginRight: theme.spacing.unit,
  },
  floatingControls: {
      position: 'absolute',
      zIndex: 10
  },
  fixedControls: {
    zIndex: 10
  },
  controlArea: {},
  center: {
      position: 'absolute',
      top: '50%',
      left: '30%',
      right: '30%',
      backgroundColor: 'white',
      marginLeft: 'auto',
      marginRight: 'auto',
      padding: '20px',
      zIndex: 100
  },
  plotMain: {
      height: 'calc(100vh - 6em)',
      width: '100%'
  },
  plot: {
      height: '100%',
      width: 'calc(100% - 5px)'
  },
  toggleProjection: {
      position: 'absolute',
      width: '30%',
      bottom: '40px',
      zIndex: 200
  },
  footer: {
      position: 'absolute',
      width: '100%',
      marginBottom: '0px',
      padding: '0px',
      bottom: '0px',
      fontSize: '0.8em',
      zIndex: 10
  },
  headerPanel: {
      whiteSpace: 'nowrap',
      margin: '0px',
      minHeight: '48px',
      padding: '0px',
      align: 'left',
      top: '0px',
      fontSize: '0.8em',
      zIndex: 10
  },
  footerPanel: {
      whiteSpace: 'nowrap',
      marginBottom: '0px',
      padding: '0px',
      align: 'center',
      bottom: '0px',
      fontSize: '0.8em',
      zIndex: 10
  },
  chooseDatasetControl: {
      position: 'absolute',
      top: 1,
      marginLeft: 1,
      zIndex: 300
  },
  searchControl: {
      margin: theme.spacing.unit,
      minWidth: 100,
      maxWidth: 300
  },
  paintControl: {
      margin: theme.spacing.unit,
      minWidth: 100,
      maxWidth: 300
  },
  paintGenesetControl: {
      margin: theme.spacing.unit,
      minWidth: 100,
      maxWidth: 300
  }
});

function ScDataViewer({ classes, settings, initialDatasetIdx }) {
    const splits = 3;
    const [dataPanels, loadDataPanel] = useDataPanels({ datasets: settings.datasets, initialDatasetIdx, nPanels: splits });
    const [maximizedPanelIdx, setMaximizedPanelIdx] = useState(0);

    function toggleMax(nv) {
        if (nv == null) return;
        setMaximizedPanelIdx(nv);
    }

    return (
        <Grid container className={classes.root} spacing={0}>
            {maximizedPanelIdx === -1
                ? Array(splits).fill(0).map((element, idx) =>
                    <Grid key={`scpanel_${idx}`} item xs={12 / splits}>
                        {dataPanels[idx].data === null
                            ? 'Loading'
                            : <ScPanel
                                classes={classes}
                                settings={settings}
                                panelDataState={dataPanels[idx]}
                                panelIdx={idx}
                                loadDataPanel={loadDataPanel}
                                maximized={false}
                                toggleMaximize={() => toggleMax(idx)}
                            />}
                    </Grid>)
                : <Grid item xs={12}>
                    {dataPanels[maximizedPanelIdx].data === null
                        ? 'Loading'
                        : <ScPanel
                            settings={settings}
                            classes={classes}
                            panelDataState={dataPanels[maximizedPanelIdx]}
                            panelIdx={maximizedPanelIdx}
                            loadDataPanel={loadDataPanel}
                            maximized={true}
                            toggleMaximize={() => toggleMax(-1)}
                        />}
                </Grid>
            }
        </Grid>
    )
}

export default withStyles(styles)(ScDataViewer);