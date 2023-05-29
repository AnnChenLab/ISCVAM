import React, { useState, useEffect } from 'react';
import { withStyles } from '@material-ui/core';

import ScDataViewer from './components/ScDataViewer'
import ChooseData from './ChooseData';
import {versions} from './versions';

const title = 'ISCVA';
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

const MainPage = (props) => {
  const {classes} = props;
  const [state, setState] = useState({
      settings: null,
      dragged: false,
      openChooseData: true,
      pristineChooseData: true,
      initialDatasetIdx: null
  })

  useEffect(() => {
    fetch(settingsURL)
      .then(response => response.json())
      .then(settings => setState(prevState => ({...prevState, settings})));
  }, []);

  const closeChooseData = (initialDatasetIdx) => {
    setState(prevState => ({...prevState, openChooseData: false, initialDatasetIdx, pristineChooseData: initialDatasetIdx === null}));
  }

  const openChooseData = () => {
    setState(prevState => ({...prevState, openChooseData: true}));
  }

  return (
    <div className="App">
        { state.settings !== null && (
          <ChooseData
            onClose={closeChooseData}
            datasets={state.settings.datasets}
            open={state.openChooseData}
            pristine={state.pristineChooseData}
            title={title}
            versions={versions}
          /> 
        )}
        { state.initialDatasetIdx !== null && (
          <ScDataViewer
            settings={state.settings}
            initialDatasetIdx={state.initialDatasetIdx}
          />
        )}
    </div>
  )
}

export default withStyles(styles)(MainPage);
