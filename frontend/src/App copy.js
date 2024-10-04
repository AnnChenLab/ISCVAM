import React from 'react';
import './App.css';

import { MuiThemeProvider, createMuiTheme } from '@material-ui/core/styles';
import indigo from '@material-ui/core/colors/indigo';
import pink from '@material-ui/core/colors/pink';
import red from '@material-ui/core/colors/red';

import MainPage from './MainPage';

const theme = createMuiTheme({
  palette: {
    primary: indigo,
    secondary: pink,
    error: red,
    contrastThreshold: 3,
    tonalOffset: 0.2,
  },
  typography: {
    useNextVariants: true
  }
});

const App = () => {
  return (
    <MuiThemeProvider theme={theme}>
      <MainPage/>
    </MuiThemeProvider>
  );
}

export default App;