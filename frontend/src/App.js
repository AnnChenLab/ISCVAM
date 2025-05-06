import React, { useState } from 'react';
import './App.css';
import { ThemeProvider, createTheme } from '@mui/material/styles';
import indigo from '@mui/material/colors/indigo';
import pink from '@mui/material/colors/pink';
import red from '@mui/material/colors/red';

import { BrowserRouter as Router, Routes, Route } from 'react-router-dom';

import Header from './Header';
import MainPage from './MainPage';
import Tutorial from './Tutorial';
import About from './About';
import Contact from './Contact';

const theme = createTheme({
  palette: {
    primary: indigo,
    secondary: pink,
    error: red,
    contrastThreshold: 3,
    tonalOffset: 0.2,
  },
  typography: {
    useNextVariants: true,
  },
});

const App = () => {
  const [showHeader, setShowHeader] = useState(true);

  const handleDatasetClick = () => {
    setShowHeader(false);
  };

  const handleShowHeader = () => {
    setShowHeader(true);
  };

  return (
    <ThemeProvider theme={theme}>
      <Router basename="/iscvam">
        {showHeader && <Header />}
        <Routes>
          <Route 
            index 
            element={
              <MainPage
                onDatasetClick={handleDatasetClick}
                onShowHeader={handleShowHeader}
              />
            }
          />
          <Route path="tutorial" element={<Tutorial />} />
          <Route path="about"    element={<About    />} />
          <Route path="contact"  element={<Contact  />} />
        </Routes>
      </Router>
    </ThemeProvider>
  );
  
};

export default App;
