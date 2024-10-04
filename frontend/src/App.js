import React, { useState } from 'react';
import './App.css';
import { ThemeProvider, createTheme } from '@mui/material/styles';
import indigo from '@mui/material/colors/indigo';
import pink from '@mui/material/colors/pink';
import red from '@mui/material/colors/red';

import { BrowserRouter as Router, Route, Routes } from 'react-router-dom';

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
    setShowHeader(false); // Hide the header when a dataset is clicked
  };

  const handleShowHeader = () => {
    setShowHeader(true); // Show the header when navigating back
  };

  return (
    <ThemeProvider theme={theme}>
      <Router>
        {showHeader && <Header />} {/* Conditionally render the header */}
        <Routes>
          <Route
            path="/"
            element={<MainPage onDatasetClick={handleDatasetClick} onShowHeader={handleShowHeader} />}
          />
          <Route
            path="/iscvam"
            element={<MainPage onDatasetClick={handleDatasetClick} onShowHeader={handleShowHeader} />}
          />
          <Route path="/iscvam/tutorial" element={<Tutorial />} />
          <Route path="/iscvam/about" element={<About />} />
          <Route path="/iscvam/contact" element={<Contact />} />
        </Routes>
      </Router>
    </ThemeProvider>
  );
};

export default App;
