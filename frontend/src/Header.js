import React from 'react';
import { AppBar, Toolbar, Typography, Tabs, Tab } from '@material-ui/core';
import { Link, useLocation } from 'react-router-dom';

const PUBLIC = process.env.PUBLIC_URL;

const Header = () => {
  const { pathname } = useLocation();

  // Strip off any trailing slash for matching
  const cleanPath = pathname.endsWith('/') && pathname !== '/' 
    ? pathname.slice(0, -1) 
    : pathname;

  // Map the current (basename‑stripped) pathname to a tab index
  const tabIndex = (() => {
    switch (cleanPath) {
      case '/':
        return 0;
      case '/tutorial':
        return 1;
      case '/about':
        return 2;
      case '/contact':
        return 3;
      default:
        return 0;
    }
  })();

  return (
    <AppBar 
      position="static" 
      style={{ 
        backgroundColor: 'rgba(42, 77, 143, 1)', 
        minHeight: '120px',  
        position: 'relative',  
      }}
    >
      {/* Background image added using ::before */}
      
      <div 
        style={{
          content: '""',  // Pseudo-element-like behavior
          position: 'absolute',
          top: '0',
          right:'25vw',  // Adjust this value to change the horizontal position
          width: '150px',  // Set the width of the image container
          height: '120px',  // Set the height of the image container
          backgroundImage: `url(${PUBLIC}/image.png)`,  // Add your image path here
          backgroundSize: 'contain',  // Make sure the image fits the container
          backgroundRepeat: 'no-repeat',  // Prevent the image from repeating
          opacity: 0.7,  // Set transparency for the image
          transform: 'rotate(-10deg)',  // Rotate the image by 15 degrees
          zIndex: '-1',  // Ensure the image stays behind the header content
        }}
      />
      <div 
        style={{
          content: '""',  // Pseudo-element-like behavior
          position: 'absolute',
          top: '0',
          right:'23vw',  // Adjust this value to change the horizontal position
          width: '150px',  // Set the width of the image container
          height: '120px',  // Set the height of the image container
          backgroundImage: `url(${PUBLIC}/image.png)`,  // Add your image path here
          backgroundSize: 'contain',  // Make sure the image fits the container
          backgroundRepeat: 'no-repeat',  // Prevent the image from repeating
          opacity: 0.6,  // Set transparency for the image
          transform: 'rotate(30deg)',  // Rotate the image by 15 degrees
          zIndex: '-1',  // Ensure the image stays behind the header content
        }}
      />
      <div 
        style={{
          content: '""',  // Pseudo-element-like behavior
          position: 'absolute',
          top: '0',
          left:'24.5vw',  // Adjust this value to change the horizontal position
          width: '120px',  // Set the width of the image container
          height: '120px',  // Set the height of the image container
          backgroundImage: `url(${PUBLIC}/heatmap_capture.png)`,  // Add your image path here
          backgroundSize: 'contain',  // Make sure the image fits the container
          backgroundRepeat: 'no-repeat',  // Prevent the image from repeating
          opacity: 0.8,  // Set transparency for the image
          zIndex: '-1',  // Ensure the image stays behind the header content
        }}
      />

      {/* (your pseudo-element divs here) */}
      <Toolbar style={{ justifyContent: 'center', flexDirection: 'column' }}>
        <Typography 
          variant="h6" 
          style={{ color: 'white', textAlign: 'center', marginBottom: 10, marginTop: 30 }}
        >
          ISCVAM: Interactive Single Cell Visual Analytics for Multiomics
        </Typography>
        <Tabs value={tabIndex} centered style={{ width: '100%' }}>
          <Tab 
            label="Home" 
            component={Link} 
            to="/" 
            style={{ color: 'white' }} 
          />
          <Tab 
            label="Tutorial" 
            component={Link} 
            to="/tutorial" 
            style={{ color: 'white' }} 
          />
          <Tab 
            label="About" 
            component={Link} 
            to="/about" 
            style={{ color: 'white' }} 
          />
          <Tab 
            label="Contact" 
            component={Link} 
            to="/contact" 
            style={{ color: 'white' }} 
          />
        </Tabs>
      </Toolbar>
    </AppBar>
  );
};

export default Header;
