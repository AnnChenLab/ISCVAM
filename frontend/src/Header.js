import React from 'react';
import { AppBar, Toolbar, Typography, Tabs, Tab } from '@material-ui/core';
import { Link, useLocation } from 'react-router-dom';

const Header = () => {
  const location = useLocation();
  const currentPath = location.pathname;

  // Determine which tab should be selected based on the current path
  const value = () => {
    switch (currentPath) {
      case '/iscvam/tutorial':
        return 1;
      case '/iscvam/about':
        return 2;
      case '/iscvam/contact':
        return 3;
      default:
        return 0;
    }
  };

  return (
    <AppBar 
      position="static" 
      style={{ 
        backgroundColor: 'rgba(42, 77, 143, 1)', 
        minHeight: '120px',  // Increase the height of the header
        position: 'relative',  // Ensure the header is the reference for absolute positioning
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
          backgroundImage: 'url(/image.png)',  // Add your image path here
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
          backgroundImage: 'url(/image.png)',  // Add your image path here
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
          backgroundImage: 'url(/heatmap_capture.png)',  // Add your image path here
          backgroundSize: 'contain',  // Make sure the image fits the container
          backgroundRepeat: 'no-repeat',  // Prevent the image from repeating
          opacity: 0.8,  // Set transparency for the image
          zIndex: '-1',  // Ensure the image stays behind the header content
        }}
      />
      <Toolbar style={{ justifyContent: 'center', flexDirection: 'column' }}>
        <Typography variant="h6" style={{ color: 'white', textAlign: 'center', marginBottom: '10px', marginTop: '30px' }}>
          ISCVAM: Interactive Single Cell Visual Analytics for Multiomics
        </Typography>
        <Tabs value={value()} centered style={{ width: '100%' }}>
          <Tab label="Home" component={Link} to="/iscvam" style={{ color: 'white' }} />
          <Tab label="Tutorial" component={Link} to="/iscvam/tutorial" style={{ color: 'white' }} />
          <Tab label="About" component={Link} to="/iscvam/about" style={{ color: 'white' }} />
          <Tab label="Contact" component={Link} to="/iscvam/contact" style={{ color: 'white' }} />
        </Tabs>
      </Toolbar>
    </AppBar>
  );
};

export default Header;
