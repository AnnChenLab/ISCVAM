import React from 'react';

const Tutorial = () => {
  return (
    <div style={{ textAlign: 'center', padding: '20px' }}>
      <iframe
        src="/ISCVAM_tutorial_2024-09-24.pdf" // Path to your PDF
        width="100%" // Adjust width to fill the container
        height="1000px" // You can increase or decrease this value
        style={{
          border: 'none', // Remove the border around the iframe
          background: 'none', // Ensure no background is applied
        }}
        title="Tutorial PDF"
      >
        Your browser does not support iframes. Download the PDF to view it: 
        <a href="/ISCVAM_tutorial_2024-09-24.pdf">Download PDF</a>.
      </iframe>
    </div>
  );
};

export default Tutorial;
