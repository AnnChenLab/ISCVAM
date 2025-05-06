import React from 'react';

const Tutorial = () => {
  // Build‐time PUBLIC_URL is "/iscvam" so this resolves to "/iscvam/ISCVAM_tutorial_2024-09-24.pdf"
  const pdfUrl = `${process.env.PUBLIC_URL}/ISCVAM_tutorial_2024-09-24.pdf`;

  return (
    <div style={{ textAlign: 'center', padding: '20px' }}>
      <iframe
        src={pdfUrl}
        width="100%"
        height="1000px"
        style={{
          border: 'none',
          background: 'none',
        }}
        title="Tutorial PDF"
      >
        Your browser doesn’t support iframes.{" "}
        <a href={pdfUrl} target="_blank" rel="noopener noreferrer">
          Download the PDF here
        </a>.
      </iframe>
    </div>
  );
};

export default Tutorial;
