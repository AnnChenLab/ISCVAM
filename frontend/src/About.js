import React from 'react';
import { Container, Typography, Link, Box, List, ListItem, ListItemText } from '@mui/material';

const About = () => {
  return (
    <Container maxWidth="md">
      <Box my={4}>
        <Typography variant="body1" paragraph>
          <em>
            ISCVAM: <b>I</b>nteractive <b>S</b>ingle <b>C</b>ell <b>V</b>isual <b>A</b>nalytics for <b>M</b>ultiomics
          </em>, is an interactive tool for analyzing and visualizing single cell multimodal data with multi-resolutions. It is deployed at{' '}
          <Link href="https://iscvam.moffitt.org/">https://iscvam.moffitt.org/</Link>, and <Link href="https://chenlab.utah.edu/iscvam/">https://chenlab.utah.edu/iscvam/</Link>.
          More than 6.3+ million single cells from ~200 cancer studies with are processes available for investigation at <Link href="https://chenlab.utah.edu/iscvam/">https://chenlab.utah.edu/iscvam/</Link>.
        </Typography>

        <Typography variant="body1" paragraph>
        <em>ISCVAM</em> is expanded from <em>ISCVA</em>, previously developed for analyzing scRNA-seq data. The information was described in this paper.
        </Typography>

        <Typography variant="body1" paragraph>
          Smalley J, Chen Z, Phadke M, Li J, Yu X, Wyatt C, Evernden B, Messina JL, Sarnaik AA, Sondak VK, Zhang C, Law V,
          Tran N, Etame A, Macaulay RJB, Eroglu Z, Forsyth PA, Rodriguez PC, Chen YA, Smalley KSM.{' '}
          <Link href="https://pubmed.ncbi.nlm.nih.gov/34035069/">
            Single-Cell Characterization of the Immune Microenvironment of Melanoma Brain and Leptomeningeal Metastases.
          </Link>{' '}
          Clin Cancer Res. 2021 Jul 15;27(14):4109-4125.
        </Typography>

        <Typography variant="h5" component="h2" gutterBottom>
          Major Components
        </Typography>

        <List>
          <ListItem>
            <ListItemText
              primary="1. A collection of Bash and R scripts using widely used algorithms in the single-cell community:"
              secondary={
                <ul>
                  <li>
                    <Link href="https://pubmed.ncbi.nlm.nih.gov/29608179/">Seurat</Link> for general processing
                  </li>
                  <li>
                    <Link href="https://pubmed.ncbi.nlm.nih.gov/30643263/">SingleR</Link> for cell-type recognition
                  </li>
                  <li>
                    <Link href="https://pubmed.ncbi.nlm.nih.gov/31294801/">Single-cell signature explorer</Link> for gene set
                    signature scoring
                  </li>
                </ul>
              }
            />
          </ListItem>
          <ListItem>
            <ListItemText
              primary="2. A web-based component using modern technologies like React.js, TensorFlow.js, and Plotly.js."
              secondary="This component allows real-time interactive exploration, backed by an HPC Node.js backend.(A note on web-browser: ISCVAM is currently available on all web browsers, except for Safari)"
            />
          </ListItem>
        </List>

        <Typography variant="h5" component="h2" gutterBottom>
          Publications with data analyzed using ISCVA/ISCVAM:
        </Typography>

        <List>
          <ListItem>
            <ListItemText
              primary={
                <Link href="https://pubmed.ncbi.nlm.nih.gov/38056899/">
                  Differential requirements for CD4+ T cells in the efficacy of the anti-PD-1+LAG-3 and anti-PD-1+CTLA-4 combinations in melanoma flank and brain metastasis models.
                </Link>
              }
            />
          </ListItem>
          <ListItem>
            <ListItemText
              primary={
                <Link href="https://pubmed.ncbi.nlm.nih.gov/33653716/">
                  Targeted Therapy Given after Anti-PD-1 Leads to Prolonged Responses in Mouse Melanoma Models through Sustained Antitumor Immunity.
                </Link>
              }
              
            />
          </ListItem>
          <ListItem>
            <ListItemText
              primary={
                <Link href="https://pubmed.ncbi.nlm.nih.gov/34035069/">
                  Single-Cell Characterization of the Immune Microenvironment of Melanoma Brain and Leptomeningeal Metastases.
                </Link>
              }
              
            />
          </ListItem>
          <ListItem>
            <ListItemText
              primary={
                <Link href="https://pubmed.ncbi.nlm.nih.gov/35247927/">
                  Single cell characterization of the cellular landscape of acral melanoma identifies novel targets for immunotherapy.
                </Link>
              }
              
            />
          </ListItem>
          <ListItem>
            <ListItemText
              primary={
                <Link href="https://pubmed.ncbi.nlm.nih.gov/35213727/">
                  A preclinical model of patient-derived cerebrospinal fluid circulating tumor cells for experimental therapeutics in leptomeningeal disease from melanoma.
                </Link>
              }
              
            />
          </ListItem>
        </List>

        {/* Removed Contact Us section */}

        <Box mt={2} mb={2}>
          <Typography variant="body2" color="textSecondary">
            <Link
              href="https://medicine.utah.edu/internal-medicine/epidemiology/chen-lab"
              target="_blank"
              rel="noopener noreferrer"
              color="inherit"
              underline="hover"
            >
              Created by Chen Lab
            </Link>
          </Typography>
        </Box>
      </Box>
    </Container>
  );
};

export default About;
