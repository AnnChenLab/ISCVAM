import React from 'react';
import { Container, Typography, Link, Box, List, ListItem, ListItemText } from '@mui/material';

const About = () => {
  return (
    <Container maxWidth="md">
      <Box my={4}>
        <Typography variant="body1" paragraph>
          <em>
            ISCVAM: <b>I</b>nteractive <b>S</b>ingle <b>C</b>ell <b>V</b>isual <b>A</b>nalytics for <b>M</b>ultiomics
          </em>, is a JavaScript-based tool for scRNA-seq data analysis. Originally developed at 
          <Link href="https://iscvam.moffitt.org/"> https://iscvam.moffitt.org/</Link>, ISCVAM has since been enhanced with 
          extended capabilities, including the 190 <Link href="http://tisch.comp-genomics.org">TISCH datasets</Link>, and is 
          now hosted at <Link href="https://chenlab.utah.edu/iscvam/"> https://chenlab.utah.edu/iscvam/</Link>.
        </Typography>

        <Typography variant="body1" paragraph>
          The computational tool <em>ISCVAM</em> is described in this paper:
        </Typography>

        <Typography variant="body1" paragraph>
          Smalley J, Chen Z, Phadke M, Li J, Yu X, Wyatt C, Evernden B, Messina JL, Sarnaik AA, Sondak VK, Zhang C, Law V, Tran N, Etame A, 
          Macaulay RJB, Eroglu Z, Forsyth PA, Rodriguez PC, Chen YA, Smalley KSM. 
          <Link href="https://pubmed.ncbi.nlm.nih.gov/34035069/">
            Single-Cell Characterization of the Immune Microenvironment of Melanoma Brain and Leptomeningeal Metastases.
          </Link> Clin Cancer Res. 2021 Jul 15;27(14):4109-4125.
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
                  <li><Link href="https://pubmed.ncbi.nlm.nih.gov/29608179/">Seurat</Link> for general processing</li>
                  <li><Link href="https://pubmed.ncbi.nlm.nih.gov/30643263/">SingleR</Link> for cell-type recognition</li>
                  <li><Link href="https://pubmed.ncbi.nlm.nih.gov/31294801/">Single-cell signature explorer</Link> for gene set signature scoring</li>
                </ul>
              }
            />
          </ListItem>
          <ListItem>
            <ListItemText
              primary="2. A web-based component using modern technologies like React.js, TensorFlow.js, and Plotly.js."
              secondary="This component allows real-time interactive exploration, backed by an HPC Node.js backend. (Currently available on all web browsers except Safari)"
            />
          </ListItem>
        </List>

        <Typography variant="h5" component="h2" gutterBottom>
          Publications Using ISCVAM
        </Typography>

        <List>
          <ListItem>
            <ListItemText
              primary={
                <Link href="https://pubmed.ncbi.nlm.nih.gov/38056899/">
                  Differential requirements for CD4+ T cells in the efficacy of the anti-PD-1+LAG-3 and anti-PD-1+CTLA-4 combinations in melanoma flank and brain metastasis models.
                </Link>
              }
              secondary="Phadke MS, Li J, Chen Z, Rodriguez PC, Mandula JK, Karapetyan L, Forsyth PA, Chen YA, Smalley KSM. J Immunother Cancer. 2023"
            />
          </ListItem>
          <ListItem>
            <ListItemText
              primary={
                <Link href="https://pubmed.ncbi.nlm.nih.gov/33653716/">
                  Targeted Therapy Given after Anti-PD-1 Leads to Prolonged Responses in Mouse Melanoma Models through Sustained Antitumor Immunity.
                </Link>
              }
              secondary="Phadke MS, Chen Z, Li J, Mohamed E, Davies MA, et al. Cancer Immunol Res. 2021"
            />
          </ListItem>
          <ListItem>
            <ListItemText
              primary={
                <Link href="https://pubmed.ncbi.nlm.nih.gov/34035069/">
                  Single-Cell Characterization of the Immune Microenvironment of Melanoma Brain and Leptomeningeal Metastases.
                </Link>
              }
              secondary="Smalley J, Chen Z, Phadke M, Li J, et al. Clin Cancer Res. 2021"
            />
          </ListItem>
          <ListItem>
            <ListItemText
              primary={
                <Link href="https://pubmed.ncbi.nlm.nih.gov/35247927/">
                  Single cell characterization of the cellular landscape of acral melanoma identifies novel targets for immunotherapy.
                </Link>
              }
              secondary="Li J, Smalley J, Chen Z, Wu JY, et al. Clin Cancer Res. 2022"
            />
          </ListItem>
        </List>

        <Box mt={4}>
          <Typography variant="h6" gutterBottom>
            Contact Us
          </Typography>
          <Typography>
            Zhihua Chen <Link href="mailto:Zhihua.Chen@moffitt.org">Zhihua.Chen@moffitt.org</Link>
          </Typography>
          <Typography>
            Ann Chen <Link href="mailto:ann.chen@hci.utah.edu">Ann.Chen@hci.utah.edu</Link>
          </Typography>
          <Typography>
            Chloe Tran <Link href="mailto:chloe.tran@hci.utah.edu">Chloe.Tran@hci.utah.edu</Link>
          </Typography>
        </Box>

        <Box mt={4} mb={2}>
          <Typography variant="body2" color="textSecondary">
            Chen Lab
          </Typography>
        </Box>
      </Box>
    </Container>
  );
}

export default About;
