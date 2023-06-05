FROM timoast/signac AS base
WORKDIR pipeline
COPY scripts_for_docker scripts
COPY iscvamR iscvamR
RUN Rscript -e 'install.packages("./iscvamR", repos=NULL, type="source")'
