FROM node:16-alpine AS base
RUN apk add g++ make python3 hdf5-dev hdf5 yarn openssl=3.0.8-r4 ncurses=6.3_p20221119-r1
RUN apk upgrade

WORKDIR app/
COPY app.js package.json datasets.json yarn.lock ./
COPY data data/
COPY hdf5.node hdf5.node/
COPY genesets genesets/

FROM base AS test
RUN yarn install

ENTRYPOINT node app.js
#RUN yarn test
