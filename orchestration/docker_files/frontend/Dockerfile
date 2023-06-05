FROM node:18-alpine AS base
RUN apk add g++ make python3 yarn 
RUN apk upgrade

FROM base AS build

WORKDIR frontend/
COPY iscva-frontend ./

RUN  NODE_TLS_REJECT_UNAUTHORIZED=0 yarn build 

FROM nginx:1.23.1-alpine
EXPOSE 80
COPY ./docker/nginx/conf.d/default.conf /etc/nginx/conf.d/default.conf
COPY --from=build /frontend/build /usr/share/nginx/html
