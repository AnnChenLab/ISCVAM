const cors = require('cors');
const express = require('express')
const app = express()
const port = 8002;
const host = '::';
const compression = require('compression');
const path = require('path');

app.use(cors());
app.use(compression());

// gene sets
var fs  = require("fs");
var genesets={}
function loadGeneSets() {
  var path = "./genesets/msigdb.v6.2.symbols.gmt";
  var lines = fs.readFileSync(path).toString().split('\n');
  lines.forEach( line => {
    var tokens = line.split('\t');
    genesets[tokens[0]] = {genes:tokens.slice(2),url:tokens[1]};
  });
}
loadGeneSets();
console.log(`${Object.keys(genesets).length} gene sets loaded.`);
  //{name: 'tirosh', filename: 'data/tirosh/test.h5' },

// var datasets = [
//     {name: 'discovery_dataset', filename: '/uufs/chpc.utah.edu/common/home/u6057318/Documents/ISCVAM_thanh_update/data/demo/test_discover_42.h5'},
//   {name: 'PBMC_multiome', filename: '/uufs/chpc.utah.edu/common/home/u6057318/Documents/ISCVAM_thanh_update/data/demo/PBMC_multiome_2.h5'},
//   {name: 'internal_validation', filename: '/uufs/chpc.utah.edu/common/home/u6057318/Documents/ISCVAM_thanh_update/data/demo/internal_validation_multiome_5.h5'},
  
// ]
// Read the JSON file
const jsonFilePath = path.join(__dirname, 'datasets.json');
let datasets;

try {
    const jsonData = fs.readFileSync(jsonFilePath, 'utf8');
    datasets = JSON.parse(jsonData);
} catch (err) {
    console.error('Error reading or parsing the JSON file:', err);
    datasets = [];
}

console.log(datasets);
//console.log('Datasets:', datasets);

var hdf5 = require('hdf5.node').hdf5;
var h5lt = require('hdf5.node').h5lt;
var h5tb = require('hdf5.node').h5tb;
var Access = require('hdf5.node/lib/globals').Access;
var H5Type = require('hdf5.node/lib/globals').H5Type;

/*
function u64ArrayToU32(buffer, len) {
  var ret = new Uint32Array(len);
  for(let i = 0 ; i < len; i++) {
    ret[i]=buffer.readUInt32LE(i*8)
  }
  return ret;
}*/
function U8ArrayToU32(buffer) {
  return Uint32Array.from(buffer)
}

function u64ArrayToU32(buffer, len) {
  for(let i = 0 ; i < len; i++) {
    buffer.writeUInt32LE(buffer.readUInt32LE(i*8), i*4)
  }
  return new Uint32Array(buffer.buffer.slice(0, len*4));
}
function min(arr) {
  if ( ! arr || !arr.length || arr.length < 1 ) {
    return undefined;
  }
  var ret = Number.MAX_VALUE;
  arr.forEach( v => {
    if(v < ret) {
      ret = v;
    }
  })
  return ret;
}
function max(arr) {
  if ( ! arr || !arr.length || arr.length < 1 ) {
    return undefined;
  }
  var ret = Number.MIN_VALUE;
  arr.forEach( v => {
    if(v > ret) {
      ret = v;
    }
  })
  return ret;
}

function readAssayData(assayGroup) {
  var children = assayGroup.getMemberNames();
  var group;
  if(children.includes('matrix') ) { 
    group = assayGroup.openGroup('matrix');
  } else if (children.includes('GRCh38') ) {
    group = assayGroup.openGroup('GRCh38');
  }
  var barcodes = h5lt.readDataset(group.id, 'barcodes');
  var data = h5lt.readDataset(group.id, 'data');
  var indices = h5lt.readDataset(group.id, 'indices');
  var barcodesLen = group.getDatasetDimensions('barcodes');
  var dataLen = group.getDatasetDimensions('data');
  var indicesLen = group.getDatasetDimensions('indices');

  var indptrLen = group.getDatasetDimensions('indptr');
 
  var indicesByteLen = indices.length / indicesLen;
  var uint32size = 4;
  console.assert(barcodesLen < Math.pow(2, uint32size*8));
  console.assert(dataLen < Math.pow(2, uint32size*8));
  // force uint64 to uint32
  //indices = u64ArrayToU32(indices, indicesLen);
  var indptr = h5lt.readDataset(group.id, 'indptr');

  if(indptrLen != indptr.length) {
    // 10xgenomics u64 array
    indptr = u64ArrayToU32(indptr, indptrLen);
  }
  var shape = h5lt.readDataset(group.id, 'shape');

  var features;
  if(group.getMemberNames().includes('gene_names') ) { 
    features = h5lt.readDataset(group.id, 'gene_names' );
  } else {
    // assuming a features group
    var featuresGroup = group.openGroup('features')
    features = h5lt.readDataset(featuresGroup.id, 'name' );
    featuresGroup.close();
  }
  group.close();
  // var featureMap = new Map(features.map( (n,i) => [n.toUpperCase(), i]));
  var featureMap = new Map(features.map( (n,i) => [n, i]));
  return {barcodes, data, indices, indicesLen, indptr, shape, featureMap};
}
function loadDataFromDisk(ds) {
  // old structure: matrix (or grch38), artifacts
  // new structure: assays - rna, atac - matrix (or grch38), artifacts

  console.log('loading from disk');
  console.time("initial lookup from disk 1");
  var file = new hdf5.File(ds.filename, Access.ACC_RDONLY);
  var children = file.getMemberNames();
  var ret = {};
  if(children.includes('assaysData')) {
    // new multiple assay format
    var assaysGroup = file.openGroup('assaysData');
    ret.assays = assaysGroup.getMemberNames();
    ret.assaysData = {};
    ret.assays.forEach(assay => {
      var assayGroup = assaysGroup.openGroup(assay);
      ret.assaysData[assay]=readAssayData(assayGroup);
      assayGroup.close();
    });
    assaysGroup.close();
  } else {
    console.log('readingAssayData');
    ret.assays = ['RNA'];
    ret.assaysData = {RNA: readAssayData(file)}
    console.log('done readingAssayData');
  }
  // read visualization artifacts
  ret.artifacts = {};
  if(file.getMemberNames().includes('artifacts') ) {
    console.log("reading artifacts");
    var artiGroup = file.openGroup('artifacts');
    var artiSubGroups = artiGroup.getMemberNames();
    var pcs, covs, discreteCovs, continuousCovs;
    artiSubGroups.forEach( ag => {
      console.log(`reading artifact ${ag}`);
      var artiSub = file.openGroup('artifacts/'+ag);
      //pcs = h5tb.readTable(artiSub.id, 'pcs' );
      covs = h5tb.readTable(artiSub.id, 'covs' );
      discreteCovs = h5lt.readDataset(artiSub.id, 'discreteCovs' );
      continuousCovs = h5lt.readDataset(artiSub.id, 'continuousCovs' );
      // ret.artifacts[ag]={pcs,covs,discreteCovs,continuousCovs};
      ret.artifacts[ag]={covs,discreteCovs,continuousCovs};
      console.log(`done reading artifact ${ag}`);
      artiSub.close();
    })
    artiGroup.close();
  }
  file.close();
  return ret;
  // return {barcodes, data, indices, indicesLen, indptr, shape, geneMap, artifacts };
}

var cache = (function () {
  var _cache = {};
  var limit = 10;
  var _cacheKeysInOrder = [];

  function get(key) {
    if(key in _cache) {
      return _cache[key]
    } else {
	throw 'key not found';
    }
  }
  function put(key, val)  {
    _cache[key] = val
    if ( ! key in _cacheKeysInOrder )  {
      _cacheKeysInOrder.push(key);
    }
    if(Object.keys(cache).length > limit) {
      // cache is full
      // eject oldest
      console.log('ejecting old cache');
      delete _cache[_cacheKeysInOrder[0]];
      _cacheKeysInOrder = _cacheKeysInOrder.slice(1);
    }
  }
  return {get, put};
})();

function loadDataset(ds) {
  console.log('loaddataset');
  console.log(ds);
  let ret;
  try {
    ret = cache.get(ds.name);
  } catch(e) {
    console.log(e);
    if (e == 'key not found') {
      console.log('key not found in cache');
      ret = loadDataFromDisk(ds);
      cache.put(ds.name, ret)
    }
  }
  return ret;
}

function lookupGenesetFromDisk(ds, genesetId) {
  const dataCache = loadDataset(ds);
  console.assert('assays' in dataCache);
  console.assert('RNA' in dataCache.assaysData);
  let  {barcodes, data, indices, indicesLen, indptr, shape, featureMap}  = dataCache.assaysData.RNA;
 
  var geneIds = new Set(genesets[genesetId]["genes"].map( g => featureMap.get(g)).filter(a=>a!==undefined));
  console.log(`${geneIds.size} genes in geneset ${genesetId}`);

  const ret = {};
  console.time('going through the barcodes');
  if(indices.length==indicesLen) {
  barcodes.forEach( (barcode, i) => {
    for(let cursor = indptr[i]; cursor < indptr[i+1]; cursor++) {
      let idx;
      idx = indices[cursor];
      if(geneIds.has(idx)) {
        if(! (barcode in ret ) ) {
          ret[barcode] = data[cursor];
        } else {
          ret[barcode] += data[cursor];
        }
      }
    }

    if(barcode in ret) {
      // normalize to total reads in barcode
      let total = 0;
      for(let cursor = indptr[i]; cursor < indptr[i+1]; cursor++) {
        total += data[cursor]
      }
      ret[barcode] = ret[barcode]/total*100
    } 
  })


  } else {

  barcodes.forEach( (barcode, i) => {
//    console.log(`number of genes: ${indptr[i+1]-indptr[i]+1}`);
    for(let cursor = indptr[i]; cursor < indptr[i+1]; cursor++) {
      let idx;
      idx = indices.readUInt32LE(cursor*8);
      //idx = indicesLen==indices.length?indices[cursor]:indices.readUInt32LE(cursor*8);
      if(geneIds.has(idx)) {
        if(! (barcode in ret ) ) {
          ret[barcode] = data[cursor];
        } else {
          ret[barcode] += data[cursor];
        }
      }
    }

    if(barcode in ret) {
      // normalize to total reads in barcode
      let total = 0;
      for(let cursor = indptr[i]; cursor < indptr[i+1]; cursor++) {
        total += data[cursor]
      }
      ret[barcode] = ret[barcode]/total*100
    } 
  })
}
  console.timeEnd('going through the barcodes');
  return ret;
}

function lookupFromDisk(ds, gene) {

  const dataCache = loadDataset(ds);
  console.assert('assays' in dataCache, 'assays in datacache');
  console.assert('RNA' in dataCache.assaysData,   'rna in assaysdata');
  const ret = {};
  // if(dataCache.assaysData.RNA.featureMap.has(gene.toUpperCase())) {
    if(dataCache.assaysData.RNA.featureMap.has(gene)) {
      ret.RNA = extractFeatureData(dataCache.assaysData.RNA, gene);
  }
  return ret.RNA;  
}

function extractFeatureData(assayData, feature) {
  console.log(`extract feature data ${feature}`);
  // console.log(assayData);
  // console.log(assayData.featureMap);
  let  {barcodes, data, indices, indicesLen, indptr, shape, featureMap}  = assayData;


  // var geneIdx = featureMap.get(feature.toUpperCase());
  var geneIdx = featureMap.get(feature);
  // var matchData = [];
    // d.indices.forEach( (ind, i) => {if(ind === geneIdx) {matchData[matchData.length]=d.data[i]}});
  
  const ret = {};
  console.time('going through the barcodes');

  if(indices.length==indicesLen) {
    barcodes.forEach( (barcode, i) => {

      let total = 0;
      for(let cursor = indptr[i]; cursor < indptr[i+1]; cursor++) {
        total += data[cursor]
      }
  //    console.log(`number of genes: ${indptr[i+1]-indptr[i]+1}`);

      for(let cursor = indptr[i]; cursor < indptr[i+1]; cursor++) {
        let idx;
        // temporary solution to support UINT8 arrays
        idx =  indices[cursor] ;
  //      console.log(idx);       
        //if(indices[cursor]==geneIdx) {
        if(idx==geneIdx) {
          ret[barcode] = Math.log1p(data[cursor]/total*10000);
          break;
        }
      }
    })  
  } else {
  
    barcodes.forEach( (barcode, i) => {
  
      let total = 0;
      for(let cursor = indptr[i]; cursor < indptr[i+1]; cursor++) {
        total += data[cursor]
      }
  //    console.log(`number of genes: ${indptr[i+1]-indptr[i]+1}`);
  
      for(let cursor = indptr[i]; cursor < indptr[i+1]; cursor++) {
        let idx;
        idx = indices.readUInt32LE(cursor*8);
        if(idx==geneIdx) {
          ret[barcode] = Math.log1p(data[cursor]/total*10000);
          break;
        }
      }
    })
  }
  console.timeEnd('going through the barcodes');
  return ret;
}

function lookupFeatureFromAssayFromDisk(ds, feature, assay) {
  const dataCache = loadDataset(ds);
  console.log(`lookupFeatureFromAssayFromDisk`);
  console.log(ds);
  console.assert('assays' in dataCache);
  // search for feature
  const ret = {};
 // if(dataCache[assay])
  if(dataCache.assays.includes(assay)) 
    if(dataCache.assaysData[assay].featureMap.has(feature))
      ret[assay]=extractFeatureData(dataCache.assaysData[assay], feature);
  return ret;
}

const getFeatureFromAssayFromDisk = (req, res) => {
  // console.log(data);
  console.time('lookup feature from disk');
  var d = datasets.filter(r=>r.name==req.params.projectId)[0];


  // prior routine: request - loadds (from cache if available) - do data look up
  // current routine: request - loadds (from cache if available) - do data look up

  var feature = req.params.featureId;
  var assay = req.params.assayId;
  var ret = lookupFeatureFromAssayFromDisk(d, feature, assay);
  console.timeEnd('lookup feature from disk');
  res.json(ret);
}


function lookupFeatureFromDisk(ds, feature) {
  const dataCache = loadDataset(ds);
  console.log(`lookupFeatureFromDisk`);
  console.log(ds);
  console.assert('assays' in dataCache);
  // search for feature
  const ret = {};
  dataCache.assays.forEach(assay => {
    // if(dataCache.assaysData[assay].featureMap.has(feature.toUpperCase())) {
      if(dataCache.assaysData[assay].featureMap.has(feature)) {
        ret[assay]=extractFeatureData(dataCache.assaysData[assay], feature);
    }
  });
  // if no match for feature, try look up nearby chromatin peaks in ATAC
  return ret;  
}

app.get('/iscvam/backend/genesets/', (req, res) => {
  res.json(Object.keys(genesets));
});
app.get('/iscvam/backend/project/:projectId/genesets/', (req, res) => {
  res.json(Object.keys(genesets));
});


const getGenesetFromDisk = (req,res) => {
  console.time('lookup gene set from disk');
  var d = datasets.filter(r=>r.name==req.params.projectId)[0];
  var genesetId = req.params.genesetId;
  if( !(genesetId in genesets)) {
    res.status(400).send(`${genesetId} not found`);
    return;
  }
  var ret = lookupGenesetFromDisk(d, genesetId);
  console.timeEnd('lookup gene set from disk');
  res.json(ret);
}
const getGeneFromDisk = (req, res) => {
  // console.log(data);
  console.time('lookup from disk');
  var d = datasets.filter(r=>r.name==req.params.projectId)[0];
  var gene = req.params.geneId;
  var ret = lookupFromDisk(d, gene);
  console.timeEnd('lookup from disk');
  res.json(ret);
}

const getFeatureFromDisk = (req, res) => {
  // console.log(data);
  console.time('lookup feature from disk');
  var d = datasets.filter(r=>r.name==req.params.projectId)[0];


  // prior routine: request - loadds (from cache if available) - do data look up
  // current routine: request - loadds (from cache if available) - do data look up

  var feature = req.params.featureId;
  var ret = lookupFeatureFromDisk(d, feature);
  console.timeEnd('lookup feature from disk');
  res.json(ret);
}


app.get('/iscvam/backend/project/:projectId/geneset/:genesetId',  getGenesetFromDisk);
app.get('/iscvam/backend/project/:projectId/gene/:geneId',  getGeneFromDisk);
app.get('/iscvam/backend/project/:projectId/feature/:featureId',  getFeatureFromDisk);

app.get('/iscvam/backend/from_disk/project/:projectId/geneset/:genesetId',  getGenesetFromDisk);
app.get('/iscvam/backend/from_disk/project/:projectId/gene/:geneId',  getGeneFromDisk);


// 
app.get('/iscvam/backend/project/:projectId/meta', (req, res) => {
  // console.log(data);
  //var d = data.get(req.params.projectId);
  var d = datasets.filter(r=>r.name==req.params.projectId)[0];
  let  {assays, assaysData}  = loadDataset(d);
  features={}
  assays.forEach(assay => features[assay] = [...assaysData[assay].featureMap.keys()]);
  res.json({assays, features});
});

app.get('/iscvam/backend/project/:projectId/:artifactId/clusterings', (req, res) => {
  // console.log(data);
  //var d = data.get(req.params.projectId);
  var ds = datasets.filter(r=>r.name==req.params.projectId)[0];
  console.log(ds);
  const projectId=req.params.projectId;
  const artifactId=req.params.artifactId;

  var ret;
  try {
    var file = new hdf5.File(ds.filename, Access.ACC_RDONLY);
    // var children = file.getMemberNames();

//    if(file.getMemberNames().includes('clusterings') ) {

      var clusteringsSub = file.openGroup('/artifacts/'+artifactId+'/clusterings');
      ret=clusteringsSub.getMemberNames();
 //   }
  } catch(e) {
    console.log('error loading '+projectId+' '+artifactId +' clusterings');
    console.log(e);
    ret=e;
  }
  res.json(ret);
  console.log(ret);
});

app.get('/iscvam/backend/project/:projectId/:artifactId/clustering/:clusteringId/:fullmatrix?', (req, res) => {
  // console.log(data);
  console.time('load clustering');
  //var d = data.get(req.params.projectId);
  var ds = datasets.filter(r=>r.name==req.params.projectId)[0];
  console.log(ds);
  const projectId=req.params.projectId;
  const artifactId=req.params.artifactId;
  const clusteringId=req.params.clusteringId;
  const returnFullMatrix=req.params.fullmatrix;

  var ret = {};
  try {
    var file = new hdf5.File(ds.filename, Access.ACC_RDONLY);
    // var children = file.getMemberNames();

//    if(file.getMemberNames().includes('clusterings') ) {

      console.log('opening the group: '+'/artifacts/'+artifactId+'/clusterings/'+clusteringId);
      var clusteringSub = file.openGroup('/artifacts/'+artifactId+'/clusterings/'+clusteringId);
      console.log('group opened');

      var assays = clusteringSub.getMemberNames();
      if(["markers","heatmap_markers_scaled", "heatmap_markers_lognorm"].some( n => assays.includes(n))) {
        ret.markers = h5tb.readTable(clusteringSub.id, `markers` ).reduce((t,c)=>{if(!c.name.includes("feature.name")) {t[c.name]=c} return t},{});
        ret.avgsContrast = h5tb.readTable(clusteringSub.id, `heatmap_markers_scaled` );
        ret.avgsLogNorm = h5tb.readTable(clusteringSub.id, `heatmap_markers_lognorm` );
        ret.heatmapMarkers = ret.avgsContrast[0];
  ret.clusterNames = ret.avgsContrast.map(r=>r.name).slice(1);
//  ret.clusterNames = ret.allAvgsContrast.map(r=>r.name).slice(1);
ret.avgsContrast = ret.avgsContrast.map(c=>c.name.includes("feature.name")?c:c.map(r=>+r));
ret.avgsLogNorm = ret.avgsLogNorm.map(c=>c.name.includes("feature.name")?c:c.map(r=>+r));
if(returnFullMatrix !=null && returnFullMatrix) {
  ret.allAvgsContrast = (h5tb.readTable(clusteringSub.id, `all_markers_scaled` )).map(c=>c.name.includes("feature.name")?c:c.map(r=>+r));
  ret.allAvgsLogNorm = (h5tb.readTable(clusteringSub.id, `all_markers_lognorm` )).map(c=>c.name.includes("feature.name")?c:c.map(r=>+r));

}
//        ret.allAvgsContrast = ret.allAvgsContrast.map(c=>c.name.includes("feature.name")?c:c.map(r=>+r));
//        ret.allAvgsLogNorm = ret.allAvgsLogNorm.map(c=>c.name.includes("feature.name")?c:c.map(r=>+r));
//  ret.markerTableColumns = ret.markers.map(r=>r.name);
//  ret.allAvgsContrast.forEach(r=> delete r.name);
//  ret.allAvgsLogNorm.forEach(r=> delete r.name);
//  ret.avgsContrast.forEach(r=> delete r.name);
//  ret.avgsLogNorm.forEach(r=> delete r.name);

} else {
assays.forEach(assay => {
    ret[assay]={};
    ret[assay].markers = h5tb.readTable(clusteringSub.id, `${assay}/markers` ).reduce((t,c)=>{if(!c.name.includes("feature.name")) {t[c.name]=c}return t},{});
    ret[assay].avgsContrast = h5tb.readTable(clusteringSub.id, `${assay}/heatmap_markers_scaled` );
    ret[assay].avgsLogNorm = h5tb.readTable(clusteringSub.id, `${assay}/heatmap_markers_lognorm` );
    ret[assay].heatmapMarkers = ret[assay].avgsContrast[0];
ret[assay].clusterNames = ret[assay].avgsContrast.map(r=>r.name).slice(1);
//      ret[assay].clusterNames = ret[assay].allAvgsContrast.map(r=>r.name).slice(1);
    ret[assay].avgsContrast =
ret[assay].avgsContrast.map(c=>c.name.includes("feature.name")?c:c.map(r=>+r));
    ret[assay].avgsLogNorm = ret[assay].avgsLogNorm.map(c=>c.name.includes("feature.name")?c:c.map(r=>+r));
if(returnFullMatrix !=null && returnFullMatrix) {
      ret[assay].allAvgsContrast =
(h5tb.readTable(clusteringSub.id, `${assay}/all_markers_scaled` )).map(c=>c.name.includes("feature.name")?c:c.map(r=>+r));
      ret[assay].allAvgsLogNorm = (h5tb.readTable(clusteringSub.id, `${assay}/all_markers_lognorm` )).map(c=>c.name.includes("feature.name")?c:c.map(r=>+r));
}

});

}
} catch(e) {
  console.log('error loading '+projectId+' '+artifactId +' '+clusteringId);
  console.log(e);
  ret=e;
}
res.json(ret);
console.timeEnd('load clustering');
});

// project  gene
app.get('/iscvam/backend/project/:projectId/:artifactId', (req, res) => {
  // console.log(data);
  console.time('load artifact');
  //var d = data.get(req.params.projectId);
  var d = datasets.filter(r=>r.name==req.params.projectId)[0];
  console.log(`load artifact`);
  console.log(d);
  let  {artifacts}  = loadDataset(d);
  // let  {pcs, covs, continuousCovs, discreteCovs}  = artifacts[req.params.artifactId];
  let  {covs, continuousCovs, discreteCovs}  = artifacts[req.params.artifactId];
//  pcs=pcs.reduce( (m,o) => {m[o.name]=Array.prototype.slice.call(o);return m}, {});
  covs=covs.reduce( (m,o) => {m[o.name]=Array.prototype.slice.call(o);return m}, {});
  // res.json({pcs,covs,discreteCovs,continuousCovs});
  res.json({covs,discreteCovs,continuousCovs});
  console.timeEnd('load artifact');
});

app.get('/project/:projectId/:assayId/:featureId',  getFeatureFromAssayFromDisk);
////// Handle the contact us
app.use(express.json());
// Define the file path where messages will be saved
const filePath = path.resolve(process.cwd(), 'user_messages/messages.json');

// Ensure the directory exists
const ensureDirectoryExists = (dirPath) => {
  if (!fs.existsSync(dirPath)) {
    fs.mkdirSync(dirPath, { recursive: true });
  }
};

// POST endpoint to save messages to a file
app.post('/iscvam/backend/api/save-message', (req, res) => {
  const { name, email, message } = req.body;

  // Basic validation of the request body
  if (!name || !email || !message) {
    return res.status(400).send({ message: 'Name, email, and message are required.' });
  }

  const newMessage = {
    name,
    email,
    message,
    date: new Date().toISOString(),
  };

  // Ensure the directory exists before saving the file
  ensureDirectoryExists(path.dirname(filePath));

  // Read the existing messages from the file
  fs.readFile(filePath, 'utf8', (err, data) => {
    if (err && err.code !== 'ENOENT') {
      return res.status(500).send({ message: 'Error reading messages', error: err });
    }

    // If the file doesn't exist or is empty, start with an empty array
    const messages = data ? JSON.parse(data) : [];

    // Add the new message to the array
    messages.push(newMessage);

    // Save the updated array back to the file
    fs.writeFile(filePath, JSON.stringify(messages, null, 2), (err) => {
      if (err) {
        return res.status(500).send({ message: 'Error saving message', error: err });
      }

      res.status(200).send({ message: 'Message saved successfully', newMessage });
    });
  });
});

// GET endpoint to retrieve all messages
app.get('/iscvam/backend/api/messages', (req, res) => {
  fs.readFile(filePath, 'utf8', (err, data) => {
    if (err && err.code !== 'ENOENT') {
      return res.status(500).send({ message: 'Error reading messages', error: err });
    }

    const messages = data ? JSON.parse(data) : [];
    res.status(200).send(messages);
  });
});



// Route for the root of your app
app.get('/iscvam/backend', (req, res) => res.send('Hello World!'));

// Route for the backend test path
app.get('/iscvam/backend', function(req, res) {
  res.send('Backend is working');
});

// Start the server on IPv6 and IPv4 (dual-stack)
app.listen(port, '::', () => {
  console.log(`Example app listening on all interfaces at port ${port}!`);
});
// Example route that should match the NGINX proxy
app.get('/iscvam/backend/project/GSE189341-acral-sc/all', function(req, res) {
  res.send('Data for project GSE189341-acral-sc all');
});

// app.get('/iscvam/', (req, res) => res.send('Hello World!'));
// app.listen(   port, host, () => console.log(`Example app listening on port ${port}!`));
// app.get('/iscvam/backend', function(req, res) {
//   res.send('Backend is working');
// });
// app.listen(8002, '0.0.0.0', () => {
//   console.log('Server is running on port 8002');
// });
// Handles any requests that don't match the above routes
// app.get('/iscvam*', (req, res) => {
//   console.log('Serving index.html');
//   res.sendFile(path.join(__dirname, '../client/build/index.html'));
// });

// app.listen(port, host, () => {
//   console.log(`Server is running on http://${host}:${port}`);
// });