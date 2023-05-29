import { useState } from 'react';

// Helper function to replace certain characters with dashes
const fixDashes = (str) => str.replace(/\u2013|\u2014|\u2212/g, "-");

export const useRemoteData = ({ dataAPI, covs }) => {
  const [remoteData, setRemoteData] = useState(null);

  const loadRemoteData = async (featureCat, feature) => {
    feature = fixDashes(feature);
    const response = await fetch(`${dataAPI}/${featureCat}/${feature}`);
    const res = await response.json();
    
    // Exit early if response data is empty
    if (Object.keys(res).length === 0) {
      setRemoteData(null);
      return null;
    }

    const assay = featureCat === 'geneset' ? 'geneset' : Object.keys(res)[0];
    const resData = featureCat === 'geneset' ? res : res[assay];

    const rd = {
      data: covs["id"].map(id => id in resData ? resData[id] : 0),
      category: assay
    };

    setRemoteData(rd);
    return rd;
  }

  return [remoteData, loadRemoteData];
}
