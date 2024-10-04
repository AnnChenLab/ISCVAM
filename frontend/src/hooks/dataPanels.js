import { useRef, useReducer, useEffect } from 'react';

export function useDataPanels({datasets, initialDatasetIdx, nPanels}) {
    const initialCache = {};
    const cache = useRef(initialCache);
    console.log("Initial Cache: ", initialCache)

    const initialState = 
        [...Array(nPanels)].map(() => ({ data: null }));
    console.log("Initial Cache: ", initialCache)   
     
    const [state, dispatch] = useReducer(reducer, initialState);

    function deleteFromCache(datasetIdx) {
        if(datasetIdx in cache.current)
            delete cache.current[datasetIdx]
    }

    function reducer(state, action) {
        const {type} = action;
        const stateToUpdate = [...state];

		switch (type) {
			case 'FETCHING':
                const currentState = stateToUpdate[action.panelIdx];
                if(currentState.data !== null) 
                    deleteFromCache(currentState.data.datasetIdx);
                currentState.data = null;
                break;

			case 'FETCHED':
                break;

			case 'FETCH_ERROR':
                break;

            case 'LOAD_PANEL':
                stateToUpdate[action.panelIdx].data = cache.current[action.datasetIdx].store;
                cache.current[action.datasetIdx].rc++;
                break;        

            case 'LOAD_ALL_PANELS':
                // stateToUpdate.forEach(p => p.data = cache.current[action.datasetIdx].store);
                // cache.current[action.datasetIdx].rc += stateToUpdate.length;
                // break;
                console.log("Action datasetIdx:", action.datasetIdx);
                console.log("Cache current:", cache.current);
            
                if (cache.current[action.datasetIdx]) {
                    stateToUpdate.forEach(p => p.data = cache.current[action.datasetIdx].store);
                    cache.current[action.datasetIdx].rc = stateToUpdate.length;
                } else {
                    console.error(`Invalid datasetIdx: ${action.datasetIdx}`);
                    console.error("Cache current:", cache.current);
                }
                break;

			default:
                throw new Error();        
        }

        return stateToUpdate;
    }

    function loadDataLayer(res, layer, data, clusterings) {
        let layerData = {};
        layerData.covs = res.covs;
        layerData.continuousCovs = res.continuousCovs;
        layerData.discreteCovs = res.discreteCovs;
        layerData.clusterings = clusterings;
        data.layersData[layer.id] = layerData;
    }

    async function loadDataset(datasetIdx) {
        if(!datasetIdx && datasetIdx !== 0) return null;

        if (!cache.current[datasetIdx]) {
            try {
                const {name, url, layers, modalities, zip} = datasets[datasetIdx];
                const data = {layers, modalities, name, url, zip, layersData:{}};

                const jobs = layers.map(async layer =>  {
                    const layerURL = url in layer ? layer.url : `${url}/${layer.id}`;
                    const response = await fetch(layerURL);
                    const j = await response.json();
                    const clusteringsURL = `${layerURL}/clusterings`;
                    const clusterings = await (await fetch(clusteringsURL)).json();
                    loadDataLayer(j, layer, data, clusterings);  
                });

                let genesets;

                jobs.push( 
                    (async () =>  {
                        const metares = await fetch(`${url}/meta`);
                        data.meta = await metares.json();
                    }) ()
                );

                jobs.push(
                    (async () => {
                        genesets = await (await fetch(`${url}/genesets`)).json();
                    })()
                );             

                await Promise.all(jobs);

                cache.current[datasetIdx] ={
                    store:{    
                        data,
                        name,
                        genesets,
                        dataAPI: url,
                        datasetIdx
                    },
                    rc: 0
                };

            } catch (error) {
                return null;
            }
        } 
    }

    async function loadDataPanel(datasetIdx, panelIdx) {
        dispatch({ type: 'FETCHING', panelIdx });
        await loadDataset(datasetIdx);
        dispatch({ type: 'LOAD_PANEL', panelIdx, datasetIdx });
    }
    
    useEffect(() => {
        const loadInitialDataset = async () => {
            await loadDataset(initialDatasetIdx);
            dispatch({ type: 'LOAD_ALL_PANELS', datasetIdx: initialDatasetIdx });
        };
        loadInitialDataset();
    }, [datasets, initialDatasetIdx, nPanels]);
    
    return [state, loadDataPanel];
};