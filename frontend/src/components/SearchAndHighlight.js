import React, {useState, useRef, useMemo} from 'react';
import {PredictiveInput, CreatablePredictiveInput} from './PredictiveInput';
import {FormControl} from '@material-ui/core';

function useSNHControls() {
    const controls = ["none", "with-select-paint"]
    const [SNHControls, setSNHControls] = useState(0)

    function toggleControls() {
        setSNHControls((SNHControls+1)%controls.length);
    }
    return [controls[SNHControls], toggleControls];
}

export default function SearchAndHighlight({classes, brushing, genesets, loadRemoteData, sNHControls, covs, discreteCovs, continuousCovs, features}) {
    const [search, setSearch] = useState([]);
    const paintControlInputEl = useRef(null);

    const searchTerms = useMemo(() => discreteCovs.map(c => ({label: c, options:
        [...new Set(covs[c])].map(cv => ({value: {group: c, term: cv}, label: cv}))})), [covs, discreteCovs]);

    const colorTerms = useMemo(() => {
        let cterms = [
            {label: 'discrete covs', options: discreteCovs.map(c => ({value: c, label: c, group: 'discrete covs'}))},
            {label: 'continuous covs', options: continuousCovs.map( c => ({value: c, label: c, group: 'continuous covs'}))}
        ];
        if(features) {
            Object.entries(features).slice().reverse().forEach(([assayName, assayFeatures]) => 
                cterms.push({label: assayName, options: assayFeatures.map(c => ({value: c, label: c, group: assayName}))})
            );
        }
        cterms.push({label: 'geneset', options: genesets.map(c => ({value: c, label: c, group: 'geneset'}))});
        return cterms;
    }, [discreteCovs, continuousCovs, features, genesets]);

    function handleSearchChange(e) {
        setSearch(e);
        brushing.setFilter(e)
    } 

    // async function changeColorBy(e) {
    //     const colorBy = e ? e.value : null;
    //     if(colorBy && !(colorBy in covs)) {
    //         let category = e.group === 'geneset' ? 'geneset' : 'feature';
    //         await loadRemoteData(category, colorBy);
    //     }
    //     brushing.setColorBy(colorBy);
    // }
    async function changeColorBy(e) {
        const colorBy = e ? e.value : null;
        if(colorBy && !(colorBy in covs)) {
            //let category = e.group === 'geneset' ? 'geneset' : 'feature';
            const category = e.group;
            await loadRemoteData(category, colorBy);
        }
        brushing.setColorBy(colorBy);
    }

    if(sNHControls === 'none') {
        return null;
    }

    return (
        <div>
            <FormControl className={classes.searchControl} variant="outlined">
                <PredictiveInput
                    classes={classes}
                    options={searchTerms}
                    value={search}
                    onChange={handleSearchChange}
                    placeholder="Select "
                    isMulti
                />
            </FormControl>
            <FormControl className={classes.paintControl} variant="outlined">
                <CreatablePredictiveInput
                    classes={classes}
                    options={colorTerms}
                    ref={paintControlInputEl}
                    formatCreateLabel={(v)=> `Look up feature: ${v}`}
                    value={brushing.colorBy ? {value: brushing.colorBy, label: brushing.colorBy} : null}
                    onChange={changeColorBy}
                    placeholder="Paint "
                    isClearable
                />
            </FormControl>
        </div>
    );
}

export {SearchAndHighlight, useSNHControls}
