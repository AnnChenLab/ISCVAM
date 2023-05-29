

import { useMemo, useState } from 'react';
import { d3 } from 'plotly.js';
import { isArray } from 'util';
import settings from '../settings';

export function usePlotting({xys, colorBy, selected, covs, discreteCovs, continuousCovs, contrast, zTransform, remoteData}) {
    function hexToRgbA(hex) {
        let c;
        if (/^#([A-Fa-f0-9]{3}){1,2}$/.test(hex)) {
            c = hex.substring(1).split('');
            if (c.length === 3) {
                c = [c[0], c[0], c[1], c[1], c[2], c[2]];
            }
            c = '0x' + c.join('');
            return 'rgba(' + [(c >> 16) & 255, (c >> 8) & 255, c & 255].join(',') + ',1)';
        }
        throw new Error('Bad Hex');
    }

    function colorModifier(color, alpha, shadeAdjust = 0) {
        if (!color || !isNaN(color)) {
            return color;
        }
        if (color.startsWith('#')) {
            color = hexToRgbA(color);
        }
        if (color.startsWith('rgba')) {
            const col = color.slice(5, -1).split(',');
            const newCol = col.slice(0, -1).map(c => Math.max(parseInt(c, 10) - shadeAdjust, 0));
            if (alpha !== null) newCol.push(alpha);
            return 'rgba(' + newCol.join() + ')';
        }
        return null;
    }

    function darker(color) {
        return colorModifier(color, null, 40);
    }

    function transparent(color, alpha) {
        return colorModifier(color, alpha);
    }

    function interpolate(a, b, r) {
        return a + (b - a) * r;
    }

    const {
        maxRowsInTooltip,
        backgroundOpacity,
        backgroundSize,
        highlightOpacity,
        highlightSize,
        highlightBorderWidth,
        backgroundBorderWidth
    } = settings.plottingDefaults;

    const brushedHighlightSize = interpolate(backgroundSize, highlightSize, contrast);

    const frames = [];
    const config = {
        // displayModeBar: null, 
        displaylogo: false,
        modeBarButtonsToRemove: ['sendDataToCloud',
        'hoverCompareCartesian', 'zoom2d',
        'pan2d', 'select2d', 'lasso2d', 'zoomIn2d',
        'zoomOut2d', 'autoScale2d', 'resetScale2d', 
        'hoverClosestCartesian', 'toggleSpikelines'
    ]};
    const useResizeHandler = true;
    const style={};
    const updateFigure = f => setLayout(f.layout);

    

    // split data according to color into traces
    // affecting: data, colors, linecolor, annos, bgcolors, opacity, size, lineopacity, linewidth
    // coordinate with selected

    const colorsAndTraces = useMemo( () => {
        const ret = {colors:null,traceIdxes:null, isContinous:false}
        if ( colorBy == null) {
            return ret;
        }
        if ( colorBy !== '' && colorBy in covs && discreteCovs.includes(colorBy)) {
            const scale = d3.scale.category20();
            // hack stable colors for TCR BCR
            const vals = [...new Set(covs[colorBy])].sort();
            vals.forEach(v => scale(v));
            const valsMap={};
            vals.forEach( (v,i)=> valsMap[v]=i);
            const legend = {};
            for(const val of vals) {
                legend[val] = scale(val);
            }
            ret.colors = covs[colorBy].map(cov => scale(cov));
            ret.traces = vals.map(val=> ({
                idxes:[],
                name: val
            }));
            covs[colorBy].forEach( (cov,i) => 
                ret.traces[valsMap[cov]].idxes.push(i)
            );
            return ret;
        } 
        if ( colorBy !== '' && colorBy in covs && continuousCovs.includes(colorBy)) {
            // if(zTransform) {
            //     const stdev = std(covs[colorBy]);
            //     if(stdev===0) ret.colors = covs[colorBy];
            //     else ret.colors = divide(subtract(covs[colorBy], mean(covs[colorBy])), stdev);
            // } else {ret.colors = covs[colorBy];}
            ret.colors=covs[colorBy];
            ret.isContinous=true;
            return ret;
        } 
        if ( colorBy !== '' && ! (colorBy in covs) && remoteData !== null ) {
            ret.colors = remoteData.data;
            ret.isContinous=true;
            // console.log(mean(remoteData.data));
            // console.log(ret.colors);
            return ret;
        }
        return ret;
    }, [colorBy, zTransform, covs,remoteData]);
    // const linecolor = colors;
    const colors = colorsAndTraces.colors;
    const linecolor = useMemo( () => colors == null ? null : 
        colors.map( c => darker(c)), [colorBy,covs,remoteData]);
    const bgcolors = useMemo( () => colors == null ? null : 
        colors.map( c => transparent(c, 0.8)), [colorBy, covs, remoteData]);


    const annos = useMemo( () => {
        if ( covs == null ) return null;
        // const covnames = Object.keys(covs);
        const covnames = Object.keys(covs).slice(0, maxRowsInTooltip);
        if ( covnames.length === 0 ) return null;
        return covs[covnames[0]].map( (x, idx) => {
            if ( selected != null && !selected[idx]) {
                return null;
            }
            let ret = '';
            if (colorBy !== '' && colorBy in covs) {
                ret += `${colorBy}: ${covs[colorBy][idx]}<br><br>`;
            } else {
                if (colorBy !== '' && colorBy !== null && ! (colorBy in covs) && remoteData != null ) {
                    ret += `${colorBy}: ${remoteData.data[idx]}<br><br>`;
                }
            }
            ret += covnames.map( cn => `${cn}: ${covs[cn][idx]}`).join('<br>')
            return ret;
        })
    }, [covs, colorBy, selected, remoteData]);

    const opacity = useMemo( () => selected === null ? 0.6 : selected.map( s => s ? highlightOpacity : backgroundOpacity),
    [selected]);
    const size =  useMemo( () => selected === null ? 6 : selected.map( s => s ? brushedHighlightSize : backgroundSize), [selected]);
    const lineopacity =  useMemo( () => selected === null ? null : selected.map( s => s ? highlightOpacity : backgroundOpacity), [selected]);
    const linewidth =  useMemo( () => selected === null ? null : selected.map( s => s ? highlightBorderWidth : backgroundBorderWidth), [selected]);


    // const [ xys, setXys ] = useState([null,null]);
    const [ layout, setLayout ] = useState({
        hovermode:'closest',
        showlegend: true,
        margin: {
            l: 0,
            r: 0,
            t: 0,
            b: 0
        },
        legend: {
            x: 1.2,
            xanchor: 'right',
            y: 0.95
          },
        xaxis: {
            showline: false,
            showgrid: false,
            zeroline: false,
            showticks: false,
            showticklabels: false
        },
        yaxis: {
            showline: false,
            showgrid: false,
            zeroline: false,
            showticks: false,
            showticklabels: false
        },
        autosize: true
    });
    const traces = colorsAndTraces.traces;
    let data;
    if(traces==null) {
        data = [{
            x: xys[0],
            y: xys[1],
            type: 'scattergl',
            mode: 'markers',
            marker: { 
                color: colors,
                colorscale: zTransform?'Picnic':null,
                colorbar: colorsAndTraces.isContinous?{title:colorBy,len:0.9}:null,
                size: size, 
                opacity: opacity, 
                line: {
                    width: linewidth,
                    opacity: lineopacity,
                    color: linecolor
            }},
            // marker: {color: this.state.colorAttr, size: 8, opacity:0.8},
            text: annos,
            hoverlabel: {bgcolor: bgcolors},
            // hoverlabel: {bgcolor: 'rgba(0,0,0,0)'},
            hoverinfo: 'text',
            name: ''
        }];  
    } else {
        data = traces.map(trace => ({
            x: trace.idxes.map(idx=>xys[0][idx]),
            y: trace.idxes.map(idx=>xys[1][idx]),
            type: 'scattergl',
            mode: 'markers',
            marker: { 
                color: trace.idxes.map(idx=>colors[idx]), 
                colorscale: zTransform?'Picnic':null,
                colorbar: colorsAndTraces.isContinous?{title:colorBy, len:0.9}:null,
                size: isArray(size)?trace.idxes.map(idx=>size[idx]):size,
                opacity: isArray(opacity)?trace.idxes.map(idx=>opacity[idx]):opacity, 
                line: {
                    width: linewidth==null?null:trace.idxes.map(idx=>linewidth[idx]),
                    opacity: lineopacity==null?null: trace.idxes.map(idx=>lineopacity[idx]),
                    color: linecolor==null?null:trace.idxes.map(idx=>linecolor[idx])
            }},
            // marker: {color: this.state.colorAttr, size: 8, opacity:0.8},
            text: trace.idxes.map(idx=>annos[idx]),
            hoverlabel: {bgcolor: bgcolors==null?null:trace.idxes.map(idx=>bgcolors[idx])},
            // hoverlabel: {bgcolor: 'rgba(0,0,0,0)'},
            hoverinfo: 'text',
            name: trace.name
        }))
    }  
    const plot = {data, frames, config, layout, useResizeHandler, style, onInitialized: updateFigure,
        onUpdate: updateFigure};
    // function changeCoords(coords) {
    //     setXys(coords === null ? [null, null] : [0,1].map(i => coords.map( c => c[i])));
    // }
    return plot;
}