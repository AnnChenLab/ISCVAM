import { useState, useMemo } from 'react';

export function useBrushing(conf) {
    const { covs } = conf;
    const [colorBy, setColorBy] = useState('');
    const [filter, setFilter] = useState([]);
    const [contrast, setContrast] = useState(1);
    const [ zTransform, setZTransform ] = useState(false);


    function selectSubpop(v) {
        setFilter(v);
    }

    const selected = useMemo(() => {
        if (filter !== null && filter.length > 0 && covs !== null && (filter[0].value.group in covs)) {
            return covs[filter[0].value.group].map((_c, i) =>
                filter.some(({ value: { group, term } }) => covs[group][i] === term)
            );
        }
        return null;
    }, [filter,covs]);
    return  {selected, setColorBy, colorBy, filter, setFilter, selectSubpop, contrast, setContrast, zTransform, setZTransform};
}
