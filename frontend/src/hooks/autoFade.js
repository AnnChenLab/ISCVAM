import {useState, useEffect} from 'react';

export const useAutoFade = (props) => {
    const timeout = props == null || props.timeout == null ? 10000: props.timeout;
    // const {timeout = 10000} = props;

    // const [active, setActive] = useState(false);
    const [active, setActive] = useState(true);
    
    function handleMouseMove() {
        console.log('handling mouse move');
        if(!active)
            setActive(true);
    }
    // useEffect(()=>{
    //     let timer;
    //     if(active) {
    //         timer = setTimeout(() => setActive(false) , timeout);
    //     }
    //     if (timer != null)
    //         return () => {console.log(`clearing timeout`);clearTimeout(timer)}
    //     else return null;
    // }, [active])
    // useEffect(() => {
    //     window.addEventListener('mousemove', handleMouseMove);
    //     return () => {
    //         window.removeEventListener(
    //           'mousemove',
    //           handleMouseMove
    //         );
    //       };
    //     }, []);
    return [active];
}
