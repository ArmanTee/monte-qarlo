// Constants 
.mc.pi:acos -1;



// Utils
.mc.utils.linspace:{[s;e;n]
    n:n-1;
    `float$+[%[e-s;n]]\[n;s]
    };
