
// Constants 
.mc.pi:acos -1;



// Utils
.mc.utils.linspace:{[s;e;n]
    n:n-1;
    `float$+[%[e-s;n]]\[n;s]
    };
/ Probability of outcome x within y and z
.mc.util.spInt:{(count where x within(y;z))%count[x]};

.mc.util.spIntInv:{[x;y;p]
   x where d=min abs d:p-.mc.util.spInt[x;y;] each x: asc x;
    };    
    
.mc.util.symInv:{ 
    n: (floor count[x]%2) cut m:.mc.norm.fn[sdev x;avg x;] x:asc x;
    if[3~ count n; n:(n[0];n[1],n[2])];
     x where m in  raze n ./:  0 1,' where each abs[d]=' min each abs d: n - y
    };

// Uniform Distribution
.mc.uni.sp:{[s;e;n]
        s+n?"f"$e-s
    };

// Inverse Sampling
.mc.inv:{[fn;n]
            fn n?1.
            };

/// Box-Muller Method
.mc.norm.bxml:{[n;s;m]
    u1:(c:ceiling[n%2])?1.;
    u2:c?1.;
    m+s*n#(sqrt[-2*log(u1)]*cos 2*.mc.pi*u2),sqrt[-2*log(u2)]*sin 2*.mc.pi*u1
    };


.mc.norm.fn:{[s;m;x]
            %[1;s*sqrt[2*.mc.pi]]*exp -0.5*%[x-m;s] xexp 2
    
    };

// REJECTION
/internal
.mc.i.rejsp:{[n;m;h;fn;s;e;x] 
    x,y where(fn y:s+n?"f"$e-s)>n?m*h
    };

.mc.rej.sph:{[fn;n;s;e;h;m]
    // rejection sampling
    // fn : target function
    // n: number of samples
    // s: lower limit
    // e: upper limit
    // m: envelop function multiple
    $[not h;
        h:max fn .mc.utils.linspace[s;e;10000];
        h:1%(e-s)
        ];
    if[not m;
        m:1.001
        ];
    / n2 = estimate of total iterations needed to get n samples,
    // multiplied by 1.1 to reduce chance of rerun.
    n2:ceiling 1.1* n*(e-s)*m*h;
    n#.mc.i.rejsp[n2;m;h;fn;s;e]/[n > count@;()]
    };
    
.mc.rej.sp:.mc.rej.sph[;;;;0b;0b];
// Importance Sampling
/Internal
.mc.i.impCalcNorm:{[x;f;g;h]
    x:(h[x]*w:f[x]%g[x]);
    :(avg[x]%avg w;sdev[x]%sqrt[count x])
    };

.mc.imp.spn:{[f;h;n;o]
    /f - target density
    /h - function on x
    /n - number of samples
    /o - options dictionary `std`g`s!(1b/0b; custom g[x]; samples from g[x])
    if[0b~o;o:()!()];
    o:(``std`g`s!raze(::;3#0b)),o;
    if[(0> type o[`s])<>100>type o[`g]; 0N!"ERROR - Incorrect g/o arguments supplied";:()];
    if[0b~o`g;
            o[`s]:.mc.norm.bxml[n;1;0];
            o[`g]:.mc.norm.fn[1;0;]
        ];
    o[`avg`sderr]:.mc.i.impCalcNorm[o`s;f;o`g;h];
    if[o`std;
        o[`std]: sqrt first[.mc.i.impCalcNorm[o`s;f;o`g;{y[x] xexp 2}[;h]]] - xexp[o`avg; 2 ]
        ];
    `s`g _ o
    };
