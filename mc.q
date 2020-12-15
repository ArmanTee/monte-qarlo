system"pwd"
\l mcutils.q
\l mcplot.q
// Inverse Sampling
.mc.inv:{[fn;n]
            fn n?1.
            };

/// Box-Muller Method
.mc.norm.bxml:{[n;s;m]
    u1:(c:ceiling[n%2])?1.;
    u2:c?1.;
    m+s*n#(sqrt[-2*log(u1)]*cos[2*.mc.pi*u2]),sqrt[-2*log(u2)]*sin[2*.mc.pi*u1]
    };


.mc.norm.fn:{[s;m;x]
            %[1;s*sqrt[2*.mc.pi]]*exp -0.5*%[x-m;s] xexp 2
    
    };

// REJECTION
/internal
.mc.i.rejsp:{[n;m;h;fn;s;e;x] 
    x,y where(fn y:s+n?"f"$e)>n?m*h
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
    if[(0> type o.s)<>100>type o.g; 0N!"ERROR - Incorrect g/o arguments supplied";:()];
    if[0b~o`g;
            o.s:.mc.norm.bxml[n;1;0];
            o.g:.mc.norm.fn[1;0;]
        ];
    o[`avg`sderr]:.mc.i.impCalcNorm[o.s;f;o.g;h];
    if[o`std;
        o[`std]: sqrt first[.mc.i.impCalcNorm[o.s;f;o.g;{y[x] xexp 2}[;h]]] - xexp[o`avg; 2 ]
        ];
    o    
    }
    
// Script

/ 
// Commands
expF:{(1%y)*exp neg[x]%y};
expF2:expF[;2];
invExpF:{neg[y]*log[1-x]};
invExpF2:invExpF[;2];

.mc.plot.fitHist.go[.mc.inv[invExpF2;10000]; expF2;0.1;1b;`small]
.mc.plot.fitHist.go[.mc.rej.sph[fn;n;0;50;0b;m]; expF2;0.1;1b;`small]
.mc.rej.sph[expF2;n;0;500;0b;0b]
n:120000
lin:{(1%40)*3+2*x};
Finv:{0.5*(-3+sqrt[9+(160*x)])};
avg .mc.inv[Finv;100000000];

avg a*lin a:.mc.utils.linspace[0;5;10000000]

/.mc.plot.fitHist.go[.mc.inv[Finv;n]; lin;0.1;1b;`small]

/.mc.plot.fitHist.go[.mc.rej.sph[fn;n;s;e;0b;m]; lin;0.1;1b;`small]
\