// Constants 
.sp.pi:3.14159265359;



// Utils
.sp.utils.linspace:{[s;e;n]
    n:n-1;
    `float$+[%[e-s;n]]\[n;s]
    };




// plot 
    
.sp.plot.fitHistInv:{[fn1;fn2;s;e;n;k]   
       // Arguments
       / fn1,  density function proper
       / fn2, inverse of the cumulitive distribution or integral
       / s, lower bound
       / e, upper bound
       / n, linspace sample
       / k, samples 
    
        t: ([] x: .sp.utils.linspace[s;e;n]);
        t: update fx: fn1 x from t;
        t2:([] y:.sp.inv[fn2;k]);
    

        .qp.go[500;500] (
            .qp.stack(
                .qp.histogram[t2;`y; .qp.s.aggr[.st.a.custom[`count__; `y; {[x;t2] count[x]%sum[t2`y]}[;t2]] ]];
                .qp.line[t;`x;`fx;::]
                ))
    };

.sp.plot.fitHistNorm:{[fn1;fn2;s;e;n;k]   
       // Arguments
       / fn1,  density function proper
       / fn2, inverse of the cumulitive distribution or integral
       / s, lower bound
       / e, upper bound
       / n, linspace sample
       / k, samples 
    
        t: ([] x: .sp.utils.linspace[s;e;n]);
        t: update fx: fn1 x from t;
        t2:([] y:fn2 k);
    

        .qp.go[500;500] (
            .qp.stack(
                .qp.histogram[t2;`y; .qp.s.aggr[.st.a.custom[`count__; `y; {[x;t2] count[x]%sum[t2`y]}[;t2]] ]];
                .qp.line[t;`x;`fx;::]
                ))
    };


// functions

.sp.inv:{[fn;n]
            fn n?1.
            };


.sp.norm.bxml:{[n;s;m]
    u1:(c:ceiling[n%2])?1.;
    u2:c?1.;
    m+ s*n#(sqrt[-2*log(u1)] * cos[2*.sp.pi*u2]), sqrt[-2*log(u2)] * sin[2*.sp.pi*u1]
    };


.sp.norm.fn:{[s;m;x]
            %[1;s*sqrt[2*.sp.pi]]*exp -0.5*%[x-m;s] xexp 2
    };

.sp.norm.fn[1;0;]
           
         
exp -0.5
\t:1000000 .sp.bxml[1;2;1]

// Script

// Exponential Distribution

expF:{(1%y)*exp(neg[x]%y)};
expF2:expF[;2];
invExpF:{neg[y]*log[1-x]};
invExpF2:invExpF[;2];
.sp.plot.fitHist[expF2;invExpF2;0;50;100;1000000];



.sp.plot.fitHistNorm[.sp.norm.fn[1;0]; .sp.norm.bxml[;1;0];-5;5;100;100000]

s: -5
e: 5
n: 10000
        t: ([] x: .sp.utils.linspace[s;e;n]);
        t: update fx: .sp.norm.fn[1;0;] x from t;
        .qp.go[500;500] 
        .qp.line[t;`x;`fx;::]
