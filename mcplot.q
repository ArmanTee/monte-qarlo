// plot 
.mc.plot.fitHistNorm.fn:{[file;x;w;ti;size;filename]
    /x list of normal variables
    /w width of histogram bins
    /t include title 1b or 0b
    s:sdev x;
    m:avg x;
    k:count x;
    t: ([] x: .mc.utils.linspace[m-3*s;m+3*s;1000]);
    t: update fx: .mc.norm.fn[s;m;] x from t;
    t2:([] y: x);
    ($[file;.qp.png[filename];.qp.go] . $[size~`large;(1000;1000);size~`medium; (750;750);(500;500)]) 
          (
          $[ti;.qp.title["Histogram with fit of density function: ",string[m], " Standard Deviation: ",string[s]] ;(),]
            .qp.stack(
                    .qp.histogram[t2;`y;]
                         .qp.s.aggr[.st.a.custom[`count__; `y; {[x;k;w] count[x]%w*k}[;k;w]]],
                         .qp.s.labels[`x`y!("X";"Density|f(x)")],
                         .qp.s.binx[`w; w; 0],
                         .qp.s.geom[`size`colour`strokewidth!(3;`steelblue;15)];
                    .qp.line[t;`x;`fx;]
                        .qp.s.labels[`x`y!("X";"Density|f(x)")],
                        .qp.s.geom[`size`colour`strokewidth!(3;`steelblue;15)]
           
                )
        )
    };
.mc.plot.fitHistNorm.go :.mc.plot.fitHistNorm.fn[0b;;;;0b];
.mc.plot.fitHistNorm.png:.mc.plot.fitHistNorm.fn[1b];
//
.mc.plot.fitHist.fn:{[file;x;fn;w;ti;size;filename]
    /x list of normal variables
    /w width of histogram bins
    /t include title 1b or 0b
    s:sdev x;
    m:avg x;
    k:count x;
    t: ([] x: .mc.utils.linspace[min x;max x;1000]);
    t: update fx: fn x from t;
    t2:([] y: x);
    ($[file;.qp.png[filename];.qp.go] . $[size~`large;(1000;1000);size~`medium; (750;750);(500;500)]) 
          (
          $[ti;.qp.title["Gaussian Fit with Mean: ",string[m], " Standard Deviation: ",string[s]] ;(),]
            .qp.stack(
                    .qp.histogram[t2;`y;]
                         .qp.s.aggr[.st.a.custom[`count__; `y; {[x;k;w] count[x]%w*k}[;k;w]]],
                         .qp.s.labels[`x`y!("X";"Density|f(x)")],
                         .qp.s.binx[`w; w; 0],
                         .qp.s.geom[`size`colour`strokewidth!(3;`steelblue;15)];
                    .qp.line[t;`x;`fx;]
                        .qp.s.labels[`x`y!("X";"Density|f(x)")],
                        .qp.s.geom[`size`colour`strokewidth!(3;`steelblue;15)]
           
                )
        )
    };

.mc.plot.fitHist.go :.mc.plot.fitHist.fn[0b;;;;;0b];
.mc.plot.fitHist.png:.mc.plot.fitHist.fn[1b];