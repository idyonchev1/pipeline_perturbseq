rule aggr:
    input:expand("{sample}/outs/web_summary.html",sample=samples)
    output: "aggr/outs/web_summary.html"
    log: "log/cellranger_aggr.log"
    resources:
        mem_mb=128000,
        disk_mb=200000,
        time="24:00:00"
    threads: 16
    run:
        cmd1="rm -r aggr"
        cmd2="cellranger aggr --id=aggr --normalize=none --csv=aggr.csv &> {log}"
        shell(cmd1)
        shell(cmd2)
