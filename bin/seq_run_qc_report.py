#!/usr/bin/env python
import argparse
import pandas as pd
import glob
import os
import time
import numpy as np

def main():
    parser = argparse.ArgumentParser(description="Derive a qc report")
    parser.add_argument("--adapter_trimming", type=str)
    parser.add_argument("--quality_trimming", type=str)
    args = parser.parse_args()
    adapter_trimming = args.adapter_trimming
    quality_trimming = args.quality_trimming

    timestr = time.strftime("%Y%m%d-%H%M%S")

    summary_dict = {}

    #if os.path.isfile("*raw_NanoStats.txt"):
    for raw_read_out in glob.glob("*raw_NanoStats.txt"):
        raw_reads = ()
        line_number = 0
        sample = (os.path.basename(raw_read_out).replace('_raw_NanoStats.txt', ''))
        with open(raw_read_out, 'r') as f:
            for line in f:
                string_to_search = ("number_of_reads")
                line_number += 1
                if string_to_search in line:
                    elements = line.rsplit('\t')
                    raw_reads = int(elements[1].strip())
                    print(raw_reads)

        summary_dict[sample] = [raw_reads]
        print(summary_dict)
        f.close()

    for qt_read_out in glob.glob("*filtered_NanoStats.txt"):
        qt_reads = ()
        line_number = 0 
        sample = (os.path.basename(qt_read_out).replace('_filtered_NanoStats.txt', ''))
        with open(qt_read_out, 'r') as f:
            for line in f:
                string_to_search = ("number_of_reads")
                line_number += 1
                if string_to_search in line:
                    elements = line.rsplit('\t')
                    qt_reads = int(elements[1].strip())
                    print(qt_reads)
            #first_line = next(f)
            #qt_reads = int(first_line[0].strip())
        f.close()
        summary_dict[sample].append(qt_reads)
    
 
    
    else:
        if adapter_trimming == "true" or quality_trimming == "true":
            run_data_df = pd.DataFrame([([k] + v) for k, v in summary_dict.items()], columns=['Sample','raw_reads','quality_filtered_reads'])
            run_data_df['percent_quality_filtered'] = run_data_df['quality_filtered_reads'] / run_data_df['raw_reads'] * 100
            run_data_df['percent_quality_filtered'] = run_data_df['percent_quality_filtered'].apply(lambda x: float("{:.2f}".format(x)))
            run_data_df = run_data_df.sort_values("Sample")
            run_data_df['raw_reads_flag'] = np.where((run_data_df['raw_reads'] < 5000), "Less than 5000 raw reads", "")
            run_data_df['qfiltered_flag'] = np.where((run_data_df['quality_filtered_reads'] < 1000), "Less than 1000 processed reads", "")
            run_data_df["QC_FLAG"] = np.where(
                (run_data_df['quality_filtered_reads'] < 1000),
                "RED",
                np.where(
                    ((run_data_df['raw_reads'] < 5000) &
                    (run_data_df['quality_filtered_reads'] > 1000)),
                    "ORANGE",
                    np.where(
                        ((run_data_df['raw_reads'] > 5000) &
                        (run_data_df['quality_filtered_reads'] > 1000)),
                        "GREEN",
                        ""
                    )
                )
            )
    
    run_data_df.to_csv("run_qc_report_" + timestr + ".txt", index = None, sep="\t")

    summary_table = run_data_df.to_html(index=False).replace('<table border="1" class="dataframe">','<table class="table table-striped">') # use bootstrap styling
    html_string = '''
    <html>
        <head>
            <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
            <style>body{ margin:0 100; background:whitesmoke; }</style>
        </head>
        <body>
            <h1>Quality check report</h1>
            <blockquote>
            <p><b> Definitions:</b></p>
            <p><b> raw_reads </b>= total raw reads sequenced </p>
            <p><b> quality_filtered_reads </b>= total reads left-over after adapter trimming and/or quality filtering </p>
            <p><b> percent_quality_filtered </b>= (quality_filtered_reads/raw_reads x 100)</p>
            </blockquote>
            <!-- *** Section 1 *** --->
            ''' + summary_table + '''
        </body>
    </html>'''

    report = open("run_qc_report_" + timestr + ".html", "w")
    report.write(html_string)
    report.close()

if __name__ == '__main__':
    main()
