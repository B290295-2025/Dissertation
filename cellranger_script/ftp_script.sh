# remotely download fast fils
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ERP146275&result=read_run&fields=fastq_ftp&format=tsv" \
| tail -n +2 | tr ';' '\n' > download_list.txt

cat download_list.txt

nohup wget -i download_list.txt -c > download.log 2>&1 &

tail -f download.log
