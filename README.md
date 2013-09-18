HCVtoolbox
==========

Tools and scripts used in various HCV projects


An example of the command used:
````
cat sample2.fastq | perl 454blastParseCleanSW.pl -inblast sample2vsref.blastnout -ref HCVref.fasta -qual testing/QualityConvertTable.txt -qcutoff 25 -out test -minlength 187 -bcfile barcodeAndPrimerS.txt -mismatches 2 -prefix testing/ -fastq sample2.fastq
````
The input reads in fastq format (sample2.fastq) will be split according to the barcodes 
specified in barcodeAndPrimerS.txt allowing for 2 mismatches in the reverse barcode, 
a minimum length of the sequence of 187 and a quality threshold cut-off of 25 (phred score).
sample2vsref.blastnout are the results of a blast against the reference set HCVref.fasta:
blastall -p blastn -d HCVref.fasta -i sample2.fa -o sample2vsref.blastnout -v 1 -b 1 &

The output files would be created in the "testing" directory that needs to be created before running the script. 