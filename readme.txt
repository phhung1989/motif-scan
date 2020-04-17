The motif scan in writtern by the python version 3.8.1


usage: python motif_scan_server.py <promoter sequence file> <TF cutoff file> <PWM folder> <output file name - species/strain name> <output path> <cds length> 
example : python motif_scan_server.py test_promoter.txt cutoffs.txt PWM SC ./ 300


<promoter sequence file>: The promoter sequence file should be fasta foramt. The example file is test_protmoter.txt.
<TF cutoff file>: The TF cutoff (cutoffs.txt) was downloaded for the ScerTF database. If you want to use the personalized cutoff. Please follow the same format as ScerTF database.
<PWM folder>: The PWMs were downloaded for the ScerTF database. If you want to use the personalized PWM. Please follow the same format as ScerTF database.
<output file name - species/strain name>: The strain or species name you want to write in the output file. The output file will be named as <output file name - species/strain name>_motifscan.txt.
<output path> : output path, the output file will be put in the assigned folder
<cds length> : Since budding yeast has samll 5'UTR, the promoter region here is defined by the upstream 1K from the translation starts site(ATG) and also we consider the 300bp downstram of the translation start site.
               Thus, the cds length is 300bp. If you do not include the coding region, the value will be 0.
               
               
The output file contains the following information: Genename, TF, start site, stop site, motif seq, score, strand.	

The motif hit start site and stop site are the relative position to the translation start site.
If it is a negative value, the motif is locate in the region upstream of the translation start site. (ex: -1 is the first nucleotide upstream of ATG)
If it is a positive value, the motif is locate in the region downstream of the translation start site. (ex: +1 is the first nucleotide in the coding region. Thus, theoretically will be "A", due to the start condon, ATG.)
