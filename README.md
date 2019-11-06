# Iwate Biotechnology Research Center (IBRC)
Python stuff from my current (circa 2019) postdoc work into analyzing rice GWAS data sets. 

The notebooks should be fairly self-explanatory. I tried to comment them where appropriate.

Performing genome-wide analyses on 3000+ recombinant inbred lines (RILs) derived from crossing Hitomebore (common parent) with one of 20 founder lines. F7 generation offspring. I use 1:1 Mendelian segregation as the expected and ignore heterozygous loci. Chi-square tests warn that f_obs should be >5 and in an F7 population you're only epxecting 1.5625% of sites to remain heterozygous (assuming Mendelian segregation), so 49.21875 % : 1.5625 % : 49.21875 %.

Mostly statistical inference analysis in these scripts. I tried using scipy.stats and statsmodels since conversion of VCF to pandas dataframes allows me to ssh tunnel to my workstation and run everything in jupyter notebooks.

Within the notebooks markdown cells are included that go over various aspects of what the analysis is trying to accomplish. Where appropriate, I've tried to add some images (basically slides exported as `png` files from Google slides).
