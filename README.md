# strand-seq
Switch point caller on strand-seq data

# Requirements
DMD v2.065

# Install
```
git clone --recursive https://github.com/CuppenResearch/strand-seq.git
```

# How to run
```
rdmd -IBioD switch_point_caller_trio.d -p <patient> -f <father> -m <mother>
```
