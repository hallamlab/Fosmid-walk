# Fosmid Walk
Estimates population size of fosmid pools

```
singularity build foswalk.sif docker://quay.io/hallam_lab/fosmid-walk

singularity run -B <workspace>:/ws foswalk.sif
foswalk -r /ws/<fosmid_read_file>.fq -o /ws/<output_folder>

```

[Singularity docs](https://singularity-userdoc.readthedocs.io/en/latest/quick_start.html)
