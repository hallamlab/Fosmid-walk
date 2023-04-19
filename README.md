# Fosmid Walk
Estimates population size of fosmid pools

Credit to [Avery Noonan](https://scholar.google.ca/citations?user=3mTbEjQAAAAJ&hl=en) and [Kateryna Ievdokymenko](https://scholar.google.ca/citations?user=_WzHYD4AAAAJ&hl=en).

```
singularity build foswalk.sif docker://quay.io/hallam_lab/fosmid-walk

singularity run -B <workspace>:/ws foswalk.sif foswalk -r /ws/<fosmid_read_file>.fq -o /ws/<output_folder>

```

[Singularity docs](https://singularity-userdoc.readthedocs.io/en/latest/quick_start.html)
