import os, sys
import argparse
from pathlib import Path
from .procedure import estimate

LINE = "################################################################"
class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\n%s: error: %s\n' % (self.prog, message))

def main():
    parser = ArgumentParser(
        prog = 'Fosmid Walk',
        description = "v1.0",
        # epilog = 'Text at the bottom of help',
    )

    # parser.add_argument('-1', metavar='FASTQ', help="paried-end fastq reads 1", required=True)
    # parser.add_argument('-t', metavar='INT', type=int, help="threads, default:16", default=16)
    parser.add_argument('-r', metavar='PATH', help="path to fosmid raw read fastq", required=True)
    parser.add_argument('-o', metavar='PATH', help="output folder", required=True)
    parser.add_argument('-t', metavar='INT', help="threads", default=False)
    parser.add_argument(
        '-b', metavar='PATH',
        help="path to vector backbone fasta, defaults to pCC1 sequence",
        default=False,
    )
    parser.add_argument('--temp-dir', metavar='PATH', help="temporary workspace", default=False)

    def check_file(p, ext_wl: list[str]=list()):
        if not os.path.exists(p):
            print(f"path [{p}] does not exist"); return False
        p = Path(p)
        if not p.is_file():
            print(f"path [{p}] is not a file"); return False
        if len(ext_wl)>0 and not any(str(p).endswith(ext) for ext in ext_wl):
            print(f"read file [{p}] doesn't end in one of {ext_wl}"); return False
        return p

    args = parser.parse_args(sys.argv[1:])
    reads = check_file(args.r, ".fastq, .fq".split(', '))
    if not reads: return
    out_dir = args.o
    out_dir = Path(out_dir)

    args_to_pass = dict(sample=reads, out_dir=out_dir)

    backbone = False
    if args.b: backbone = check_file(args.b, ".fasta, .fa".split(', '))
    if backbone: args_to_pass["backbone"] = backbone

    if args.o: args_to_pass["temp_dir"] = Path(args.temp_dir)
    if args.t: args_to_pass["threads"] = args.t

    estimate(**args_to_pass)


    # if not os.path.exists(args.r) or os.path.isfile:
    #     return
    # os.makedirs(args.o, exist_ok=True)
    # ref_path = os.path.abspath(args.r)
    # out_path = os.path.abspath(args.o)

    # subsample = min(1, args.s)

    # cmd = f"""\
    #     snakemake -s /app/main.smk --configfile /app/config.yaml -d /ws\
    #         {"-n" if args.mock else ""} \
    #         --config sample={args.i} ss={subsample} \
    #         --keep-going --keep-incomplete --cores {args.t}"""

    # print_header()
    # if os.path.exists(f'{ref_path}/{IMAGE}.sif'):
    #     success = 0==shell_singularity(ref_path, out_path, cmd)
    # else:
    #     success = 0==shell_docker(ref_path, out_path, cmd)
    # if success and not args.mock:
    #     print(f"{LINE}\nsubsample rate: {args.s}\nfinal output at\n{out_path}/{args.i}/diamond/")
