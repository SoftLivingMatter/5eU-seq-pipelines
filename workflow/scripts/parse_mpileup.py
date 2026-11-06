import argparse
import gzip
import re


def parse_args():
    parser = argparse.ArgumentParser(
        description='Parse samtools mpileup output for G→T transversion frequencies'
    )
    parser.add_argument('--input', required=True, help='Input pileup file (gzipped)')
    parser.add_argument('--output', required=True, help='Output TSV file')
    return parser.parse_args()


INDEL_RE = re.compile(r'[+-](\d+)')


def clean_bases(bases):
    """Remove read-start markers (^<mapq>), read-end markers ($) and indels
    ([+-]<len><seq>).  The indel sequence is exactly <len> characters, so it
    is skipped by length rather than by regex to avoid consuming the base
    of the next read."""
    out = []
    i = 0
    n = len(bases)
    while i < n:
        c = bases[i]
        if c == '^':
            i += 2
        elif c == '$':
            i += 1
        elif c in '+-':
            m = INDEL_RE.match(bases, i)
            i = m.end() + int(m.group(1)) if m else i + 1
        else:
            out.append(c)
            i += 1
    return ''.join(out)


def parse_pileup(pileup_path, output_path):
    opener = gzip.open if pileup_path.endswith('.gz') else open

    with opener(pileup_path, 'rt') as fin, open(output_path, 'w') as fout:
        fout.write('\t'.join(['chr', 'pos', 'ref', 'depth', 'T_count', 'G_to_T_freq']) + '\n')

        for line in fin:
            fields = line.strip().split('\t')
            if len(fields) < 6:
                continue
            chrom, pos, ref, depth, bases, quals = fields[:6]

            if ref.upper() != 'G':
                continue

            depth = int(depth)

            bases = clean_bases(bases).upper()

            matches = bases.count('.') + bases.count(',')
            T = bases.count('T')
            total = matches + T + bases.count('A') + bases.count('C')

            if total == 0:
                continue

            fout.write(f'{chrom}\t{pos}\t{ref}\t{depth}\t{T}\t{T/total:.4f}\n')


def main():
    args = parse_args()
    parse_pileup(args.input, args.output)


if __name__ == '__main__':
    main()
