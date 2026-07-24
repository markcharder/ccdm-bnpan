#!/usr/bin/env python3

import argparse
import cyvcf2
import pysam
import re
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--chromids", required=False, default=None)
    parser.add_argument("--paf", required=True)
    parser.add_argument("--bcf", required=True)
    parser.add_argument("--fasta", required=True, help="Indexed reference FASTA")
    parser.add_argument(
        "--flank", type=int,
        default=10_000, help="Flanking region around gene"
    )
    return parser.parse_args()

def iter_primary_alignments(paf, chromids=None):

    with open(paf) as f:

        for line in f:

            if line.startswith("#"):
                continue

            parts = line.strip().split("\t")
            tags = {t[:4]: t for t in parts[12:]}

            if tags.get("tp:A") != "tp:A:P":
                continue

            if chromids:

                if parts[5] not in chromids:
                    continue
                
                else:
                    parts[5] = chromids[parts[5]]

            yield {
                "query": parts[0],
                "strand": parts[4],
                "chrom": parts[5],
                "target_start": int(parts[7]),
                "target_end": int(parts[8]),
                "cs": tags.get("cs:Z") or tags.get("cs:z"),
            }

def parse_cs_to_exons(rec):
    cs = rec.get("cs", "")

    if not cs:
        return []

    cs_str = cs.split(":", 1)[-1]

    exons, pos = [], rec["target_start"]

    for op, val in re.findall(r"([=:*+-])([a-zA-Z0-9]+)", cs_str):
        length = len(val) if op == "=" else (int(val) if op in (":", "+", "-") else 1)

        if op in ("=", ":", "*"):
            exons.append((pos, pos + length))
            pos += length

        elif op == "-":
            pos += length

    return exons

def classify_pos(pos, gene_start, gene_end, exons):

    if pos < gene_start:
        return "upstream"

    if pos >= gene_end:
        return "downstream"

    for ex_s, ex_e in exons:

        if ex_s <= pos < ex_e:
            return "CDS/exon"

    return "intron"

def build_alignment(ref_seq, region_start, all_hap_vars):
    ref_at = {}
    max_ins = {}

    for hap_vars in all_hap_vars:

        for pos, ref, alt in hap_vars:
            ref_at[pos] = ref
            ins = max(0, len(alt) - len(ref))
            max_ins[pos] = max(max_ins.get(pos, 0), ins)

    var_positions = sorted(ref_at.keys())

    def _ref_seq():
        result, cur = [], 0

        for pos in var_positions:
            offset = pos - region_start
            ref_allele = ref_at[pos]
            result.append(ref_seq[cur:offset])
            result.append(ref_allele.upper())
            result.append("-" * max_ins.get(pos, 0))
            cur = offset + len(ref_allele)

        result.append(ref_seq[cur:])
        return "".join(result)

    def _hap_seq(hap_vars):
        var_lookup = {pos: (ref, alt) for pos, ref, alt in hap_vars}
        result, cur = [], 0

        for pos in var_positions:
            offset = pos - region_start
            ref_allele = ref_at[pos]
            ins_pad = max_ins.get(pos, 0)
            result.append(ref_seq[cur:offset])

            if pos in var_lookup:
                _, alt = var_lookup[pos]
                allele_str = alt.lower() if alt != ref_allele else ref_allele

            else:
                allele_str = ref_allele

            result.append(allele_str)
            del_gap  = max(0, len(ref_allele) - len(allele_str))
            this_ins = max(0, len(allele_str) - len(ref_allele))
            ins_gap  = ins_pad - this_ins
            result.append("-" * (del_gap + ins_gap))
            cur = offset + len(ref_allele)

        result.append(ref_seq[cur:])
        return "".join(result)
    return [_ref_seq()] + [_hap_seq(hv) for hv in all_hap_vars]

def main():
    args = parse_args()
    fasta = pysam.FastaFile(args.fasta)
    vcf = cyvcf2.VCF(args.bcf, gts012=True)
    samples = list(vcf.samples)
    chromids = None

    if args.chromids:
        chromids = {}
        chromids_rev = {}

        with open(args.chromids) as f:

            for line in f:
            
                if line.startswith("#"):
                    continue

                parts = line.strip().split("\t")

                if len(parts) < 2:
                    continue

                chromids[parts[0]] = parts[1]
                chromids_rev[parts[1]] = parts[0]

    for rec in iter_primary_alignments(args.paf, chromids=chromids):
        query_id = rec["query"]
        chrom = rec["chrom"]
        gene_start = rec["target_start"]
        gene_end = rec["target_end"]
        exons = parse_cs_to_exons(rec)
        region_start = max(0, gene_start - args.flank)
        region_end = gene_end + args.flank

        print(f"\n## {query_id}")
        print(f"# Gene   : {chrom}:{gene_start+1}-{gene_end}  strand={rec['strand']}")
        print(f"# Region : {chrom}:{region_start+1}-{region_end}")

        if exons:
            print(f"# Exons  : {exons}")
        print()

        if chromids:
            chrom_fasta = chromids_rev.get(chrom, chrom)

        else:
            chrom_fasta = chrom

        ref_seq = fasta.fetch(chrom_fasta, region_start, region_end)

        # Reset per-gene accumulators.
        hap_vars = [[[], []] for _ in samples]
        var_summary = []
        region_str = f"{chrom}:{region_start}-{region_end}"

        for var in vcf(region_str):
            pos = var.POS - 1
            ref = var.REF
            alts = var.ALT or ["."]
            loc = classify_pos(pos, gene_start, gene_end, exons)
            var_summary.append((pos, ref, alts, loc))

            for si, sample in enumerate(samples):
                gt = var.genotypes[si]

                for hi, allele_idx in enumerate(gt[:2]):

                    if allele_idx is None or allele_idx < 0:
                        alt = ref

                    elif allele_idx == 0:
                        alt = ref

                    else:
                        alt = alts[allele_idx - 1]

                    hap_vars[si][hi].append((pos, ref, alt))

        print("## Variants")
        print(f"{'POS(1-based)':>13}  {'REF':<8}  {'ALT':<20}  LOC")
        print("-" * 58)

        for pos, ref, alts, loc in var_summary:
            print(f"{pos+1:>13}  {ref[:8]:<8}  {','.join(alts)[:20]:<20}  {loc}")
        print()

        out_fa = f"haplotypes_{query_id}.fa"
        print("## Haplotype sequences  (non-ref bases in lowercase)")
        print(f"  NOTE: region is {len(ref_seq)} bp — writing full sequences to {out_fa}")
        print()

        all_hap_vars = [hap_vars[si][hi] for si in range(len(samples)) for hi in [0, 1]]
        aligned = build_alignment(ref_seq, region_start, all_hap_vars)
        labels  = [f"REF_{chrom}_{region_start+1}_{region_end}"] + [
            f"{query_id}__{sample}_hap{hi+1}"
            for si, sample in enumerate(samples)
            for hi in [0, 1]
        ]

        with open(out_fa, "w") as out:
            for label, seq in zip(labels, aligned):
                out.write(f">{label}\n{seq}\n")

        if var_summary:
            positions = [p for p, *_ in var_summary]
            header = " ".join(f"{p+1:>8}" for p in positions)
            print(f"{'':>20}  {header}")
            ref_row = " ".join(f"{r:>8}" for _, r, _, _ in var_summary)
            print(f"{'REF':>20}  {ref_row}")

            for si, sample in enumerate(samples):

                for hi in [0, 1]:
                    vd  = {p: a for p, r, a in hap_vars[si][hi]}
                    row = " ".join(f"{vd.get(p, r):>8}" for p, r, _, _ in var_summary)
                    print(f"{sample+'_h'+str(hi+1):>20}  {row}")

        print(f"\n# Full sequences written to {out_fa}")

    vcf.close()
    fasta.close()

if __name__ == "__main__":
    main()