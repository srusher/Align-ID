#!/usr/bin/env python3
"""
trim_fastq_prefixes.py

Trims reads in a FASTQ file (and the corresponding PHRED scores). Two
trimming modes are supported:

  1. Sequence-match mode (--sequences FILE):
     Trims a read if it begins with one of a list of sequences supplied
     in a separate text file.
       - Exact prefix match only (no mismatches allowed).
       - Each read is trimmed at most once, using the FIRST sequence in
         the list that matches the start of the read (list order =
         priority).
       - Duplicate sequences in the list are ignored (first occurrence
         wins) with a warning.

  2. Fixed-length mode (--headcrop N and/or --tailcrop N):
     Unconditionally trims N bases/scores from the start (--headcrop)
     and/or end (--tailcrop) of every processed read. Both may be given
     together to trim from both ends in the same pass. At least one of
     --headcrop/--tailcrop must be > 0 if fixed-length mode is used.

Exactly one trimming mode must be selected: either --sequences, or one
or both of --headcrop/--tailcrop.

You can also restrict processing to forward reads, reverse reads, both,
or single-end data, since paired FASTQ is commonly provided as two
separate files (R1 = forward, R2 = reverse) rather than interleaved.

Usage examples:

    # Sequence-match trimming on a forward-read file only
    python trim_fastq_prefixes.py \\
        --fastq R1.fastq --read-type forward \\
        --sequences adapters.txt --output R1.trimmed.fastq

    # Headcrop 8bp, reverse reads
    python trim_fastq_prefixes.py \\
        --fastq R2.fastq --read-type reverse \\
        --headcrop 8 --output R2.trimmed.fastq

    # Headcrop 8bp and tailcrop 5bp together, no direction filtering
    python trim_fastq_prefixes.py \\
        --fastq reads.fastq --read-type both \\
        --headcrop 8 --tailcrop 5 --output reads.trimmed.fastq

    # Single-end data, sequence-match trimming
    python trim_fastq_prefixes.py \\
        --fastq se_reads.fastq --read-type single \\
        --sequences adapters.txt --output se_reads.trimmed.fastq

    # Single-end data, headcrop only
    python trim_fastq_prefixes.py \\
        --fastq se_reads.fastq --read-type single \\
        --headcrop 8 --output se_reads.trimmed.fastq

Supports gzipped input/output automatically based on the .gz extension.
"""

import argparse
import gzip
import sys


def open_maybe_gzip(path, mode):
    """Open a file transparently handling gzip based on extension."""
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t" if "b" not in mode else mode)
    return open(path, mode)


def load_sequences(path):
    """
    Load trim sequences from a text file, one per line.
    Preserves file order (used as match priority).
    Skips blank lines and lines starting with '#'.
    """
    sequences = []
    seen = set()
    with open(path, "r") as f:
        for line_num, line in enumerate(f, start=1):
            seq = line.strip()
            if not seq or seq.startswith("#"):
                continue
            seq_upper = seq.upper()
            if seq_upper in seen:
                print(
                    f"Warning: duplicate sequence '{seq_upper}' found on line "
                    f"{line_num} of {path} — ignoring duplicate entry.",
                    file=sys.stderr,
                )
                continue
            seen.add(seq_upper)
            sequences.append(seq_upper)

    if not sequences:
        raise ValueError(f"No valid sequences found in {path}")

    return sequences


def find_matching_prefix(read_seq, trim_sequences):
    """
    Return the first sequence (in list order) that the read starts with,
    or None if no match is found.
    """
    read_upper = read_seq.upper()
    for trim_seq in trim_sequences:
        if read_upper.startswith(trim_seq):
            return trim_seq
    return None


def fastq_records(fh):
    """
    Generator that yields (header, seq, plus_line, qual) tuples from an
    open FASTQ file handle. Raises ValueError on malformed records.
    """
    while True:
        header = fh.readline()
        if not header:
            break  # EOF
        seq = fh.readline()
        plus_line = fh.readline()
        qual = fh.readline()

        if not seq or not plus_line or not qual:
            raise ValueError("Malformed FASTQ: incomplete record at end of file")

        header = header.rstrip("\n")
        seq = seq.rstrip("\n")
        plus_line = plus_line.rstrip("\n")
        qual = qual.rstrip("\n")

        if not header.startswith("@"):
            raise ValueError(f"Malformed FASTQ: expected '@' header, got: {header}")
        if not plus_line.startswith("+"):
            raise ValueError(f"Malformed FASTQ: expected '+' line, got: {plus_line}")
        if len(seq) != len(qual):
            raise ValueError(
                f"Malformed FASTQ: SEQ/QUAL length mismatch for read '{header}' "
                f"({len(seq)} vs {len(qual)})"
            )

        yield header, seq, plus_line, qual


def detect_read_direction(header):
    """
    Best-effort detection of forward (1) vs reverse (2) read from a FASTQ
    header, supporting common Illumina conventions:
        - Old style:  .../1 or .../2 at end of header
        - New style:  contains " 1:" or " 2:" (e.g. " 1:N:0:ATCG")
    Returns '1', '2', or None if undetermined.
    """
    header = header.strip()

    # New Illumina style: "@INSTR:RUN:FLOWCELL:LANE:TILE:X:Y 1:N:0:INDEX"
    if " 1:" in header:
        return "1"
    if " 2:" in header:
        return "2"

    # Old style: "@READID/1" or "@READID/2"
    if header.endswith("/1"):
        return "1"
    if header.endswith("/2"):
        return "2"

    return None


def apply_crop(seq, qual, headcrop, tailcrop):
    """
    Trim headcrop bases from the start and tailcrop bases from the end
    of seq/qual. Returns (new_seq, new_qual, was_trimmed).
    Guards against over-cropping past the end of the read.
    """
    read_len = len(seq)
    start = min(headcrop, read_len)
    end = max(read_len - tailcrop, start)  # never let end < start

    new_seq = seq[start:end]
    new_qual = qual[start:end]
    was_trimmed = (start > 0) or (end < read_len)
    return new_seq, new_qual, was_trimmed


def trim_fastq(fastq_path, output_path, read_type, trim_sequences=None,
               headcrop=0, tailcrop=0):
    total_reads = 0
    processed_reads = 0
    trimmed_reads = 0
    skipped_wrong_direction = 0
    undetermined_direction = 0
    trim_counts = {}
    if trim_sequences:
        trim_counts = {seq: 0 for seq in trim_sequences}

    with open_maybe_gzip(fastq_path, "r") as infile, \
         open_maybe_gzip(output_path, "w") as outfile:

        for header, seq, plus_line, qual in fastq_records(infile):
            total_reads += 1

            do_process = True
            if read_type not in ("both", "single"):
                direction = detect_read_direction(header)
                if direction is None:
                    undetermined_direction += 1
                    # Can't confirm direction — process it, but warn once
                    # per file via the summary rather than spamming stderr.
                elif read_type == "forward" and direction != "1":
                    do_process = False
                    skipped_wrong_direction += 1
                elif read_type == "reverse" and direction != "2":
                    do_process = False
                    skipped_wrong_direction += 1

            if do_process:
                processed_reads += 1

                if trim_sequences is not None:
                    match = find_matching_prefix(seq, trim_sequences)
                    if match is not None:
                        trim_len = len(match)
                        seq = seq[trim_len:]
                        qual = qual[trim_len:]
                        trimmed_reads += 1
                        trim_counts[match] += 1
                else:
                    # Fixed-length headcrop/tailcrop mode
                    seq, qual, was_trimmed = apply_crop(seq, qual, headcrop, tailcrop)
                    if was_trimmed:
                        trimmed_reads += 1

            outfile.write(f"{header}\n{seq}\n{plus_line}\n{qual}\n")

    return {
        "total_reads": total_reads,
        "processed_reads": processed_reads,
        "trimmed_reads": trimmed_reads,
        "skipped_wrong_direction": skipped_wrong_direction,
        "undetermined_direction": undetermined_direction,
        "trim_counts": trim_counts,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Trim reads (sequence-match or headcrop/tailcrop) from FASTQ reads.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--fastq", required=True, help="Input FASTQ file (.fastq or .fastq.gz)")
    parser.add_argument("--output", required=True, help="Output FASTQ file path (.fastq or .fastq.gz)")
    parser.add_argument(
        "--read-type", choices=["forward", "reverse", "both", "single"], default="both",
        help=(
            "Which reads to process. 'forward' = read 1 only (by header "
            "orientation), 'reverse' = read 2 only, 'both' = process every "
            "read in a paired file regardless of detected orientation "
            "(default), 'single' = single-end data — skips orientation "
            "detection entirely and processes every read in the file. "
            "Reads whose direction can't be determined from the header "
            "under 'forward'/'reverse' are still processed, with a note "
            "in the summary."
        ),
    )

    trim_group = parser.add_mutually_exclusive_group(required=True)
    trim_group.add_argument(
        "--sequences",
        help="Text file listing sequences to trim as exact prefixes, one per line"
    )
    trim_group.add_argument(
        "--headcrop", type=int, default=0,
        help="Unconditionally trim this many bases/scores from the START of every processed read"
    )

    # tailcrop is separate (not part of the mutually exclusive group) so it
    # can be combined with --headcrop in the same invocation.
    parser.add_argument(
        "--tailcrop", type=int, default=0,
        help="Unconditionally trim this many bases/scores from the END of every processed read"
    )

    args = parser.parse_args()

    if args.headcrop < 0:
        parser.error("--headcrop must be a non-negative integer")
    if args.tailcrop < 0:
        parser.error("--tailcrop must be a non-negative integer")
    if args.sequences and (args.headcrop > 0 or args.tailcrop > 0):
        parser.error("--sequences cannot be combined with --headcrop/--tailcrop; choose one trimming mode")
    if not args.sequences and args.headcrop == 0 and args.tailcrop == 0:
        parser.error("Must specify --sequences, or at least one of --headcrop/--tailcrop greater than 0")

    trim_sequences = None
    if args.sequences:
        trim_sequences = load_sequences(args.sequences)
        print(f"Loaded {len(trim_sequences)} unique sequence(s) to trim.", file=sys.stderr)
    else:
        print(
            f"Fixed-length trim mode: headcrop={args.headcrop}, tailcrop={args.tailcrop}.",
            file=sys.stderr,
        )

    print(f"Read-type filter: {args.read_type}", file=sys.stderr)

    stats = trim_fastq(
        args.fastq,
        args.output,
        read_type=args.read_type,
        trim_sequences=trim_sequences,
        headcrop=args.headcrop,
        tailcrop=args.tailcrop,
    )

    print(f"\nTotal reads in file: {stats['total_reads']}", file=sys.stderr)
    if args.read_type not in ("both", "single"):
        print(f"Reads skipped (wrong direction): {stats['skipped_wrong_direction']}", file=sys.stderr)
        if stats["undetermined_direction"] > 0:
            print(
                f"Warning: {stats['undetermined_direction']} read(s) had no "
                f"detectable /1 or /2 orientation in the header and were "
                f"processed anyway.",
                file=sys.stderr,
            )
    print(f"Reads processed (eligible for trimming): {stats['processed_reads']}", file=sys.stderr)
    denom = stats["processed_reads"] if stats["processed_reads"] > 0 else 1
    print(
        f"Reads trimmed: {stats['trimmed_reads']} "
        f"({stats['trimmed_reads']/denom*100:.2f}% of processed reads)",
        file=sys.stderr,
    )

    if trim_sequences:
        print("\nPer-sequence trim counts:", file=sys.stderr)
        for seq, count in stats["trim_counts"].items():
            if count > 0:
                print(f"  {seq}: {count}", file=sys.stderr)


if __name__ == "__main__":
    main()