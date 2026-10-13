BEGIN {
    OFS = "\t"
    if (manifest == "" || samples == "" || outdir == "") {
        print "split_mpileup_chunks.awk: manifest, samples, and outdir are required" > "/dev/stderr"
        failed = 1
        exit 2
    }

    n_samples = split(samples, sample_order, ",")
    if (n_samples < 1) {
        print "split_mpileup_chunks.awk: no sample names supplied" > "/dev/stderr"
        failed = 1
        exit 2
    }

    n_manifest = 0
    while ((getline manifest_line < manifest) > 0) {
        if (manifest_line ~ /^chunk_id\tchrom\tlength$/) {
            continue
        }
        split(manifest_line, field, "\t")
        if (field[1] == "" || field[2] == "" || field[3] !~ /^[0-9]+$/) {
            print "split_mpileup_chunks.awk: invalid manifest record: " manifest_line > "/dev/stderr"
            failed = 1
            exit 2
        }
        if (field[2] in chrom_to_chunk) {
            print "split_mpileup_chunks.awk: duplicate sequence in manifest: " field[2] > "/dev/stderr"
            failed = 1
            exit 2
        }
        chrom_to_chunk[field[2]] = field[1]
        expected_chunks[field[1]] = 1
        n_manifest++
    }
    close(manifest)

    if (n_manifest == 0) {
        print "split_mpileup_chunks.awk: manifest contains no sequences" > "/dev/stderr"
        failed = 1
        exit 2
    }

    expected_fields = 5 * n_samples
}

NF != expected_fields {
    print "split_mpileup_chunks.awk: expected " expected_fields " pasted fields, got " NF > "/dev/stderr"
    failed = 1
    exit 2
}

!($1 in chrom_to_chunk) {
    print "split_mpileup_chunks.awk: sequence is absent from manifest: " $1 > "/dev/stderr"
    failed = 1
    exit 2
}

{
    # Each sample TSV has chr, pos, ref_base, strand, value. paste therefore lays
    # out five fields per sample. Check the key columns rather than blindly
    # concatenating rows, so a sample ordering or -aa coverage mismatch fails early.
    for (i = 5; i <= NF; i += 5) {
        if ($(i-4) != $1 || $(i-3) != $2 || $(i-2) != $3 || $(i-1) != $4 || $i == "") {
            print "split_mpileup_chunks.awk: sample-key mismatch at " $1 ":" $2 > "/dev/stderr"
            failed = 1
            exit 2
        }
    }

    chunk_id = chrom_to_chunk[$1]
    if (!(chunk_id in seen_chunk)) {
        out_path = outdir "/" chunk_id ".tsv"
        header = "chr\tpos\tref_base\tstrand"
        for (i = 1; i <= n_samples; i++) {
            header = header "\t" sample_order[i]
        }
        print header > out_path
        seen_chunk[chunk_id] = 1
    }

    out_path = outdir "/" chunk_id ".tsv"
    rebuilt = $1 OFS $2 OFS $3 OFS $4
    for (i = 5; i <= NF; i += 5) {
        rebuilt = rebuilt OFS $i
    }
    print rebuilt > out_path
    seen_sequence[$1] = 1
    n_rows++
}

END {
    if (failed) {
        exit 2
    }
    if (n_rows == 0) {
        print "split_mpileup_chunks.awk: no input rows" > "/dev/stderr"
        exit 2
    }
    for (sequence in chrom_to_chunk) {
        if (!(sequence in seen_sequence)) {
            print "split_mpileup_chunks.awk: expected sequence has no rows: " sequence > "/dev/stderr"
            failed = 1
        }
    }
    for (chunk_id in expected_chunks) {
        if (!(chunk_id in seen_chunk)) {
            print "split_mpileup_chunks.awk: expected chunk has no rows: " chunk_id > "/dev/stderr"
            failed = 1
        }
    }
    if (failed) {
        exit 2
    }
}
