// geno2mdb -- convert a GEMMA BIMBAM file into the lmdb genotype store pangemma
// reads with `-g X.mdb`.
//
// The Odin counterpart of bin/geno2mdb.py, which it has to match *byte for
// byte*: the store layout is what pangemma reads, so a different packing is a
// different answer. Why a second implementation at all: the python one needs an
// interpreter plus numpy and python-lmdb at run time, and in this deployment
// that meant exporting PYTHONPATH by hand on the API host and carrying three
// python packages in the ravanan job's manifest for one converter. This is one
// binary and liblmdb.
//
// The format, unchanged from the python and from geno2mdb.rb before it:
//
//   the lmdb environment is a single file (MDB_NOSUBDIR)
//   "geno"    key   chr, pos, line as "S>L>L>" (uint16, uint32, uint32 1-based)
//             value Gf: one little-endian float32 per sample, missing = NaN
//                   Gb: one byte per sample, round(g*127) clipped to 0..254,
//                       missing = 255
//   "marker"  key   the same key -> "name<TAB>allele1<TAB>allele0"
//   "info"    "format" -> Gf|Gb, "numsamples"/"nummarkers" (uint64 LE), "meta"
//
// The chromosome is a number, or the ASCII code of the letter for X/Y/M
// (pangemma's CHR_X is 'X', not 23), and it is part of the key because that is
// how LOCO walks the store: skipping the chromosome in the key. Keys are
// therefore written in chromosome and position order, not input order.
package geno2mdb

import "core:fmt"
import "core:math"
import "core:mem"
import "core:os"
import "core:slice"
import "core:strconv"
import "core:strings"

CHRPOS_NAME :: "S>L>L>"

Row :: struct {
	chr:     u16,
	pos:     u32,
	name:    string,
	allele1: string,
	allele0: string,
	values:  []f32,
	missing: int,
}

Snp :: struct {
	chr: u16,
	pos: u32,
}

Format :: enum {
	Gf,
	Gb,
}

chr_code :: proc(name: string) -> (u16, bool) {
	if code, ok := strconv.parse_u64(name); ok && code <= 65535 {
		return u16(code), true
	}
	upper := strings.to_upper(name)
	switch upper {
	case "X":
		return 88, true // 'X'
	case "Y":
		return 89, true // 'Y'
	case "M", "MT":
		return 77, true // 'M'
	}
	return 0, false
}

// fields_of splits a line on commas *and* whitespace, in place: the annotation
// comes either way ("marker,pos,chr" from R, "marker pos chr" from BXD_snps.txt).
fields_of :: proc(line: string, out: []string) -> int {
	is_sep :: proc(c: u8) -> bool {
		return c == ',' || c == ' ' || c == '\t' || c == '\r'
	}
	n, i := 0, 0
	for i < len(line) {
		for i < len(line) && is_sep(line[i]) {
			i += 1
		}
		start := i
		for i < len(line) && !is_sep(line[i]) {
			i += 1
		}
		if i > start {
			if n >= len(out) {
				break
			}
			out[n] = line[start:i]
			n += 1
		}
	}
	return n
}

// Lines walks a file in fixed-size chunks so a multi-GB BIMBAM is never held in
// memory; the returned slice is valid until the next call.
Lines :: struct {
	file:  ^os.File,
	buf:   []u8,
	start: int,
	end:   int,
	eof:   bool,
}

lines_open :: proc(path: string, size: int) -> (Lines, bool) {
	file, err := os.open(path)
	if err != nil {
		return {}, false
	}
	return Lines{file = file, buf = make([]u8, size)}, true
}

lines_close :: proc(l: ^Lines) {
	os.close(l.file)
	delete(l.buf)
	l.buf = nil
}

lines_next :: proc(l: ^Lines) -> (line: string, ok: bool) {
	for {
		// Look for a newline in what is buffered.
		for i in l.start ..< l.end {
			if l.buf[i] == '\n' {
				line = string(l.buf[l.start:i])
				l.start = i + 1
				// Strip a trailing CR so DOS-terminated files parse the same.
				if len(line) > 0 && line[len(line) - 1] == '\r' {
					line = line[:len(line) - 1]
				}
				return line, true
			}
		}
		if l.eof {
			if l.start < l.end {
				line = string(l.buf[l.start:l.end])
				l.start = l.end
				return line, true
			}
			return "", false
		}
		// No newline: keep the remainder, read more. `remaining` can be 0 when
		// the buffer ended exactly on a line boundary, and &buf[len] is out of
		// bounds even for a zero-length copy (measured: bounds-check panic on a
		// 1 MiB buffer that filled exactly).
		remaining := l.end - l.start
		if l.start > 0 {
			if remaining > 0 {
				mem.copy(&l.buf[0], &l.buf[l.start], remaining)
			}
			l.start = 0
			l.end = remaining
		}
		if l.end == len(l.buf) {
			return "", false // a line longer than the buffer; the caller reports it
		}
		n, err := os.read(l.file, l.buf[l.end:])
		if err != nil || n == 0 {
			l.eof = true
		} else {
			l.end += n
		}
	}
}

read_snps :: proc(path: string) -> (map[string]Snp, int, int) {
	lines, ok := lines_open(path, 1 << 20)
	if !ok {
		fmt.eprintfln("geno2mdb: cannot read %s", path)
		os.exit(1)
	}
	defer lines_close(&lines)

	snps := make(map[string]Snp, 1 << 18)
	fields := make([]string, 8)
	defer delete(fields)
	unknown := 0
	lineno := 0
	for {
		line, more := lines_next(&lines)
		if !more {
			break
		}
		lineno += 1
		n := fields_of(line, fields)
		if n < 3 {
			continue
		}
		pos_float, parsed := strconv.parse_f64(fields[1])
		if !parsed {
			fmt.eprintfln("geno2mdb: %s:%d: expected 'rs pos chr', got %q", path, lineno, line)
			os.exit(1)
		}
		code, known := chr_code(fields[2])
		if !known {
			unknown += 1
			continue
		}
		snps[strings.clone(fields[0])] = Snp{chr = code, pos = u32(math.round(pos_float))}
	}
	if len(snps) == 0 {
		fmt.eprintfln("geno2mdb: %s: no usable markers -- expected 'rs pos chr' per line", path)
		os.exit(1)
	}
	if unknown > 0 {
		fmt.eprintfln(
			"WARNING: %d markers in %s have a chromosome this store cannot key (only numbers and X/Y/M)",
			unknown,
			path,
		)
	}
	return snps, len(snps), unknown
}

// round_half_even matches numpy's rint, which the python converter uses: half
// way cases go to the even neighbour, not away from zero.
round_half_even :: proc(x: f64) -> f64 {
	f := math.floor(x)
	diff := x - f
	if diff < 0.5 {
		return f
	}
	if diff > 0.5 {
		return f + 1
	}
	return f if i64(f) % 2 == 0 else f + 1
}

put_key :: proc(key: ^[10]u8, chr: u16, pos: u32, line: u32) {
	key[0] = u8(chr >> 8)
	key[1] = u8(chr & 0xff)
	key[2] = u8(pos >> 24)
	key[3] = u8((pos >> 16) & 0xff)
	key[4] = u8((pos >> 8) & 0xff)
	key[5] = u8(pos & 0xff)
	key[6] = u8(line >> 24)
	key[7] = u8((line >> 16) & 0xff)
	key[8] = u8((line >> 8) & 0xff)
	key[9] = u8(line & 0xff)
}

usage :: proc() -> ! {
	fmt.eprintln(
		"usage: geno2mdb --snps FILE [--format Gf|Gb] [--out FILE] [--no-sort] BIMBAM...\n" +
		"\nConverts a BIMBAM mean-genotype file to the lmdb store pangemma reads (-g X.mdb).\n" +
		"--snps is the annotation, 'rs pos chr' per line. --out defaults to <input>.mdb.",
	)
	os.exit(2)
}

main :: proc() {
	files: [dynamic]string
	snps_path: string
	out_path: string
	format := Format.Gf
	no_sort := false

	args := os.args
	i := 1
	for i < len(args) {
		arg := args[i]
		switch {
		case arg == "--snps":
			i += 1
			if i >= len(args) { usage() }
			snps_path = args[i]
		case strings.has_prefix(arg, "--snps="):
			snps_path = arg[len("--snps="):]
		case arg == "--format":
			i += 1
			if i >= len(args) { usage() }
			switch args[i] {
			case "Gf":
				format = .Gf
			case "Gb":
				format = .Gb
			case:
				fmt.eprintfln("geno2mdb: --format must be Gf or Gb")
				os.exit(2)
			}
		case arg == "--out":
			i += 1
			if i >= len(args) { usage() }
			out_path = args[i]
		case arg == "--no-sort":
			no_sort = true
		case arg == "--geno-json":
			fmt.eprintln(
				"geno2mdb: --geno-json is not implemented in this build; the store itself does not need it",
			)
			os.exit(2)
		case arg == "-h" || arg == "--help":
			usage()
		case strings.has_prefix(arg, "--"):
			fmt.eprintfln("geno2mdb: unknown option %s", arg)
			os.exit(2)
		case:
			append(&files, arg)
		}
		i += 1
	}
	if snps_path == "" || len(files) == 0 {
		usage()
	}

	snps, _, _ := read_snps(snps_path)
	defer delete(snps)

	for file in files {
		out := out_path
		if out == "" {
			out = strings.concatenate({file, ".mdb"})
		}
		convert(file, out, snps, format, !no_sort)
		if out_path != "" {
			break // --out names one store
		}
	}
}

convert :: proc(file, out: string, snps: map[string]Snp, format: Format, sort_rows: bool) {
	if strings.has_suffix(file, ".gz") {
		// The python version read gzip through zlib; this build does not link
		// zlib yet, and reading a gzip stream as text would produce garbage
		// rather than fail, so say so.
		fmt.eprintfln(
			"geno2mdb: %s is gzipped, which this build does not read -- gunzip it first (the .geno input path will not need this)",
			file,
		)
		os.exit(2)
	}
	lines, ok := lines_open(file, 8 << 20)
	if !ok {
		fmt.eprintfln("geno2mdb: cannot read %s", file)
		os.exit(1)
	}
	defer lines_close(&lines)

	rows: [dynamic]Row
	defer {
		for row in rows {
			delete(row.name)
			delete(row.allele1)
			delete(row.allele0)
			delete(row.values)
		}
		delete(rows)
	}

	fields := make([]string, 0)
	defer delete(fields)
	nsamples := -1
	not_in_annotation := 0
	line_no := 0
	for {
		line, more := lines_next(&lines)
		if !more {
			break
		}
		line_no += 1
		if len(strings.trim_space(line)) == 0 {
			continue
		}
		// Grow the field slice to whatever this line needs; BIMBAM lines are
		// the same width, so this happens at most twice.
		needed := strings.count(line, ",") + strings.count(line, " ") + 4
		if needed > len(fields) {
			delete(fields)
			fields = make([]string, needed)
		}
		n := fields_of(line, fields)
		if n < 4 {
			continue
		}
		count := n - 3
		if nsamples < 0 {
			nsamples = count
		} else if count != nsamples {
			fmt.eprintfln(
				"geno2mdb: %s:%d: %d genotypes, expected %d (line 1)",
				file,
				line_no,
				count,
				nsamples,
			)
			os.exit(1)
		}
		name := fields[0]
		snp, known := snps[name]
		if !known {
			not_in_annotation += 1
			continue
		}
		values := make([]f32, nsamples)
		missing := 0
		for field, k in fields[3:n] {
			if field == "NA" || field == "nan" {
				values[k] = math.nan_f32()
				missing += 1
				continue
			}
			v, parsed := strconv.parse_f32(field)
			if !parsed {
				fmt.eprintfln("geno2mdb: %s:%d: %q is not a genotype", file, line_no, field)
				os.exit(1)
			}
			values[k] = v
		}
		append(
			&rows,
			Row {
				chr = snp.chr,
				pos = snp.pos,
				name = strings.clone(name),
				allele1 = strings.clone(fields[1]),
				allele0 = strings.clone(fields[2]),
				values = values,
				missing = missing,
			},
		)
	}
	if len(rows) == 0 {
		fmt.eprintfln("geno2mdb: %s: no usable markers", file)
		os.exit(1)
	}
	if not_in_annotation > 0 {
		fmt.eprintfln(
			"WARNING: %d markers in %s are not in the annotation, skipped (a marker with no chromosome would be invisible to LOCO)",
			not_in_annotation,
			file,
		)
	}
	if sort_rows {
		slice.stable_sort_by(rows[:], proc(a, b: Row) -> bool {
			if a.chr != b.chr {
				return a.chr < b.chr
			}
			return a.pos < b.pos
		})
	}

	total_missing := 0
	for row in rows {
		total_missing += row.missing
	}

	env: ^MDB_env
	if rc := env_create(&env); rc != 0 { die("mdb_env_create", rc) }
	if rc := env_set_maxdbs(env, 4); rc != 0 { die("mdb_env_set_maxdbs", rc) }
	if rc := env_set_mapsize(env, 8 << 30); rc != 0 { die("mdb_env_set_mapsize", rc) }
	if rc := env_open(env, fmt.ctprint(out), ENV_NOSUBDIR, 0o644); rc != 0 {
		die(fmt.tprintf("mdb_env_open %s", out), rc)
	}
	txn: ^MDB_txn
	if rc := txn_begin(env, nil, 0, &txn); rc != 0 { die("mdb_txn_begin", rc) }
	geno_db, marker_db, info_db: MDB_dbi
	if rc := dbi_open(txn, "geno", DBI_CREATE, &geno_db); rc != 0 { die("mdb_dbi_open geno", rc) }
	if rc := dbi_open(txn, "marker", DBI_CREATE, &marker_db); rc != 0 {
		die("mdb_dbi_open marker", rc)
	}
	if rc := dbi_open(txn, "info", DBI_CREATE, &info_db); rc != 0 { die("mdb_dbi_open info", rc) }

	key: [10]u8
	quantised := make([]u8, nsamples)
	value: []u8
	defer delete(quantised)
	for row, index in rows {
		// Odin wants declarations before statements in a block, hence `value`
		// above the loop.
		put_key(&key, row.chr, row.pos, u32(index + 1))
		key_val := MDB_val{size = len(key), data = raw_data(key[:])}
		if format == .Gb {
			for v, k in row.values {
				if math.is_nan(v) {
					quantised[k] = 255
					continue
				}
				q := round_half_even(f64(v) * 127.0)
				q = math.clamp(q, 0, 254)
				quantised[k] = u8(q)
			}
			value = quantised
		} else {
			value = mem.slice_to_bytes(row.values)
		}
		data_val := MDB_val{size = len(value), data = raw_data(value)}
		if rc := put(txn, geno_db, &key_val, &data_val, 0); rc != 0 { die("mdb_put geno", rc) }

		marker_text := strings.concatenate({row.name, "\t", row.allele1, "\t", row.allele0})
		defer delete(marker_text)
		marker_val := MDB_val{size = len(marker_text), data = raw_byte_data(marker_text)}
		if rc := put(txn, marker_db, &key_val, &marker_val, 0); rc != 0 {
			die("mdb_put marker", rc)
		}
	}

	put_info :: proc(txn: ^MDB_txn, db: MDB_dbi, name: string, bytes: []u8) {
		key_val := MDB_val{size = len(name), data = raw_byte_data(name)}
		data_val := MDB_val{size = len(bytes), data = raw_data(bytes)}
		if rc := put(txn, db, &key_val, &data_val, 0); rc != 0 { die("mdb_put info", rc) }
	}

	format_name := format == .Gb ? "Gb" : "Gf"
	format_bytes := transmute([]u8)format_name
	put_info(txn, info_db, "format", format_bytes)

	u64_bytes :: proc(n: u64) -> (out: [8]u8) {
		for i in 0 ..< 8 {
			out[i] = u8(n >> (8 * u32(i)))
		}
		return
	}
	ns := u64_bytes(u64(nsamples))
	nm := u64_bytes(u64(len(rows)))
	put_info(txn, info_db, "numsamples", ns[:])
	put_info(txn, info_db, "nummarkers", nm[:])

	rec := format == .Gb ? "C*" : "f*"
	// Built by concatenation, not printf: Odin's format strings treat { and } as
	// placeholders, and this is a JSON document (measured: "%!(MISSING CLOSE
	// BRACE)type": ..." was the only difference between this store and the
	// python one).
	meta := strings.concatenate(
		{
			`{"type": "gemma-geno", "format": "`,
			format_name,
			`", "version": 1.0, "key-format": "`,
			CHRPOS_NAME,
			`", "rec-format": "`,
			rec,
			`"}`,
		},
	)
	defer delete(meta)
	put_info(txn, info_db, "meta", transmute([]u8)meta)

	if rc := txn_commit(txn); rc != 0 { die("mdb_txn_commit", rc) }
	if rc := env_sync(env, 1); rc != 0 { die("mdb_env_sync", rc) }
	env_close(env)

	pct := 0.0
	if nsamples > 0 && len(rows) > 0 {
		pct = 100.0 * f64(total_missing) / f64(len(rows) * nsamples)
	}
	fmt.printfln(
		"%d markers x %d samples, format %s, %d missing values (%.1f%%) -> %s",
		len(rows),
		nsamples,
		format_name,
		total_missing,
		pct,
		out,
	)
}

// raw_byte_data gives an immutable byte view of a string's storage for lmdb.
raw_byte_data :: proc(s: string) -> rawptr {
	return rawptr(uintptr(raw_data(s)))
}
