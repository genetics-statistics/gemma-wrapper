// The parts of liblmdb this tool needs, hand-bound.
//
// lmdb's C API is small and stable (0.9.x); a package would have to be packaged
// for Guix before it could be used inside a ravanan job, and this way the only
// runtime dependency is the same liblmdb the rest of the deployment already
// links. The struct shapes below follow lmdb.h -- MDB_val is the only one whose
// layout matters on our side.
package geno2mdb

import "core:fmt"
import "core:os"

foreign import lmdb "system:lmdb"

MDB_env :: struct {}
MDB_txn :: struct {}
MDB_dbi :: distinct u32

MDB_val :: struct {
	size: uint,
	data: rawptr,
}

// lmdb.h: MDB_NOSUBDIR 0x4000, MDB_RDONLY 0x20000, MDB_CREATE 0x40000.
ENV_NOSUBDIR :: u32(0x4000)
DBI_CREATE :: u32(0x40000)

foreign lmdb {
	@(link_name = "mdb_env_create")
	env_create :: proc(env: ^^MDB_env) -> i32 ---

	@(link_name = "mdb_env_set_maxdbs")
	env_set_maxdbs :: proc(env: ^MDB_env, dbs: MDB_dbi) -> i32 ---

	@(link_name = "mdb_env_set_mapsize")
	env_set_mapsize :: proc(env: ^MDB_env, size: uint) -> i32 ---

	@(link_name = "mdb_env_open")
	env_open :: proc(env: ^MDB_env, path: cstring, flags: u32, mode: u32) -> i32 ---

	@(link_name = "mdb_txn_begin")
	txn_begin :: proc(env: ^MDB_env, parent: ^MDB_txn, flags: u32, txn: ^^MDB_txn) -> i32 ---

	@(link_name = "mdb_txn_commit")
	txn_commit :: proc(txn: ^MDB_txn) -> i32 ---

	@(link_name = "mdb_txn_abort")
	txn_abort :: proc(txn: ^MDB_txn) ---

	@(link_name = "mdb_dbi_open")
	dbi_open :: proc(txn: ^MDB_txn, name: cstring, flags: u32, dbi: ^MDB_dbi) -> i32 ---

	@(link_name = "mdb_put")
	put :: proc(txn: ^MDB_txn, dbi: MDB_dbi, key: ^MDB_val, data: ^MDB_val, flags: u32) -> i32 ---

	@(link_name = "mdb_env_sync")
	env_sync :: proc(env: ^MDB_env, force: i32) -> i32 ---

	@(link_name = "mdb_env_close")
	env_close :: proc(env: ^MDB_env) ---

	@(link_name = "mdb_strerror")
	strerror :: proc(err: i32) -> cstring ---
}

// die prints lmdb's own explanation and exits: far more useful than the
// numeric code, which says nothing on its own.
die :: proc(what: string, code: i32) -> ! {
	fmt.eprintfln("geno2mdb: %s: %s (mdb error %d)", what, string(strerror(code)), code)
	os.exit(1)
}
