#!/bin/bash
# Source this to set up the likegt build environment using Guix libraries.
#
#   source ./env.sh
#   cargo build
#   cargo build --release
#   cargo install --path .
#   cargo test
#
# Why: This system (Debian Buster) has glibc 2.28 and cmake 3.13, but likegt
# needs a modern Rust/C++ build stack and the lib_wfa2 dependency needs C++17.
# Guix provides the newer toolchain and libraries. To avoid mixing glibc
# versions, we run cargo and rustc under Guix's dynamic linker so everything
# uses glibc 2.35 consistently.
#
# Prerequisites (one-time):
#   guix install jemalloc

_likegt_env_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

# --- Guix store paths ---
_guix_gcc="/gnu/store/rfm800pq3q2midj29a4xlikdzjp1ps2l-gcc-toolchain-12.3.0"
_guix_gcc_lib="/gnu/store/gmv6n5vy5qcsn71pkapg2hnknyn1p7g3-gcc-12.3.0-lib/lib"
_guix_glibc="/gnu/store/gsjczqir1wbz8p770zndrpw4rnppmxi3-glibc-2.35/lib"
_guix_cmake="/gnu/store/r0g7ygnvmgmlv13375fkv2dn6r694k11-cmake-3.25.1/bin"
_guix_pkgconfig="/gnu/store/jz5dwdxq4di29cd0rjjzkw356dhkzjil-pkg-config-0.29.2/bin"
_guix_profile="/gnu/store/2r5ryq2ibvy44jkz9diar0fvf5cm7q4c-profile"
_guix_ldso="$_guix_glibc/ld-linux-x86-64.so.2"

# --- Rust toolchain ---
_rust_toolchain="$HOME/.rustup/toolchains/stable-x86_64-unknown-linux-gnu"
_lib_path="$_guix_glibc:$_guix_gcc_lib:$_guix_profile/lib:$_rust_toolchain/lib"

# --- Create cargo/rustc wrappers that run under Guix's glibc ---
_wrapper_dir="$_likegt_env_dir/.guix/bin"
mkdir -p "$_wrapper_dir"

cat > "$_wrapper_dir/cargo" << SCRIPT
#!/bin/bash
exec $_guix_ldso --library-path "$_lib_path" $_rust_toolchain/bin/cargo "\$@"
SCRIPT

cat > "$_wrapper_dir/rustc" << SCRIPT
#!/bin/bash
exec $_guix_ldso --library-path "$_lib_path" $_rust_toolchain/bin/rustc "\$@"
SCRIPT

cat > "$_wrapper_dir/rustdoc" << SCRIPT
#!/bin/bash
exec $_guix_ldso --library-path "$_lib_path" $_rust_toolchain/bin/rustdoc "\$@"
SCRIPT

chmod +x "$_wrapper_dir/cargo" "$_wrapper_dir/rustc" "$_wrapper_dir/rustdoc"

# --- Environment ---
export PATH="$_wrapper_dir:$_guix_gcc/bin:$_guix_cmake:$_guix_pkgconfig:$_guix_profile/bin:$PATH"
export CC="$_guix_gcc/bin/gcc"
export CXX="$_guix_gcc/bin/g++"
export LIBRARY_PATH="$_likegt_env_dir/.guix-stubs:$_guix_gcc/lib:$_guix_profile/lib"
export C_INCLUDE_PATH="$_guix_profile/include"
export CPLUS_INCLUDE_PATH="$_guix_profile/include"
export PKG_CONFIG_PATH="$_guix_profile/lib/pkgconfig"
export LIBCLANG_PATH="$_guix_profile/lib"
export CXXFLAGS="-include limits -include cstdint"
# Suppress locale warnings from /bin/bash spawned by build scripts:
# system bash (glibc 2.28) can't read Guix's 2.35 locale data.
export LC_ALL=C
export LANG=C
unset LC_CTYPE
# Use /scratch for compiler temporaries so we don't fill /tmp on the root fs.
export TMPDIR=/scratch
export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER="$_guix_gcc/bin/gcc"
export RUSTC="$_wrapper_dir/rustc"
export RUSTDOC="$_wrapper_dir/rustdoc"

echo "likegt build env ready - Guix GCC 12.3 + glibc 2.35 (cmake $(cmake --version 2>/dev/null | head -1 | awk '{print $3}'))"

# cleanup temp vars
unset _likegt_env_dir _guix_gcc _guix_gcc_lib _guix_glibc _guix_cmake
unset _guix_pkgconfig _guix_profile _guix_ldso _rust_toolchain _lib_path _wrapper_dir

# librt stub for linking
mkdir -p /moosefs/raid5/salehi/crispr2/crispr-progress/crispr-bwt/.guix-stubs
echo "GROUP ( /gnu/store/10krix03rl5hqjv2c0qmj44ic9bgd8rc-gcc-toolchain-13.3.0/lib/librt.so.1 )" > /moosefs/raid5/salehi/crispr2/crispr-progress/crispr-bwt/.guix-stubs/librt.so
export LIBRARY_PATH="/moosefs/raid5/salehi/crispr2/crispr-progress/crispr-bwt/.guix-stubs:/gnu/store/gsjczqir1wbz8p770zndrpw4rnppmxi3-glibc-2.35/lib"
export CARGO_TARGET_DIR="./target"
