# skses-rust
## Dependencies
- [Rust](https://www.rust-lang.org/tools/install)
- [Fortanix EDP](https://edp.fortanix.com/docs/installation/guide/)
## Build
To build in the SGX Enclave mode:
```bash
make
```
To build in the non-SGX Enclave mode:
```bash
make ENCLAVE_MODE=0
```
## Run
To run in the SGX Enclave mode:
```bash
make run DATA_DIR=<directory_of_data_files>
```
To run in the non-SGX Enclave mode:
```bash
make run DATA_DIR=<directory_of_data_files> ENCLAVE_MODE=0
```
