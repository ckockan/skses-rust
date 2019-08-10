FTXSGX_RUNNER:=ftxsgx-runner
FTXSGX_EFL2SGXS:=ftxsgx-elf2sgxs
ENCLAVE_CARGO:=cargo
ENCLAVE_CARGO_PARAMS:=build --release

ENCLAVE_FILE:=build/enclave.sgxs
LAUNCHER=launcher/target/release/launcher

CONFIG_FILE?=config
ENCLAVE_MODE?=1
SMALL_TABLE?=0

ifeq ($(SMALL_TABLE), 1)
    ENCLAVE_CARGO_PARAMS+=--features small-table
endif

TARGET=$(LAUNCHER)
LAUNCH=$(LAUNCHER) $(DATA_DIR)
ifeq ($(ENCLAVE_MODE), 1)
    TARGET+=$(ENCLAVE_FILE)
    ENCLAVE_EFL:=enclave/target/x86_64-fortanix-unknown-sgx/release/enclave
    ENCLAVE_CARGO+=+nightly
    ENCLAVE_CARGO_PARAMS+=--target=x86_64-fortanix-unknown-sgx
    include $(CONFIG_FILE) 
    LAUNCH+=$(FTXSGX_RUNNER) $(ENCLAVE_FILE)
else
    ENCLAVE_EFL:=enclave/target/release/enclave
    TARGET+=$(ENCLAVE_EFL)
    LAUNCH+=$(ENCLAVE_EFL)
endif

.PHONY: build run clean $(ENCLAVE_EFL) $(LAUNCHER)

build: $(TARGET) 
	cd enclave; \
	$(ENCLAVE_CARGO) $(ENCLAVE_CARGO_PARAMS) \

run: build $(DATA_DIR)
	$(if $(DATA_DIR),, \
	    $(error ERROR: DATA_DIR is not set))
	$(LAUNCH)

$(ENCLAVE_FILE): $(ENCLAVE_EFL) $(CONFIG_FILE)
	mkdir -p build
	$(FTXSGX_EFL2SGXS) $(ENCLAVE_EFL) --heap-size $(HEAP) --stack-size $(STACK) \
	    -t $(THREADS) -o $(ENCLAVE_FILE) 

$(ENCLAVE_EFL):
	cd enclave; \
	$(ENCLAVE_CARGO) $(ENCLAVE_CARGO_PARAMS) 

$(LAUNCHER):
	cd launcher; \
	cargo build --release

clean:
	rm -rf build
	cd enclave; cargo clean
	cd launcher; cargo clean
