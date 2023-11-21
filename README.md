# CoACD for WebAssembly (WASM)

This project removes a bunch of dependencies, allowing CoACD to be compiled easily for WebAssembly.

* Manifold preprocessing is not supported.
* Multi-threading is not supported.

## Compile from source (Windows)

### (1) Clone the code

```
git clone --recurse-submodules https://github.com/kermado/CoACD.git
```

### (2) Dependencies
Install dependencies: [emscripten](https://emscripten.org/), [MinGW 64](https://winlibs.com/)

### (3) Compile

```
cd CoACD \
&& mkdir build \
&& cd build \
&& emcmake cmake .. -DCMAKE_BUILD_TYPE=Release -DWASM=1 -DPARALLEL=1 \
&& make -j
```
